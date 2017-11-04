package edu.ufl.cise.bsmock.graph.ksp;

import edu.ufl.cise.bsmock.graph.*;
import edu.ufl.cise.bsmock.graph.util.*;
import java.util.*;

/**
 * This is "simplified Eppstein's algorithm", a version of Eppstein's algorithm for computing the K shortest paths
 * between two nodes in a graph, which has been simplified by Brandon Smock.
 *
 * This simplified version eliminates the pre-processing and additional storage required by Eppstein's
 * algorithm. It is faster for smaller values of K due to the reduction in pre-processing, but it incurs more
 * computation per path, so for large values of K, Eppstein's algorithm is faster.
 *
 * It is primarily intended to be a useful first step for understanding how Eppstein's algorithm works and how it can be
 * implemented.
 *
 * Copyright (C) 2016  Brandon Smock (dr.brandon.smock@gmail.com, GitHub: bsmock)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created by Brandon Smock on April 6, 2016.
 * Last updated by Brandon Smock on April 7, 2016.
 */
public class SimpleEppstein implements KSPAlgorithm {

    public boolean isLoopless() {
        return false;
    }

    public SimpleEppstein() {};

    /**
     * Computes the K shortest paths (allowing cycles) in a graph from node s to node t in graph G using a simple
     * version of Eppstein's algorithm. ("Finding the k Shortest Paths", Eppstein)
     *
     * @param graph         the graph on which to compute the K shortest paths from s to t
     * @param sourceLabel   the starting node for all of the paths
     * @param targetLabel   the ending node for all of the paths
     * @param K             the number of shortest paths to compute
     */
    public List<Path> ksp(Graph graph, String sourceLabel, String targetLabel, int K) {
        return kspCutoff(graph, sourceLabel, targetLabel, K, Double.MAX_VALUE);
    }

    /**
     * Computes the K shortest paths (allowing cycles) in a graph from node s to node t in graph G using a simplified
     * version of Eppstein's algorithm. ("Finding the k Shortest Paths", Eppstein)
     *
     * See Eppstein.java for some explanatory notes about how Eppstein's algorithm works.
     *
     * - In this simplified version, like Eppstein's algorithm we represent each path using its unique sequence of
     * sidetrack edges.
     * - Like Eppstein's algorithm, there is a min heap which partially orders these paths, and we can use this fact
     * to efficiently search the space of paths to find the K shortest.
     * - Like Eppstein's algorithm, we never fully generate this heap, and we can use a second heap to store the
     * portions of the first heap that we have traversed/generated.
     *
     * DIFFERENCES IN THIS SIMPLIFIED VERSION:
     * - The first heap has O(|E|) children per node, and one of the big contributions of Eppstein was a scheme for
     * re-organizing this heap to have at most 4 children per node
     * - This requires some not insignificant computation at the beginning of the procedure, but then requires much less
     * computation per path, since the search for paths expands at a much smaller rate, which is independent of the size
     * of the graph.
     * - In this simplified version, we do not re-organize the path heap, nor keep track of children in any way.
     * - This saves space and eliminates a lot of pre-processing computation but does require more computation per path.
     * - Thus for small to moderate values of K, this simplified algorithm ought to be faster than Eppstein's algorithm,
     * but for large values of K, Eppstein's algorithm will eventually overtake this algorithm in efficiency.
     *
     * @param graph         the graph on which to compute the K shortest paths from s to t
     * @param sourceLabel   the starting node for all of the paths
     * @param targetLabel   the ending node for all of the paths
     * @param K             the number of shortest paths to compute
     * @param threshold     the maximum cost allowed for a path
     * @return              a list of the K shortest paths from s to t, ordered from shortest to longest
     */
    public List<Path> kspCutoff(Graph graph, String sourceLabel, String targetLabel, int K, double threshold) {
        /* Compute the shortest path tree, T, for the target node (the shortest path from every node in the graph to the
            target) */
        ShortestPathTree tree;
        try {
            tree = Dijkstra.shortestPathTree(graph.transpose(), targetLabel);
        } catch (Exception e) {
            tree = new ShortestPathTree(targetLabel);
        }

        // Compute the set of sidetrack edge costs
        HashMap<String,Double> sidetrackEdgeCostMap = computeSidetrackEdgeCosts(graph, tree);

        // Initialize the containers for the candidate k shortest paths and the actual found k shortest paths
        ArrayList<Path> ksp = new ArrayList<Path>();
        PriorityQueue<ImplicitPath> pathPQ = new PriorityQueue<ImplicitPath>();

        // Place the shortest path in the candidate-path priority queue
        pathPQ.add(new ImplicitPath(new Edge(null,sourceLabel,0),-1, tree.getNodes().get(sourceLabel).getDist()));

        /* Pop k times from the candidate-path priority queue to determine the k shortest paths */
        for (int k = 0; k < K && pathPQ.size() > 0; k++) {
            /* Get the next shortest path, which is implicitly represented as:
                1) A parent, shorter path, p, from s (source) to t (target)
                2) A sidetrack edge which branches off of path p at node u, and points to node v
                3) The shortest path (in the shortest path tree) from node v to t */
            ImplicitPath kpathImplicit = pathPQ.poll();

            // Convert from the implicit path representation to the explicit path representation
            Path kpath = kpathImplicit.explicitPath(ksp, tree);

            /* Optional/added step:
                Stop if this path is above the cost/length threshold (if a threshold exists) */
            if (kpath.getTotalCost() > threshold)
                return ksp;

            // Add explicit path to the list of K shortest paths
            ksp.add(kpath);

            // Push the O(|E|) children of this path within the path heap onto the priority queue as new candidates
            addChildrenToQueue(graph, sidetrackEdgeCostMap, kpathImplicit, k, pathPQ);
        }

        // Return the set of k shortest paths
        return ksp;
    }

    /**
     * Compute the set of sidetrack edge costs.
     *
     * Each sidetrack edge (u,v) is an edge in graph G that does not appear in the shortest path tree, T.
     * For every sidetrack edge (u,v), compute S(u,v) = w(u,v) + d(v) - d(u), where w(u,v) is the cost of edge (u,v);
     * and d(v) is the cost of the shortest path from node v to the target.
     *
     * @param graph     the graph on which to compute the K shortest paths from s to t
     * @param tree      the shortest path tree, T, rooted at the target node, t
     * @return
     */
    protected static HashMap<String,Double> computeSidetrackEdgeCosts(Graph graph, ShortestPathTree tree) {
        HashMap<String, Double> sidetrackEdgeCostMap = new HashMap<String, Double>();
        List<Edge> edgeList = graph.getEdgeList();
        for (Edge edge : edgeList) {
            /* Check to see if the target node is reachable from the outgoing vertex of the current edge,
                and check to see if the current edge is a sidetrack edge. If so, calculate its sidetrack cost.*/
            String tp = tree.getParentOf(edge.getFromNode());
            if (tp == null || !tp.equals(edge.getToNode())) {
                double sidetrackEdgeCost = edge.getWeight() + tree.getNodes().get(edge.getToNode()).getDist() - tree.getNodes().get(edge.getFromNode()).getDist();
                sidetrackEdgeCostMap.put(edge.getFromNode() + "," + edge.getToNode(), sidetrackEdgeCost);
            }
        }

        return sidetrackEdgeCostMap;
    }

    /**
     * Push the children of the given (kth) path, onto the priority queue.
     *
     * @param kpathImplicit     implicit representation of the (kth) path
     * @param sidetrackMap      map container with all of the costs associated with sidetrack edges
     * @param kpathImplicit     implicit representation of the previous/parent (kth) shortest path
     * @param k                 k, the index of the previous/parent shortest path
     * @param pathPQ            priority queue of candidate paths
     */
    protected static void addChildrenToQueue(Graph graph, HashMap<String,Double> sidetrackMap, ImplicitPath kpathImplicit, int k, PriorityQueue<ImplicitPath> pathPQ) {
        double kpathCost = kpathImplicit.getCost();

        /* Each path is represented as a sequence of sidetrack edges.
            Each path's children are all of the paths with the same sequence of sidetrack edges followed by one
            additional sidetrack.
            Therefore, starting from the last sidetrack edge of the parent (kth) path, its children correspond to each
            sidetrack edge reachable by traversing non-sidetrack edges only in the graph.
            These can be found using a depth-first search. */

        // Initialize the stack for the DFS
        Stack<Edge> edgeStack = new Stack<Edge>();

        // Add the neighbors of the last sidetrack edge in the graph to the stack for DFS
        for (Edge outgoingEdge : graph.getNode(kpathImplicit.getSidetrackEdge().getToNode()).getEdges()) {
            edgeStack.push(outgoingEdge);
        }

        // Iterate/execute the DFS
        while (!edgeStack.empty()) {
            Edge poppedEdge = edgeStack.pop();
            String edgeString = poppedEdge.getFromNode() + "," + poppedEdge.getToNode();

            if (sidetrackMap.containsKey(edgeString)) {
                // Base case for DFS: sidetrack edge found
                // Add the child/candidate path represented by the sidetrack edge to the priority queue
                ImplicitPath candidate = new ImplicitPath(poppedEdge, k, kpathCost + sidetrackMap.get(edgeString));
                pathPQ.add(candidate);
            }
            else {
                // Recursive case for DFS: sidetrack edge not reached, keep going (if current node has outgoing
                // edges)
                for (Edge outgoingEdge : graph.getNode(poppedEdge.getToNode()).getEdges()) {
                    edgeStack.push(outgoingEdge);
                }
            }
        }
    }
}

/**
 * Data structure for representing a source-target path implicitly inside the priority queue of candidate k shortest
 * paths during the execution of the simplified version of Eppstein's algorithm.
 */
class ImplicitPath implements Comparable<ImplicitPath> {
    Edge sidetrackEdge; // last sidetrack edge in this candidate path
    int parentPath; // index of the shorter path that this path sidetracks from
    Double cost; // the total cost of the path

    public ImplicitPath(Edge sidetrackEdge, int prefPath, Double cost) {
        this.sidetrackEdge = sidetrackEdge;
        this.parentPath = prefPath;
        this.cost = cost;
    }

    public Edge getSidetrackEdge() {
        return sidetrackEdge;
    }

    public void setSidetrackEdge(Edge sidetrackEdge) {
        this.sidetrackEdge = sidetrackEdge;
    }

    public int getPrefPath() {
        return parentPath;
    }

    public void setPrefPath(int prefPath) {
        this.parentPath = prefPath;
    }

    public Double getCost() {
        return cost;
    }

    public void setCost(Double cost) {
        this.cost = cost;
    }

    // Convert from the implicit representation of the path to an explicit listing of all of the edges in the path
    // There are potentially three pieces to the path:
    // 1) the path from node s (source) to node u in the parent path
    // 2) the sidetrack edge (u,v)
    // 3) the shortest path (in the shortest path tree) from node v to node t (target)
    public Path explicitPath(List<Path> ksp, ShortestPathTree tree) {
        Path explicitPath = new Path();

        // If path is not the shortest path in the graph...
        if (parentPath >= 0) {
            // Get the explicit representation of the shorter parent path that this path sidetracks from
            Path explicitPrefPath = ksp.get(parentPath);

            // 1a) Identify the s-u portion of the path
            // Identify and add the segment of the parent path up until the point where the current path sidetracks off
            // of it.
            // In other words, if (u,v) is the sidetrack edge of the current path off of the parent path, look for the
            // last instance of node u in the parent path.
            LinkedList<Edge> edges = explicitPrefPath.getEdges();
            int lastEdgeNum = -1;
            for (int i = edges.size()-1; i >= 0; i--) {
                Edge currentEdge = edges.get(i);
                if (currentEdge.getToNode().equals(sidetrackEdge.getFromNode())) {
                    lastEdgeNum = i;
                    break;
                }
            }

            // 1b) Add the s-u portion of the path
            // Copy the explicit parent path up to the identified point where the current/child path sidetracks
            explicitPath = new Path();
            for (int i = 0; i <= lastEdgeNum; i++) {
                explicitPath.add(edges.get(i));
            }

            // 2) Add the (u,v) portion of the path
            // Add the last sidetrack edge to the explicit path representation
            explicitPath.add(sidetrackEdge);
        }

        // 3) Add the v-t portion of the path
        // Add the shortest path from v (either the source node, or the incoming node of the sidetrack edge associated
        // with the current path) to the explicit path representation
        String current = sidetrackEdge.getToNode();
        while (!current.equals(tree.getRoot())) {
            String next = tree.getParentOf(current);
            Double edgeWeight = tree.getNodes().get(current).getDist() - tree.getNodes().get(next).getDist();
            explicitPath.add(new Edge(current, next, edgeWeight));
            current = next;
        }

        return explicitPath;
    }

    public int compareTo(ImplicitPath comparedNode) {
        double cost1 = this.cost;
        double cost2 = comparedNode.getCost();
        if (cost1 == cost2)
            return 0;
        if (cost1 > cost2)
            return 1;
        return -1;
    }

}