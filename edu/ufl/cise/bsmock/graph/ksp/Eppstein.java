package edu.ufl.cise.bsmock.graph.ksp;

import edu.ufl.cise.bsmock.graph.*;
import edu.ufl.cise.bsmock.graph.util.*;
import java.util.*;

/**
 * Eppstein's algorithm for computing the K shortest paths between two nodes in a graph.
 *
 * Copyright (C) 2015  Brandon Smock (dr.brandon.smock@gmail.com, GitHub: bsmock)
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
 * Created by Brandon Smock on October 5, 2015.
 * Last updated by Brandon Smock on December 24, 2015.
 */
public class Eppstein implements KSPAlgorithm {

    public boolean isLoopless() {
        return false;
    }

    public Eppstein() {};

    /**
     * Computes the K shortest paths (allowing cycles) in a graph from node s to node t in graph G using Eppstein's
     * algorithm. ("Finding the k Shortest Paths", Eppstein)
     *
     * Some explanatory notes about how Eppstein's algorithm works:
     * - Start with the shortest path in the graph from s to t, which can be calculated using Dijkstra's algorithm.
     * - The second shortest path must diverge away from this shortest path somewhere along the path at node u, with
     * edge (u,v) known as a "sidetrack" edge. It then follows the shortest path from v to t.
     * - In general, let T represent the shortest path tree rooted at target node t, containing for each node v, an edge
     * (v,w) to the parent of node v (w) along its shortest path to node t.
     * - All other edges (u,v) in the graph are "sidetrack" edges. In other words, this includes any edge (u,v) for
     * which node v is not on the shortest path from u to t.
     * - All paths from s to t can be uniquely represented by the sequence of sidetrack edges that appear in the path.
     * - All non-sidetrack edges in the path are represented implicitly in the shortest path tree, T.
     * - All paths from s to t can be represented in a tree, H, where the root node has no sidetrack edges, its children
     * are all of the possible one sidetrack paths, the children of each of these sidetracks are the paths with a
     * second sidetrack edge, and so on.
     * - Each path from s to t corresponds to a path from the root of H to some descendant node in H.
     * - The tree is a heap, as each additional sidetrack creates a new path whose cost is either the same as or greater
     * than the cost of its parent path.
     * - One possible way to find the k shortest paths, then, is to maintain a priority queue of candidate paths where
     * once a path is pulled from the priority queue, its children in the heap are added to the priority queue.
     * - Two observations about this heap:
     *   - Each node in H has at most O(E) children, where E is the number of edges in the graph.
     *   - If G has looped paths, then this heap is infinite.
     * - Eppstein's algorithm is basically a scheme for re-organizing this heap so that each heap node has at most 4
     * children, which is O(1) instead of O(E) - performing potentially much less work for each path found.
     * - This re-organized heap R has the following form, in terms of the original heap H:
     *   - For each heap node in H, select its best child, B.
     *   - This child is known as a "cross edge" child.
     *   - Child B has N siblings in H.
     *   - These N siblings are removed from its parent in heap H and are made instead to be descendants of B in heap R.
     *   - The best child, B, has up to three of its siblings placed as direct children (along with its own best "cross
     *   edge" child in the original heap H, this yields a total of four possible children in R).
     *   - Among its sibling-children, one of its children in R is the root of a new binary heap containing all of the
     *   siblings of B in H that are sidetrack edges directed away from the same node in G.
     *   - All of the remaining siblings of B in H are placed in a binary heap with B as the root in R.
     *   - Because this is a binary heap, B has up to two children of this type.
     *   - Thus, B has at most four children in R. Any node in R that is a "cross-edge" child has up to four children,
     *   whereas other nodes fall inside binary heaps and are limited to having at most two children.
     *
     * @param graph         the graph on which to compute the K shortest paths from s to t
     * @param sourceLabel   the starting node for all of the paths
     * @param targetLabel   the ending node for all of the paths
     * @param K             the number of shortest paths to compute
     * @return              a list of the K shortest paths from s to t, ordered from shortest to longest
     */
    public List<Path> ksp(Graph graph, String sourceLabel, String targetLabel, int K) {
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

        /* Make indexes to give fast access to these heaps later */
        // Heap H_out(v) for every node v
        HashMap<String,EppsteinHeap> nodeHeaps = new HashMap<String, EppsteinHeap>(graph.numNodes());
        HashMap<String,EppsteinHeap> edgeHeaps = new HashMap<String, EppsteinHeap>(graph.numEdges());
        // Heap H_T(v) for every node v
        HashMap<String,EppsteinHeap> outrootHeaps = new HashMap<String, EppsteinHeap>();

        /* COMPUTE EPPSTEIN HEAP, Part 1: Compute sub-heap H_out(v) for each node v.
            -- H_out(v) is a heap of all of the outgoing sidetrack edges of v. */
        for (String nodeLabel : graph.getNodes().keySet()) {
            computeOutHeap(nodeLabel, graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps);
        }

        /* COMPUTE EPPSTEIN HEAP, Part 2: Compute sub-heap H_T(v) for each node v.
            -- H_T(v) is a heap of all of the "best" sidetrack edges for each node on the shortest path from v to T.
            -- H_T(v) is computed by adding the lowest cost sidetrack edge of v to heap H_T(nextT(v)),
            where nextT(v) is the parent node of node v in the shortest path tree rooted at the target node T.
            -- Therefore, can compute H_T(v) recursively. But instead of a top-down recursion, we will compute each
            H_T(v) bottom-up, starting with the root node of the tree, T.
            -- To facilitate bottom-up computation, reverse the edges of the shortest path tree so each node points to
            its children instead of its parent. */
        Graph reversedSPT = new Graph();
        for (DijkstraNode node: tree.getNodes().values()) {
            reversedSPT.addEdge(node.getParent(),node.getLabel(),graph.getNode(node.getLabel()).getNeighbors().get(node.getParent()));
        }

        /* Use a depth-first search from node T to perform the bottom-up computation, computing each H_T(v) given
            H_T(nextT(v)). */
        // Create the initial (empty) heap for the root node T to build from.
        EppsteinArrayHeap rootArrayHeap = new EppsteinArrayHeap();
        // Perform the DFS (recursively initiating additional depth-first searches for each child)
        recursiveOutrootHeaps(targetLabel, rootArrayHeap, nodeHeaps, outrootHeaps, reversedSPT);

        // Create a virtual/dummy heap that is the root of the overall Eppstein heap. It represents the best path from
        // the source node to the target node, which does not involve any sidetrack edges.
        EppsteinHeap hg = new EppsteinHeap(new Edge(sourceLabel,sourceLabel,0));

        // Initialize the containers for the candidate k shortest paths and the actual found k shortest paths
        ArrayList<Path> ksp = new ArrayList<Path>();
        PriorityQueue<EppsteinPath> pathPQ = new PriorityQueue<EppsteinPath>();

        // Place root heap in priority queue
        pathPQ.add(new EppsteinPath(hg, -1, tree.getNodes().get(sourceLabel).getDist()));

        /* Pop k times from the priority queue to determine the k shortest paths */
        for (int i = 0; i < K && pathPQ.size() > 0; i++) {
            /* Get the next shortest path, which is implicitly represented as:
                1) Some shorter path, p, from s (source) to t (target)
                2) A sidetrack edge which branches off of path p at node u, and points to node v
                3) The shortest path in the shortest path tree from node v to t */
            EppsteinPath kpathImplicit = pathPQ.poll();

            // Convert from the implicit path representation to the explicit path representation
            Path kpath = kpathImplicit.explicitPath(ksp, tree);

            // Add explicit path to the list of K shortest paths
            ksp.add(kpath);

            // Push the (up to 3) children of this path within the Eppstein heap onto the priority queue
            addExplicitChildrenToQueue(kpathImplicit, ksp, pathPQ);

            /* Check for the existence of a potential fourth child, known as a "cross edge", to push onto the queue.
                This heap edge/child does not need to be explicitly represented in the Eppstein heap because it is easy
                to check for its existence. */
            addCrossEdgeChildToQueue(outrootHeaps, kpathImplicit, i, ksp, pathPQ);
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
            // Check to see if the target node is reachable from the outgoing vertex of the current edge,
            // and check to see if the current edge is a sidetrack edge. If so, calculate its sidetrack cost.
            String tp = tree.getParentOf(edge.getFromNode());
            if (tp == null || !tp.equals(edge.getToNode())) {
                double sidetrackEdgeCost = edge.getWeight() + tree.getNodes().get(edge.getToNode()).getDist() - tree.getNodes().get(edge.getFromNode()).getDist();
                sidetrackEdgeCostMap.put(edge.getFromNode() + "," + edge.getToNode(), sidetrackEdgeCost);
            }
        }

        return sidetrackEdgeCostMap;
    }

    /**
     * Compute sub-heap H_out(v) for node v.
     *
     * @param nodeLabel             node v
     * @param graph                 the graph, G, on which to compute the K shortest paths from s to t
     * @param sidetrackEdgeCostMap  the cost of each sidetrack edge in G
     * @param nodeHeaps             an index/hash table of heap H_out(v), for each node v in the graph
     * @param edgeHeaps             an index/hash table of heaps H_out(v), but indexed by sidetrack edge
     */
    protected static void computeOutHeap(String nodeLabel, Graph graph, HashMap<String,Double> sidetrackEdgeCostMap, HashMap<String,EppsteinHeap> nodeHeaps, HashMap<String,EppsteinHeap> edgeHeaps) {
        Node node = graph.getNode(nodeLabel);
        // This list holds the 2nd through last sidetrack edges, ordered by sidetrack cost
        ArrayList<Edge> sidetrackEdges = new ArrayList<Edge>();
        Edge bestSidetrack = null;
        double minSidetrackCost = Double.MAX_VALUE;
        // Iterate over the outgoing edges of v
        for (String neighbor : node.getAdjacencyList()) {
            String edgeLabel = nodeLabel+","+neighbor;
            // Check to see if the current edge is a sidetrack edge
            if (sidetrackEdgeCostMap.containsKey(edgeLabel)) {
                double sidetrackEdgeCost = sidetrackEdgeCostMap.get(edgeLabel);
                // Check to see if the current sidetrack edge has the lowest cost discovered so far for node v
                if (sidetrackEdgeCost < minSidetrackCost) {
                    // If there was a previously-known best sidetrack edge, add it to the list of non-best
                    // sidetrack edges
                    if (bestSidetrack != null) {
                        sidetrackEdges.add(bestSidetrack);
                    }
                    // Set the new best (lowest cost) sidetrack edge to be the current one
                    bestSidetrack = new Edge(nodeLabel, neighbor, node.getNeighbors().get(neighbor));
                    minSidetrackCost = sidetrackEdgeCost;
                }
                // If current sidetrack edge is not the one with the lowest cost, add it to the list of non-best
                // sidetrack edges
                else {
                    sidetrackEdges.add(new Edge(nodeLabel, neighbor, node.getNeighbors().get(neighbor)));
                }
            }
        }
        // If v was found to have at least one outgoing sidetrack edge...
        if (bestSidetrack != null) {
            // ...make a heap of the outgoing sidetrack edges of v, with the lowest-cost sidetrack edge put as the
            // root
            EppsteinHeap bestSidetrackHeap = new EppsteinHeap(bestSidetrack,sidetrackEdgeCostMap.get(bestSidetrack.getFromNode()+","+bestSidetrack.getToNode()));

            // Make another heap (a binary heap) out of the rest of the sidetrack edges of v
            EppsteinArrayHeap arrayHeap = new EppsteinArrayHeap();
            if (sidetrackEdges.size() > 0) {
                bestSidetrackHeap.setNumOtherSidetracks(bestSidetrackHeap.getNumOtherSidetracks()+1);
                for (Edge edge : sidetrackEdges) {
                    EppsteinHeap sidetrackHeap = new EppsteinHeap(edge,sidetrackEdgeCostMap.get(edge.getFromNode()+","+edge.getToNode()));
                    edgeHeaps.put(edge.getFromNode()+","+edge.getToNode(), sidetrackHeap);
                    arrayHeap.add(sidetrackHeap);
                }

                // Add the binary heap of 2nd-through-last lowest cost sidetrack edges as a child (the only child)
                // of the lowest-cost sidetrack edge, forming the overall heap H_out(v)
                bestSidetrackHeap.addChild(arrayHeap.toEppsteinHeap());
            }

            // Index H_out(v) by node v, for easy access later
            nodeHeaps.put(nodeLabel, bestSidetrackHeap);
            // Index H_out(v) by its lowest cost sidetrack edge, for easy access later
            edgeHeaps.put(bestSidetrack.getFromNode()+","+bestSidetrack.getToNode(), bestSidetrackHeap);
        }
    }

    /**
     * Push the (up to 3) children (within the Eppstein heap) of the given (kth) path, onto the priority queue.
     *
     * @param kpathImplicit     implicit representation of the (kth) path
     * @param ksp               list of shortest paths found so far
     * @param pathPQ            priority queue of candidate paths
     */
    protected static void addExplicitChildrenToQueue(EppsteinPath kpathImplicit, ArrayList<Path> ksp, PriorityQueue<EppsteinPath> pathPQ) {
        double kpathCost = kpathImplicit.getCost();
        for(EppsteinHeap childHeap : kpathImplicit.getHeap().getChildren()) {
            // Get the index of the previous shorter path off of which this candidate sidetracks/branches
            int prefPath = kpathImplicit.getPrefPath();

            // Calculate the path cost of the new child/candidate
            Double candidateCost = ksp.get(prefPath).getTotalCost() + childHeap.getSidetrackCost();

            // Add the child/candidate to the priority queue
            EppsteinPath candidate = new EppsteinPath(childHeap, prefPath, candidateCost);
            pathPQ.add(candidate);
        }
    }

    /**
     *
     * @param outrootHeaps      an index of heaps H_T(v) for each node v
     * @param kpathImplicit     implicit representation of the (kth) path
     * @param prefPath          the index k of the path off which this cross-edge child sidetracks
     * @param ksp               list of shortest paths found so far
     * @param pathPQ            priority queue of candidate paths
     */
    protected static void addCrossEdgeChildToQueue(HashMap<String,EppsteinHeap> outrootHeaps, EppsteinPath kpathImplicit, int prefPath, ArrayList<Path> ksp, PriorityQueue<EppsteinPath> pathPQ) {
        if (outrootHeaps.containsKey(kpathImplicit.getHeap().getSidetrack().getToNode())) {
            EppsteinHeap childHeap = outrootHeaps.get(kpathImplicit.getHeap().getSidetrack().getToNode());

            // Calculate the path cost of the new child/candidate
            Double candidateCost = ksp.get(prefPath).getTotalCost() + childHeap.getSidetrackCost();

            // Add the child/candidate to the priority queue
            EppsteinPath candidate = new EppsteinPath(childHeap, prefPath, candidateCost);
            pathPQ.add(candidate);
        }
    }

    /**
     * Generate sub-heap H_T(v) for node v and its children in the shortest path tree T, using a recursive depth-first
     * search over the transpose graph, T', of the shortest path tree, T.
     *
     * The transpose graph is necessary because tree T is represented with pointers to parents instead of with pointers
     * to children.
     *
     * @param nodeLabel         node v
     * @param currentArrayHeap  the heap of v's parent in the shortest path tree; H_T(nextT(v))
     * @param nodeHeaps         an index/hash table of heap H_out(v), for each node v in the graph
     * @param outrootHeaps      an index/hash table of heap H_T(v), for each node v in the graph
     * @param reversedSPT       the transpose graph, T', of the shortest path tree, T
     */
    protected static void recursiveOutrootHeaps(String nodeLabel, EppsteinArrayHeap currentArrayHeap, HashMap<String,EppsteinHeap> nodeHeaps, HashMap<String,EppsteinHeap> outrootHeaps, Graph reversedSPT) {
        // Get H_out(v)
        EppsteinHeap sidetrackHeap = nodeHeaps.get(nodeLabel);

        // Check to see if node v (nodeLabel) has a sidetrack edge
        if (sidetrackHeap != null) {
            // If so, need to add its best sidetrack edge to H_T(nextT(v))
            // H_T(nextT(v)) is in variable currentArrayHeap, which was passed onto v when this function was called by
            // nextT(v)

            // The goal of Eppstein's algorithm is to re-use heap structures where possible
            // When adding the best sidetrack edge of v to heap H_T(nextT(v)) to form H_T(v), this means that when
            // the best sidetrack of v is added to H_T(nextT(v)), we need to make a new copy of all of the heap nodes
            // on the path from the position where the sidetrack is added in the heap, to the root of the heap. All
            // other heap nodes have the same structure in H_T(v) as in H_T(nextT(v)) and can simply be pointed to.

            // Give H_T(v) a new set of pointers to the sub-heap nodes in H_T(nextT(v))
            currentArrayHeap = currentArrayHeap.clone();

            // Add the best sidetrack edge of v to heap H_T(v) (place the sidetrack edge in the next unoccupied space in
            // the heap and bubble it up the tree to its rightful position) and make all new copies of the sub-heap
            // nodes that fall on the path where the new edge is added in the heap, to the root of the heap
            currentArrayHeap.addOutroot(sidetrackHeap);
        }

        // Convert from an array representation of the heap (EppsteinArrayHeap), which is convenient for accessing and
        // manipulating a binary heap, to a pointer representation (EppsteinHeap), which is consistent with the overall
        // heap, which is not strictly a binary heap since some nodes can have up to 4 children.
        EppsteinHeap currentHeap = currentArrayHeap.toEppsteinHeap2();

        // If v has any children (so H_T(v) exists), index heap H_T(v) in a list of heaps for fast access later
        if (currentHeap != null) {
            outrootHeaps.put(nodeLabel, currentHeap);
        }

        // Continue the depth-first search (recursively initiating additional depth-first searches for each child)
        for (String neighbor : reversedSPT.getNode(nodeLabel).getNeighbors().keySet()) {
            recursiveOutrootHeaps(neighbor, currentArrayHeap, nodeHeaps, outrootHeaps, reversedSPT);
        }
    }
}

/**
 * A pointer representation of an N-ary heap with data structures that aid in representing the heap constructed by
 * Eppstein's algorithm.
 * Does not contain functions for adding/removing elements while maintaining the heap property.
 * The children added to the root are sub-heaps whose elements are guaranteed to have a greater cost than the root
 * element of the heap.
 */
class EppsteinHeap {
    private Edge sidetrack; // the sidetrack edge (u,v) associated with the root of this heap or sub-heap
    private double sidetrackCost = 0.0;
    private ArrayList<EppsteinHeap> children; // supports N children but Eppstein is limited to 4
    private int numOtherSidetracks = 0; // number of elements of H_out(u) - 1

    public EppsteinHeap(Edge sidetrack) {
        this.sidetrack = sidetrack;
        this.children = new ArrayList<EppsteinHeap>();
    }

    public EppsteinHeap(Edge sidetrack, Double sidetrackCost) {
        this.sidetrack = sidetrack;
        this.sidetrackCost = sidetrackCost;
        this.children = new ArrayList<EppsteinHeap>();
    }

    public EppsteinHeap(Edge sidetrack, double sidetrackCost, ArrayList<EppsteinHeap> children, int numOtherSidetracks) { //, boolean bestChild, int copy) {
        this.sidetrack = sidetrack;
        this.sidetrackCost = sidetrackCost;
        this.children = children;
        this.numOtherSidetracks = numOtherSidetracks;
    }

    public Edge getSidetrack() {
        return sidetrack;
    }

    public void setSidetrack(Edge sidetrack) {
        this.sidetrack = sidetrack;
    }

    public double getSidetrackCost() {
        return sidetrackCost;
    }

    public void setSidetrackCost(double sidetrackCost) {
        this.sidetrackCost = sidetrackCost;
    }

    public ArrayList<EppsteinHeap> getChildren() {
        return children;
    }

    public void setChildren(ArrayList<EppsteinHeap> children) {
        this.children = children;
    }

    public void addChild(EppsteinHeap child) {
        this.children.add(child);
    }

    public int getNumOtherSidetracks() {
        return numOtherSidetracks;
    }

    public void setNumOtherSidetracks(int numOtherSidetracks) {
        this.numOtherSidetracks = numOtherSidetracks;
    }

    public EppsteinHeap clone() {
        ArrayList<EppsteinHeap> children_clone = new ArrayList<EppsteinHeap>(children.size());
        for (EppsteinHeap eh: children) {
            children_clone.add(eh);
        }

        return new EppsteinHeap(sidetrack, sidetrackCost, children_clone, numOtherSidetracks);
    }
}

/**
 * An array representation of a binary heap with additional functions specialized for implementing Eppstein's algorithm.
 */
class EppsteinArrayHeap {
    private ArrayList<EppsteinHeap> arrayHeap;

    public EppsteinArrayHeap() {
        arrayHeap = new ArrayList<EppsteinHeap>(0);
    }

    public ArrayList<EppsteinHeap> getArrayHeap() {
        return arrayHeap;
    }

    public void setArrayHeap(ArrayList<EppsteinHeap> arrayHeap) {
        this.arrayHeap = arrayHeap;
    }

    public int getParentIndex(int i) {
        return (i-1)/2;
    }

    public void add(EppsteinHeap h) {
        arrayHeap.add(h);
        bubbleUp(arrayHeap.size()-1);
    }

    public void addOutroot(EppsteinHeap h) {
        int current = arrayHeap.size();

        while (current > 0) {
            int parent = getParentIndex(current);
            EppsteinHeap newHeap = arrayHeap.get(parent).clone();
            arrayHeap.set(parent,newHeap);
            current = parent;
        }
        arrayHeap.add(h);
        bubbleUp(arrayHeap.size() - 1);
    }

    private void bubbleUp(int current) {
        if (current == 0)
            return;

        int parent = getParentIndex(current);
        if (arrayHeap.get(current).getSidetrackCost() >= arrayHeap.get(parent).getSidetrackCost())
            return;

        EppsteinHeap temp = arrayHeap.get(current);
        arrayHeap.set(current, arrayHeap.get(parent));
        arrayHeap.set(parent, temp);
        bubbleUp(parent);
    }

    // Convert from an array representation of a binary heap to a pointer representation of a binary heap, which can fit
    // consistently within an overall N-ary heap
    public EppsteinHeap toEppsteinHeap() {
        int heapsize = arrayHeap.size();
        if (arrayHeap.size() == 0)
            return null;

        EppsteinHeap eh = arrayHeap.get(0);
        for (int i = 1; i < arrayHeap.size(); i++) {
            EppsteinHeap h = arrayHeap.get(i);
            arrayHeap.get(getParentIndex(i)).addChild(h);
        }

        return eh;
    }

    // Convert from an array representation of a binary heap to a pointer representation of a binary heap, which can fit
    // consistently within an overall non-binary heap.
    public EppsteinHeap toEppsteinHeap2() {
        int current = arrayHeap.size()-1;
        if (current == -1)
            return null;

        while (current >= 0) {
            EppsteinHeap childHeap = arrayHeap.get(current);
            while (childHeap.getChildren().size() > childHeap.getNumOtherSidetracks()) {
                childHeap.getChildren().remove(childHeap.getChildren().size()-1);
            }

            int child1 = current * 2 + 1;
            int child2 = current * 2 + 2;

            if (child1 < arrayHeap.size()) {
                arrayHeap.get(current).addChild(arrayHeap.get(child1));
            }
            if (child2 < arrayHeap.size()) {
                arrayHeap.get(current).addChild(arrayHeap.get(child2));
            }
            if (current > 0) {
                current = getParentIndex(current);
            }
            else {
                current = -1;
            }
        }

        return arrayHeap.get(0);
    }

    public EppsteinArrayHeap clone() {
        EppsteinArrayHeap clonedArrayHeap = new EppsteinArrayHeap();
        for (EppsteinHeap heap: arrayHeap) {
            clonedArrayHeap.add(heap);
        }

        return clonedArrayHeap;
    }
}

/**
 * Data structure for representing a source-target path implicitly inside the priority queue of candidate k shortest
 * paths during the execution of Eppstein's algorithm.
 */
class EppsteinPath implements Comparable<EppsteinPath> {
    EppsteinHeap heap; // pointer to the heap node and last sidetrack edge in this candidate path
    int prefPath; // index of the shorter path that this path sidetracks from
    Double cost; // the total cost of the path

    public EppsteinPath(EppsteinHeap heap, int prefPath, Double cost) {
        this.heap = heap;
        this.prefPath = prefPath;
        this.cost = cost;
    }

    public int getPrefPath() {
        return prefPath;
    }

    public void setPrefPath(int prefPath) {
        this.prefPath = prefPath;
    }

    public EppsteinHeap getHeap() {
        return heap;
    }

    public void setHeap(EppsteinHeap heap) {
        this.heap = heap;
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
        if (prefPath >= 0) {
            // Get the explicit representation of the shorter parent path that this path sidetracks from
            Path explicitPrefPath = ksp.get(prefPath);

            // 1a) Identify the s-u portion of the path
            // Identify and add the segment of the parent path up until the point where the current path sidetracks off
            // of it.
            // In other words, if (u,v) is the sidetrack edge of the current path off of the parent path, look for the
            // last instance of node u in the parent path.
            LinkedList<Edge> edges = explicitPrefPath.getEdges();
            int lastEdgeNum = -1;
            Edge heapSidetrack = heap.getSidetrack();
            for (int i = edges.size()-1; i >= 0; i--) {
                Edge currentEdge = edges.get(i);
                if (currentEdge.getToNode().equals(heapSidetrack.getFromNode())) {
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
            explicitPath.add(heap.getSidetrack());
        }

        // 3) Add the v-t portion of the path
        // Add the shortest path from v (either the source node, or the incoming node of the sidetrack edge associated
        // with the current path) to the explicit path representation
        String current = heap.getSidetrack().getToNode();
        while (!current.equals(tree.getRoot())) {
            String next = tree.getParentOf(current);
            Double edgeWeight = tree.getNodes().get(current).getDist() - tree.getNodes().get(next).getDist();
            explicitPath.add(new Edge(current, next, edgeWeight));
            current = next;
        }

        return explicitPath;
    }

    public int compareTo(EppsteinPath comparedNode) {
        double cost1 = this.cost;
        double cost2 = comparedNode.getCost();
        if (cost1 == cost2)
            return 0;
        if (cost1 > cost2)
            return 1;
        return -1;
    }

}
