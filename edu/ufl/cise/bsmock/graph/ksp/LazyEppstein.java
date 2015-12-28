package edu.ufl.cise.bsmock.graph.ksp;

import edu.ufl.cise.bsmock.graph.*;
import edu.ufl.cise.bsmock.graph.util.*;
import java.util.*;

/**
 * Lazy version of Eppstein's algorithm (by Jimenez and Marzal) for computing the K shortest paths
 * between two nodes in a graph.
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
public final class LazyEppstein extends Eppstein implements KSPAlgorithm {

    public LazyEppstein() {};

    /**
     * Computes the K shortest paths (allowing cycles) in a graph from node s to node t in graph G using the lazy
     * version of Eppstein's algorithm. ("A lazy version of Eppstein's K shortest paths algorithm", Jimenez & Marzal)
     *
     * See Eppstein.java for explanatory notes about Eppstein's algorithm.
     *
     * Some explanatory notes about how lazy Eppstein's algorithm works, given Eppstein's algorithm:
     * - Main idea: Instead of building the entire Eppstein heap before finding any of the K shortest paths, build the
     * heap as it is traversed, meaning build it as candidate paths are pulled from the priority queue.
     * - When candidate paths are pulled from the priority queue, their children need to be placed in the queue - this
     * requires that their children be in the heap as children, or be known about some other way.
     * - Each heap node is associated with a sidetrack edge (u,v), and has at most 3 explicit children in the heap, as
     * well as a fourth, implicit "cross-edge" child associated with node v.
     * - For the 3 explicit children to be in the heap, their explicit children must be in the heap, as well, which
     * requires some amount of heap-building in advance.
     * - The lazy heap-building strategy arises from the fact that since the implicit "cross-edge" child is always known
     * about without being built into the heap, nothing below it in the heap needs to be built until it is reached -
     * meaning, until its parent is pulled from the priority queue.
     * - When the parent of a cross-edge child is pulled from the priority queue and the cross-edge child is added to
     * the priority queue, the explicit children of the cross-edge child need to be computed (built into the heap).
     * - This means, if the heaps H_out(v) and H_T(v) associated with the cross-edge child have not been built, build
     * them (this is where all of the heap-building actually happens).
     * - Building H_T(v) requires H_T(nextT(v)) to be built (see Eppstein.java for more info), so this may initiate a
     * recursive build of multiple H_T heaps until a heap is reached that has been built. Thus some portions of the heap
     * that are never actually reached may be built in order to ensure H_T(v) can be built for node v.
     * - This entire lazy heap-building process can be initiated by simply adding the shortest path in the graph to the
     * priority queue - the shortest path has only a "cross-edge" child, so the first parts of the heap will be built
     * when the shortest path is removed from the priority queue.
     *
     * @param graph       the graph on which to compute the K shortest paths from s to t
     * @param sourceLabel the starting node for all of the paths
     * @param targetLabel the ending node for all of the paths
     * @param K           the number of shortest paths to compute
     * @return a list of the K shortest paths from s to t, ordered from shortest to longest
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

        /* Make indexes to give (fast) access to these heaps later */
        // Heap H_out(v) for every node v
        HashMap<String,EppsteinHeap> nodeHeaps = new HashMap<String, EppsteinHeap>(graph.numNodes());
        HashMap<String,EppsteinHeap> edgeHeaps = new HashMap<String, EppsteinHeap>(graph.numEdges());
        // Heap H_T(v) for every node v
        HashMap<String,EppsteinHeap> outrootHeaps = new HashMap<String, EppsteinHeap>();
        HashMap<String, EppsteinArrayHeap> arrayHeaps = new HashMap<String, EppsteinArrayHeap>(graph.numNodes());

        // Create a virtual/dummy heap that is the root of the overall Eppstein heap. It represents the best path from
        // the source node to the target node, which does not involve any sidetrack edges.
        EppsteinHeap hg = new EppsteinHeap(new Edge(sourceLabel, sourceLabel, 0));

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
            ksp.add(kpath);

            // Push the (up to 3) children of this path within the Eppstein heap onto the priority queue
            addExplicitChildrenToQueue(kpathImplicit, ksp, pathPQ);

            /* Check for the existence of a potential fourth child, known as a "cross edge", to push onto the queue.
                This heap edge/child does not need to be explicitly represented in the Eppstein heap because it is easy
                to check for its existence. */
            // If the cross-edge child does not exist, try building its part of the heap, then check again.
            // Note: this is where all of the heap building happens in the lazy version of Eppstein's algorithm.
            if (!outrootHeaps.containsKey(kpathImplicit.getHeap().getSidetrack().getToNode())) {
                buildHeap(kpathImplicit.getHeap().getSidetrack().getToNode(), graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps, outrootHeaps, arrayHeaps, tree);
            }
            addCrossEdgeChildToQueue(outrootHeaps, kpathImplicit, i, ksp, pathPQ);
        }

        return ksp;
    }

    /**
     * Build H_out(v) and H_T(v) for node v, and recursively build H_T(nextT(v)) if it does not exist.
     *
     * @param nodeLabel             node v
     * @param graph                 the graph, G, on which to compute the K shortest paths from s to t
     * @param sidetrackEdgeCostMap  the cost of each sidetrack edge in G
     * @param nodeHeaps             an index/hash table of heap H_out(v), for each node v in the graph
     * @param edgeHeaps             an index/hash table of heaps H_out(v), but indexed by sidetrack edge
     * @param outrootHeaps          an index of heaps H_T(v) for each node v
     * @param arrayHeaps            an index of the array representation of heaps H_T(v) for each node v, as they are
     *                              being built
     * @param tree                  the shortest path tree, T, rooted at the target node, t
     */
    private static void buildHeap(String nodeLabel, Graph graph, HashMap<String, Double> sidetrackEdgeCostMap, HashMap<String, EppsteinHeap> nodeHeaps, HashMap<String, EppsteinHeap> edgeHeaps, HashMap<String, EppsteinHeap> outrootHeaps, HashMap<String, EppsteinArrayHeap> arrayHeaps, ShortestPathTree tree) {
        /* Part 1: Compute sub-heap H_out(v). */
        computeOutHeap(nodeLabel, graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps);

        /* PART 2: Compute sub-heap H_T(v) */

        /* Need to get H_T(nextT(v)) before H_T(v) can be built. If H_T(nextT(v)) has not been built, build it. */
        EppsteinArrayHeap currentArrayHeap;
        if (nodeLabel.equals(tree.getRoot())) {
            // Case 0 (recursive base case): v = the target node, t, so there is no H_T(nextT(v))
            currentArrayHeap = new EppsteinArrayHeap();
        }
        else {
            if (!outrootHeaps.containsKey(tree.getParentOf(nodeLabel))) {
                // Case 1: H_T(nextT(v)) has not been built, so build it
                buildHeap(tree.getParentOf(nodeLabel), graph, sidetrackEdgeCostMap, nodeHeaps, edgeHeaps, outrootHeaps, arrayHeaps, tree);
            }
            // Case 2: H_T(nextT(v)) has been built
            currentArrayHeap = arrayHeaps.get(tree.getParentOf(nodeLabel));
        }

        /* Create H_T(v) from H_T(nextT(v) */
        EppsteinHeap sidetrackHeap = nodeHeaps.get(nodeLabel);
        if (sidetrackHeap != null) {
            currentArrayHeap = currentArrayHeap.clone();
            currentArrayHeap.addOutroot(sidetrackHeap);
        }

        // Index the array representation of H_T(v), which is in a modifiable form and is used only by buildHeap() to
        // create further H_T heaps for which the current heap is a dependency.
        arrayHeaps.put(nodeLabel,currentArrayHeap);

        // Index the static, pointer representation of H_T(v), which is traversed by Eppstein's algorithm to add
        // candidate paths to the priority queue.
        EppsteinHeap currentHeap = currentArrayHeap.toEppsteinHeap2();
        if (currentHeap != null) {
            outrootHeaps.put(nodeLabel, currentHeap);
        }
    }
}