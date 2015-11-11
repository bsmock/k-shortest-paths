package edu.ufl.cise.bsmock.graph.util;

/**
 * Created by brandonsmock on 6/1/15.
 */
import java.util.*;
import edu.ufl.cise.bsmock.graph.*;

public final class Dijkstra {

    private Dijkstra() {}

    public static ShortestPathTree shortestPathTree(Graph graph, String sourceLabel) throws Exception {
        HashMap<String,Node> nodes = graph.getNodes();
        if (!nodes.containsKey(sourceLabel))
            throw new Exception("Source node not found in graph.");
        ShortestPathTree predecessorTree = new ShortestPathTree(sourceLabel);
        Set<DijkstraNode> visited = new HashSet<DijkstraNode>();
        PriorityQueue<DijkstraNode> pq = new PriorityQueue<DijkstraNode>();
        for (String nodeLabel:nodes.keySet()) {
            DijkstraNode newNode = new DijkstraNode(nodeLabel);
            newNode.setDist(Double.MAX_VALUE);
            newNode.setDepth(Integer.MAX_VALUE);
            predecessorTree.add(newNode);
        }
        DijkstraNode sourceNode = predecessorTree.getNodes().get(predecessorTree.getRoot());
        sourceNode.setDist(0);
        sourceNode.setDepth(0);
        pq.add(sourceNode);

        int count = 0;
        while (!pq.isEmpty()) {
            DijkstraNode current = pq.poll();
            String currLabel = current.getLabel();
            visited.add(current);
            count++;
            HashMap<String, Double> neighbors = nodes.get(currLabel).getNeighbors();
            for (String currNeighborLabel:neighbors.keySet()) {
                DijkstraNode neighborNode = predecessorTree.getNodes().get(currNeighborLabel);
                Double currDistance = neighborNode.getDist();
                Double newDistance = current.getDist() + nodes.get(currLabel).getNeighbors().get(currNeighborLabel);
                if (newDistance < currDistance) {
                    DijkstraNode neighbor = predecessorTree.getNodes().get(currNeighborLabel);

                    pq.remove(neighbor);
                    neighbor.setDist(newDistance);
                    neighbor.setDepth(current.getDepth() + 1);
                    neighbor.setParent(currLabel);
                    pq.add(neighbor);
                }
            }
        }

        return predecessorTree;
    }

    public static Path shortestPath(Graph graph, String sourceLabel, String targetLabel) throws Exception {
        //if (!nodes.containsKey(sourceLabel))
        //    throw new Exception("Source node not found in graph.");
        HashMap<String,Node> nodes = graph.getNodes();
        ShortestPathTree predecessorTree = new ShortestPathTree(sourceLabel);
        PriorityQueue<DijkstraNode> pq = new PriorityQueue<DijkstraNode>();
        for (String nodeLabel:nodes.keySet()) {
            DijkstraNode newNode = new DijkstraNode(nodeLabel);
            newNode.setDist(Double.MAX_VALUE);
            newNode.setDepth(Integer.MAX_VALUE);
            predecessorTree.add(newNode);
        }
        DijkstraNode sourceNode = predecessorTree.getNodes().get(predecessorTree.getRoot());
        sourceNode.setDist(0);
        sourceNode.setDepth(0);
        pq.add(sourceNode);

        int count = 0;
        while (!pq.isEmpty()) {
            DijkstraNode current = pq.poll();
            String currLabel = current.getLabel();
            if (currLabel.equals(targetLabel)) {
                Path shortestPath = new Path();
                String currentN = targetLabel;
                String parentN = predecessorTree.getParentOf(currentN);
                while (parentN != null) {
                    shortestPath.addFirst(new Edge(parentN,currentN,nodes.get(parentN).getNeighbors().get(currentN)));
                    currentN = parentN;
                    parentN = predecessorTree.getParentOf(currentN);
                }
                return shortestPath;
            }
            count++;
            HashMap<String, Double> neighbors = nodes.get(currLabel).getNeighbors();
            for (String currNeighborLabel:neighbors.keySet()) {
                DijkstraNode neighborNode = predecessorTree.getNodes().get(currNeighborLabel);
                Double currDistance = neighborNode.getDist();
                Double newDistance = current.getDist() + nodes.get(currLabel).getNeighbors().get(currNeighborLabel);
                if (newDistance < currDistance) {
                    DijkstraNode neighbor = predecessorTree.getNodes().get(currNeighborLabel);

                    pq.remove(neighbor);
                    neighbor.setDist(newDistance);
                    neighbor.setDepth(current.getDepth() + 1);
                    neighbor.setParent(currLabel);
                    pq.add(neighbor);
                }
            }
        }

        return null;
    }
}
