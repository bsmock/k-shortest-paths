package edu.ufl.cise.bsmock.graph;

/**
 * The Graph class implements a weighted, directed graph using an adjacency list representation.
 *
 * Created by brandonsmock on 6/1/15.
 */

import java.io.*;
import java.util.*;

public class Graph {
    private HashMap<String,Node> nodes;

    public Graph() {
        nodes = new HashMap<String,Node>();
    }

    public Graph(String filename) {
        this();
        readFromFile(filename);
    }

    public Graph(HashMap<String,Node> nodes) {
        this.nodes = nodes;
    }

    public int numNodes() {
        return nodes.size();
    }

    public int numEdges() {
        int edgeCount = 0;
        for (Node node : nodes.values()) {
            edgeCount += node.getEdges().size();
        }
        return edgeCount;
    }

    public void addNode(String label) {
        if (!nodes.containsKey(label))
            nodes.put(label,new Node(label));
    }

    public void addNode(Node node) {
        String label = node.getLabel();
        if (!nodes.containsKey(label))
            nodes.put(label,node);
    }

    public void addEdge(String label1, String label2, Double weight) {
        if (!nodes.containsKey(label1))
            addNode(label1);
        if (!nodes.containsKey(label2))
            addNode(label2);
        nodes.get(label1).addEdge(label2,weight);
    }

    public void addEdge(Edge edge) {
        addEdge(edge.getFromNode(),edge.getToNode(),edge.getWeight());
    }

    public void addEdges(List<Edge> edges) {
        for (Edge edge : edges) {
            addEdge(edge);
        }
    }

    public Edge removeEdge(String label1, String label2) {
        if (nodes.containsKey(label1)) {
            double weight = nodes.get(label1).removeEdge(label2);
            if (weight != Double.MAX_VALUE) {
                return new Edge(label1, label2, weight);
            }
        }

        return null;
    }

    public double getEdgeWeight(String label1, String label2) {
        if (nodes.containsKey(label1)) {
            Node node1 = nodes.get(label1);
            if (node1.getNeighbors().containsKey(label2)) {
                return node1.getNeighbors().get(label2);
            }
        }

        return Double.MAX_VALUE;
    }

    public HashMap<String,Node> getNodes() {
        return nodes;
    }

    public List<Edge> getEdgeList() {
        List<Edge> edgeList = new LinkedList<Edge>();

        for (Node node : nodes.values()) {
            edgeList.addAll(node.getEdges());
        }

        return edgeList;
    }

    public Set<String> getNodeLabels() {
        return nodes.keySet();
    }

    public Node getNode(String label) {
        return nodes.get(label);
    }

    public List<Edge> removeNode(String label) {
        LinkedList<Edge> edges = new LinkedList<Edge>();
        if (nodes.containsKey(label)) {
            Node node = nodes.remove(label);
            edges.addAll(node.getEdges());
            edges.addAll(removeEdgesToNode(label));
        }

        return edges;
    }

    public List<Edge> removeEdgesToNode(String label) {
        List<Edge> edges = new LinkedList<Edge>();
        for (Node node : nodes.values()) {
            if (node.getAdjacencyList().contains(label)) {
                double weight = node.removeEdge(label);
                edges.add(new Edge(node.getLabel(),label,weight));
            }
        }
        return edges;
    }



    public Graph transpose() {
        HashMap<String,Node> newNodes = new HashMap<String, Node>();

        Iterator<String> it = nodes.keySet().iterator();
        while (it.hasNext()) {
            String nodeLabel = it.next();
            newNodes.put(nodeLabel,new Node(nodeLabel));
        }

        it = nodes.keySet().iterator();
        while (it.hasNext()) {
            String nodeLabel = it.next();
            Node node = nodes.get(nodeLabel);
            Set<String> adjacencyList = node.getAdjacencyList();
            Iterator<String> alIt = adjacencyList.iterator();
            HashMap<String, Double> neighbors = node.getNeighbors();
            while (alIt.hasNext()) {
                String neighborLabel = alIt.next();
                newNodes.get(neighborLabel).addEdge(nodeLabel,neighbors.get(neighborLabel));
            }
        }

        return new Graph(newNodes);
    }

    public void clear() {
        nodes = new HashMap<String,Node>();
    }

    public void readFromFile(String fileName) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(fileName));

            String line = in.readLine();

            while (line != null) {
                String[] edgeDescription = line.split("\\s");
                if (edgeDescription.length == 3) {
                    addEdge(edgeDescription[0],edgeDescription[1],Double.parseDouble(edgeDescription[2]));
                    //addEdge(edgeDescription[1],edgeDescription[0],Double.parseDouble(edgeDescription[2]));
                }
                line = in.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String toString() {
        StringBuilder graphStringB = new StringBuilder();
        Iterator<String> it = nodes.keySet().iterator();
        while (it.hasNext()) {
            String nodeLabel = it.next();
            graphStringB.append(nodeLabel.toString());
            graphStringB.append(": {");
            Node node = nodes.get(nodeLabel);
            Set<String> adjacencyList = node.getAdjacencyList();
            Iterator<String> alIt = adjacencyList.iterator();
            HashMap<String, Double> neighbors = node.getNeighbors();
            while (alIt.hasNext()) {
                String neighborLabel = alIt.next();
                graphStringB.append(neighborLabel.toString());
                graphStringB.append(": ");
                graphStringB.append(neighbors.get(neighborLabel));
                if (alIt.hasNext())
                    graphStringB.append(", ");
            }
            graphStringB.append("}");
            graphStringB.append("\n");
        }

        return graphStringB.toString();
    }

    public void graphToFile(String filename) {
        BufferedWriter writer = null;
        try {
            File subgraphFile = new File(filename);

            // This will output the full path where the file will be written to...
            System.out.println(subgraphFile.getCanonicalPath());

            writer = new BufferedWriter(new FileWriter(subgraphFile));
            writer.write(Integer.toString(nodes.size()) + "\n\n");

            Iterator<Node> it = nodes.values().iterator();
            while (it.hasNext()) {
                Node node = it.next();
                String nodeLabel = node.getLabel();
                if (nodes.containsKey(nodeLabel)) {
                    HashMap<String,Double> neighbors = node.getNeighbors();
                    Iterator<String> it2 = neighbors.keySet().iterator();
                    while (it2.hasNext()) {
                        String nodeLabel2 = it2.next();
                        if (nodes.containsKey(nodeLabel2)) {
                            writer.write(nodeLabel + " " + nodeLabel2 + " " + neighbors.get(nodeLabel2) + "\n");
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                // Close the writer regardless of what happens...
                writer.close();
            } catch (Exception e) {
            }
        }
    }
}
