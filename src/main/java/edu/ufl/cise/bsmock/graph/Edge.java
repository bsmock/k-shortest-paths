package edu.ufl.cise.bsmock.graph;

/**
 * The Edge class implements standard properties and methods for a weighted edge in a directed graph.
 *
 * Created by Brandon Smock on 6/19/15.
 */
public class Edge implements Cloneable {
    private String fromNode;
    private String toNode;
    private double weight;

    public Edge() {
        this.fromNode = null;
        this.toNode = null;
        this.weight = Double.MAX_VALUE;
    }

    public Edge(String fromNode, String toNode, double weight) {
        this.fromNode = fromNode;
        this.toNode = toNode;
        this.weight = weight;
    }

    public String getFromNode() {
        return fromNode;
    }

    public void setFromNode(String fromNode) {
        this.fromNode = fromNode;
    }

    public String getToNode() {
        return toNode;
    }

    public void setToNode(String toNode) {
        this.toNode = toNode;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public Edge clone() {
        return new Edge(fromNode, toNode, weight);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        sb.append(fromNode);
        sb.append(",");
        sb.append(toNode);
        sb.append("){");
        sb.append(weight);
        sb.append("}");

        return sb.toString();
    }

    public boolean equals(Edge edge2) {
        if (hasSameEndpoints(edge2) && weight == edge2.getWeight())
            return true;

        return false;
    }

    public boolean hasSameEndpoints(Edge edge2) {
        if (fromNode.equals(edge2.getFromNode()) && toNode.equals(edge2.getToNode()))
            return true;

        return false;
    }
}
