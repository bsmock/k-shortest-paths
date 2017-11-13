package edu.ufl.cise.bsmock.graph;

/**
 * The WeightedEdge class implements standard properties and methods for a weighted edge in a directed graph.
 *
 * Created by Brandon Smock on 6/6/15.
 */
public class WeightedEdge implements Comparable<WeightedEdge> {
    private String sourceLabel;
    private String targetLabel;
    private double edgeWeight = 0.0;

    public WeightedEdge(String targetLabel, double edgeWeight) {
        this.targetLabel = targetLabel;
        this.edgeWeight = edgeWeight;
    }

    public WeightedEdge(String sourceLabel, String targetLabel, double edgeWeight) {
        this.sourceLabel = sourceLabel;
        this.targetLabel = targetLabel;
        this.edgeWeight = edgeWeight;
    }

    public String getSourceLabel() {
        return sourceLabel;
    }

    public void setSourceLabel(String sourceLabel) {
        this.sourceLabel = sourceLabel;
    }

    public String getTargetLabel() {
        return targetLabel;
    }

    public void setTargetLabel(String targetLabel) {
        this.targetLabel = targetLabel;
    }

    public double getEdgeWeight() {
        return edgeWeight;
    }

    public void setEdgeWeight(double edgeWeight) {
        this.edgeWeight = edgeWeight;
    }

    public int compareTo(WeightedEdge comparedObject) {
        double weight1 = this.getEdgeWeight();
        double weight2 = comparedObject.getEdgeWeight();

        if (weight1 == weight2)
            return 0;
        if (weight1 > weight2)
            return 1;
        return -1;
    }
}
