package edu.ufl.cise.bsmock.graph.ksp;

import edu.ufl.cise.bsmock.graph.Graph;
import edu.ufl.cise.bsmock.graph.util.Path;
import java.util.List;

/**
 * Created by brandonsmock on 12/24/15.
 */
public interface KSPAlgorithm {
    public boolean isLoopless();

    public List<Path> ksp(Graph graph, String sourceLabel, String targetLabel, int K);
}
