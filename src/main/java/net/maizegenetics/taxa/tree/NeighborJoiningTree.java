// NeighborJoiningTree.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
// computational complexity O(numSeqs^3)
package net.maizegenetics.taxa.tree;

import net.maizegenetics.taxa.distance.DistanceMatrix;

/**
 * constructs a neighbor-joining tree from pairwise distances
 * <br><br>
 * Saitou, N., and Nei, M., (1987) The neighbor-joining method: A new method for
 * reconstructing phylogenetic trees. <i> Mol. Biol. Evol,</i> 4(4):406-425,
 * <br>
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class NeighborJoiningTree extends SimpleTree {
    
    /**
     * construct NJ tree
     *
     * @param m distance matrix
     */
    public NeighborJoiningTree(DistanceMatrix m) {
        if (m.getSize() < 3) {
            throw new IllegalArgumentException("NeighborJoiningTree: Less than 3 taxa in distance matrix.");
        }
        if (!m.isSymmetric()) {
            throw new IllegalArgumentException("NeighborJoiningTree: Unsymmetrix Distance Matrix: Probably due to taxa with large proportion of missing sites.");
        }

        init(m);

        //while (numClusters > 3)
        while (true) {
            findNextPair();
            newBranchLengths();
            if (numClusters == 3) {
                break;
            }
            newCluster();
        }

        finish();
    }

    private int numClusters;
    private int besti, abi;
    private int bestj;
    private int[] alias;
    private double[][] distance;
    private double[] r;
    private double scale;

    private double getDist(int a, int b) {
        return distance[alias[a]][alias[b]];
    }

    private void init(DistanceMatrix m) {
        numClusters = m.getSize();

        distance = m.getClonedDistances();

        for (int i = 0; i < numClusters; i++) {
            Node tmp = NodeFactory.createNode();
            tmp.setIdentifier(m.getTaxon(i));
            getRoot().addChild(tmp);
        }

        alias = new int[numClusters];
        for (int i = 0; i < numClusters; i++) {
            alias[i] = i;
        }

        r = new double[numClusters];
    }

    private void finish() {
        if (besti != 0 && bestj != 0) {
            getRoot().getChild(0).setBranchLength(updatedDistance(besti, bestj, 0));
        } else if (besti != 1 && bestj != 1) {
            getRoot().getChild(1).setBranchLength(updatedDistance(besti, bestj, 1));
        } else {
            getRoot().getChild(2).setBranchLength(updatedDistance(besti, bestj, 2));
        }
        distance = null;

        // make node heights available also
        NodeUtils.lengths2Heights(getRoot());
    }

    private void findNextPair() {
        for (int i = 0; i < numClusters; i++) {
            r[i] = 0;
            for (int j = 0; j < numClusters; j++) {
                r[i] += getDist(i, j);
            }
        }

        besti = 0;
        bestj = 1;
        double smax = -1.0;
        scale = 1.0 / (numClusters - 2);
        for (int i = 0; i < numClusters - 1; i++) {
            for (int j = i + 1; j < numClusters; j++) {
                double sij = (r[i] + r[j]) * scale - getDist(i, j);

                if (sij > smax) {
                    smax = sij;
                    besti = i;
                    bestj = j;
                }
            }
        }
        abi = alias[besti];
    }

    private void newBranchLengths() {
        double dij = getDist(besti, bestj);
        double li = (dij + (r[besti] - r[bestj]) * scale) * 0.5;
        double lj = dij - li; // = (dij + (r[bestj]-r[besti])*scale)*0.5

        getRoot().getChild(besti).setBranchLength(li);
        getRoot().getChild(bestj).setBranchLength(lj);
    }

    private void newCluster() {
        // Update distances
        for (int k = 0; k < numClusters; k++) {
            if (k != besti && k != bestj) {
                int ak = alias[k];
                distance[ak][abi] = distance[abi][ak] = updatedDistance(besti, bestj, k);
            }
        }
        distance[abi][abi] = 0.0;

        // Replace besti with new cluster
        NodeUtils.joinChilds(getRoot(), besti, bestj);

        // Update alias
        for (int i = bestj; i < numClusters - 1; i++) {
            alias[i] = alias[i + 1];
        }

        numClusters--;
    }

    /**
     * compute updated distance between the new cluster (i,j) to any other
     * cluster k
     */
    private double updatedDistance(int i, int j, int k) {
        return (getDist(k, i) + getDist(k, j) - getDist(i, j)) * 0.5;
    }

}
