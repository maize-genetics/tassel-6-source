// Tree.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.taxa.tree;

import java.io.Serializable;
import net.maizegenetics.taxa.Taxon;

/**
 * Interface for a phylogenetic or genealogical tree.
 *
 * @author Alexei Drummond
 */
public interface Tree extends Serializable, UnitsProvider {

    /**
     * @return the root node of this tree.
     */
    Node getRoot();

    /**
     * This method constructs a tree from the given root node.
     *
     * @param root the root node of the tree to construct.
     */
    void setRoot(Node root);

    /**
     * @return a count of the number of external nodes (tips) in this tree.
     */
    int getExternalNodeCount();

    /**
     * @return a count of the number of internal nodes (and hence clades) in
     * this tree.
     */
    int getInternalNodeCount();

    /**
     * @return the ith external node in the tree.
     */
    Node getExternalNode(int i);

    /**
     * @return the ith internal node in the tree.
     */
    Node getInternalNode(int i);

    /**
     * This method is called to ensure that the calls to other methods in this
     * interface are valid.
     */
    void createNodeList();

    /**
     * Gets the units that this tree's branch lengths and node heights are
     * expressed in.
     */
    int getUnits();
    
    public int whichIdNumber(Taxon t);

    /**
     * Sets an named attribute for a given node.
     *
     * @param node the node whose attribute is being set.
     * @param name the name of the attribute.
     * @param value the new value of the attribute.
     */
    void setAttribute(Node node, String name, Object value);

    /**
     * @return an object representing the named attributed for the numbered
     * node.
     * @param node the node being interrogated.
     * @param name the name of the attribute of interest.
     */
    Object getAttribute(Node node, String name);

    /**
     * @return a clone of this tree
     */
    public Tree getCopy();

}
