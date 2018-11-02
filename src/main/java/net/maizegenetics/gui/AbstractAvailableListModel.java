/*
 * AbstractAvailableListModel
 */
package net.maizegenetics.gui;

import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import javax.swing.*;
import java.util.regex.Pattern;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractAvailableListModel extends AbstractListModel<String> {

    private final BitSet myShownIndices;
    private boolean myShowAll = true;
    private final int myRealSize;

    public AbstractAvailableListModel(int realSize) {
        myRealSize = realSize;
        myShownIndices = new OpenBitSet(myRealSize);
        showAll();
    }

    @Override
    public int getSize() {
        return (int) myShownIndices.cardinality();
    }

    @Override
    public String getElementAt(int index) {
        if (myShowAll) {
            return getRealElementAt(index);
        } else {
            return getRealElementAt(myShownIndices.indexOfNthSetBit(index + 1));
        }
    }

    public void setShown(String search) {

        if ((search == null) || (search.length() == 0)) {
            if (showAll()) {
                fireContentsChanged(this, 0, getSize());
            }
            return;
        }

        try {

            if (search.indexOf('*', search.length() - 1) == -1) {
                search = search + "*";
            }

            search = search.replaceAll("\\.", "\\\\.");
            search = search.replaceAll("\\*", "(.*)");

            showNone();
            Pattern pattern = Pattern.compile(search);
            for (int i = 0, n = getRealSize(); i < n; i++) {
                if (pattern.matcher(getRealElementAt(i)).matches()) {
                    myShownIndices.fastSet(i);
                }
            }

            if (myShownIndices.cardinality() == getRealSize()) {
                myShowAll = true;
            } else {
                myShowAll = false;
            }

        } catch (Exception e) {
            showNone();
        }

        fireContentsChanged(this, 0, getSize());

    }

    public int[] translateToRealIndices(int[] shownIndices) {
        for (int i = 0; i < shownIndices.length; i++) {
            shownIndices[i] = myShownIndices.indexOfNthSetBit(shownIndices[i] + 1);
        }
        return shownIndices;
    }

    /**
     * Sets bits to show all elements in list.
     *
     * @return whether bits changed
     */
    private boolean showAll() {
        myShowAll = true;
        if (myShownIndices.cardinality() == getRealSize()) {
            return false;
        } else {
            myShownIndices.set(0, getRealSize());
            return true;
        }
    }

    private void showNone() {
        myShowAll = false;
        myShownIndices.clear(0, getRealSize());
    }

    public int getRealSize() {
        return myRealSize;
    }

    abstract public String getRealElementAt(int index);
}
