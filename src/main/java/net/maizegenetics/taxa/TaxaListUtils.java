/*
 *  TaxaListUtils
 */
package net.maizegenetics.taxa;

import net.maizegenetics.util.Tuple;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Terry Casstevens
 */
public class TaxaListUtils {

    private TaxaListUtils() {
        // utility class
    }

    /**
     * Intersect joins the specified taxa.
     *
     * @param group1 an TaxaList
     * @param group2 another TaxaList
     *
     * @return the taxa in the intersection of groups 1 and 2, sorted in ascending order
     */
    public static TaxaList getCommonTaxa(TaxaList group1, TaxaList group2) {
        return getCommonTaxa(new TaxaList[]{group1, group2});
    }

    /**
     * Intersect joins the specified taxa.
     *
     * @param groups groups to join.
     *
     * @return The taxa from the intersect join sorted alphabetically
     */
    public static TaxaList getCommonTaxa(TaxaList[] groups) {
        return getCommonTaxa(groups, true);
    }

    /**
     * Intersect joins the specified taxa.
     *
     * @param groups groups to join.
     * @param sorted whether to sort taxa alphabetically
     *
     * @return The taxa from the intersect join
     */
    public static TaxaList getCommonTaxa(TaxaList[] groups, boolean sorted) {

        if ((groups == null) || (groups.length == 0)) {
            return null;
        } else if (groups.length == 1 && !sorted) {
            return groups[0];
        }

        Set<Taxon> intersectIds = new LinkedHashSet<>();
        for (Taxon current : groups[0]) {
            intersectIds.add(current);
        }
        for (int i = 1; i < groups.length; i++) {
            List<Taxon> temp = new ArrayList<>();
            for (Taxon current : groups[i]) {
                temp.add(current);
            }
            intersectIds.retainAll(temp);
        }

        TaxaListBuilder builder = new TaxaListBuilder();
        builder.addAll(intersectIds);
        if (sorted) {
            builder.sortTaxaAlphabetically();
        }
        return builder.build();

    }

    /**
     * Union joins the specified taxa.
     *
     * @param group1 an TaxaList
     * @param group2 another TaxaList
     *
     * @return the taxa in the union of taxa 1 and 2, sorted in ascending order
     */
    public static TaxaList getAllTaxa(TaxaList group1, TaxaList group2) {
        return getAllTaxa(new TaxaList[]{group1, group2});
    }

    public static TaxaList getAllTaxa(TaxaList[] groups) {
        return getAllTaxa(groups, true);
    }

    /**
     * Union joins the specified taxa.
     *
     * @param lists taxa to join.
     * @param sorted whether to sort taxa alphabetically
     *
     * @return The taxa from the union join
     */
    public static TaxaList getAllTaxa(TaxaList[] lists, boolean sorted) {

        if ((lists == null) || (lists.length == 0)) {
            return null;
        } else if (lists.length == 1) {
            return lists[0];
        }

        List<Taxon> allIds = new ArrayList<>();

        for (Taxon current : lists[0]) {
            allIds.add(current);
        }

        for (int i = 1; i < lists.length; i++) {
            int insertIndex = 0;
            for (Taxon current : lists[i]) {
                int index = allIds.indexOf(current);
                if (index == -1) {
                    allIds.add(insertIndex++, current);
                } else {
                    insertIndex = Math.max(index + 1, insertIndex);
                }
            }
        }

        TaxaListBuilder builder = new TaxaListBuilder();
        builder.addAll(allIds);
        if (sorted) {
            builder.sortTaxaAlphabetically();
        }
        return builder.build();

    }

    /**
     * Return whether all taxa same between lists
     *
     * @param lists taxa lists
     *
     * @return whether all taxa same between lists
     */
    public static boolean areAllTaxaSame(TaxaList[] lists) {

        if ((lists == null) || (lists.length <= 1)) {
            return true;
        }

        TaxaList first = lists[0];
        for (int i = 1; i < lists.length; i++) {
            if (lists[i].size() != first.size()) {
                return false;
            }
            if (!first.containsAll(lists[i])) {
                return false;
            }
        }

        return true;

    }

    public static Tuple<List<Taxon>, Set<Taxon>> getUnionAndIntersectTaxa(TaxaList[] lists) {

        if ((lists == null) || (lists.length <= 1)) {
            return null;
        }

        List<Taxon> allIds = new ArrayList<>();
        Set<Taxon> intersectIds = new LinkedHashSet<>();

        for (Taxon current : lists[0]) {
            allIds.add(current);
            intersectIds.add(current);
        }

        for (int i = 1; i < lists.length; i++) {

            List<Taxon> temp = new ArrayList<>();
            int insertIndex = 0;
            for (Taxon current : lists[i]) {
                temp.add(current);
                int index = allIds.indexOf(current);
                if (index == -1) {
                    allIds.add(insertIndex++, current);
                } else {
                    insertIndex = Math.max(index + 1, insertIndex);
                }
            }
            intersectIds.retainAll(temp);

        }

        return new Tuple<>(allIds, intersectIds);

    }

}
