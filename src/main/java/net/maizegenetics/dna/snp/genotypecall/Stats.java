/*
 *  Stats
 *
 *  Created on Jan 6, 2017
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

import static net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache.HETEROZYGOUS_COUNT;

/**
 * @author Terry Casstevens
 */
public class Stats {

    private final int[][] myAlleleCounts;
    private final int[] myStats;
    private final int myNumIndices;
    private final int myIndex;

    private Stats(int[][] alleleCounts, int[] stats, int numIndices, int index) {
        myAlleleCounts = alleleCounts;
        myStats = stats;
        myIndex = index;
        myNumIndices = numIndices;
    }

    public static Stats getInstance(int[][] alleleCounts, int[] stats, int numIndices, int index) {
        return new Stats(alleleCounts, stats, numIndices, index);
    }

    public static Stats getInstance(Stats orig, int newIndex) {
        return new Stats(orig.myAlleleCounts, orig.myStats, orig.myNumIndices, newIndex);
    }

    public byte majorAllele() {
        if (myAlleleCounts[0].length >= 1) {
            return (byte) myAlleleCounts[0][0];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public double majorAlleleFrequency() {

        int numAlleles = myAlleleCounts[0].length;
        if (numAlleles >= 1) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing += myAlleleCounts[1][i];
            }
            return (double) myAlleleCounts[1][0] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    public byte minorAllele() {
        if (myAlleleCounts[0].length >= 2) {
            return (byte) myAlleleCounts[0][1];
        } else {
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    public double minorAlleleFrequency() {

        int numAlleles = myAlleleCounts[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing += myAlleleCounts[1][i];
            }
            return (double) myAlleleCounts[1][1] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    public int numAllMinorAlleles() {

        int result = 0;
        for (int i = 1; i < myAlleleCounts[0].length; i++) {
            result += myAlleleCounts[1][i];
        }
        return result;

    }

    public int totalGametesNonMissingForSite() {
        int numAlleles = myAlleleCounts[0].length;
        int result = 0;
        for (int i = 0; i < numAlleles; i++) {
            result += myAlleleCounts[1][i];
        }
        return result;
    }

    public double proportionHeterozygous() {
        return (double) myStats[HETEROZYGOUS_COUNT] / (double) myNumIndices;
    }

    public double percentNotMissing() {
        return (double) (myNumIndices - myStats[AlleleFreqCache.UNKNOWN_COUNT]) / (double) myNumIndices;
    }

    public boolean hasIndel() {
        int numAlleles = myAlleleCounts[0].length;
        for (int i = 0; i < numAlleles; i++) {
            if (myAlleleCounts[0][i] == NucleotideAlignmentConstants.GAP_ALLELE || myAlleleCounts[0][i] == NucleotideAlignmentConstants.INSERT_ALLELE) {
                return true;
            }
        }
        return false;
    }

    public int index() {
        return myIndex;
    }

}
