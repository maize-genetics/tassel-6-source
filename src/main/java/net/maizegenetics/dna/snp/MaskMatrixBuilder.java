/*
 *  MaskMatrixBuilder
 *
 *  Created on Dec 10, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.Stats;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

import java.util.Random;
import java.util.function.BiPredicate;
import java.util.function.Predicate;

/**
 * @author Terry Casstevens
 */
public class MaskMatrixBuilder {

    private BitSet[] myBitSets;
    private final boolean myIsSiteOptimized;
    private final int myNumTaxa;
    private final int myNumSites;
    private int myNextSite = 0;
    private int myNextTaxon = 0;
    private long myNextCount = 0;

    private MaskMatrixBuilder(int numTaxa, int numSites, boolean isSiteOptimized) {
        myNumTaxa = numTaxa;
        myNumSites = numSites;
        myIsSiteOptimized = isSiteOptimized;
        if (myIsSiteOptimized) {
            myBitSets = new BitSet[myNumSites];
            for (int s = 0; s < myNumSites; s++) {
                myBitSets[s] = new OpenBitSet(myNumTaxa);
            }
        } else {
            myBitSets = new BitSet[myNumTaxa];
            for (int t = 0; t < myNumTaxa; t++) {
                myBitSets[t] = new OpenBitSet(myNumSites);
            }
        }
    }

    public static MaskMatrixBuilder getInstance(int numTaxa, int numSites, boolean isSiteOptimized) {
        return new MaskMatrixBuilder(numTaxa, numSites, isSiteOptimized);
    }

    public static MaskMatrixBuilder getInstance(MaskMatrix orig) {

        if (orig == null) {
            throw new IllegalArgumentException("MaskMatrixBuilder: getInstance: must specific orig");
        }

        if (orig instanceof MaskSiteMatrix) {
            MaskSiteMatrix temp = (MaskSiteMatrix) orig;
            MaskMatrixBuilder result = new MaskMatrixBuilder(orig.numTaxa(), orig.numSites(), true);
            for (int s = 0; s < temp.numSites(); s++) {
                result.myBitSets[s].or(temp.maskForSite(s));
            }
            return result;
        } else if (orig instanceof MaskTaxaMatrix) {
            MaskTaxaMatrix temp = (MaskTaxaMatrix) orig;
            MaskMatrixBuilder result = new MaskMatrixBuilder(orig.numTaxa(), orig.numSites(), false);
            for (int t = 0; t < temp.numTaxa(); t++) {
                result.myBitSets[t].or(temp.maskForTaxon(t));
            }
            return result;
        } else {
            throw new IllegalArgumentException("MaskMatrixBuilder: getInstance: don't know type: " + orig.getClass().getName());
        }

    }

    public static MaskMatrix getInstanceRemoveMinorSNPs(GenotypeCallTable genotype) {
        return getInstance(genotype, (Byte t, Stats u) -> {
            byte major = u.majorAllele();
            byte minor = u.minorAllele();
            if (((t & 0xF) != major) && ((t & 0xF) != minor)) {
                return true;
            } else {
                return ((t >>> 4) != major) && ((t >>> 4) != minor);
            }
        });
    }

    public static MaskMatrix getInstanceRemoveMinCountSNPs(GenotypeCallTable genotype, int minCountOfMinors) {
        return getInstance(genotype, new BiPredicate<Byte, Stats>() {
            @Override
            public boolean test(Byte aByte, Stats stats) {
                if (aByte != GenotypeTable.UNKNOWN_DIPLOID_ALLELE && stats.numAllMinorAlleles() <= minCountOfMinors) {
                    byte major = stats.majorAllele();
                    if ((aByte & 0xF) != major || (aByte >>> 4) != major) {
                        return true;
                    }
                }
                return false;
            }
        });
    }

    public static MaskMatrix getInstanceRemoveHeterozygous(GenotypeCallTable genotype) {
        return getInstance(genotype, (Byte t) -> (t & 0xF) != (t >>> 4));
    }

    public static MaskMatrix getInstanceRemoveHomozygous(GenotypeCallTable genotype) {
        return getInstance(genotype, (Byte t) -> (t & 0xF) == (t >>> 4));
    }

    public static MaskMatrix getInstanceRemoveIndels(GenotypeCallTable genotype) {
        Predicate<Byte> predicate = diploid -> (((diploid >>> 4) & 0xf) == NucleotideAlignmentConstants.GAP_ALLELE ||
                ((diploid >>> 4) & 0xf) == NucleotideAlignmentConstants.INSERT_ALLELE ||
                (diploid & 0xf) == NucleotideAlignmentConstants.GAP_ALLELE ||
                (diploid & 0xf) == NucleotideAlignmentConstants.INSERT_ALLELE);
        return getInstance(genotype, predicate);
    }

    public static MaskMatrix getInstance(GenotypeCallTable genotype, Predicate<Byte> predicate) {
        return new MaskGenotypeMatrix(genotype, predicate);
    }

    public static MaskMatrix getInstance(GenotypeCallTable genotype, BiPredicate<Byte, Stats> predicate) {
        return new MaskGenotypeStatsMatrix(genotype, predicate);
    }

    public boolean get(int taxon, int site) {
        if (myIsSiteOptimized) {
            return myBitSets[site].fastGet(taxon);
        } else {
            return myBitSets[taxon].fastGet(site);
        }
    }

    public void set(int taxon, int site) {
        if (myIsSiteOptimized) {
            myBitSets[site].fastSet(taxon);
        } else {
            myBitSets[taxon].fastSet(site);
        }
    }

    public void setNext(boolean value) {
        if (myIsSiteOptimized) {
            if (value) {
                myBitSets[myNextTaxon].fastSet(myNextSite);
            }
            myNextCount++;
            myNextTaxon = (int) (myNextCount % myNumTaxa);
            myNextSite = (int) (myNextCount / myNumTaxa);
        } else {
            if (value) {
                myBitSets[myNextTaxon].fastSet(myNextSite);
            }
            myNextCount++;
            myNextTaxon = (int) (myNextCount / myNumSites);
            myNextSite = (int) (myNextCount % myNumSites);
        }
    }

    public long reduceMaskTo(double percent) {
        if (myIsSiteOptimized) {
            return siteReduceMaskTo(percent);
        } else {
            return taxaReduceMaskTo(percent);
        }
    }

    private long taxaReduceMaskTo(double percent) {

        Random random = new Random();
        double remainder = 0.0;
        long totalNumMasked = 0;

        for (BitSet current : myBitSets) {
            long numMasksThisTaxon = current.cardinality();
            double percentOfMasksToKeep = numMasksThisTaxon * percent;
            int numOfMasksToKeep = (int) Math.floor(percentOfMasksToKeep);
            remainder += (percentOfMasksToKeep - (double) numOfMasksToKeep);
            if (remainder > 1.0) {
                numOfMasksToKeep++;
                remainder -= 1.0;
            }

            if (numOfMasksToKeep != 0) {

                BitSet copy = new OpenBitSet(current);
                int numCleared = 0;
                while (true) {
                    int site = random.nextInt(myNumSites);
                    if (copy.getAndClear(site)) {
                        numCleared++;
                        if (numCleared >= numOfMasksToKeep) {
                            break;
                        }
                    }
                }
                current.xor(copy);
                totalNumMasked += numOfMasksToKeep;

            }
        }

        return totalNumMasked;

    }

    private long siteReduceMaskTo(double percent) {

        Random random = new Random();
        double remainder = 0.0;
        long totalNumMasked = 0;

        for (BitSet current : myBitSets) {
            long numMasksThisSite = current.cardinality();
            double percentOfMasksToKeep = numMasksThisSite * percent;
            int numOfMasksToKeep = (int) Math.floor(percentOfMasksToKeep);
            remainder += (percentOfMasksToKeep - (double) numOfMasksToKeep);
            if (remainder > 1.0) {
                numOfMasksToKeep++;
                remainder -= 1.0;
            }

            if (numOfMasksToKeep != 0) {

                BitSet copy = new OpenBitSet(current);
                int numCleared = 0;
                while (true) {
                    int taxon = random.nextInt(myNumTaxa);
                    if (copy.getAndClear(taxon)) {
                        numCleared++;
                        if (numCleared >= numOfMasksToKeep) {
                            break;
                        }
                    }
                }
                current.xor(copy);
                totalNumMasked += numOfMasksToKeep;

            }
        }

        return totalNumMasked;

    }

    public MaskMatrix build() {

        BitSet[] temp = myBitSets;
        myBitSets = null;

        if (myIsSiteOptimized) {
            return new MaskSiteMatrix(temp, myNumTaxa, myNumSites);
        } else {
            return new MaskTaxaMatrix(temp, myNumTaxa, myNumSites);
        }

    }

}
