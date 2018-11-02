package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;
import org.apache.commons.math3.distribution.BinomialDistribution;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;

/**
 * Defines the methods for merging the calls from two taxa.  The merge rules need to be defined at the level of
 * genotypic calls and for read depth.  In general if depth is available, it will be used to merge.
 *
 * @author Ed Buckler
 * @author Jeff Glaubitz
 */
public class BasicGenotypeMergeRule implements GenotypeMergeRule {
    private final double errorRate;
    private final int maxCountAtGeno=100; // can return to 500 once we fix the BinomialDistribution below, which exceeded 30%
    private final int[] likelihoodRatioThreshAlleleCnt;  // index = sample size; value = min count of less tagged allele for likelihood ratio > 1
        // if less tagged allele has counts < likelihoodRatioThreshAlleleCnt[totalCount], call it a homozygote
        // where likelihood ratio = (binomial likelihood het) / (binomial likelihood all less tagged alleles are errors)


    private void setLikelihoodThresh() {   // initialize the likelihood ratio cutoffs for quantitative SNP calling
//        System.out.println("\n\nInitializing the cutoffs for quantitative SNP calling likelihood ratio (pHet/pErr) >1\n");
//        System.out.println("totalReadsForSNPInIndiv\tminLessTaggedAlleleCountForHet");
        for (int trials = 0; trials < 2; ++trials) {
            likelihoodRatioThreshAlleleCnt[trials] = 1;
        }
        int lastThresh = 1;       
        for (int trials = 2; trials < likelihoodRatioThreshAlleleCnt.length; ++trials) {
            BinomialDistribution binomHet = new BinomialDistribution(trials, 0.5);
            BinomialDistribution binomErr = new BinomialDistribution(trials, errorRate);
            double LikeRatio;
            try {
                LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                while (LikeRatio <= 1.0) {
                    ++lastThresh;
                    LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                }
                likelihoodRatioThreshAlleleCnt[trials] = lastThresh;
               // double resultDouble = (double)lastThresh/trials; 
               // System.out.println(trials + "\t" + lastThresh + "\t" +  LikeRatio + "\t" + resultDouble);
               // System.out.println(trials + "\t" + lastThresh); // if just want to see total depth and min num for het call
                
            } catch (Exception e) {
                System.err.println("Error in the TagsAtLocus.BinomialDistributionImpl");
            }
        }
        System.out.println("\n");
    }

    public BasicGenotypeMergeRule(double errorRate) {
        this.errorRate=errorRate;
        likelihoodRatioThreshAlleleCnt = new int[maxCountAtGeno];
        setLikelihoodThresh();
    }


    @Override
    public boolean isMergePossible() {
        return true;
    }

    @Override
    public byte mergeCalls(byte geno1, byte geno2) {
        if(geno1==geno2) return geno1;
        if(geno1==UNKNOWN_DIPLOID_ALLELE) return geno2;
        if(geno2==UNKNOWN_DIPLOID_ALLELE) return geno1;
        return UNKNOWN_DIPLOID_ALLELE;
    }

    @Override
    public byte[] mergeWithDepth(byte[] geno1depths, byte[] geno2depths) {
        if(geno1depths.length!=geno2depths.length) throw new IllegalStateException("Depth arrays must be same length");
        byte[] result=new byte[geno1depths.length];
        for (int i=0; i<result.length; i++) {
            result[i]=AlleleDepthUtil.addByteDepths(geno1depths[i],geno2depths[i]);
        }
        return result;
    }

    @Override
    public byte callBasedOnDepth(byte[] genoDepths) {
        return resolveHetGeno(AlleleDepthUtil.depthByteToInt(genoDepths));
    }

    @Override
    public byte callBasedOnDepth(int[] genoDepths) {
        return resolveHetGeno(genoDepths);
    }

    private byte resolveHetGeno(int[] depths) {
        int max = 0;
        byte maxAllele = GenotypeTable.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = GenotypeTable.UNKNOWN_ALLELE;
        for (int a = 0; a < depths.length; a++) {
            if (depths[a] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = depths[a];
                maxAllele = (byte) a;
            } else if (depths[a] > nextMax) {
                nextMax = depths[a];
                nextMaxAllele = (byte) a;
            }
        }
        // use the Glaubitz/Buckler LR method (if binomialPHet/binomialPErr > 1, call it a het)
        int totCount = max + nextMax;
        if (totCount < maxCountAtGeno) {
            if (nextMax < likelihoodRatioThreshAlleleCnt[totCount]) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        } else {
            // Must cast to double.  In java, int1/int2 where int1 < int2 always returns 0
            // as the conversion to float happens AFTER the division.
            double resultDouble = (double)nextMax/totCount;
            if (resultDouble < 0.1) {               
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        }
    }
}
