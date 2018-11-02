package net.maizegenetics.dna.snp.genotypecall;

/**
 * Defines the methods for merging the calls from two taxa.  The merge rules need to be defined at the level of
 * genotypic calls and for read depth.  In general if depth is available, it will be used to merge.
 *
 * @author Ed Buckler
 */
public interface GenotypeMergeRule {
    /**
     * Whether merge is even possible
     * @return whether merging is allowed.
     */
    boolean isMergePossible();

    /**
     * Merges diploid genotypic calls into one
     * @param geno1 genotype call of taxa 1
     * @param geno2 genotype call of taxa 2
     * @return merged genotype call
     */
    byte mergeCalls(byte geno1, byte geno2);

    /**
     * Merges sequencing depths of two taxa
     * @param geno1depths allele depth of taxa 1
     * @param geno2depths allele depths of taxa 2
     * @return merged depths
     */
    byte[] mergeWithDepth(byte[] geno1depths, byte[] geno2depths);

    /**
     * Makes a genotypic call based on allele depths
     * @param genoDepths allele depth of taxa
     * @return genotype call
     */
    byte callBasedOnDepth(byte[] genoDepths);

    /**
     * Makes a genotypic call based on allele depths
     * @param genoDepths allele depth of taxa
     * @return genotype call
     */
    byte callBasedOnDepth(int[] genoDepths);
}
