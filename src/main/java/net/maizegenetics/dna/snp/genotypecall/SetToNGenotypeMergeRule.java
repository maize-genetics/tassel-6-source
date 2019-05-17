package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_GENOTYPE;

public class SetToNGenotypeMergeRule implements GenotypeMergeRule {
    @Override
    public boolean isMergePossible() {
        return true;
    }

    @Override
    public byte mergeCalls(byte geno1, byte geno2) {
        if(geno1==geno2) return geno1;
        if(geno1== UNKNOWN_GENOTYPE) return geno2;
        if(geno2== UNKNOWN_GENOTYPE) return geno1;
        return UNKNOWN_GENOTYPE;
    }

    @Override
    public byte[] mergeWithDepth(byte[] geno1depths, byte[] geno2depths) {
        if(geno1depths.length!=geno2depths.length) throw new IllegalStateException("Depth arrays must be same length");
        byte[] result=new byte[geno1depths.length];
        for (int i=0; i<result.length; i++) {
            result[i]= AlleleDepthUtil.addByteDepths(geno1depths[i],geno2depths[i]);
        }
        return result;
    }

    @Override
    public byte callBasedOnDepth(byte[] genoDepths) {
        return callBasedOnDepth(AlleleDepthUtil.depthByteToInt(genoDepths));
    }

    @Override
    public byte callBasedOnDepth(int[] genoDepths) {
        int max = 0;
        byte maxAllele = GenotypeTable.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = GenotypeTable.UNKNOWN_ALLELE;
        for (int a = 0; a < genoDepths.length; a++) {
            if (genoDepths[a] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = genoDepths[a];
                maxAllele = (byte) a;
            }
            else if (genoDepths[a] > nextMax) {
                nextMax = genoDepths[a];
                nextMaxAllele = (byte) a;
            }
        }

        if(nextMaxAllele == GenotypeTable.UNKNOWN_ALLELE ) {
            return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
        }
        else {
            return GenotypeTable.UNKNOWN_GENOTYPE;
        }
    }
}
