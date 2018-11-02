/*
 *  HybridGenotypeCallTable
 * 
 *  Created on May 18, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.util.List;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;

/**
 *
 * @author Terry Casstevens
 */
public class HybridGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable myBase;
    private final int[] myFirstParents;
    private final int[] mySecondParents;

    public HybridGenotypeCallTable(GenotypeCallTable base, List<Integer> firstParents, List<Integer> secondParents) {
        super(firstParents.size(), base.numberOfSites(), base.isPhased(), base.alleleDefinitions());
        int numParents = firstParents.size();
        if (numParents != secondParents.size()) {
            throw new IllegalArgumentException("HybridGenotypeCallTable: init: number of first and second parents must be the same");
        }
        myFirstParents = new int[numParents];
        mySecondParents = new int[numParents];
        for (int i = 0; i < numParents; i++) {
            myFirstParents[i] = firstParents.get(i);
            mySecondParents[i] = secondParents.get(i);
        }
        myBase = base;
    }

    @Override
    public byte genotype(int taxon, int site) {
        byte first = myBase.genotype(myFirstParents[taxon], site);
        byte second = myBase.genotype(mySecondParents[taxon], site);
        if (GenotypeTableUtils.isHomozygous(first) && GenotypeTableUtils.isHomozygous(second)) {
            return (byte) ((first & 0xF0) | (second & 0xF));
        } else {
            return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myBase.diploidAsString(site, genotype(taxon, site));
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        return myBase.genotypeAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myBase.diploidAsString(site, value);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBase.transposeData(siteInnerLoop);
    }

    @Override
    public boolean isSiteOptimized() {
        return myBase.isSiteOptimized();
    }

}
