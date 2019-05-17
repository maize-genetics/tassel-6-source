/*
 *  MaskGenotypeCallTable
 *
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.MaskMatrix;
import net.maizegenetics.util.BitSet;

/**
 * @author Terry Casstevens
 */
public class MaskGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable myBase;
    private final MaskMatrix myMask;

    public MaskGenotypeCallTable(GenotypeCallTable base, MaskMatrix mask) {
        super(base.numberOfTaxa(), base.numberOfSites(), base.isPhased(), base.alleleDefinitions());
        myBase = base;
        myMask = mask;
    }

    @Override
    public byte genotype(int taxon, int site) {
        if (myMask.get(taxon, site)) {
            return GenotypeTable.UNKNOWN_GENOTYPE;
        } else {
            return myBase.genotype(taxon, site);
        }
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        BitSet mask = myMask.maskForSite(site);
        byte[] result = myBase.genotypeForAllTaxa(site);
        for (int t = 0; t < numberOfTaxa(); t++) {
            if (mask.fastGet(t)) {
                result[t] = GenotypeTable.UNKNOWN_GENOTYPE;
            }
        }
        return result;
    }

    @Override
    public byte[] genotypeForAllSites(int taxon) {
        BitSet mask = myMask.maskForTaxon(taxon);
        byte[] result = myBase.genotypeForAllSites(taxon);
        for (int s = 0; s < numberOfSites(); s++) {
            if (mask.fastGet(s)) {
                result[s] = GenotypeTable.UNKNOWN_GENOTYPE;
            }
        }
        return result;
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myBase.diploidAsString(site, genotype(taxon, site));
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(genotypeAsString(taxon, i));
        }
        return builder.toString();
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
