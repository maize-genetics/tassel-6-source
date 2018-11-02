package net.maizegenetics.dna.snp.genotypecall;

import java.util.List;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;

public class DifferenceGenotypeCallTable extends AbstractGenotypeCallTable {
    private final GenotypeCallTable myBase;
    private final int[] myHybrids;
    private final int[] myParents;

    public DifferenceGenotypeCallTable(GenotypeCallTable base, List<Integer> hybrids, List<Integer> parents) {
        super(hybrids.size(), base.numberOfSites(), base.isPhased(), base.alleleDefinitions());
        int numParents = parents.size();
        if (numParents != hybrids.size()) {
            throw new IllegalArgumentException("DifferenceGenotypeCallTable: init: number of parents and hybrids must be the same");
        }
        myHybrids = hybrids.stream().mapToInt(Integer::intValue).toArray();
        myParents = parents.stream().mapToInt(Integer::intValue).toArray();
        myBase = base;
    }

    @Override
    public byte genotype(int taxon, int site) {
        byte hybrid = myBase.genotype(myHybrids[taxon], site);
        byte parent = myBase.genotype(myParents[taxon], site);
//        if (GenotypeTableUtils.isHomozygous(first) && GenotypeTableUtils.isHomozygous(second)) {
//            return (byte) ((first & 0xF0) | (second & 0xF));
//        } else {
//            return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
//        }
        if (GenotypeTableUtils.isHomozygous(hybrid)) {
        	if (hybrid == parent) return hybrid;
        	else return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        } else {
        	if (GenotypeTableUtils.isHomozygous(parent)) {
        		byte[] hybridAlleles = GenotypeTableUtils.getDiploidValues(hybrid);
        		byte parentAllele = GenotypeTableUtils.getDiploidValues(parent)[0];
        		if (parentAllele == hybridAlleles[0]) return GenotypeTableUtils.getDiploidValue(hybridAlleles[1], hybridAlleles[1]);
        		else return GenotypeTableUtils.getDiploidValue(hybridAlleles[0], hybridAlleles[0]);
        	} else return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
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
