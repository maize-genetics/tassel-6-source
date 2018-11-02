package net.maizegenetics.dna.snp.genotypecall;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.map.DonorHaplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.NavigableSet;

/**
 * Projection genotype use defined haplotypes and breakpoints that point to a high
 * density genotypes (base GenotypeTable). These are used to efficiently store and
 * connect low density maps with imputed high density genotypes.
 * <p>
 * </p>
 * The alignment built by this builder is a CoreGenotypeTable with a
 * ProjectionGenotypeCallTable. The taxa indices come from the projection alignment file,
 * while the site indices are the same as the base GenotypeTable. TODO this
 * implement a projection interface with the getDonorHaplotypes method
 *
 * @author Ed Buckler
 */
public class ProjectionGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeTable myBaseGenoTable;  //high density marker alignment that is being projected. It was suggested that this
    //just have a pointer to a genotype, which would work, excepting for saving the file, when the base taxa names are needed.
    private ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints;

    private enum BaseMode {

        General, Site, Taxa
    };
    private BaseMode currMode = BaseMode.Taxa;
    private ArrayList<RangeMap<Integer, DonorSiteHaps>> breakMaps;
    private byte[] donorForCachedSite;
    private byte[] projForCachedTaxon;
    private int cachedSite = -1;
    private int cachedTaxon = -1;
    int[] primDSH; //startSite,endSite,parent1,parent2 array for the

  public ProjectionGenotypeCallTable(GenotypeTable hdAlign, ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints) {
        super(allBreakPoints.size(), hdAlign.numberOfSites(), false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myBaseGenoTable= hdAlign;
        this.allBreakPoints = allBreakPoints;
        breakMaps = new ArrayList<>(numberOfTaxa());
        for (NavigableSet<DonorHaplotypes> allBreakPoint : allBreakPoints) {
            RangeMap<Integer, DonorSiteHaps> tRM = TreeRangeMap.create();
            for (DonorHaplotypes dh : allBreakPoint) {
                int[] siteRange = siteRangeForDonor(dh);
                DonorSiteHaps dsh = new DonorSiteHaps(siteRange[0], siteRange[1], dh.getParent1index(), dh.getParent2index());
                tRM.put(Range.closed(siteRange[0], siteRange[1]), dsh);
                //TODO consider putting in blank range maps
            }
            breakMaps.add(tRM);
        }
        primDSH = new int[myTaxaCount * 4];
        Arrays.fill(primDSH, Integer.MIN_VALUE);
    }

    public NavigableSet<DonorHaplotypes> getDonorHaplotypes(int taxon) {
        return allBreakPoints.get(taxon);
    }

    private int[] siteRangeForDonor(DonorHaplotypes dh) {
        int start = myBaseGenoTable.siteOfPhysicalPosition(dh.getStartPosition(), dh.getChromosome());
        if (start < 0) {
            start = -(start + 1);
        }
        int end = myBaseGenoTable.siteOfPhysicalPosition(dh.getEndPosition(), dh.getChromosome());
        if (end < 0) {
            end = -(end + 1);
        }
        return new int[]{start, end};
    }

    @Override
    public byte genotype(int taxon, int site) {
        if (currMode == BaseMode.Site) {
            return getBaseSite(taxon, site);
        }
//        if(currMode==BaseMode.Taxa) {return getBaseTaxon(taxon, site);}
        return getBaseGeneral(taxon, site);
    }
    
    public int[] taxonDonors(int taxon, int site) {
        int primPos = taxon << 2;
        if ((site < primDSH[primPos++]) || (site > primDSH[primPos++])) {
            DonorSiteHaps currentDSH = breakMaps.get(taxon).get(site);
            primPos = taxon << 2;
            primDSH[primPos++] = currentDSH.getStartSite(); //NOTE:-used to be currentDSH.getStartSite(), but it threw an exception
            primDSH[primPos++] = currentDSH.getStartSite(); //NOTE:-used to be currentDSH.getEndSite(), but it threw an exception
            primDSH[primPos++] = currentDSH.getParent1index();
            primDSH[primPos++] = currentDSH.getParent2index();
            primPos = (taxon << 2) + 2;
            //TODO consider null
        }
        return new int[]{primDSH[primPos], primDSH[primPos + 1]};
    }

    /**
     * Returns the high density base genotypeTable of the projection genotypeTable.
     *
     * @return base GenotypeTable
     */
    public GenotypeTable getBaseGenotypeTable() {
        return myBaseGenoTable;
    }

    private byte getBaseGeneral(int taxon, int site) {
        DonorSiteHaps currentDSH = breakMaps.get(taxon).get(site);
        if (currentDSH == null) {
            return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        }
        byte p1 = myBaseGenoTable.genotype(currentDSH.getParent1index(), site);
        byte p2 = myBaseGenoTable.genotype(currentDSH.getParent2index(), site);
        return GenotypeTableUtils.getUnphasedDiploidValueNoHets(p1, p2);
    }

    //Currently this is no faster than general genotype.  it should be possible to make this faster.
    private byte getBaseTaxon(int taxon, int site) {
        if (taxon != cachedTaxon) {
            projForCachedTaxon = new byte[mySiteCount];
            Arrays.fill(projForCachedTaxon, GenotypeTable.RARE_DIPLOID_ALLELE);
            cachedTaxon = taxon;
        }
        byte result = projForCachedTaxon[site];
        if (result == GenotypeTable.RARE_DIPLOID_ALLELE) {
            DonorSiteHaps currentDSH = breakMaps.get(taxon).get(site);
            byte[] r = myBaseGenoTable.genotypeRange(currentDSH.getParent1index(), currentDSH.getStartSite(), currentDSH.getEndSite() + 1);
            System.arraycopy(r, 0, projForCachedTaxon, currentDSH.getStartSite(), r.length);
            result = projForCachedTaxon[site];
        }
        return result;
    }

    private byte getBaseSite(int taxon, int site) {
        //test transpose problems
        if (site != cachedSite) {
            donorForCachedSite = myBaseGenoTable.genotypeMatrix().genotypeForAllTaxa(site);
            cachedSite = site;
        }
        int primPos = taxon << 2;
        if ((site < primDSH[primPos++]) || (site > primDSH[primPos++])) {
            DonorSiteHaps currentDSH = breakMaps.get(taxon).get(site);
            primPos = taxon << 2;
            primDSH[primPos++] = currentDSH.getStartSite();
            primDSH[primPos++] = currentDSH.getEndSite();
            primDSH[primPos++] = currentDSH.getParent1index();
            primDSH[primPos] = currentDSH.getParent2index();
            primPos = (taxon << 2) + 2;
            //TODO consider null
        }
        //       if(primDSH[primPos]==primDSH[primPos+1]) return donorForCachedSite[primDSH[primPos]];
        return GenotypeTableUtils.getUnphasedDiploidValueNoHets(donorForCachedSite[primDSH[primPos]], donorForCachedSite[primDSH[primPos + 1]]);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype(taxon, site));
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBaseGenoTable.genotypeMatrix().transposeData(siteInnerLoop);
        if (siteInnerLoop) {
            currMode = BaseMode.Site;
        } else {
            currMode = BaseMode.General;
        }

    }

    @Override
    public boolean isSiteOptimized() {
        if (currMode == BaseMode.Site) {
            return false;
        } else {
            return true;
        }
    }

    private class DonorSiteHaps {

        private final int startSite;
        private final int endSite;
        private final int parent1index;
        private final int parent2index;

        private DonorSiteHaps(int startSite, int endSite, int parent1index, int parent2index) {
            this.startSite = startSite;
            this.endSite = endSite;
            this.parent1index = parent1index;
            this.parent2index = parent2index;
        }

        private int getStartSite() {
            return startSite;
        }

        private int getEndSite() {
            return endSite;
        }

        private int getParent1index() {
            return parent1index;
        }

        private int getParent2index() {
            return parent2index;
        }

        private boolean containsSite(int site) {
            if ((site < startSite) || (site > endSite)) {
                return false;
            }
            return true;
        }
    }

}
