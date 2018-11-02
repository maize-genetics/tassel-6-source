package net.maizegenetics.dna.snp;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.stream.Stream;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.MergedGenotypeCallTable;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.GeneralAnnotationStorage;

public class MergedGenotypeTable implements GenotypeTable {

    private final GenotypeTable[] myGenotypeTables;
    private final GenotypeCallTable myGenotype;
    private final GenotypeCallTable[] myGenotypes;
    private final Map<Chromosome, GenotypeTable> myChromosomes = new HashMap<>();
    private Chromosome[] myChromosomesList;
    private final TaxaList myTaxaList;
    private String[][] myAlleleStates;
    private PositionList myPositions = null;
    private boolean mergedMode = false;

    public MergedGenotypeTable(GenotypeTable[] genoTables, TaxaList taxaList, PositionList positionList) {
        myGenotypeTables = genoTables;
        myTaxaList = taxaList;
        myPositions = positionList;
        myGenotypes = new GenotypeCallTable[myGenotypeTables.length];
        for (int i = 0; i < myGenotypeTables.length; i++) {
            myGenotypes[i] = myGenotypeTables[i].genotypeMatrix();
        }
        //Create Taxon Mapping

        int[][] taxonMap = new int[taxaList.size()][genoTables.length];
        for (int i = 0; i < myTaxaList.numberOfTaxa(); i++) {
            int[] taxaIndexPerTable = new int[myGenotypeTables.length];
            for (int j = 0; j < taxaIndexPerTable.length; j++) {
                if (myGenotypeTables[j].taxa().indexOf(myTaxaList.get(i)) == -1) {
                    taxaIndexPerTable[j] = -1;
                } else {
                    taxaIndexPerTable[j] = myGenotypeTables[j].taxa().indexOf(myTaxaList.get(i));
                }
            }
            taxonMap[i] = taxaIndexPerTable;
        }

        //Create Position Mapping
        int[][] positionMap = new int[myPositions.size()][myGenotypeTables.length];
        for (int i = 0; i < positionMap.length; i++) {
            int[] posIndexPerTable = new int[positionMap[i].length];
            for (int j = 0; j < posIndexPerTable.length; j++) {
                if (myGenotypeTables[j].positions().indexOf(myPositions.get(i)) == -1) {
                    posIndexPerTable[j] = -1;
                } else {
                    posIndexPerTable[j] = myGenotypeTables[j].positions().indexOf(myPositions.get(i));
                }
            }
            positionMap[i] = posIndexPerTable;
        }

        myGenotype = MergedGenotypeCallTable.getInstance(myGenotypes, taxonMap, positionMap);

    }

    public static GenotypeTable getInstance(GenotypeTable[] genoTables, BiFunction taxaMergeRule, BiFunction positionMergeRule) {
        TaxaList txl = (TaxaList) taxaMergeRule.apply(genoTables[0].taxa(), genoTables[1].taxa());
        PositionList posList = (PositionList) positionMergeRule.apply(genoTables[0].positions(), genoTables[1].positions());
        return new MergedGenotypeTable(genoTables, txl, posList);
    }

    @Override
    public boolean hasGenotype() {
        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (myGenotypeTables[i].hasGenotype()) {
                return true;
            }
        }
        return false;
    }

    @Override
    public GenotypeCallTable genotypeMatrix() {
        return myGenotype;
    }

    @Override
    public byte genotype(int taxon, int site) {
        return myGenotype.genotype(taxon, site);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return myGenotype.genotypeArray(taxon, site);
    }

    @Override
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition) {
        return myGenotype.genotype(taxon, myPositions.siteOfPhysicalPosition(physicalPosition, chromosome));
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeRange(taxon, startSite, endSite);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        return myGenotype.genotypeAllSites(taxon);
    }

    @Override
    public byte[] genotypeAllTaxa(int site) {
        byte[] genotypes = new byte[myTaxaList.numberOfTaxa()];
        for (int i = 0; i < genotypes.length; i++) {
            genotypes[i] = myGenotype.genotype(i, site);
        }
        return genotypes;
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("CombineGenotypeTable: getAllelePresenceForAllSites: This operation isn't possible as it spans multiple GenotypeTables.");
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, WHICH_ALLELE allele, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("CombineGenotypeTable: getAllelePresenceForSitesBlock: This operation isn't possible as it spans multiple GenotypeTables.");
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, WHICH_ALLELE allele, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myGenotype.genotypeAsString(taxon, site);

    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeAsStringRange(taxon, startSite, endSite);

    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        return myGenotype.genotypeAsStringRow(taxon);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        return myGenotype.genotypeAsStringArray(taxon, site);
    }

    @Override
    public byte referenceAllele(int site) {
        // Need to combine calls and figure out which one is reference
        return GenotypeTable.UNKNOWN_ALLELE;
    }

    @Override
    public byte alternateAllele(int site) {
        // Need to combine calls and figure out which one is reference
        return GenotypeTable.UNKNOWN_ALLELE;
    }

    @Override
    public byte[] referenceAlleles(int startSite, int endSite) {
        // TODO From CombineGenotypeTable
        //int numSites = endSite - startSite;
        //byte[] result = new byte[numSites];
        //for (int i = 0; i < numSites; i++) {
        //    result[i] = referenceAllele(startSite + i);
        //}
        //return result;
        return null;
    }

    @Override
    public byte[] referenceAlleleForAllSites() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean hasReference() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public int heterozygousCount(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public String siteName(int site) {
        return myPositions.get(site).getSNPID();
    }

    @Override
    public int numberOfSites() {
        return myPositions.size();

    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int[] firstLastSiteOfChromosome(Chromosome chromosome) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaList.size();
    }

    @Override
    public PositionList positions() {
        return myPositions;

    }

    @Override
    public int chromosomalPosition(int site) {
        return myPositions.chromosomalPosition(site);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        return myPositions.siteOfPhysicalPosition(physicalPosition, chromosome);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        return myPositions.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
    }

    @Override
    public int[] physicalPositions() {
        return myPositions.physicalPositions();
    }

    @Override
    public String chromosomeName(int site) {
        return myPositions.chromosomeName(site);
    }

    @Override
    public Chromosome chromosome(int site) {
        return myPositions.chromosome(site);
    }

    @Override
    public Chromosome chromosome(String name) {
        return myPositions.chromosome(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myPositions.chromosomes();
    }

    @Override
    public int numChromosomes() {
        return myPositions.numChromosomes();
    }

    @Override
    public int[] chromosomesOffsets() {
        return myPositions.chromosomesOffsets();
    }

    @Override
    public boolean hasDepth() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean hasAlleleProbabilities() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean hasReferenceProbablity() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean hasDosage() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public AlleleProbability alleleProbability() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public float alleleProbability(int taxon, int site, SITE_SCORE_TYPE type) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public ReferenceProbability referenceProbability() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public float referenceProbability(int taxon, int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public Dosage dosage() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public byte dosage(int taxon, int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int indelSize(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public boolean isIndel(int site) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean isAllPolymorphic() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean isPolymorphic(int site) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public byte majorAllele(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public String majorAlleleAsString(int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public byte minorAllele(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public String minorAlleleAsString(int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public byte[] minorAlleles(int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public byte[] alleles(int site) {
        return myGenotype.alleles(site);
    }

    @Override
    public double minorAlleleFrequency(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double majorAlleleFrequency(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.taxaName(index);
    }

    @Override
    public String genomeVersion() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public GenotypeTable[] compositeAlignments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myGenotype.allelesSortedByFrequency(site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean isPhased() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public boolean retainsRareAlleles() {
        // TODO Auto-generated method stub
        return false;
    }

    @Override
    public String[][] alleleDefinitions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String[] alleleDefinitions(int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String diploidAsString(int site, byte value) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int maxNumAlleles() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int totalNonMissingForSite(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int minorAlleleCount(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int majorAlleleCount(int site) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public Object[][] genoCounts() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Object[][] majorMinorCounts() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public AlleleDepth depth() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, WHICH_ALLELE allele) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public BitStorage bitStorage(WHICH_ALLELE allele) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public GeneralAnnotationStorage annotations() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Stream<Byte> streamGenotype() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public Stream<Byte> streamGenotype(int taxon) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean hasSiteTranslations() {
        return false;
    }

    @Override
    public int[] siteTranslations() {
        return null;
    }

}
