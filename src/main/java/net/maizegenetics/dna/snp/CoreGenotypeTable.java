/*
 * CoreGenotypeTable
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.bit.DynamicBitStorage;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.GeneralAnnotationStorage;

import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Basic implementation of a {@link GenotypeTable}. Use the GenotypeTableBuilder
 * to construct.
 *
 * @see GenotypeTable
 * @see GenotypeTableBuilder
 *
 * @author Terry Casstevens
 */
public class CoreGenotypeTable implements GenotypeTable {

    private static final Logger myLogger = Logger.getLogger(CoreGenotypeTable.class);
    private final GenotypeCallTable myGenotype;
    private final Map<WHICH_ALLELE, BitStorage> myBitStorage = new HashMap<>();
    private final PositionList myPositionList;
    private final TaxaList myTaxaList;
    private final AlleleProbability myAlleleProbability;
    private final ReferenceProbability myReferenceProbabily;
    private final AlleleDepth myAlleleDepth;
    private final Dosage myDosage;
    private final int mySiteCount;
    private final int myTaxaCount;
    private final GeneralAnnotationStorage myAnnotations;
    private final Translate myGenotypeTranslate;

    CoreGenotypeTable(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, AlleleDepth alleleDepth, AlleleProbability alleleProbability, ReferenceProbability referenceProbability, Dosage dosage, GeneralAnnotationStorage annotations, Translate genotypeTranslate) {
        myPositionList = positionList;
        myTaxaList = taxaList;
        mySiteCount = myPositionList.numberOfSites();
        myTaxaCount = myTaxaList.numberOfTaxa();

        if ((genotype != null) && (genotype.numberOfTaxa() != myTaxaCount)) {
            throw new IllegalArgumentException("CoreGenotypeTable: init: genotype number of taxa: " + genotype.numberOfTaxa() + " doesn't match taxa list: " + myTaxaCount);
        }
        if ((genotype != null) && (genotype.numberOfSites() != mySiteCount)) {
            throw new IllegalArgumentException("CoreGenotypeTable: init: genotype number of sites: " + genotype.numberOfSites() + " doesn't match position list: " + mySiteCount);
        }

        myGenotype = genotype;
        myAlleleDepth = alleleDepth;
        myAlleleProbability = alleleProbability;
        myReferenceProbabily = referenceProbability;
        myDosage = dosage;
        myAnnotations = annotations;
        myGenotypeTranslate = genotypeTranslate;
    }

    CoreGenotypeTable(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, AlleleDepth alleleDepth, AlleleProbability alleleProbability, ReferenceProbability referenceProbability, Dosage dosage, GeneralAnnotationStorage annotations) {
        this(genotype, positionList, taxaList, alleleDepth, alleleProbability, referenceProbability, dosage, annotations, null);
    }

    CoreGenotypeTable(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList) {
        this(genotype, positionList, taxaList, null, null, null, null, null);
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
        return myGenotype.genotype(taxon, myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome));
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeRange(taxon, startSite, endSite);
    }

    @Override
    public byte[] genotypeAllTaxa(int site) {
        return myGenotype.genotypeForAllTaxa(site);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        return myGenotype.genotypeAllSites(taxon);
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, WHICH_ALLELE allele) {
        return bitStorage(allele).allelePresenceForAllSites(taxon);
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, WHICH_ALLELE allele, int startBlock, int endBlock) {
        return bitStorage(allele).allelePresenceForSitesBlock(taxon, startBlock, endBlock);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, WHICH_ALLELE allele) {
        return bitStorage(allele).haplotypeAllelePresenceForAllSites(taxon, firstParent);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, WHICH_ALLELE allele) {
        return bitStorage(allele).haplotypeAllelePresenceForAllTaxa(site, firstParent);
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, WHICH_ALLELE allele, int startBlock, int endBlock) {
        return bitStorage(allele).haplotypeAllelePresenceForSitesBlock(taxon, firstParent, startBlock, endBlock);
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
        return myPositionList.allele(WHICH_ALLELE.Reference, site);
    }

    @Override
    public byte alternateAllele(int site) {
        return myPositionList.allele(WHICH_ALLELE.Alternate, site);
    }

    @Override
    public byte[] referenceAlleles(int startSite, int endSite) {
        return myPositionList.alleles(WHICH_ALLELE.Reference, startSite, endSite);
    }

    @Override
    public byte[] referenceAlleleForAllSites() {
        return myPositionList.alleleForAllSites(WHICH_ALLELE.Reference);
    }

    @Override
    public boolean hasReference() {
        return myPositionList.hasReference();
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        return myGenotype.isHeterozygous(taxon, site);
    }

    @Override
    public int heterozygousCount(int site) {
        return myGenotype.heterozygousCount(site);
    }

    @Override
    public PositionList positions() {
        return myPositionList;
    }

    @Override
    public String siteName(int site) {
        return myPositionList.siteName(site);
    }

    @Override
    public int numberOfSites() {
        return mySiteCount;
    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        return myPositionList.chromosomeSiteCount(chromosome);
    }

    @Override
    public int[] firstLastSiteOfChromosome(Chromosome chromosome) {
        return myPositionList.startAndEndOfChromosome(chromosome);
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaCount;
    }

    @Override
    public int chromosomalPosition(int site) {
        return myPositionList.chromosomalPosition(site);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
    }

    @Override
    public int[] physicalPositions() {
        return myPositionList.physicalPositions();
    }

    @Override
    public String chromosomeName(int site) {
        return myPositionList.chromosomeName(site);
    }

    @Override
    public Chromosome chromosome(int site) {
        return myPositionList.chromosome(site);
    }

    @Override
    public Chromosome chromosome(String name) {
        return myPositionList.chromosome(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myPositionList.chromosomes();
    }

    @Override
    public int numChromosomes() {
        return myPositionList.numChromosomes();
    }

    @Override
    public int[] chromosomesOffsets() {
        return myPositionList.chromosomesOffsets();
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        Set<SITE_SCORE_TYPE> result = new HashSet<>();
        if (hasAlleleProbabilities()) {
            result.addAll(myAlleleProbability.siteScoreTypes());
        }
        if (hasDosage()) {
            result.addAll(myDosage.siteScoreTypes());
        }
        if (hasReferenceProbablity()) {
            result.addAll(myReferenceProbabily.siteScoreTypes());
        }
        return result;
    }

    @Override
    public int indelSize(int site) {
        return myPositionList.indelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return myPositionList.isIndel(site);
    }

    @Override
    public boolean isAllPolymorphic() {
        return myGenotype.isAllPolymorphic();
    }

    @Override
    public boolean isPolymorphic(int site) {
        return myGenotype.isPolymorphic(site);
    }

    @Override
    public byte majorAllele(int site) {
        return myGenotype.majorAllele(site);
    }

    @Override
    public String majorAlleleAsString(int site) {
        return myGenotype.majorAlleleAsString(site);
    }

    @Override
    public byte minorAllele(int site) {
        return myGenotype.minorAllele(site);
    }

    @Override
    public String minorAlleleAsString(int site) {
        return myGenotype.minorAlleleAsString(site);
    }

    @Override
    public byte[] minorAlleles(int site) {
        return myGenotype.minorAlleles(site);
    }

    @Override
    public byte[] alleles(int site) {
        return myGenotype.alleles(site);
    }

    @Override
    public double minorAlleleFrequency(int site) {
        return myGenotype.minorAlleleFrequency(site);
    }

    @Override
    public double majorAlleleFrequency(int site) {
        return myGenotype.majorAlleleFrequency(site);
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
        return myPositionList.genomeVersion();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myPositionList.isPositiveStrand(site);
    }

    @Override
    public GenotypeTable[] compositeAlignments() {
        return new GenotypeTable[]{this};
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myGenotype.allelesSortedByFrequency(site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        return myGenotype.genosSortedByFrequency(site);
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] alleleDefinitions() {
        return myGenotype.alleleDefinitions();
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return myGenotype.alleleDefinitions(site);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return myGenotype.genotypeAsString(site, value);
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myGenotype.diploidAsString(site, value);
    }

    @Override
    public int maxNumAlleles() {
        return myGenotype.maxNumAlleles();
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        return myGenotype.totalGametesNonMissingForSite(site);
    }

    @Override
    public int totalNonMissingForSite(int site) {
        return myGenotype.totalNonMissingForSite(site);
    }

    @Override
    public int minorAlleleCount(int site) {
        return myGenotype.minorAlleleCount(site);
    }

    @Override
    public int majorAlleleCount(int site) {
        return myGenotype.majorAlleleCount(site);
    }

    @Override
    public Object[][] genoCounts() {
        return myGenotype.genoCounts();
    }

    @Override
    public Object[][] majorMinorCounts() {
        return myGenotype.majorMinorCounts();
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        return myGenotype.totalGametesNonMissingForTaxon(taxon);
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        return myGenotype.heterozygousCountForTaxon(taxon);
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        return myGenotype.totalNonMissingForTaxon(taxon);
    }

    @Override
    public boolean hasGenotype() {
        return myGenotype != null;
    }

    @Override
    public boolean hasDepth() {
        return myAlleleDepth != null;
    }

    @Override
    public boolean hasAlleleProbabilities() {
        return myAlleleProbability != null;
    }

    @Override
    public boolean hasReferenceProbablity() {
        return myReferenceProbabily != null;
    }

    @Override
    public boolean hasDosage() {
        return myDosage != null;
    }

    @Override
    public AlleleDepth depth() {
        return myAlleleDepth;
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        return myAlleleDepth.values(taxon, site);
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        switch (scope) {
            case Frequency:
                return alleles(site);
            default:
                myLogger.warn("getAllelesByScope: Unsupported type: " + scope);
                return null;
        }
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, WHICH_ALLELE allele) {
        return bitStorage(allele).allelePresenceForAllTaxa(site);
    }

    @Override
    public BitStorage bitStorage(WHICH_ALLELE allele) {
        BitStorage result = myBitStorage.get(allele);
        if (result != null) {
            return result;
        }

        switch (allele) {
            case Major:
                result = new DynamicBitStorage(myGenotype, allele, myGenotype.majorAlleleForAllSites());
                break;
            case Minor:
                result = new DynamicBitStorage(myGenotype, allele, myGenotype.minorAlleleForAllSites());
                break;
            case Minor2:
                result = new DynamicBitStorage(myGenotype, allele, myGenotype.thirdAlleleForAllSites());
                break;
            case Unknown:
                result = DynamicBitStorage.getUnknownInstance(myGenotype);
                break;
            default:
                myLogger.warn("bitStorage: Unsupported allele: " + allele);
                return null;
        }

        myBitStorage.put(allele, result);
        return result;
    }

    @Override
    public AlleleProbability alleleProbability() {
        return myAlleleProbability;
    }

    @Override
    public float alleleProbability(int taxon, int site, SITE_SCORE_TYPE type) {
        return myAlleleProbability.value(taxon, site, type);
    }

    @Override
    public ReferenceProbability referenceProbability() {
        return myReferenceProbabily;
    }

    @Override
    public float referenceProbability(int taxon, int site) {
        return myReferenceProbabily.value(taxon, site);
    }

    @Override
    public Dosage dosage() {
        return myDosage;
    }

    @Override
    public byte dosage(int taxon, int site) {
        return myDosage.value(taxon, site);
    }

    @Override
    public GeneralAnnotationStorage annotations() {
        return myAnnotations;
    }

    @Override
    public Stream<Byte> streamGenotype() {
        return myGenotype.stream();
    }

    @Override
    public Stream<Byte> streamGenotype(int taxon) {
        return myGenotype.stream(taxon);
    }

    @Override
    public boolean hasSiteTranslations() {
        return myGenotypeTranslate.hasSiteTranslations();
    }

    @Override
    public int[] siteTranslations() {
        return myGenotypeTranslate.siteTranslations();
    }

}
