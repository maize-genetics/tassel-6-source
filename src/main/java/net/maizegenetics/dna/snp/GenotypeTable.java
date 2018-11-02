package net.maizegenetics.dna.snp;

import java.util.Set;
import java.util.stream.Stream;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.GeneralAnnotationStorage;

/**
 * A representation of the SNP and indel variation for a set of taxa and genomic
 * positions.
 * <p>
 * </p>
 * GenotypeTable always consist of a TaxaList, PositionList, and
 * GenotypeCallTable. Additionally, as needed they also can represent allele in
 * a bit form, with sequencing depth, or other scores (e.g. quality scores).
 * <p>
 * </p>
 * Use GenotypeTableBuilder to create GenotypeTable.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public interface GenotypeTable {

    /**
     * This encoding is used to lump together allele values with frequencies too
     * low to be retained as one of the maximum number of alleles.
     */
    public static byte RARE_ALLELE = 0xE;
    public static byte RARE_DIPLOID_ALLELE = (byte) 0xEE;
    public static String RARE_ALLELE_STR = "Z";
    public static byte UNKNOWN_ALLELE = 0xF;
    public static byte UNKNOWN_DIPLOID_ALLELE = (byte) 0xFF;
    public static String UNKNOWN_ALLELE_STR = "N";
    public static String UNKNOWN_DIPLOID_ALLELE_STR = "N:N";
    public static char UNKNOWN_ALLELE_CHAR = 'N';

    /**
     * Annotations
     */
    public static final String ANNOTATION_DATA_SET_NAME = "dataSetName";
    public static final String ANNOTATION_DATA_SET_DESCRIPTION = "dataSetDescription";

    public static final String[] GENOTYPE_TABLE_ANNOTATIONS = new String[]{ANNOTATION_DATA_SET_NAME, ANNOTATION_DATA_SET_DESCRIPTION};

    /**
     * Possible Components of a Genotype Table.
     */
    public static enum GENOTYPE_TABLE_COMPONENT {

        Genotype, ReferenceProbability, Dosage, Depth, AlleleProbability
    };

    /**
     * This defines the possible allele scope types.
     */
    public static enum ALLELE_SORT_TYPE {

        /**
         * This is the default where alleles are sorted by frequency. Same as
         * alleles().
         */
        Frequency,
        /**
         * This sorts alleles based on there depth value.
         */
        Depth,
        /**
         * This uses the allele frequency of a base/global Genotype table
         * determine sort order of alleles. That Genotype table is usually a
         * superset.
         */
        Global_Frequency,
        /**
         * This sorts alleles based on the reference sequence.
         */
        Reference
    };

    /**
     * Returns whether this Genotype Table has Genotype values.
     *
     * @return has genotype values.
     */
    public boolean hasGenotype();

    /**
     * Returns the immutable Genotype matrix. Taxa and Positions are not part of
     * the matrix. This method is used for copying Genotype tables, when either
     * the Taxa or Positions have changed.
     *
     * @return genotype matrix
     */
    public GenotypeCallTable genotypeMatrix();

    /**
     * Returns diploid value (genotype) for a given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return high four bits generally encode the more frequent allele and the
     * lower four bits encode the less frequent allele.
     */
    public byte genotype(int taxon, int site);

    /**
     * Returns diploid values for given taxon and site. Same values as
     * genotype(), except two values are already separated into two bytes.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return first byte (index 0) holds first allele value in right-most four
     * bits. second byte (index 1) holds second allele value in right-most four
     * bits.
     */
    public byte[] genotypeArray(int taxon, int site);

    /**
     * Returns diploid values for given taxon, chromosome, and physical
     * position. The chromosome and physical position should map to an unique
     * site.
     *
     * @param taxon taxon
     * @param chromosome chromosome
     * @param physicalPosition physical position
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition);

    /**
     * Returns sequence of diploid allele values for given taxon in specified
     * range (end site excluded). Each value in array is what would be returned
     * by genotype().
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return sequence of diploid allele values.
     */
    public byte[] genotypeRange(int taxon, int startSite, int endSite);

    /**
     * Returns sequence of diploid allele values for all sites for given taxon.
     * Each value in array is what would be returned by genotype().
     *
     * @param taxon taxon
     *
     * @return sequence of diploid allele values.
     */
    public byte[] genotypeAllSites(int taxon);

    /**
     * Returns sequence of diploid allele values for all taxa for given site.
     * Each value in array is what would be returned by genotype().
     *
     * @param site site
     *
     * @return sequence of diploid allele values.
     */
    public byte[] genotypeAllTaxa(int site);

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site matches a specific allele (based on frequency). Allele number of
     * value 0 would be the major allele. Allele number of value 1 would be the
     * minor allele. Allele number of value 2 would be the third most frequent
     * allele value and so on.
     *
     * @param taxon taxon
     * @param allele allele
     *
     * @return sequence of true/false values.
     */
    public BitSet allelePresenceForAllSites(int taxon, WHICH_ALLELE allele);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) matches a specific allele.
     *
     * @param taxon taxon
     * @param allele allele
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] allelePresenceForSitesBlock(int taxon, WHICH_ALLELE allele, int startBlock, int endBlock);

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site for given parent matches a specific allele.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param allele allele
     *
     * @return sequence of true/false values.
     */
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, WHICH_ALLELE allele);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param site site
     * @param firstParent true for first parent (false for second parent)
     * @param allele allele
     *
     * @return sequence of true/false values.
     */
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, WHICH_ALLELE allele);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) for given parent matches a specific allele (based on
     * frequency). Allele number of value 0 would be the major allele. Allele
     * number of value 1 would be the minor allele. Allele number of value 2
     * would be the third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param allele allele
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, WHICH_ALLELE allele, int startBlock, int endBlock);

    /**
     * Returns string representation of diploid values returned by genotype()
     * for given taxon and site. The default implementation separates The two
     * allele values with a colon (:) delimiter. Nucleotide data will be
     * represented by a single letter IUPAC code.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return string representation of diploid values.
     */
    public String genotypeAsString(int taxon, int site);

    /**
     * Returns string representation of diploid alleles for given taxon in
     * specified range (end site excluded). Each value in string is what would
     * be returned by genotypeAsString().
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return string representation of alleles in range
     */
    public String genotypeAsStringRange(int taxon, int startSite, int endSite);

    /**
     * Returns string representation of diploid alleles for given taxon for all
     * sites. Each value in string is what would be returned by
     * genotypeAsString().
     *
     * @param taxon taxon
     *
     * @return string representation of alleles
     */
    public String genotypeAsStringRow(int taxon);

    /**
     * Returns string representation of diploid values returned by
     * genotypeArray() for given taxon and site. Same two allele values as
     * genotypeAsString(), except already separated into two Strings.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return string representations of diploid values.
     */
    public String[] genotypeAsStringArray(int taxon, int site);

    /**
     * Return (haploid) reference allele value at given site.
     *
     * @param site site
     *
     * @return first four bits are the encoded reference value.
     */
    public byte referenceAllele(int site);

    /**
     * Return (haploid) alternate allele value at given site.
     *
     * @param site site
     *
     * @return first four bits are the encoded alternate value.
     */
    public byte alternateAllele(int site);

    /**
     * Returns (haploid) reference alleles in specified range. End site not
     * included.
     *
     * @param startSite start site
     * @param endSite end site
     *
     * @return reference allele values.
     */
    public byte[] referenceAlleles(int startSite, int endSite);

    /**
     * Returns (haploid) reference alleles for all sites.
     *
     * @return reference allele values.
     */
    public byte[] referenceAlleleForAllSites();

    /**
     * Return whether this genotype table has defined reference sequence.
     *
     * @return true if this genotype table has reference sequence.
     */
    public boolean hasReference();

    /**
     * Returns whether allele values at given taxon and site are heterozygous.
     * If two values returned by genotype() are different, this will return
     * false.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return whether heterozygous
     */
    public boolean isHeterozygous(int taxon, int site);

    /**
     * Returns number of heterozygous taxa at given site.
     *
     * @param site site
     *
     * @return number of heterozygous taxa
     */
    public int heterozygousCount(int site);

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String siteName(int site);

    /**
     * Returns total number of sites of this genotype table.
     *
     * @return number of sites
     */
    public int numberOfSites();

    /**
     * Return number of sites for given chromosome.
     *
     * @param chromosome chromosome
     *
     * @return number of sites
     */
    public int chromosomeSiteCount(Chromosome chromosome);

    /**
     * Get the first (inclusive) and last (inclusive) site of the specified
     * chromosome in this genotype table.
     *
     * @param chromosome chromosome
     *
     * @return first and last site
     */
    public int[] firstLastSiteOfChromosome(Chromosome chromosome);

    /**
     * Returns number of taxa
     *
     * @return number of taxa
     */
    public int numberOfTaxa();

    /**
     * Return the position list for the genotype table.
     *
     * @return PositionList for all sites.
     */
    public PositionList positions();

    /**
     * Returns the physical position at given site.
     *
     * @param site site
     *
     * @return physical position
     */
    public int chromosomalPosition(int site);

    /**
     * Return site of given physical position in chromosome. If the physical
     * position doesn't exist, (-(insertion point) - 1) is returned. If
     * chromosome is not found, an exception is thrown.
     *
     * @param physicalPosition physical position
     * @param chromosome chromosome. if null, the first chromosome is used.
     *
     * @return index
     */
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome);

    /**
     * Return site of given physical position / SNP ID in chromosome. If the
     * physical position doesn't exist, (-(insertion point) - 1) is returned. If
     * chromosome is not found, an exception is thrown. This is to support
     * multiple sites with the same physical position but different SNP IDs.
     *
     * @param physicalPosition physical position
     * @param chromosome chromosome. if null, the first chromosome is used.
     * @param snpName SNP ID
     *
     * @return index
     */
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName);

    /**
     * Returns all physical positions.
     *
     * @return physical positions.
     */
    public int[] physicalPositions();

    /**
     * Return Chromosome Name for given site.
     *
     * @param site site
     *
     * @return Chromosome Name
     */
    public String chromosomeName(int site);

    /**
     * Return Chromosome for given site.
     *
     * @param site site
     *
     * @return Chromosome
     */
    public Chromosome chromosome(int site);

    /**
     * Return Chromosome with matching name. First to match will be returned.
     *
     * @param name name
     *
     * @return Chromosome
     */
    public Chromosome chromosome(String name);

    /**
     * Return all chromosomes.
     *
     * @return chromosomes
     */
    public Chromosome[] chromosomes();

    /**
     * Return number of chromosomes.
     *
     * @return number of chromosomes
     */
    public int numChromosomes();

    /**
     * Returns starting site for each chromosome.
     *
     * @return starting site for each chromosome.
     */
    public int[] chromosomesOffsets();

    /**
     * Returns true if this genotype table has sequencing depth.
     *
     * @return true if this genotype table has sequencing depth.
     */
    public boolean hasDepth();

    /**
     * Returns true if this genotype table has allele probabilities.
     *
     * @return true if this genotype table has allele probabilities.
     */
    public boolean hasAlleleProbabilities();

    /**
     * Returns true if this genotype table has reference probability.
     *
     * @return true if this genotype table has reference probability.
     */
    public boolean hasReferenceProbablity();

    /**
     * Returns true if this genotype table has dosage.
     *
     * @return true if this genotype table has dosage.
     */
    public boolean hasDosage();

    /**
     * Return the site scores types that have value for this genotype table.
     *
     * @return site score types.
     */
    public Set<SITE_SCORE_TYPE> siteScoreTypes();

    /**
     * Returns allele probability object (null if not present)
     *
     * @return allele probability associated with genotypeTable
     */
    public AlleleProbability alleleProbability();

    public float alleleProbability(int taxon, int site, SITE_SCORE_TYPE type);

    /**
     * Returns reference probability object (null if not present)
     *
     * @return reference probability associated with genotypeTable
     */
    public ReferenceProbability referenceProbability();

    public float referenceProbability(int taxon, int site);

    /**
     * Returns dosage object (null if not present)
     *
     * @return dosage associated with genotypeTable
     */
    public Dosage dosage();

    public byte dosage(int taxon, int site);

    /**
     * Return size of indel at given site.
     *
     * @param site site
     *
     * @return indel size
     */
    public int indelSize(int site);

    /**
     * Returns whether give site is an indel.
     *
     * @param site site
     *
     * @return true if indel
     */
    public boolean isIndel(int site);

    /**
     * Returns whether all sites are polymorphic.
     *
     * @return true if all sites are polymorphic.
     */
    public boolean isAllPolymorphic();

    /**
     * Return whether given site is polymorphic.
     *
     * @param site site
     *
     * @return true if given site is polymorphic.
     */
    public boolean isPolymorphic(int site);

    /**
     * Return most common allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common allele
     */
    public byte majorAllele(int site);

    /**
     * Return most common allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common allele as String
     */
    public String majorAlleleAsString(int site);

    /**
     * Return most common minor allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common minor allele
     */
    public byte minorAllele(int site);

    /**
     * Return most common minor allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common minor allele as String
     */
    public String minorAlleleAsString(int site);

    /**
     * Return all minor alleles at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return all minor alleles
     */
    public byte[] minorAlleles(int site);

    /**
     * Returns all alleles at given site in order of frequency. Gap is included
     * as state. Heterozygous count one for each allele value. Homozygous counts
     * two for the allele value.
     *
     * @param site site
     *
     * @return all alleles
     */
    public byte[] alleles(int site);

    /**
     * Return frequency for most common minor allele at given site. Gap is
     * included as state. Heterozygous count one for each allele value.
     * Homozygous counts two for the allele value.
     *
     * @param site site
     *
     * @return frequency
     */
    public double minorAlleleFrequency(int site);

    /**
     * Return frequency for major allele at given site. Gap is included as
     * state. Heterozygous count one for each allele value. Homozygous counts
     * two for the allele value.
     *
     * @param site site
     *
     * @return frequency
     */
    public double majorAlleleFrequency(int site);

    /**
     * Return taxa list of this genotype table.
     *
     * @return taxa list.
     */
    public TaxaList taxa();

    /**
     * Return taxa name at given index.
     *
     * @param index
     *
     * @return taxa name
     */
    public String taxaName(int index);

    /**
     * Gets the Genome Assembly.
     *
     * @return the genome assembly.
     */
    public String genomeVersion();

    /**
     * Return whether is positive strand at given site.
     *
     * @param site site
     *
     * @return whether is positive strand.
     */
    public boolean isPositiveStrand(int site);

    /**
     * Returns individual genotype tables within this genotype table.
     *
     * @return list of genotype tables.
     */
    public GenotypeTable[] compositeAlignments();

    /**
     * Return sorted list of alleles from highest frequency to lowest at given
     * site in genotype table. Resulting double dimension array holds alleles
     * (bytes) in result[0]. And the counts are in result[1]. Counts haploid
     * values twice and diploid values once. Higher ploids are not supported.
     *
     * @param site site
     *
     * @return sorted list of alleles and counts
     */
    public int[][] allelesSortedByFrequency(int site);

    /**
     * Return sorted list of diploid vales from highest frequency to lowest at
     * given site in genotype table. Resulting double dimension array holds
     * diploids (Strings) in result[0]. And the counts are in result[1]
     * (Integers).
     *
     * @param site site
     *
     * @return sorted list of diploids and counts
     */
    public Object[][] genosSortedByFrequency(int site);

    /**
     * Returns whether this genotype table is phased.
     *
     * @return true if phased.
     */
    public boolean isPhased();

    /**
     * Returns true if this genotype table retains rare alleles. If false, rare
     * alleles are recorded as unknown.
     *
     * @return whether rare alleles are retained.
     */
    public boolean retainsRareAlleles();

    /**
     * Returns allele values as strings for all sites. The first dimension of
     * the array indexes the sites. The second dimension indexes the allele
     * values for given site. The indices for the allele values are used as the
     * codes to store data. These codes (indices) are returned by the genotype()
     * methods. If only one array of allele values is returned, that is the
     * encoding for all sites.
     *
     * @return allele values for all sites.
     */
    public String[][] alleleDefinitions();

    /**
     * Same as alleleDefinitions() for only one site.
     *
     * @param site site
     *
     * @return allele values for given site.
     */
    public String[] alleleDefinitions(int site);

    /**
     * Returns String representation of allele value at site.
     *
     * @param site site
     * @param value allele value
     *
     * @return String representation
     */
    public String genotypeAsString(int site, byte value);

    /**
     * Returns String representation of diploid allele value at site.
     *
     * @param site site
     * @param value diploid allele value
     *
     * @return String representation
     */
    public String diploidAsString(int site, byte value);

    /**
     * Return max number of alleles defined for any given site.
     *
     * @return max number of alleles.
     */
    public int maxNumAlleles();

    /**
     * Returns total number of non-missing allele values for given site. This
     * can be twice the number of taxa, as diploid values are supported.
     *
     * @param site site
     * @return number of non-missing allele values.
     */
    public int totalGametesNonMissingForSite(int site);

    /**
     * Returns total number of non-missing taxa for given site. Taxa are
     * considered missing only if both allele values are Unknown (N).
     *
     * @param site site
     *
     * @return number of non-missing taxa..
     */
    public int totalNonMissingForSite(int site);

    /**
     * Returns the minor allele count for given site.
     *
     * @param site site
     * @return minor allele count
     */
    public int minorAlleleCount(int site);

    /**
     * Returns the major allele count for given site.
     *
     * @param site site
     * @return major allele count
     */
    public int majorAlleleCount(int site);

    /**
     * Returns counts of all diploid combinations from highest frequency to
     * lowest for whole genotype table. Resulting double dimension array holds
     * diploids (Strings) in result[0]. And the counts are in result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] genoCounts();

    /**
     * Returns counts of all major/minor allele combinations from highest
     * frequency to lowest for whole genotype table. Resulting double dimension
     * array holds major/minor allele (Strings) in result[0]. And the counts are
     * in result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] majorMinorCounts();

    /**
     * Returns total number of non-missing allele values for given taxon. This
     * can be twice the number of sites, as diploid values are supported.
     *
     * @param taxon taxon
     *
     * @return number of non-missing allele values.
     */
    public int totalGametesNonMissingForTaxon(int taxon);

    /**
     * Returns number of heterozygous sites at given taxon.
     *
     * @param taxon taxon
     *
     * @return number of heterozygous sites
     */
    public int heterozygousCountForTaxon(int taxon);

    /**
     * Returns total number of non-missing sites for given taxon. Sites are
     * considered missing only if both allele values are Unknown (N).
     *
     * @param taxon taxon
     *
     * @return number of non-missing sites.
     */
    public int totalNonMissingForTaxon(int taxon);

    /**
     * Returns allele depth object (null if not present)
     *
     * @return allele depth associated with genotypeTable
     */
    public AlleleDepth depth();

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return two counts
     */
    public int[] depthForAlleles(int taxon, int site);

    /**
     * Returns all alleles at given site in order defined by scope.
     *
     * @param scope scope
     * @param site site
     *
     * @return alleles
     */
    public byte[] allelesBySortType(GenotypeTable.ALLELE_SORT_TYPE scope, int site);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele.
     *
     * @param site site
     * @param allele allele
     *
     * @return sequence of true/false values.
     */
    public BitSet allelePresenceForAllTaxa(int site, WHICH_ALLELE allele);

    /**
     * Returns BitStorage for this Genotype
     *
     * @param allele allele
     *
     * @return BitStorage
     */
    public BitStorage bitStorage(WHICH_ALLELE allele);

    /**
     * Annotations of this Genotype Table. Null if there are none.
     *
     * @return Annotations or null
     */
    public GeneralAnnotationStorage annotations();

    public Stream<Byte> streamGenotype();

    public Stream<Byte> streamGenotype(int taxon);

    public boolean hasSiteTranslations();

    public int[] siteTranslations();
}
