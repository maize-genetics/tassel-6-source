/*
 * Genotype
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.util.stream.Stream;

/**
 * Interface for genotype calls for a table of taxa and sites. GenotypeCallTable
 * only contain information on the genotype calls - the call table does not have
 * the TaxaList or PositionList. The calls are diploid genotype calls stored and
 * retrieved as bytes. Conversions to string are also provided.
 * <p>
 * </p>
 * Generally these methods are accessed through the GenotypeTable interface.
 *
 * @see net.maizegenetics.dna.snp.GenotypeTable;
 *
 * @author Terry Casstevens
 */
public interface GenotypeCallTable {

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
     * Same as genotypeAsStringArray(taxon, site), except given value is
     * converted for given site.
     *
     * @param site site
     * @param value diploid allele value
     *
     * @return String representation
     */
    public String[] genotypeAsStringArray(int site, byte value);

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
     * Returns whether this alignment is phased.
     *
     * @return true if phased.
     */
    public boolean isPhased();

    /**
     * Returns true if this Alignment retains rare alleles. If false, rare
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
     * Same as getAlleleEncodings() for only one site.
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
     * lowest for whole alignment. Resulting double dimension array holds
     * diploids (Strings) in result[0]. And the counts are in result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] genoCounts();

    /**
     * Returns counts of all major/minor allele combinations from highest
     * frequency to lowest for whole alignment. Resulting double dimension array
     * holds major/minor allele (Strings) in result[0]. And the counts are in
     * result[1] (Longs).
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
     * Return sorted list of alleles from highest frequency to lowest at given
     * site in alignment. Resulting double dimension array holds alleles (bytes)
     * in result[0]. And the counts are in result[1]. Counts haploid values
     * twice and diploid values once. Higher ploids are not supported.
     *
     * @param site site
     *
     * @return sorted list of alleles and counts
     */
    public int[][] allelesSortedByFrequency(int site);

    /**
     * Return sorted list of diploid vales from highest frequency to lowest at
     * given site in alignment. Resulting double dimension array holds diploids
     * (Strings) in result[0]. And the counts are in result[1] (Integers).
     *
     * @param site site
     *
     * @return sorted list of diploids and counts
     */
    public Object[][] genosSortedByFrequency(int site);

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
     * Returns major allele for all sites.
     *
     * @return all major alleles
     */
    public byte[] majorAlleleForAllSites();

    /**
     * Returns minor allele for all sites.
     *
     * @return all minor alleles
     */
    public byte[] minorAlleleForAllSites();

    /**
     * Return the second minor allele at given site (Third allele state). Gap is
     * included as state. Heterozygous count one for each allele value.
     * Homozygous counts two for the allele value.
     *
     * @param site site
     *
     * @return second minor allele
     */
    public byte thirdAllele(int site);

    /**
     * Returns third (second minor) allele for all sites.
     *
     * @return all third alleles
     */
    public byte[] thirdAlleleForAllSites();

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
     * Returns number of taxa (samples) in this genotype
     *
     * @return number of taxa
     */
    public int numberOfTaxa();

    /**
     * Returns total number of sites in this genotype.
     *
     * @return number of sites
     */
    public int numberOfSites();

    /**
     * Get all genotypes for given taxon.
     *
     * @param taxon taxon
     *
     * @return genotypes
     */
    public byte[] genotypeForAllSites(int taxon);

    /**
     * Get all genotypes for given taxon from start site (inclusive) to end site
     * (exclusive).
     *
     * @param taxon taxon
     * @param start start
     * @param end end
     *
     * @return genotypes
     */
    public byte[] genotypeForSiteRange(int taxon, int start, int end);

    /**
     * Get all genotypes for given site.
     *
     * @param site site
     *
     * @return genotypes
     */
    public byte[] genotypeForAllTaxa(int site);
    
    public Stats siteStats(int site);
    
    public Stats taxonStats(int taxon);

    /**
     * Tells this Genotype to transpose it's data to optimize performance for
     * given iteration nesting. If siteInnerLoop is true, performance better
     * when looping through sites inside taxa loop. If false, performance better
     * when looping through taxa inside site loop.
     *
     * @param siteInnerLoop flag for which iteration
     */
    public void transposeData(boolean siteInnerLoop);

    /**
     * This returns true if this Genotype performs better when processing whole
     * sites at a time. Return false if performance is better when processing
     * whole taxa.
     *
     * @return true if optimized for site processing
     */
    public boolean isSiteOptimized();

    /**
     * Returns a Stream over the genotype calls.
     *
     * @return Stream over the genotype calls
     */
    public Stream<Byte> stream();

    public Stream<Byte> stream(int taxon);
}
