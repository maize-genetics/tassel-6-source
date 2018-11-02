package net.maizegenetics.dna.map;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;

import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;

/**
 * List of positions in the genome. This type is used by every
 * {@link GenotypeTable}, but it can also be used list of GWAS results and other
 * genomic annotations.
 *
 * @author Terry Casstevens and Ed Buckler
 */
public interface PositionList extends List<Position> {

    /**
     * Return the (haploid) reference allele at given site.
     *
     * @param site site
     *
     * @return byte from 0-15
     */
    public byte allele(WHICH_ALLELE alleleType, int site);

    /**
     * Returns reference sequence alleles in specified range. End site not
     * included. One haploid allele for each site.
     *
     * @param startSite start site
     * @param endSite end site (not included in result)
     *
     * @return reference sequence of haploid allele values.
     */
    public byte[] alleles(WHICH_ALLELE alleleType, int startSite, int endSite);

    /**
     * Returns sequence alleles. One haploid allele for each site.
     *
     * @return reference sequence of haploid allele values.
     */
    public byte[] alleleForAllSites(WHICH_ALLELE alleleType);

    /**
     * Return whether this alignment has defined reference sequence.
     *
     * @return true if this alignment has reference sequence.
     */
    public boolean hasReference();

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String siteName(int site);

    /**
     * Returns total number of sites of this alignment.
     *
     * @return number of sites
     */
    public int numberOfSites();

    /**
     * Return number of sites for given Chromosome
     *
     * @param chromosome
     * @return number of sites
     */
    public int chromosomeSiteCount(Chromosome chromosome);

    /**
     * Get the first (inclusive) and last (inclusive) site of the specified
     * chromosome in this alignment.
     *
     * @param chromosome chromosome
     *
     * @return first and last site
     */
    public int[] startAndEndOfChromosome(Chromosome chromosome);

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
     * Returns PositionList Collector that validates order of Positions.
     *
     * @return collector
     */
    public static Collector<Position, ?, PositionList> collectValidateOrder() {
        return new PositionListCollector(true);
    }

    /**
     * Returns PositionList Collector that reorders position if necessary.
     *
     * @return collector
     */
    public static Collector<Position, ?, PositionList> collectReorder() {
        return new PositionListCollector(false);
    }

    public static class PositionListCollector implements Collector<Position, PositionListBuilder, PositionList> {

        private final boolean myValidateOrder;

        public PositionListCollector(boolean validateOrder) {
            myValidateOrder = validateOrder;
        }

        @Override
        public Supplier<PositionListBuilder> supplier() {
            return PositionListBuilder::new;
        }

        @Override
        public BiConsumer<PositionListBuilder, Position> accumulator() {
            return PositionListBuilder::add;
        }

        @Override
        public BinaryOperator<PositionListBuilder> combiner() {
            return (left, right) -> {
                left.addAll(right);
                return left;
            };
        }

        @Override
        public Function<PositionListBuilder, PositionList> finisher() {
            return (result) -> {
                if (myValidateOrder && !result.validateOrdering()) {
                    throw new IllegalStateException("PositionList: PositionListCollector: Postions are not in order.");
                } else {
                    return result.build();
                }
            };
        }

        @Override
        public Set<Collector.Characteristics> characteristics() {
            if (myValidateOrder) {
                return Collections.EMPTY_SET;
            } else {
                return Collections.unmodifiableSet(EnumSet.of(Collector.Characteristics.UNORDERED));
            }
        }

    }

}
