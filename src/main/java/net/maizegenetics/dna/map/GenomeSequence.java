/**
 * Interface for Genome Sequence
 */
package net.maizegenetics.dna.map;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Tuple;


/**
 * Defines the genome sequence of a chromosome
 *
 * @author Lynn Johnson
 *
 */
public interface GenomeSequence {

    /**
     * Returns a list of chromosomes whose sequences have been
     * stored in the chromsomeSequence map of the class implementing
     * this interface.  Return empty set if empty.
     *
     * @return a Set of Chromosome objects
     */
    public Set<Chromosome> chromosomes();

    /**
     * Takes a Chromosome object and returns the stored byte array representing
     * the genomic sequence for the specified chromosome.
     *
     * @param chrom; a Chromosome object representing the chromosome whose
     * 				sequence will be returned
     * @return A byte array containing the chromosome alleles in NucleotideAlignmentConstant
     * 			form packed in half bytes
     */
    public byte[] chromosomeSequence(Chromosome chrom);

    /**
     * Returns the partial genomic sequence for a  chromosome, from the specified start
     * position to the specified end position.  THe start/end positions are inclusive and
     * the request is 1-based (though the alleles are stored in a 0-based byte array).
     *
     * @param chrom:  the chromosome whose partial sequence will be returned.
     * @param startSite:  the 1-based position in the sequence to start the pull.
     * @param endSite:  the 1-based position in the sequence that will be the last allele in the pull
     * @return A byte array of alleles in NucleotideAlignmentConstant form that is packed into
     * 			half bytes.
     */
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int endSite);

    /**
     * Returns the partial genomic sequence from the specified start
     * position to the specified end position.  THe start/end positions are inclusive and
     * the request is 0-based (though the alleles are stored in a 0-based byte array).  Note the difference with
     * chromosomes, which start with 1.  Can only return 2.1 billion sites per call.
     *
     * @param startSite:  the 0-based position in the sequence to start the pull.
     * @param lastSite:  the 0-based position in the sequence that will be the last allele in the pull
     * @return A byte array of alleles in NucleotideAlignmentConstant
     */
    public byte[] genomeSequence(long startSite, long lastSite);

    /**
     * Returns the partial genomic sequence from the specified start
     * position to the specified end position.  The start/end positions are inclusive and
     * the request is 0-based (though the alleles are stored in a 0-based byte array).  Note the difference with
     * chromosomes, which start with 1.  Can only return 2.1 billion sites per call.
     *
     * @param startSite:  the 0-based position in the sequence to start the pull.
     * @param lastSite:  the 0-based position in the sequence that will be the last allele in the pull
     * @return A String of the sequence
     */
    default String genomeSequenceAsString(long startSite, long lastSite) {
        return NucleotideAlignmentConstants.nucleotideBytetoString(genomeSequence(startSite, lastSite));
    }

    /**
     * Takes a list of coordinates from the full genome sequence and for each returns
     * the corresponding chromosome and a coordinate relative to the start of
     * that chromosome.  The request is 0-based, as are the arrays where the
     * alleles are stored, and the results. 
     *
     * @param coordinates:  list of coordinates from the full genome sequence to be like mapped
     * @return A map of <Long, Tuple<Chromosome/coordinate>> for each global coordinate
     *          passed where "Long" is the global ref, and Tuple<> is the chrom/chrom-position
     *          to which this relates.
     */
    public Map<Long, Tuple<Chromosome, Integer>> fullRefCoordinateToChromCoordinate(ArrayList<Long> coordinates);

    /**
     * Returns the length of the current chromsome
     * @param chromosome
     * @return
     */
    public int chromosomeSize(Chromosome chromosome);

    /**
     * Returns the length of the entire genome
     */
    public long genomeSize();

    /**
     * Returns the number of chromosomes
     */
    public int numberOfChromosomes();

    /**
     * Returns the allele value in a byte for the specified PHYSICAL position (1-based)
     * @param chrom  Chromosome object we wish to query
     * @param position Position on the chromosome whose value will be returned
     * @return
     */
    byte genotype(Chromosome chrom, int position);

    /**
     * Returns the TASSEL encoding for the allele value in a byte for the specified PHYSICAL position
     * @param chrom
     * @param positionObject  Position object from which the physical will be obtained.
     * @return
     */
    byte genotype(Chromosome chrom, Position positionObject);


    /**
     * Returns the haplotype allele value in a String for the specified physical position on the specified chromosome
     * @param chrom
     * @param position
     * @return
     */
    String genotypeAsString(Chromosome chrom, int position);

    /**
     * Returns the haplotype allele value in a String for the specified physical position on the specified chromosome
     * @param chrom
     * @param positionObject
     * @return
     */
    String genotypeAsString(Chromosome chrom, Position positionObject);

    /**
     * Returns a string of haplotype allele values for the specified physical start and
     * physical end positions on the specified chromosome.
     * @param chrom
     * @param startSite
     * @param endSite
     * @return
     */
    String genotypeAsString(Chromosome chrom, int startSite, int endSite);


}
