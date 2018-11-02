/**
 * 
 */
package net.maizegenetics.dna.tag;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.repgen.AlignmentInfo;
import net.maizegenetics.analysis.gbs.repgen.RefTagData;
import net.maizegenetics.analysis.gbs.repgen.TagCorrelationInfo;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;

/**
 * @author lcj34
 *
 */
public interface RepGenData {
    /**
     * Returns the TaxaDistribution for a given taxa.  If distribution data is not available, null is returned.
     */
    TaxaDistribution getTaxaDistribution(Tag tag);

    /**
     * Provides all SNP allele calls associated with a given tag.
     * @param tag used for query
     * @return Map of position with allele state.  The positions are the positions of the actual SNP
     * (not the start of the tag).  The Byte value is defined by {@link net.maizegenetics.dna.snp.NucleotideAlignmentConstants}.
     * If no allele are available for a tag, the map will be empty.
     */
    Set<Allele> getAlleles(Tag tag);


    /**
     * Provides an iterator that goes through all tags with non-empty alleleMaps.
     * @return an iterator of alleleMaps (Map of position with allele state.  The positions are the positions of the actual SNP
     * (not the start of the tag).  The Byte value is defined by {@link net.maizegenetics.dna.snp.NucleotideAlignmentConstants}.)
     */
    Multimap<Tag,Allele> getAlleleMap();

    /**
     * Return all tags that call a particular SNP position allele.  Empty list is returned if there are no tags.
     */
    Set<Tag> getTagsForAllele(Position position, byte allele);

    /**
     * Return all tags that call a particular SNP position allele.  Empty list is returned if there are no tags.
     */
    Set<Tag> getTagsForAllele(Allele allele);

    /**
     * Set of all tags
     */
    Set<Tag> getTags();

    /**
     * Map of all tags with associated names
     */
    Map<Tag,String> getTagsNameMap();

    /**
     * Create the unique list of SNPs (variants) based on all alleles.
     * @return PositionList of all SNPs
     */
    PositionList getSNPPositions();

    /**
     * Create the unique list of SNPs (variants) based on all alleles.
     * Only alleles with support values greater than
     * the minSupportValue are considered.
     * @return PositionList of all SNPs
     */
    PositionList getSNPPositions(int minSupportValue);
    /**
     * Create the unique list of SNPs (variants) based on all alleles.
     * Only positions with minimum quality score equal to or greater than
     * the u are considered.
     * @return PositionList of all SNPs
     */
    PositionList getSNPPositions(double minQualityScore);
    
    /**
     * Get SNPs for specified chromosomes
     * @param starting chromosome number
     * @param ending chromosome number
     * @return  PostionList of SNP positions for the requested chromosomes
     */
    PositionList getSNPPositionsForChromosomes(Integer startChr,Integer endChr);

    /**
     * Returns the list of tags present of a taxon.  Note this could be a very compute intensive request.
     */
    Set<Tag> getTagsForTaxon (Taxon taxon);


    /*
     * Return a map of Tags (keys) and depth (value) for a given taxon and SNP position.
     */
    Map<Tag, Integer> getTagDepth(Taxon taxon, Position position);

    /*
     * Return a map of Tags (keys) and depth (value) greater or equal to a minimum depth
     */
    Map<Tag, Integer> getTagsWithDepth(int minimumDepth);


    /**
     * Get all of the genomic cut positions associated with all tags.  Positions in the position list are annotated with
     * cigarAlignment, isBest, alignmentApproach, and supportValue.
     * @return iterator of annotated positions
     */
    PositionList getTagCutPositions(boolean onlyBest);

    /**
     * Get the unique list of position that tags to for a specific chromosome within a range.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @param onlyBest only return positions flagged as best
     * @return List of positions
     */
    PositionList getTagCutPositions(Chromosome chromosome, int firstPosition, int lastPosition, boolean onlyBest);

    /**
     * Get the cut position associated with a list of tags.  Return a map of Tag/CutPosition
     * The Position object is returned to facilitate grabbing both physical position 
     * and strand information.
     * @param list of tags for which a cut position is requested.  
     * @return map of tag, Tuple<Strand,Cutposition>
     */
     Map<Tag, Position> getTagCutPosition(Set<Tag> tag);
     
    /*
    Map of alignment approaches names (key) and their protocols (value)
     */
    Map<String,String> getTagAlignmentApproaches();

    /**
     * Map of positions and with associated map of Tags and their taxa distribution and their alignment direction.
     * Warning:  This can be a very large data structure for entire chromosomes. Only the best positions are returned.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @return  Map of maps for position (key) to Map of Tag(key) to the Tuple(Direction,TaxaDistribution)(Value)
     */
    Map<Position, Map<Tag, Tuple<Boolean,TaxaDistribution>>> getCutPositionTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition);
    /**
     * For a given snp position,
     * Returns a Map of Allele with its associated Tag/TaxaDistribution 
     * Differs from getCutPositionTagTaxaMap() in that it now specifies on which
     * strand the returned tags must appear.
     * @return Map or null if not available.
     */
    Map<Position, Map<Tag, TaxaDistribution>> getCutPosForStrandTagTaxaMap(Chromosome chromosome, int firstPosition, int lastPosition, boolean strand);


    /**
     * Returns the taxa list associated with taxa distribution
     * @return TaxaList or null if not available.
     */
    TaxaList getTaxaList();
    /**
     * For a given map of chromosome/posiiton,
     * Returns a map of chromsome (string) and Tuple (position,qualityScore)
     * @return Map or null if not available.
     */
    ListMultimap<String, Tuple<Integer, Float>> getSNPPositionQS(HashMultimap<String, Integer> myMap);

    /**
     * For a given snp position,
     * Returns a Map of Allele with its associated Tag/TaxaDistribution 
     * @return Map or null if not available.
     */
    Multimap<Allele, Map<Tag, TaxaDistribution>> getAllelesTagTaxaDistForSNP(
                        Position position);     
        
    /**
     * Return all chromosomes stored in database.
     * Returns a list of distinct chromsome objects from the cutPosition table
     * @return List<Chromosome> or null if none found
     */
     List<Chromosome> getChromosomesFromCutPositions();

     /**
      * Return all tag/taxadistribution stored in database.
      * Returns a Hashmap of tag.taxadistributions
      * @return Map<Tag, TaxaDistribution>  - map is empty if db table is empty
      */
     Map<Tag, TaxaDistribution> getAllTagsTaxaMap();

     /**
      * Get all of the physical map positions associated with the ref tags. 
      * @return iterator of annotated positions
      */
    PositionList getPhysicalMapPositions();

    /**
     * Get the unique list of positions for ref tags for a specific chromosome within a range.
     * @param chromosome chromosome
     * @param firstPosition first physical position in genome (value <0 will return physical >=0)
     * @param lastPosition inclusive last physical position (value <0 will assume Integer.MAX)
     * @return List of positions
     */
    PositionList getPhysicalMapPositions(Chromosome chromosome, int firstPosition, int lastPosition);
    
    /**
     * Set of all reference tags
     */
    Set<RefTagData> getRefTags();
    

    /**
     * Method returns a map of all tag (non-ref) alignments for specified non-ref tags.
     * All tags on the list must be non-ref tags
     * as the tagID is pulled from  the tagTagIDMap.
     *
     * @return
     */
    Multimap<Tag,AlignmentInfo> getTagAlignmentsForTags(List<Tag> tags, int minscore);
    
    /**
     * Map of all reftag alignments for specified non-ref tags.
     * All tags on the list must be non-ref tags
     * as the tagID is pulled from  the tagTagIDMap.
     * 
     * @param tags  list of tags for which we want the ref tag alignment
     * @param minscore  minimum score value for alignment to be included
     * @return
     */
    Multimap<Tag,AlignmentInfo> getRefAlignmentsForTags(List<Tag> tags,int minscore);
    
    /**
     * Map of all non-reftag alignments for specified ref tag.
     * 
     * @param refTag  A ref tag for which we want the non-ref tag alignment
     * @param minscore  minimum score value for alignment to be included
     * @return
     */
    Multimap<Tag,AlignmentInfo> getTagAlignmentsForRefTag(RefTagData refTag,int minscore);
    
    /**
     * Map of all non-ref tags and their alignments
     * This method calls getTagAlignmentsForTags().
     * All tag-tag alignments are return
     * 
     * @param minscore  the minimum score required to pull the alignment
     *                  minscore=0 means get them all
     * @return
     */
    Multimap<Tag,AlignmentInfo> getAllNonRefTagAlignments(int minscore);
    
    /**
     * Map of all alignments for reference tags
     * This returns reftag to reftag alignments.
     * 
     * @param minscore the minimum score required to pull the alignment
     * 
     * @return
     */
    Multimap<RefTagData,AlignmentInfo> getAllRefTagAlignments(int minscore);

    /**
     * Map of all alignments for specified ref tags having
     * minimum specified score.  minscore = 0 means get all tags.
     *
     * This returns refTag to refTag alignments from the reftag_reftag_Alignments table.
     * It does not include refTag to non-refTag alignments. 
     * 
     * @param minscore the minimum score required to pull the alignment
     * @return
     */
    Multimap<RefTagData, AlignmentInfo> getRefAlignmentsForRefTags(List<RefTagData> refTags, int minscore);  
    
    /**
     * Grab correlations for specified tags
     * @param  List<Tags> List of tags whose correlations will be pulled
     */

    Multimap<Tag,TagCorrelationInfo> getCorrelationsForTags(List<Tag> tags);
}
