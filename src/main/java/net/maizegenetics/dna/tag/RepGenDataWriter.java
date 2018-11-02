/**
 * 
 */
package net.maizegenetics.dna.tag;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Multimap;

import net.maizegenetics.analysis.gbs.repgen.AlignmentInfo;
import net.maizegenetics.analysis.gbs.repgen.RefTagData;
import net.maizegenetics.analysis.gbs.repgen.TagCorrelationInfo;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.Tuple;

/**
 * @author lcj34
 *
 */
public interface RepGenDataWriter extends RepGenData {
    /**
     * Add a tag to list of known tags
     * @return true if this set did not already contain the specified element
     * @param tags
     * @param tagInstanceAverageQS map showing number of tags, and average quality score for each tag
     */
    boolean putAllTag(Set<Tag> tags,Map<Tag,Tuple<Integer,String>> tagInstanceAverageQS);

    /**
     * Add a tag to list of known tags with their associated names (good for named contigs)
     * @return true if this set did not already contain the specified element
     * @param tagNameMap
     */
    boolean putAllNamesTag(Map<Tag,String> tagNameMap);
    
    /**
     * Add a reference tag to list of known reference tags
     * @return true if this set did not already contain the specified element
     * @param refTagPositionap  map containing refTag sequence and all chrom/positions where it occurs
     * @param referenceGenome  name of reference genome
     */
    boolean putAllRefTag(Multimap<Tag,Position> refTagPositionMap, String refGenome);

    /**
     * Associates a map full of the specified Tag (key) with the specified TaxaDistribution (value).
     * If there was a prior association is it replaced, as all pairs are unique.
     */
    void putTaxaDistribution(Map<Tag, TaxaDistribution> tagTaxaDistributionMap);


    /**
     * Associates the specified Reference Tag with the specified site Position (value). 
     * @param tagAnnotatedPositionMap Map of specific tag with chrom/positions specified.
     *       
     */
    void putRefTagMapping(Multimap<Tag, Position> tagAnnotatedPositionMap, String refGenome);
    
    /**
     * Stores the Smith Waterman score from2 tag alignments. 
     * tag2 chrom/pos comes from the AlignmentInfo object.  tag1 chrom/pos are separate parameters
     * @param tagAlignInfoMap Map of specific tag to tag2 alignment data
     * 
     */
    void putTagTagAlignments(Multimap<Tag,AlignmentInfo> tagAlignInfoMap);
    
    /**
     * Stores the Smith Waterman score from2 tag alignments. 
     * tag2 chrom/pos comes from the AlignmentInfo object.  tag1 chrom/pos are separate parameters
     * @param tagAlignInfoMap Map of specific tag to tag2 alignment data
     * @param refGeome    String used to determine refernce genome if one of the tags is a reference
     */
    void putTagRefTagAlignments(Multimap<Tag,AlignmentInfo> tagAlignInfoMap, String refGenome);
    
    /**
     * Adds entries to the tagAlignments table for ref-ref alignments
     * @param tagAlignInfoMap holds alignment info for each reftag-reftag pair
     * @param refGenome  name of the reference genome
     * @throws SQLException
     */
    void putRefRefAlignments(Multimap<RefTagData, AlignmentInfo> tagAlignInfoMap, String refGenome);
    
    /*
    Set the specified Tag and Position combination to best, and all others were set to false.
     */
    void setTagAlignmentBest(Tag tag, Position position, boolean isBest);

    /*
    Associates a specific Tag with an Allele (a specified SNP position and allele call (plus optional support value)).
    Prior associations at the same Tag-Allele combination are replaced.
     */
    boolean putTagAlleles(Multimap<Tag, Allele> tagAlleleMap);

    /*
    Adds a new Alignment approach name to a detailed protocol.
     */
    boolean putTagAlignmentApproach(String tagAlignmentName, String protocol);

    /*Sets the taxaList for given set of Taxa, this is the order in which the taxa distribution is recorded*/
    void putTaxaList(TaxaList taxaList);
    
    /**
     * Stores a quality position in the snpposition table for each chromosome/position
     * @param qsMap
     */
    void putSNPPositionQS(PositionList qsPosL);
    
    /**
     * Removes all data from the tagtaxadistribution table. 
     * This should be called from GBSSeqToTagDBPlugin when a user requests an "append"
     * option (ie, db exists, and user opted NOT to clear it).  AFter we grab the
     * existing data it is cleared.  It will be re-entered when GBSSeqToTagDBPlugin completes.
     * 
     */
    void clearTagTaxaDistributionData();
    
    /**
     * Removes all data from the DB that was added from SAMToGBSDBPlugin call.
     * The tables cleared are  CutPosition and TagCutPosition 
     */
    void clearAlignmentData();
    
    /**
     * Removes all data from the DB that was added from the DiscoverySNPCallerPluginV2
     * The tables cleared are  Allele, TagAllele and SNPPosition 
     */
    void clearDiscoveryData();
    
    /**
     * Removes all data from the snpQuality table
     */
    void clearSNPQualityData();

    /**
     * Adds a mapping approach strategy to the mappingApproach table
     * @param name
     */
    void addMappingApproach(String name);

    /**
     * Adds a mapping approach strategy to the mappingApproach table
     * @param name
     * @throws SQLException
     */
    void addReferenceGenome(String refGenome);
    
    /**
     * Adds entries to the allelepair table
     * @param name
     * @throws SQLExceptionƒ
     */
    void putAllelePairs(Multimap<Tag,Tuple<Tag,Integer>> tagTagAlignMap);


    
    /**
     * Adds entries to the tagCorrelations table for tag-tag correlation data 
     * @param tagCorrelationMap holds correltaions info for each tag/taxa-depth to tag/taxa-depth vector pair
     * @throws SQLException
     */
    void putTagTagCorrelationMatrix(Multimap<Tag,TagCorrelationInfo> tagCorrelationMap);
}
