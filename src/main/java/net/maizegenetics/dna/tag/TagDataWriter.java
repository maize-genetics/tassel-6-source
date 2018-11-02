package net.maizegenetics.dna.tag;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.Allele;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaTissueDist;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.Tuple;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Defines writer to modify the GBS tag data store
 *
 * @author Ed Buckler
 */
public interface TagDataWriter extends TagData {

    /**
     * Add a tag to list of known tags
     * @return true if this set did not already contain the specified element
     * @param tags
     */
    boolean putAllTag(Set<Tag> tags);

    /**
     * Add a tag to list of known tags with their associated names (good for named contigs)
     * @return true if this set did not already contain the specified element
     * @param tagNameMap
     */
    boolean putAllNamesTag(Map<Tag,String> tagNameMap);

    /**
     * Associates a map full of the specified Tag (key) with the specified TaxaDistribution (value).
     * If there was a prior association is it replaced, as all pairs are unique.
     */
    void putTaxaDistribution(Map<Tag, TaxaDistribution> tagTaxaDistributionMap);


    /**
     * Associates the specified Tag (key) with the specified cut site Position (value).  Multiple associations are allowed, as
     * Tags can map to multiple locations.  Each tag should only have one best annotation.
     * @param tagAnnotatedPositionMap Map of specific tag with Annotated Position of the tag cut site.
     *                                Annotations should be cigarAlignment, isBest, alignmentApproach, forward, and supportValue
     */
    void putTagAlignments(Multimap<Tag, Position> tagAnnotatedPositionMap);

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
     * Adds list of tissues to the tissue table
     * @param masterTissueList
     * @return
     */
    boolean putAllTissue(ArrayList<String> masterTissueList);
    
    /**
     * Set distribution of taxa/tissue for each tag
     * @param taxon
     * @param tags
     * @param counts
     * @return
     */
    boolean putTaxaTissueDistribution(String taxon, String tissue, List<Tag> tags, int[] counts);
}
