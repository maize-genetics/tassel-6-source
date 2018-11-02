package net.maizegenetics.dna.map;

import com.google.common.collect.BoundType;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import net.maizegenetics.util.DirectedGraph;
import net.maizegenetics.util.Utils;
import org.apache.commons.lang.StringUtils;
import org.json.simple.JSONObject;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by jgw87 on 7/2/14.
 * A map to hold genome features for lookup by name and by location. The features themselves are hierarchical and so
 * can be traced up and down the tree.
 * <p/>
 * As methods are added to return based on different filters, try to unify the results. That is, all filter-like methods
 * should return the same sort of data structure, such as a HashSet of GenomeFeatures. It may be worthwhile to create
 * a custom GenomeFeatureSet class to be able to string such operations together (mygenes = mygenes.ofType("exon").onChrom(1).inRange(1, 10000);),
 * although whether such filters are needed for the bulk of this class's purpose (matching SNPs to genome annotations) is yet
 * to be seen.
 * <p/>
 * This class shouldn't be created directly, but should instead use the GenomeFeatureMapBuilder class .build() method
 */
//TODO: Add functionality to write out a full constructed map to a file of some sort
public class GenomeFeatureMap {


    //Graph of all the genome features, rooted at the genome itself
    DirectedGraph<GenomeFeature> featureTree = null;

    //Lookups to identify GenomeFeatures by their name, type, and location
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private Multimap<String, GenomeFeature> typeLookup = null;
    private HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup = null;

    /**
     * Default constructor for creating a GenomeFeatureMap from pre-made data structures. This should ONLY be called from
     * the {@link GenomeFeatureMapBuilder} class
     *
     * @param nameLookup     Lookup table of unique IDs -> {@link GenomeFeature} objects
     * @param locationLookup Lookup table to retrieve {@link GenomeFeature}s by their genomic location
     * @param typeLookup     Lookup table to retrieve {@link GenomeFeature}s by their type (exon, UTR, etc)
     * @param featureTree    The graph of genomic features. Should be a directed graph, with the first two levels being genome and chromosome
     */
    GenomeFeatureMap(HashMap<String, GenomeFeature> nameLookup, Multimap<String, GenomeFeature> typeLookup,
                     HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup, DirectedGraph<GenomeFeature> featureTree) {
        this.typeLookup = typeLookup;
        this.nameLookup = nameLookup;
        this.locationLookup = locationLookup;
        this.featureTree = featureTree;
    }

    public GenomeFeature getFeatureFromId(String id) {
        return nameLookup.get(id);
    }

    /**
     * Get a {@link HashSet} of {@link GenomeFeature}s at a specified genome location. Takes chromsome as a String
     * for ones like "Pt", "scaffold487", etc.
     * @param chrom The chromosome name
     * @param start Beginning physical position
     * @param end End physical position
     * @return A {@link HashSet} of GenomeFeatures
     */
    public HashSet<GenomeFeature> getFeaturesInRange(String chrom, int start, int end) {
        Range myrange = Range.closed(start, end); //'Closed' = inclusive, so closed(1,3) = 1,2,3 and closed(1,1) = 1
        Map<Range<Integer>, HashSet<GenomeFeature>> chromMap = locationLookup.get(chrom).subRangeMap(myrange).asMapOfRanges();
        HashSet<GenomeFeature> featureSet = new HashSet<>();
        for(Range r: chromMap.keySet()){
            featureSet.addAll(chromMap.get(r));
        }
        return featureSet;
    }

    /**
     * Get a {@link HashSet} of {@link GenomeFeature}s at a specified genome location
     * @param chrom Chromosome number (should be the same as its name)
     * @param position Physical position (base pair)
     * @return A HashSet of GenomeFeatures
     */
    public HashSet<GenomeFeature> getFeaturesAtLocation(int chrom, int position) {
        return getFeaturesInRange(chrom, position, position);
    }

    /**
     * Get a {@link HashSet} of {@link GenomeFeature}s at a specified genome location Takes chromsome as a String
     * for ones like "Pt", "scaffold487", etc.
     * @param chrom Chromosome name
     * @param position Physical position (base pair)
     * @return A HashSet of GenomeFeatures
     */
    public HashSet<GenomeFeature> getFeaturesAtLocation(String chrom, int position) {
        return getFeaturesInRange(chrom, position, position);
    }

    /**
     * Get a {@link HashSet} of {@link GenomeFeature}s at a specified genome location. Takes chromsome as a String
     * for ones like "Pt", "scaffold487", etc.
     * @param chrom The chromosome number (should be the same as its name)
     * @param start Beginning physical position
     * @param end End physical position
     * @return A HashSet of GenomeFeatures
     */
    public HashSet<GenomeFeature> getFeaturesInRange(int chrom, int start, int end) {
        return getFeaturesInRange("" + chrom, start, end);
    }

    /** Get all {@link GenomeFeature}s of a specified type
     * @param type The type of feature to get
     * @return A {@link HashSet} of GenomeFeatures
     */
    public HashSet<GenomeFeature> getFeaturesOfType(String type){
        HashSet<GenomeFeature> featureSet = new HashSet<>();
        featureSet.addAll(typeLookup.get(type));
        return featureSet;
    }

    /**
     * Write just the location lookup to a tab-delimited file. This is mostly to check that your locations loaded properly,
     * as there is no way to read them back in.
     *
     * @param filename The output file to be written to
     */
    public void writeLocationLookupToFile(String filename) {
        try {
            BufferedWriter writer = Utils.getBufferedWriter(filename);
            writer.append("chrom\tstart\tstop\tfeatures\n");
            String[] chroms = locationLookup.keySet().toArray(new String[0]);
            Arrays.sort(chroms);
            for (String chrom : chroms) {
                Map<Range<Integer>, HashSet<GenomeFeature>> itermap = locationLookup.get(chrom).asMapOfRanges();
                for (Range<Integer> r : itermap.keySet()) {
                    int start = r.lowerEndpoint();
                    int stop = r.upperEndpoint();
                    //Adjust start-stop positions if ranges are open on that end (= up to but not including that number)
                    if (r.lowerBoundType() == BoundType.OPEN) {
                        start++;
                    }
                    if (r.upperBoundType() == BoundType.OPEN) {
                        stop--;
                    }
                    writer.append(chrom + "\t" + start + "\t" + stop + "\t");
                    //List of features
                    HashSet<GenomeFeature> features = itermap.get(r);
                    for (GenomeFeature f : features) {
                        writer.append(f.id() + ";");
                    }
                    writer.append("\n");
                }
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Write the map data as a JSON file (which can be read in by {@link GenomeFeatureMapBuilder}). Core attributes
     * (unique ID, chromosome, start, stop, and parent ID) are output to all features. Any additional attributes are
     * output only for those features that have them. Attributes are output in alphabetical order. Since the attribute
     * name has to be output for every feature, this can waste space if all your features have the same attributes. In
     * that case a tab-delimited flat file or precompiled binary file is probably a better choice.
     *
     * @param filename The output file to be written to
     */
    public void writeMapAsJsonFile(String filename) {
        try {
            BufferedWriter writer = Utils.getBufferedWriter(filename);
            writer.append("[\n");
            String[] sortedNames = nameLookup.keySet().toArray(new String[0]);  //Sorted array so ordering is consistent
            Arrays.sort(sortedNames);

            //Go through with a for loop instead of foreach to know when hit last element
            for (int i=0; i<sortedNames.length; i++){
                GenomeFeature feature = nameLookup.get(sortedNames[i]);
                JSONObject json = new JSONObject();
                json.putAll(feature.annotations());

                //Skip the terminal comma on the last item
                if(i < sortedNames.length - 1){
                    writer.append(json.toJSONString() + ",\n");
                }else{
                    writer.append(json.toJSONString() + "\n");
                }
            }
            writer.append("]");
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     * Write the map data as a flat, tab-delimited file (which can be read in by {@link GenomeFeatureMapBuilder}). The
     * writer first compiles a set of all attribute types across the map and sets these as the columns (in alphabetical
     * order). Attributes that don't apply to a given feature are output as "NA". This can end up being very wasteful
     * if you have some attributes that only apply to a small percentage of features; in that case, a JSON or binary
     * file is probably a better choice.
     *
     * @param filename The output file to be written to
     */
    public void writeMapAsFlatFile(String filename) {
        //Get list of all attributes
        HashSet<String> attributes = new HashSet<>();
        for(GenomeFeature feature: nameLookup.values()){
            attributes.addAll(feature.annotations().keySet());
        }

        //Put column heads in alphabetical order (except for "id"; that goes first and _should_ always be included)
        ArrayList<String> header = new ArrayList<>();
        if(attributes.contains("id")){
            attributes.remove("id");
            header.add("id");
        }
        String[] tempNames = attributes.toArray(new String[0]);
        Arrays.sort(tempNames);
        for(String s: tempNames){
            header.add(s);
        }

        //Write out file
        try {
            BufferedWriter writer = Utils.getBufferedWriter(filename);

            //Header line
            String firstline = StringUtils.join(header, "\t") + "\n";
            writer.append(firstline);

            //Bulk of data
            String[] names = nameLookup.keySet().toArray(new String[0]);
            Arrays.sort(names);
            for(String id: names){
                GenomeFeature feature = nameLookup.get(id);
                ArrayList<String> data = new ArrayList<>();
                for(String column: header){
                    data.add(feature.getAnnotation(column));
                }
                writer.append(StringUtils.join(data, "\t") + "\n");
            }
            writer.close();
        }catch(IOException e){
            e.printStackTrace();
        }
    }

    /**
     * Write the map data as a precompiled binary that can be read back in by {@link GenomeFeatureMapBuilder}). Unlike
     * JSON or flatfile format, this file is not human-readable and is simply a a representation of the Java data
     * structures saved onto a disk. This makes it very compact and fast to read back in. This is the preferred format
     * for long-term storage of a {@link GenomeFeatureMap} since there is (almost) no risk of someone accidentally modifying
     * the file.
     *
     * @param filename The output file to be written to
     */
    //TODO: Test and use this
    //TODO: Implement Serializable interface in order to write all these things out
    public void writeMapAsBinaryFile(String filename) {

    }
}
