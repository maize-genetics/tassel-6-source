package net.maizegenetics.dna.map;

import com.google.common.collect.*;
import net.maizegenetics.util.DirectedGraph;
import net.maizegenetics.util.Graph;
import net.maizegenetics.util.GraphBuilder;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * Created by jgw87 on 7/2/14.
 * A Builder class to create a GenomeFeatureMap to identify genomic features. Can build piecemeal or read in from a file
 * For now, not implementing functionality to make a new builder from an existing map.
 */
//TODO: Add functionality for reading in a precompiled map
public class GenomeFeatureMapBuilder {

    private static final Logger myLogger = Logger.getLogger(GenomeFeatureMapBuilder.class);

    //Graph of all the genome features connecting to each other
    DirectedGraph<GenomeFeature> featureTree = null;

    //Lookups to identify GenomeFeatures by their name, location, and type (exon, gene, etc)
    private HashMap<String, GenomeFeature> nameLookup = new HashMap<>();
    private HashMap<String, RangeMap<Integer, HashSet<GenomeFeature>>> locationLookup = new HashMap<>(); // Complex.
    private Multimap<String, GenomeFeature> typeLookup = HashMultimap.create();

    //Helper variables used to store information to build the feature map
    HashSet<String> chromosomes = new HashSet<>();


    public GenomeFeatureMap build(){
        buildGenomeTree();
        buildLocationLookup();
        return new GenomeFeatureMap(nameLookup, typeLookup, locationLookup, featureTree);
    }

    /**
     * Build the genome feature tree, linking everything to its parent. If something has no parent, it's assigned to the
     * chromosome (and if no chromosome, then the genome, but that's not very useful). It also checks that the genome itself
     * and chromosomes exist as their own {@link GenomeFeature}s, and creates them if needed.
     */
    private void buildGenomeTree(){
        addRootFeaturesIfNeeded();

        //Go through each feature and add parent IDs for those that don't have them (link to chromosome or genome)
        GenomeFeature genome = getGenome();
        ArrayList<GenomeFeature> toReplace = new ArrayList();
        for(GenomeFeature feature: nameLookup.values()){
            if(feature.id() == getGenome().id()){   //Don't assign the genome as its own parent
                continue;
            }
            if(feature.parentId().equals("NA") || !nameLookup.containsKey(feature.parentId())){
                GenomeFeatureBuilder replacement = new GenomeFeatureBuilder(feature);
                if(feature.chromosome().equals("NA")){
                    replacement.parentId(genome.id());
                }else{
                    replacement.parentId(feature.chromosome());
                }
                toReplace.add(replacement.build());
            }
        }
        for(GenomeFeature f:toReplace){
            replaceFeature(f);
        }

        //Build the actual genome graph based on parent information
        GraphBuilder builder = new GraphBuilder(Graph.GraphType.DIRECTED);
        for(GenomeFeature feature: nameLookup.values()){
            if(feature.id() == getGenome().id()){   //Don't assign the genome as its own parent
                continue;
            }
            if(!nameLookup.containsKey(feature.parentId())){
                myLogger.warn("WARNING! Unable to find parent " + feature.parentId() + " for feature " + feature.id());
                continue;
            }
            GenomeFeature parent = nameLookup.get(feature.parentId());
            builder.addEdge(parent, feature);
        }
        featureTree = (DirectedGraph<GenomeFeature>) builder.build();
    }

    /**
     * Add the genome itself and each chromosome as an overarching GenomeFeatures if they haven't been added already
     */
    private void addRootFeaturesIfNeeded(){

        //Check if a feature names "genome" has been added. If so, confirm there's only one. If not, create it.
        GenomeFeature genome;
        if(typeLookup.keySet().contains("genome")){
            genome = getGenome();
        }else{
            genome = new GenomeFeatureBuilder().id("GENOME").type("genome").build();
            addFeature(genome);
        }

        //Go through and check/create chromosomes, assigning their parent to the genome if need be.
        ArrayList<GenomeFeature> toReplace = new ArrayList();
        ArrayList<GenomeFeature> toAdd = new ArrayList();
        for(String c: chromosomes){
            //If chromosome already entered, check that parent is assigned to genome
            if(nameLookup.containsKey(c)){
                GenomeFeature mychrom = nameLookup.get(c);
                //No parent assigned -> Assign to the genome
                if(mychrom.parentId().equals("NA")){
                    GenomeFeature newChrom = new GenomeFeatureBuilder(mychrom).parentId(genome.id()).build();
                    toReplace.add(newChrom);
                //parent assigned, but to different feature -> throw error
                }else if(!mychrom.parentId().equals(genome.id())){
                    throw new UnsupportedOperationException("Error: Chromosome " + mychrom.id() + " assigned with parent "
                            + mychrom.parentId() + " that does not match the recognized genome " + genome.id());
                }
            }else{ //Chromosome not entered yet -> make a new one
                System.out.println("Making chromosome " + c + " from scratch");
                GenomeFeature newChrom = new GenomeFeatureBuilder().id(c).parentId(genome.id()).chromosome(c).type("chromosome").build();
                toAdd.add(newChrom);
            }
        }

        //Actually add & replace the relevant GenomeFeatures; Done here to avoid ConcurrentModificationExceptions
        for(GenomeFeature f: toAdd){
            addFeature(f);
        }
        for(GenomeFeature f: toReplace){
            replaceFeature(f);
        }
    }

    private GenomeFeature getGenome(){
        if(typeLookup.keySet().contains("genome")) {
            Collection<GenomeFeature> genomes = typeLookup.get("genome");
            //If somehow can't pull feature out, throw an error
            if (genomes.size() == 0) {
                throw new UnsupportedOperationException("Error: Cannot find feature with type 'genome' even though supposedly loaded");
            }

            //If more than one feature annotated as "genome", throw an error
            if (genomes.size() > 1) {
                StringBuilder multiples = new StringBuilder();
                for (GenomeFeature g : genomes) {
                    multiples.append("\n\t" + g.id());
                }
                throw new UnsupportedOperationException("Error: Attempt to build a GenomeFeatureMap with more than one feature annotated as type 'genome':" + multiples);
            }
            return(genomes.iterator().next()); //Take first (and by this point, only) feature
        }else{
            myLogger.warn("Warning! Attempting to retrieve the root genome when it hasn't been created yet");
            return null;
        }
    }

    private void buildLocationLookup(){
        //Initialize each chromosome as a new Rangemap
        for(String c: chromosomes){
            RangeMap<Integer, HashSet<GenomeFeature>> newmap = TreeRangeMap.create();
            locationLookup.put(c, newmap);
        }

        for(GenomeFeature feature: nameLookup.values()){
            RangeMap mymap = locationLookup.get(feature.chromosome());
            addFeatureToRangemap(mymap, feature);
        }
    }

    /**
     * Add a GenomeFeature to a RangeMap, stacking it on top of any existing features instead of overwriting them
     * (how RangeMap behaves natively)
     * @param masterMap Rangemap object to be modified
     * @param feature GenomeFeature to be added. Only start and stop position are checked
     */
    public static void addFeatureToRangemap(
            RangeMap<Integer, HashSet<GenomeFeature>>  masterMap, GenomeFeature feature){

        //Check that range is valid (all >=0). Features with invalid ranges (e.g., default chromosomes) are ignored
        if(feature.start() <0 || feature.stop() < 0){
            return;
        }

        //Make ranges closed-open b/c otherwise results in a lot of 1-bp redundant ranges where adjacent features were added
        Range<Integer> featureRange = Range.closedOpen(feature.start(), feature.stop() + 1);

        //First, extract out subrange and save any annotations already there. (Requires making a copy so later modification of
        // masterMap doesn't change things
        ArrayList<Range<Integer>> rangeList = new ArrayList<>();
        ArrayList<HashSet<GenomeFeature>> hashList = new ArrayList<>();
        Map<Range<Integer>, HashSet<GenomeFeature>> subranges = masterMap.subRangeMap(featureRange).asMapOfRanges();
        for(Range<Integer> r: subranges.keySet()){
            rangeList.add(Range.closedOpen(r.lowerEndpoint(), r.upperEndpoint()));
            hashList.add(new HashSet<GenomeFeature>(subranges.get(r)));
        }

        //Next, set entire range equal to the new feature to cover any areas not covered by the existing ranges
        HashSet<GenomeFeature> newset = new HashSet<>();
        newset.add(feature);
        masterMap.put(featureRange, newset);

        //Take the extracted, existing annotations and add back in, adding this feature to the list since they overwrite
        // everything already there
        for(int i=0; i<rangeList.size(); i++){
            HashSet<GenomeFeature> tempset = hashList.get(i);
            Range<Integer> temprange = rangeList.get(i);
            tempset.add(feature); //Add the feature on top of existing ones
            masterMap.put(temprange, tempset);
        }
    }

    /**
     * Adds a GenomeFeature to the map. If the feature's unique ID has already been loaded, it throws an UnsupportedOperationException
     * @param feature The GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureMapBuilder addFeature(GenomeFeature feature){
        String ID=feature.id();
        if(nameLookup.containsKey(ID)){
            throw new UnsupportedOperationException("Error: Attempt to add a GenomeFeature whose unique ID is already loaded: "
            + ID);
        }
        return this.addOrReplaceFeature(feature);   //Done to avoid code duplication
    }

    /**
     * Replaces the GenomeFeature with a specified ID with a new one. If unique ID isn't already in the map, it throws
     * an UnsupportedOperationException.
     * @param feature The GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureMapBuilder replaceFeature(GenomeFeature feature){
        String ID=feature.id();
        if(!nameLookup.containsKey(ID)){
            throw new UnsupportedOperationException("Error: Attempt to replace a GenomeFeature whose unique ID has not been loaded yet: "
                    + ID);
        }
        return this.addOrReplaceFeature(feature);   //Done to avoid code duplication
    }

    /**
     * Adds a GenomeFeature to the map, regardless of whether or not it's already been added. This method throws no
     * warnings if you'll overwrite existing data, so use it with caution.
     * @param feature The GenomeFeature to be added
     * @return This builder
     */
    public GenomeFeatureMapBuilder addOrReplaceFeature(GenomeFeature feature){
        nameLookup.put(feature.id(), feature);
        typeLookup.put(feature.type(), feature);
        if(!feature.chromosome().equals("NA")){
            chromosomes.add(feature.chromosome());
        }
        return this;
    }

    /**
     * Load in data from a GFF (Gene Feature Format) file. Since GFF files have only a loose standard for the 9th column,
     * this involves several ad-hoc heuristics about what things to look for. As such, it is not the preferred way to
     * read in annotations. (That is JSON or tab-delimited format.)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files)
     *
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromGffFile(String filename){
        myLogger.warn("GenomeFeatureMapBuilder - Loading genome annotations from GFF file. Will try to parse annotations " +
                "field as best as possible. (JSON or tab-delimited formats are preferred.)");
        try {
            BufferedReader reader = Utils.getBufferedReader(filename);
            String line=reader.readLine();
            while(line != null ){
                if(line.startsWith("#")){    //Skip comment lines
                    line=reader.readLine();
                    continue;
                }
                GenomeFeature newFeature = new GenomeFeatureBuilder().parseGffLine(line).build();
                addFeature(newFeature);
                line=reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return this;
    }


    /**
     * Load in data from a JSON-formatted file. JSON format is defined at http://www.json.org/, and consists of structured
     * key-value pairs. For genome features, the key is the name of an attribute and the value is (obviously) its value.
     * (For example: "chromosome":1). Note that if you have more than one feature per file (the normal case), all but the
     * last closing brace ('}') should be followed by a comma, and the whole group should be within square braces ('[...]'
     * That is, the first character of the file should be '[' and the last should be ']'). This makes it a properly-formatted
     * JSON array.
     *
     * Common attributes (keys) are listed below. Although only the "id" attribute is required, a feature is pretty
     * useless without some sort of positional information (chromosome, start/stop, etc.).
     *   "id":       Unique identifier for this feature. Repeated identifiers throw an error. (Also accepts "name".) Required.
     *   "chrom":    Which chromosome it occurs on (Also accepts "chr" or "chromosome")
     *   "start":    Start position on the chromosome
     *   "stop":     Stop position on the chromosome (Also accepts "end")
     *   "position": Postion on chromosome (in place of "start" and "stop" for features that are a single nucleotide)
     *   "parent_id":   What other named feature this descends from (eg, Gene -> Transcript -> Exon). If none given, this
     *                   will default to the chromosome (or the genome if chromosome isn't supplied)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files).
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromJsonFile(String filename) {
        JSONParser jparse = new JSONParser();
        BufferedReader reader = Utils.getBufferedReader(filename);
        try {
            JSONArray jarray = (JSONArray) jparse.parse(reader);
            Iterator iter = jarray.iterator();
            while(iter.hasNext()){
                JSONObject json = (JSONObject) iter.next();
                GenomeFeature newFeature = new GenomeFeatureBuilder().parseJsonObject(json).build();
                addFeature(newFeature);
            }
            reader.close();
        } catch (IOException e) {
            myLogger.error("Error loading data from JSON file " + filename);
            e.printStackTrace();
        } catch (ParseException e) {
            myLogger.error("Error parsing information in JSON file " + filename);
            e.printStackTrace();
        }

        return this;
    }

    //Alternate implementation with manually reading in the JSON data; not recommended
    /*public GenomeFeatureMapBuilder addFromJsonFile(String filename) {
        JSONParser jparse = new JSONParser();
        BufferedReader reader = Utils.getBufferedReader(filename);
        try {
            String inline = reader.readLine();
            String tempJson = "";
            while(inline != null){
                tempJson += inline;
                //If has closing brace (mark of end of object), then parse
                while(tempJson.contains("}")){  //While loop in case multiple objects on one line
                    //Subtract out JSON object
                    int start=tempJson.indexOf('{');
                    int end=tempJson.indexOf('}');
                    String json = tempJson.substring(start, end+1);
                    System.out.println("jsonAsString:" + tempJson);
                    System.out.println("\ttake out:" + json);
                    tempJson = tempJson.substring(end + 1, tempJson.length());    //Remove everything
                    System.out.println("\tleft:" + tempJson);

                    //Add data as feature
                    JSONObject featureData = (JSONObject) jparse.parse(json);
                    GenomeFeature newFeature = new GenomeFeatureBuilder().parseJsonObject(featureData).build();
                    addFeature(newFeature);
                }
                inline = reader.readLine();
            }
        } catch (IOException e) {
            myLogger.error("Error loading data from JSON file " + filename);
            e.printStackTrace();
        } catch (ParseException e) {
            myLogger.error("Error parsing information in JSON file " + filename);
            e.printStackTrace();
        }

        return this;
    }*/

    /**
     * Load in data from a flat, tab-delimited text file. The first row should be a header identifying what attribute is
     * in each column, and each subsequent row should correspond to a single feature. Columns that don't apply to a
     * given feature should use "NA" or left empty. Common attributes (columns) are listed below. Although only the "id"
     * attribute is required, a feature is pretty useless without some sort of positional information (chromosome,
     * start/stop, etc.).
     *   "id":       Unique identifier for this feature. Repeated identifiers throw an error (Also accepts "name".) Required.
     *   "chrom":    Which chromosome it occurs on (Also accepts "chr" or "chromosome")
     *   "start":    Start position on the chromosome
     *   "stop":     Stop position on the chromosome (Also accepts "end")
     *   "position": Postion on chromosome (in place of "start" and "stop" for features that are a single nucleotide)
     *   "parent_id":   What other named feature this descends from (eg, Gene -> Transcript -> Exon). If none given, this
     *                   will default to the chromosome (or the genome if chromosome isn't supplied)
     *
     * This method does not build the map, so you can string multiple calls together (if, for example, you have
     * different annotations in different files).
     * @param filename
     */
    public GenomeFeatureMapBuilder addFromFlatFile(String filename) {
        try {
            BufferedReader reader = Utils.getBufferedReader(filename);
            String line=reader.readLine();
            String[] header = null;
            int n=1;    //Keep track of line numbers
            while(line != null ){
                n++;
                if(line.startsWith("#")){    //Skip comment lines
                    line=reader.readLine();
                    continue;
                }
                String[] tokens = line.split("\t");

                if(header == null){ //Save header data
                    header=tokens;
                    line=reader.readLine();
                    continue;
                }

                //Check that number of fields is correct
                if(tokens.length != header.length){
                    myLogger.error("Error: line " + n + " has a different number of fields (" + tokens.length + ") than the header (" + header.length + ")");
                }

                //Load everything into a hapmap
                HashMap<String, String> data = new HashMap<>();
                for(int i=0; i<tokens.length; i++){
                    data.put(header[i], tokens[i]);
                }

                //Add to map
                GenomeFeature newFeature = new GenomeFeatureBuilder().loadAll(data).build();
                addFeature(newFeature);
                line=reader.readLine();
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return this;
    }

}
