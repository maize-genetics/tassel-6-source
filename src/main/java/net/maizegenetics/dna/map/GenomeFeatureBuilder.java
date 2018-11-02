package net.maizegenetics.dna.map;

import org.apache.log4j.Logger;

import org.json.simple.JSONObject;

import java.util.HashMap;
import java.util.Locale;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by jgw87 on 7/2/14.
 * Builder class to create a GenomeFeature. All annoations are stored in a HashMap. Any annotation can be added through
 * the addAnnotation() method, but the more common fields have their own convenience methods. Only the feature's own ID
 * is required; all other annotations are optional
 */
public class GenomeFeatureBuilder {

    private static final Logger myLogger = Logger.getLogger(GenomeFeatureBuilder.class);

    //Variables to store the information on the feature
    private HashMap<String, String> myannotations = null;

    /**
     * Generic constructor which does nothing special
     */
    public GenomeFeatureBuilder() {
        myannotations = new HashMap<String, String>();
    }

    /**
     * Constructor to build a new feature off of an existing one.
     *
     * @param feature The genome feature to copy
     */
    public GenomeFeatureBuilder(GenomeFeature feature) {
        this.myannotations = feature.annotations();    //the annotations() method returns a shallow copy, so should be safe
    }

    /**
     * Public accessor method to get a new GenomeFeatureBuilder based off an existing GenomeFeature
     * TODO: Get this matching other getInstance methods better; doesn't seem to fit, so commented out for now
     * //@param feature
     * //@return
     */
    /*public static GenomeFeatureBuilder getInstance(GenomeFeature feature){
        return new GenomeFeatureBuilder(feature);
    }*/
    public GenomeFeature build() {
        validateData();
        return new GenomeFeature(myannotations);
    }

    private void validateData() {

        //Test that feature has the required fields
        if ( (!myannotations.containsKey("id")) || myannotations.get("id") == "NA") {
            throw new UnsupportedOperationException("GenomeFeatureBuilder: Cannot build a feature without a personal identifier (field 'id')");
        }


        //Test if start or stop is negative
        if(myannotations.containsKey("start") && myannotations.containsKey("stop")) {
            int mystart = Integer.parseInt(myannotations.get("start"));
            int mystop = Integer.parseInt(myannotations.get("stop"));
            if (mystart < 0) {
                throw new UnsupportedOperationException("GenomeFeatureBuilder: Start coordinate is negative for " +
                        myannotations.get("id") + ": " + mystart + " (possibly unassigned?)");
            }
            if (mystop < 0) {
                throw new UnsupportedOperationException("GenomeFeatureBuilder: Stop coordinate is negative for " +
                        myannotations.get("id") + ": " + mystop + " (possibly unassigned?)");
            }

            //Test that start is less than stop
            if (mystart > mystop) {
                throw new UnsupportedOperationException("GenomeFeatureBuilder: Start coordinate is greater than stop " +
                        "coordinate for " + myannotations.get("id") + ": " + mystart + " vs " + mystop);
            }
        }
    }

    public GenomeFeatureBuilder id(String id) {
        return addAnnotation("id", id);
    }

    public GenomeFeatureBuilder type(String type) {
        return addAnnotation("type", type);
    }

    public GenomeFeatureBuilder parentId(String parentId) {
        return addAnnotation("parent_id", parentId);
    }

    public GenomeFeatureBuilder chromosome(Chromosome chr) {
        return addAnnotation("chromosome", chr.getName());
    }

    public GenomeFeatureBuilder chromosome(String chr) {
        return addAnnotation("chromosome", chr);
    }

    public GenomeFeatureBuilder chromosome(int chr) {
        return addAnnotation("chromosome", "" + chr);
    }

    public GenomeFeatureBuilder start(int start) {
        return addAnnotation("start", "" + start);
    }

    public GenomeFeatureBuilder start(String start) {
        return addAnnotation("start", start);
    }

    public GenomeFeatureBuilder stop(int stop) {
        return addAnnotation("stop", "" + stop);
    }

    public GenomeFeatureBuilder stop(String stop) {
        return addAnnotation("stop", stop);
    }

    public GenomeFeatureBuilder position(String position) {
        return addAnnotation("position", position);
    }

    public GenomeFeatureBuilder position(int position) {
        return addAnnotation("position", "" + position);
    }

    public GenomeFeatureBuilder addAnnotation(String key, String value) {
        key = synonymizeKeys(key);
        if("".equals(value)){   //Convert empty strings to 'NA'
            value="NA";
        }
        myannotations.put(key, value);  //All annotations kept in the hash

        //If the key is "position", convert to identical start-stop coordinates as well.
        if (key == "position") {
            this.start(value);
            this.stop(value);
        }
        return this;
    }

    /**
     * Method that takes common synonyms of annotation types and standardizes them according to the following rules:
     * (1) Make lowercase
     * (2) Standardize according to following rules. (Any not on this list are returned as just lowercased)
     * name, id -> id
     * chr, chrom, chromosome -> chromosome
     * stop, end -> stop
     * parentid, parent_id, parent -> parent_id
     * pos, position -> position
     *
     * @param key The key to standardize
     * @return
     */
    public static String synonymizeKeys(String key) {
        key = key.toLowerCase(Locale.ENGLISH);
        switch (key) {
            case "name":
            case "id":
                return "id";

            case "chr":
            case "chrom":
            case "chromosome":
                return "chromosome";

            case "end":
            case "stop":
                return "stop";

            case "parent":
            case "parentid":
            case "parent_id":
                return "parent_id";

            case "pos":
            case "position":
                return "position";

            default:
                return key;
        }
    }

    /**
     * Load all annotations from a hashmap. Keys become the annotations, and values the annotation value. Each key-value
     * pair is added individually (instead of using a putAll() method) to allow for key standardization, etc.
     * @return This builder, with the new data loaded from the hashmap
     */
    public GenomeFeatureBuilder loadAll(HashMap<String, String> newAnnotations){
        for(String key: newAnnotations.keySet()){
            addAnnotation(key, newAnnotations.get(key));
        }
        return this;
    }

    /**
     * Create a GenomeFeature from a line of GFF file. This method is modified from the BioJava source code for the same
     * purpose, in biojava3-genome/src/main/java/org/biojava3/genome/GFF3Reader.java
     *
     * @param line A single line from a GFF file as a string
     * @return This builder, with the data loaded from the line
     */
    public GenomeFeatureBuilder parseGffLine(String line) {
        //Field identifiers for GFF format
        int gffSeq = 0, gffSource = 1, gffFeatureType = 2, gffStart = 3, gffStop = 4, gffScore = 5, gffStrand = 6, gffFrame = 7, gffAttributes = 8;

        //Get all the easy data stored in its own fields
        String[] tokens = line.split("\t");
        this.chromosome(tokens[gffSeq].trim());
        this.type(tokens[gffFeatureType].trim());
        this.start(tokens[gffStart]);
        this.stop(tokens[gffStop]);
        addAnnotation("strand", tokens[gffStrand].trim());

        //Extract the parent from the attributes field
        String parentID = getParentFromGffAttributes(tokens[gffAttributes]);
        this.parentId(parentID);

        //Extract the unique identifier for this feature. If none, build one from available info
        String myID = getFeatureIdFromGffAttributes(tokens[gffAttributes]);
        if (myID == null) {
            myID = myannotations.get("type") + "_" + myannotations.get("chromosome") + "_" + myannotations.get("start") + "_" + myannotations.get("stop");
        }
        this.id(myID);

        return this;
    }

    /**
     * Parse a GFF attribute field to identify the parent of the current GenomeFeature. Tricky b/c of the different ways it
     * can be represented. There's a hierarchy of accepted answers, with 'parent_id', 'Parent=', 'transcript_id', and
     * 'gene_id' taken in that order. If nothing is found, returns an empty string ("")
     *
     * @param attributes The string from the attribute field of the GFF file
     * @return The parent's ID string
     */
    public static String getParentFromGffAttributes(String attributes) {
        //Match the appropriate string with regular expressions
        Matcher matcher;

        //Pattern for 'parent_id "GRMZM2G005232"'
        matcher = Pattern.compile("parent_id \"(\\w+)\"").matcher(attributes);
        if (matcher.find()) {
            return matcher.group(1);
        }

        //Pattern for 'Parent=GRMZM2G005232' and  'Parent=gene:GRMZM2G005232'
        matcher = Pattern.compile("Parent=(\\w+:){0,1}(\\w+)").matcher(attributes);
        if (matcher.find()) {
            return matcher.group(2);
        }

        //Pattern for 'gene_id "GRMZM2G005232"'
        matcher = Pattern.compile("transcript_id \"(\\w+)\"").matcher(attributes);
        if (matcher.find()) {
            return matcher.group(1);
        }

        //Pattern for 'transcript_id "GRMZM2G005232_T01"'
        matcher = Pattern.compile("gene_id \"(\\w+)\"").matcher(attributes);
        if (matcher.find()) {
            return matcher.group(1);
        }

        return "";
    }

    /**
     * Parse a GFF attribute field to identify the name of the current GenomeFeature. Looks for 'ID=' and 'Name=' fields
     *
     * @param attributes The string from the attribute field of the GFF file. REturns null if not found
     * @return The feature's ID string
     */
    public static String getFeatureIdFromGffAttributes(String attributes) {
        //Match the appropriate string with regular expressions
        Matcher matcher;

        //Pattern for 'ID=GRMZM2G005232' and  'ID=gene:GRMZM2G005232'
        matcher = Pattern.compile("(Name|ID)=(\\w+:){0,1}(\\w+)").matcher(attributes);
        if (matcher.find()) {
            return matcher.group(3);
        }

        return null;
    }

    public GenomeFeatureBuilder parseJsonObject (JSONObject featureData){
        HashMap<String, String> jsonHash = new HashMap<>();
        for (String key : (Set<String>) featureData.keySet()) {
            addAnnotation(key, featureData.get(key).toString());
        }
        return this;
    }

}
