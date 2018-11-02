package net.maizegenetics.dna.map;

import java.util.HashMap;

/**
 * Created by Jason Wallace on 7/2/14.
 * This class stores a generic "feature" on a genome, such as a gene, transcript, exon, miRNA, etc. The intent is for
 * this information to be read in from an external file (such as an Ensembl GFF/GTF file) and collated into a
 * GenomeFeatureMap, but it could be used in other ways.
 */
public class GenomeFeature {

    private int start, stop;    //Location on the chromosome (start and stop should be inclusive)
    private HashMap<String, String> annotations;    //Hashmap of all annotations, stored as strings.

    GenomeFeature(HashMap<String, String> myannotations){
       this.annotations=myannotations;

       //Assign position values based on annotations. Lookup is 100-1000x faster this way than having to convert from String each time
       this.start = assignPosition("start");
       this.stop = assignPosition("stop");
    }

    /**
     * Parse the stored annotation on a position into an int. Returns -1 if not found.
     * @param key
     * @return An integer value of the position, or -1 if parsing it fails
     */
    private int assignPosition(String key){
        String value = getAnnotation(key);
        try{
            return Integer.parseInt(value);
        }catch(NumberFormatException e){
            return -1;
        }
    }

    //Various convenience methods to get the most common annotations
    public String id(){
        return getAnnotation("id");
    }

    public String type(){
        return getAnnotation("type");
    }

    public String parentId(){
        return getAnnotation("parent_id");
    }

    public String chromosome(){
        return getAnnotation("chromosome");
    }

    public int start(){
        return this.start;
    }

    public String startAsString(){
        return getAnnotation("start");
    }

    public int stop(){
        return this.stop;
    }

    public String stopAsString(){
        return getAnnotation("stop");
    }

    /**
     * Get any annotation based on its key. If this feature lacks that annotation, it returns 'NA'
     * @param key The name of the annotation to look for
     * @return The value of that annotation, or 'NA' if not found
     */
    public String getAnnotation(String key){
        if(annotations.containsKey(key)){
            return annotations.get(key);
        }else{
            return "NA";
        }
    }

    /**
     * Returns a (shallow) copy of the Hashmap that keeps all annotations for this feature. Since the hashmap just
     * stores Strings, a shallow copy is still safe to modify b/c it won't be reflected in the original.
     * @return A copy of the Hashmap storing this feature's annotations.
     */
    public HashMap<String, String> annotations(){
        HashMap<String, String> annotationsCopy = new HashMap<>(annotations);
        return annotationsCopy;
    }

    @Override
    public String toString() {
        return "GenomeFeature{" +
                "chr=" + chromosome() +
                ", id=" + id() +
                ", parentid=" + parentId() +
                ", start=" + start +
                ", stop=" + stop +
                '}';
    }
}
