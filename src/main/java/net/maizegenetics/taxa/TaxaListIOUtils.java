package net.maizegenetics.taxa;

import com.google.common.base.Splitter;
import com.google.common.collect.*;

import net.maizegenetics.util.Utils;
import net.maizegenetics.util.TableReportUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

import org.apache.log4j.Logger;

/**
 * Utilities for reading and writing IdGroup and PedigreeIdGroups.
 *
 * @author Ed Buckler
 */
public class TaxaListIOUtils {

    private static final Logger myLogger = Logger.getLogger(TaxaListIOUtils.class);

    private static final String DELIMITER = "\t";

    private TaxaListIOUtils() {
    }

    /**
     * Create a Multimap of all the taxa associated with a particular annotation
     * value.
     *
     * @param taxaList input taxa list with annotation associated with
     * @param annotation annotation key used to create the multimap, the values
     * of these keys become the key of the resulting Multimap
     * @return Map of AnnotationValues -> Taxon
     */
    public static Multimap<String, Taxon> getMapOfTaxonByAnnotation(TaxaList taxaList, String annotation) {
        ImmutableMultimap.Builder<String, Taxon> annoMap = new ImmutableMultimap.Builder<String, Taxon>().orderKeysBy(Ordering.natural());
        for (Taxon taxon : taxaList) {
            for (String value : taxon.getAnnotation().getTextAnnotation(annotation)) {
                annoMap.put(value, taxon);
            }
        }
        return annoMap.build();
    }

    /**
     * Create a Map of all the taxa associated with a particular annotation
     * value.  If there would be a duplicate mapping, then an Optional.empty() is returned.
     *
     * @param taxaList input taxa list with annotation associated with
     * @param annotation annotation key used to create the map, the values
     * of these keys become the key of the resulting map
     * @return Map of AnnotationValues -> Taxon
     */
    public static Optional<Map<String, Taxon>> getUniqueMapOfTaxonByAnnotation(TaxaList taxaList, String annotation) {
        Map<String, Taxon> annoMap = new TreeMap<>();
        for (Taxon taxon : taxaList) {
            for (String value : taxon.getAnnotation().getTextAnnotation(annotation)) {
                if(annoMap.containsKey(value)) return Optional.empty();
                annoMap.put(value, taxon);
            }
        }
        return Optional.of(ImmutableSortedMap.copyOf(annoMap));
    }

    /**
     * Returns a subsetted taxa list based on annotation value. For example,
     * return all taxa where {@literal GermType=Inbred}.
     *
     * @param baseTaxaList base annotated taxa list
     * @param annotation annotation name (key)
     * @param annoValue annotation value being tested for
     * @return TaxaList equal to the annotation value
     */
    public static TaxaList subsetTaxaListByAnnotation(TaxaList baseTaxaList, String annotation, String annoValue) {
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (Taxon taxon : baseTaxaList) {
            for (String value : taxon.getAnnotation().getTextAnnotation(annotation)) {
                if (value.equals(annoValue)) {
                    tlb.add(taxon);
                    break;
                }
            }
        }
        return tlb.build();
    }

    /**
     * Creates a new taxa list with the taxa only retaining annotations within a
     * specified list. All taxa are retained, only the annotations are changed.
     *
     * @param baseTaxaList
     * @param annotationsToKeep the retained keys annotation
     * @return new TaxaList with a subset of the annotations
     */
    public static TaxaList retainSpecificAnnotations(TaxaList baseTaxaList, String[] annotationsToKeep) {
        Set<String> keepers = new ImmutableSet.Builder<String>().addAll(Arrays.asList(annotationsToKeep)).build();
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (Taxon taxon : baseTaxaList) {
            Taxon.Builder tb = new Taxon.Builder(taxon.getName());
            for (Map.Entry<String, String> entry : taxon.getAnnotation().getAllAnnotationEntries()) {
                if (keepers.contains(entry.getKey())) {
                    tb.addAnno(entry.getKey(), entry.getValue());
                }
            }
            tlb.add(tb.build());
        }
        return tlb.build();
    }

    /**
     * Creates a new taxa list with the taxa retaining annotations EXCEPT those
     * specified by the list. All taxa are retained, only the annotations are
     * changed.
     *
     * @param baseTaxaList
     * @param annotationsToRemove the retained keys annotation
     * @return new TaxaList with a subset of the annotations
     */
    public static TaxaList removeSpecificAnnotations(TaxaList baseTaxaList, String[] annotationsToRemove) {
        Set<String> keepers = new ImmutableSet.Builder<String>().addAll(Arrays.asList(annotationsToRemove)).build();
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (Taxon taxon : baseTaxaList) {
            Taxon.Builder tb = new Taxon.Builder(taxon.getName());
            for (Map.Entry<String, String> entry : taxon.getAnnotation().getAllAnnotationEntries()) {
                if (!keepers.contains(entry.getKey())) {
                    tb.addAnno(entry.getKey(), entry.getValue());
                }
            }
            tlb.add(tb.build());
        }
        return tlb.build();
    }

    /**
     * Provides the set of all annotation key found in any of taxa
     *
     * @param baseTaxaList
     * @return
     */
    public static Set<String> allAnnotationKeys(TaxaList baseTaxaList) {
        ImmutableSet.Builder<String> keepers = new ImmutableSet.Builder<String>();
        for (Taxon taxon : baseTaxaList) {
            for (Map.Entry<String, String> entry : taxon.getAnnotation().getAllAnnotationEntries()) {
                keepers.add(entry.getKey());
            }
        }
        return keepers.build();
    }

    public static void exportAnnotatedTaxaListTable(TaxaList taxa, String filename) {
        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            writer.append("<TaxaList>\n");
            TableReportUtils.saveDelimitedTableReport(new TaxaListTableReport(taxa), DELIMITER, writer, true);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TaxaListIOUtils: exportAnnotatedTaxaListTable: problem saving file: " + filename);
        }
    }

    public static TaxaList importAnnotatedTaxaList(String filename) {
        TaxaListBuilder builder = new TaxaListBuilder();
        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            String header = reader.readLine().trim();
            if (!header.equalsIgnoreCase("<TaxaList>")) {
                throw new IllegalArgumentException("TaxaListIOUtils: importAnnotatedTaxaList: This file doesn't start with <TaxaList>: " + filename);
            }
            header = reader.readLine().trim();
            String[] columns = header.split(DELIMITER);
            for (int i = 0; i < columns.length; i++) {
                columns[i] = columns[i].trim();
            }
            if (!columns[0].equalsIgnoreCase("Taxa")) {
                throw new IllegalArgumentException("TaxaListIOUtils: importAnnotatedTaxaList: First column should be Taxa: " + filename);
            }
            int numColumns = columns.length;
            int lineNum = 2;
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                lineNum++;
                String[] annotations = line.split(DELIMITER);
                if (numColumns != annotations.length) {
                    throw new IllegalStateException("TaxaListIOUtils: importAnnotatedTaxaList: number of annotations doesn't match number of columns line: " + lineNum + " taxon: " + annotations[0].trim());
                }
                Taxon.Builder currentTaxon = new Taxon.Builder(annotations[0].trim());
                for (int i = 1; i < numColumns; i++) {
                    String value = annotations[i].trim();
                    if (!value.isEmpty()) {
                        currentTaxon.addAnno(columns[i], value);
                    }
                }
                builder.add(currentTaxon.build());
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("TaxaListIOUtils: importAnnotatedTaxaList: Problem reading file: " + filename + "\n" + e.getMessage());
        }

        return builder.build();
    }

    /**
     * Returns an annotated TaxaList from a text annotation file in matrix
     * format. This is a tab delimited file. First row in the file with the
     * field {@literal taxaNameField} is the header row.
     * {@literal taxaNameField} indicated the taxon name, all other fields are
     * user defined. The fields become the keys for the taxa annotation.
     * Quantitative fields should be tagged with "#" sign, e.g.
     * {@literal <#INBREEDF>}. Multiple values are supported per key, and
     * additional values can be either described with an additional column or
     * ";" to delimit values with the same key.
     * <p>
     * </p>
     * Filters are a map of filters to be applied. Key are the fields, and value
     * are what are tested for equality. Only taxa rows true for filters are
     * retained.
     * <p>
     * </p> {@literal <Name>	<GermType>	<Set>	<#InbreedF>	<Set>}<br></br>
     * {@literal B73	Inbred	Goodman282	0.98    ISU;IBMFounder}<br></br>
     * {@literal MO17	Inbred	Goodman282	0.98    UMC;IBMFounder}<br></br>
     * <p>
     * </p>
     * Produces:<br></br>
     * {@literal B73<GermType=Inbred,Set=Goodman282,Set=ISU,Set=IBMFounder,f=0.98>}<br></br>
     * {@literal MO17<GermType=Inbred,Set=Goodman282,Set=UMC,Set=IBMFounder,f=0.98>}<br></br>
     * The standardized keys are described in the
     * {@link net.maizegenetics.taxa.Taxon}, and these constant fields are all
     * upper case.
     *
     * @param fileName with complete path
     * @param taxaNameField field name with the taxon name
     * @param filters Map of filter to determine which rows to retain as the
     * file is processed.
     * @return TaxaList with annotations
     */
    
    // New version calls readTaxaAnnotationFileAL, then creates TaxaList from ArrayList
    public static TaxaList readTaxaAnnotationFile(String fileName, String taxaNameField, Map<String, String> filters, boolean mergeSameNames) {    
    	// create list
    	ArrayList<Taxon> taxaAL = readTaxaAnnotationFileAL( fileName,  taxaNameField, filters);
    	if (taxaAL == null) return null;
        TaxaListBuilder tlb = new TaxaListBuilder();
    	taxaAL.stream().forEach(taxa -> {
    		if (mergeSameNames) {
    			tlb.addOrMerge(taxa);
    		} else {
    			tlb.add(taxa);
    		}
    	});
    	return tlb.sortTaxaAlphabetically().build();   	
    }
    
    //  Version of readTaxaAnnotationFile that returns an ArrayList.  This is called
    // from places that wish to allow duplicate taxa.  The TaxaList does not allow for
    // duplicate entries.  An ArrayList does.
    public static ArrayList<Taxon> readTaxaAnnotationFileAL(String fileName, String taxaNameField, Map<String, String> filters) {
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1 << 16);
            String line = fileIn.readLine();
            ArrayList<Taxon> taxaAL = new ArrayList<Taxon>();
            int indexOfName = 0;
            //parse headers
            List<String> headers = new ArrayList<>();
            List<Boolean> isQuant = new ArrayList<>();
            if (line.contains(taxaNameField)) {
                int i = 0;
                for (String header : line.split("\\t")) {
                    if (header.equals(taxaNameField)) {
                        indexOfName = i;
                    }
                    isQuant.add(header.startsWith("#") || header.startsWith("<#"));
                    headers.add(header.replace(">", "").replace("<", "").replace("#", ""));
                    i++;
                }
            } else {
                fileIn.reset();
            }
            //parse taxa rows
            while ((line = fileIn.readLine()) != null) {
                String[] s = line.split("\\t");
                Taxon.Builder anID = new Taxon.Builder(s[indexOfName]);
                for (int i = 0; i < s.length; i++) {
                    if (i == indexOfName) {
                        continue;
                    }
                    String[] cs = s[i].split(";");
                    for (String ta : cs) {
                        if (ta == null || ta.isEmpty()) {
                            continue;
                        }
                        if (isQuant.get(i)) {
                            if (ta.equals("NA")) {
                                anID.addAnno(headers.get(i), Double.NaN);
                            } else {
                                anID.addAnno(headers.get(i), Double.parseDouble(ta));
                            }
                        } else {
                            anID.addAnno(headers.get(i), ta);
                        }
                    }
                }
                Taxon t = anID.build();
                if (doesTaxonHaveAllAnnotations(t, filters)) {           	
                	taxaAL.add(t);
                }
            }
            // Sort alphabetically based on name.  This is to remain consistent
            // with the readTaxaAnnotationFile(), which alphabetizes the taxaList
            Collections.sort(taxaAL, new Comparator<Taxon>() {
                public int compare(Taxon taxa1, Taxon taxa2) {
                    return taxa1.getName().compareTo(taxa2.getName());
                }
            });
            return taxaAL;
        } catch (Exception e) {
            System.err.println("Error in Reading Annotated Taxon File:" + fileName);
            e.printStackTrace();
        }
        return null;
    }

    /**
     * Returns an annotated TaxaList from a text annotation file in matrix
     * format. This is a tab delimited file. First row in the file with the
     * field {@literal taxaNameField} is the header row.
     * {@literal taxaNameField} indicated the taxon name, all other fields are
     * user defined. The fields become the keys for the taxa annotation.
     * Quantitative fields should be tagged with "#" sign, e.g.
     * {@literal <#INBREEDF>}. Multiple values are supported per key, and
     * additional values can be either described with an additional column or
     * ";" to delimit values with the same key.
     * <p>
     * </p> {@literal <Name>	<GermType>	<Set>	<#InbreedF>	<Set>}<br></br>
     * {@literal B73	Inbred	Goodman282	0.98    ISU;IBMFounder}<br></br>
     * {@literal MO17	Inbred	Goodman282	0.98    UMC;IBMFounder}<br></br>
     * <p>
     * </p>
     * Produces:<br></br>
     * {@literal B73<GermType=Inbred,Set=Goodman282,Set=ISU,Set=IBMFounder,f=0.98>}<br></br>
     * {@literal MO17<GermType=Inbred,Set=Goodman282,Set=UMC,Set=IBMFounder,f=0.98>}<br></br>
     * The standardized keys are described in the
     * {@link net.maizegenetics.taxa.Taxon}, and these constant fields are all
     * upper case.
     *
     * @param fileName with complete path
     * @param taxaNameField field name with the taxon name
     * @return TaxaList with annotations
     */
    public static TaxaList readTaxaAnnotationFile(String fileName, String taxaNameField) {
        return readTaxaAnnotationFile(fileName, taxaNameField, new HashMap<String, String>(), false);
    }

    /**
     * Tests whether a taxon has annotation values in the map
     *
     * @param taxon
     * @param filters
     * @return true if all present, false is otherwise
     */
    public static boolean doesTaxonHaveAllAnnotations(Taxon taxon, Map<String, String> filters) {
        SetMultimap<String, String> taxonAnno = taxon.getAnnotation().getAnnotationAsMap();
        boolean keep = true;
        for (Map.Entry<String, String> entry : filters.entrySet()) {
            keep = false;
            for (String s1 : taxonAnno.get(entry.getKey())) {
                if (s1.equals(entry.getValue())) {
                    keep = true;
                } 
            }
            if (keep == false) {
                break;
            }
        }
        return keep;
    }

    /**
     * Parses a VCF header with the taxa names and annotations into a multimap.
     * The taxa name is return as the "ID" key, as used by the VCF format.
     *
     * @param s
     * @return
     */
    public static SetMultimap<String, String> parseVCFHeadersIntoMap(String s) {
        if (s == null) {
            return null;
        }
        if (!(s.startsWith("<") && s.endsWith(">"))) {
            return null;
        }
        String value = s.substring(1, s.length() - 1);
        ImmutableSetMultimap.Builder<String, String> im = new ImmutableSetMultimap.Builder<String, String>()
                .orderKeysBy(Ordering.natural()).orderValuesBy(Ordering.natural());
        for (String s1 : Splitter.on(",").trimResults().split(value)) {
            String[] ssEntry = s1.split("=", 2);
            im.put(ssEntry[0], ssEntry[1]);
        }
        return im.build();
    }
    
    /**
     * This method takes a key file and creates a SortedSet<String> that
     * contains a set of the tissue values.  The set will be null if
     * no tissues are present
     * @param fileName - name of Keyfile containing Tissue header
     * @param tissueNameField - field name
     * @return
     */
    public static List<String> readTissueAnnotationFile(String fileName, String tissueNameField) {
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1 << 16);
            String line = fileIn.readLine();
            List<String> tissues = new ArrayList<String>();
            int indexOfTissue = -1;
            //parse headers
            if (line.contains(tissueNameField)) {
                int idx = 0;
                for (String header : line.split("\\t")) {
                    if (header.equals(tissueNameField)) {
                        indexOfTissue = idx;
                        break;
                    }
                    idx++;
                }
                if (indexOfTissue == -1) {
                    // Tissue header not found - return null
                    return null;
                }
            } else {
                fileIn.reset();
            }
            if (indexOfTissue == -1) return null;
            // Found tissue header, read values - no duplicates, into set
            while ((line = fileIn.readLine()) != null) {
                String[] items = line.split("\\t");
                for (int idx=0; idx < items.length; idx++) {
                    if (idx == indexOfTissue) {
                        if (!(tissues.contains(items[idx]))) {
                            tissues.add(items[idx]);
                        }
                        continue;
                    }
                }  
            }
            Collections.sort(tissues);
                
            return tissues;
        } catch (Exception e) {
            System.err.println("Error in Reading Annotated Tissue File:" + fileName);
            e.printStackTrace();
        }
        return null;
    }
}
