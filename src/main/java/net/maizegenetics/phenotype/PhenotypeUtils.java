/*
 *  PhenotypeUtils
 *
 *  Created on Oct 27, 2014
 */
package net.maizegenetics.phenotype;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Terry Casstevens
 */
public class PhenotypeUtils {

    private static final Logger myLogger = LogManager.getLogger(PhenotypeUtils.class);

    private static final String DELIMITER = "\t";

    private PhenotypeUtils() {
        // utility
    }

    public static void write(Phenotype phenotype, String filename) {

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {

            writer.write("<Phenotype>\n");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if (i != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(phenotype.attributeType(i).name());
            }
            writer.write("\n");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if (i != 0) {
                    writer.write(DELIMITER);
                }
                writer.write(phenotype.attributeName(i));
            }
            writer.write("\n");

            TableReportUtils.saveDelimitedTableReport(phenotype, DELIMITER, writer, false);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("PhenotypeUtils: write: problem saving file: " + filename);
        }

    }

    public static void writePlink(Phenotype phenotype, String filename) {

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {

            writer.write("FID");
            writer.write(DELIMITER);
            writer.write("IID");

            for (int i = 0; i < phenotype.numberOfAttributes(); i++) {
                if ((phenotype.attributeType(i) == Phenotype.ATTRIBUTE_TYPE.data)
                        || (phenotype.attributeType(i) == Phenotype.ATTRIBUTE_TYPE.covariate)) {
                    writer.write(DELIMITER);
                    writer.write(phenotype.attributeName(i));
                }
            }
            writer.write("\n");

            int numObservations = phenotype.numberOfObservations();
            for (int i = 0; i < numObservations; i++) {
                String taxonName = phenotype.value(i, 0).toString();
                writer.write(taxonName);
                writer.write(DELIMITER);
                writer.write(taxonName);
                writer.write(DELIMITER);
                for (int j = 1; j < phenotype.numberOfAttributes(); j++) {
                    if ((phenotype.attributeType(j) == Phenotype.ATTRIBUTE_TYPE.data)
                            || (phenotype.attributeType(j) == Phenotype.ATTRIBUTE_TYPE.covariate)) {
                        writer.write(DELIMITER);
                        String value = phenotype.value(i, j).toString();
                        if (value.equalsIgnoreCase("NaN")) {
                            writer.write("NA");
                        } else {
                            writer.write(value);
                        }
                    }
                }
                writer.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("PhenotypeUtils: writePlink: problem saving file: " + filename);
        }

    }

    public static Phenotype createPhenotypeObjectFromDB(ArrayList<ArrayList<String>> phenos) {
        //Sort the taxa names and get the unique ones out
        ArrayList<String> taxaNames = (ArrayList<String>) phenos.get(1).stream().distinct().collect(Collectors.toList());
        Collections.sort(taxaNames);
        //Add these to a TaxaList

        //Loop through each row 
        //Associate the row with the correct phenotype
        //use the following structure:
        //HashMap<TaxaName,Multimap<VariableName,Value>>
        HashMap<String, ListMultimap<String, String>> phenotypeMapping = new HashMap<>();//ArrayListMultimap.create();
        //put a new multimap object for each taxa
        for (String taxa : taxaNames) {
            phenotypeMapping.put(taxa, ArrayListMultimap.create());
        }

        for (int i = 0; i < phenos.get(0).size(); i++) {
            phenotypeMapping.get(phenos.get(1).get(i)).put(phenos.get(3).get(i), phenos.get(4).get(i));
        }

        //Get the unique Variable names:
        ArrayList<String> variableNames = (ArrayList<String>) phenos.get(3).stream().distinct().collect(Collectors.toList());
        Collections.sort(variableNames);

        //CorePhenotype(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types, String name)
        ArrayList<Taxon> taxaList = new ArrayList<Taxon>();

        for (String taxaName : taxaNames) {
            taxaList.add(Taxon.instance(taxaName));

        }

        HashMap<String, double[]> attributeMap = new HashMap<String, double[]>();
        //Loop through the variables
        for (int i = 0; i < variableNames.size(); i++) {
            double[] currentAttribute = new double[taxaNames.size()];
            //loop through the taxa
            for (int j = 0; j < taxaNames.size(); j++) {
                currentAttribute[j] = Double.parseDouble(phenotypeMapping.get(taxaNames.get(j)).get(variableNames.get(i)).get(0));
            }

            attributeMap.put(variableNames.get(i), currentAttribute);
        }


        ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(variableNames.size() + 1);
        ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(variableNames.size() + 1);

        attributes.add(new TaxaAttribute(taxaList));
        types.add(ATTRIBUTE_TYPE.taxa);

        for (String variable : variableNames) {
            //NumericAttribute(String name, double[] doubleValues)
            attributes.add(new NumericAttribute(variable, attributeMap.get(variable)));
            types.add(ATTRIBUTE_TYPE.data);
        }

        return new CorePhenotype(attributes, types, "B4R_Phenotype");

    }

    public static Phenotype createPhenotypeObjectFromDB2(ArrayList<ArrayList<String>> phenos) {
        ArrayList<Taxon> taxaList = new ArrayList<Taxon>();

        HashMap<Integer, String> plotNoToTaxaNameMap = new HashMap<Integer, String>();
        ArrayList<Integer> plotNoMapping = new ArrayList<Integer>();
        for (int i = 0; i < phenos.get(0).size(); i++) {
            // System.out.println(phenos.get(1).get(i));
            plotNoToTaxaNameMap.put(Integer.parseInt(phenos.get(2).get(i)), phenos.get(1).get(i));

        }
        System.out.println("Size of map: " + plotNoToTaxaNameMap.size());

        //Sort the taxa names and get the unique ones out
        ArrayList<String> taxaNames = new ArrayList<String>();//(ArrayList<String>)phenos.get(1);

        for (Integer key : plotNoToTaxaNameMap.keySet()) {
            plotNoMapping.add(key);
            taxaNames.add(plotNoToTaxaNameMap.get(key));
            taxaList.add(Taxon.instance(plotNoToTaxaNameMap.get(key)));
        }

        //Get the unique Variable names:
        ArrayList<String> variableNames = (ArrayList<String>) phenos.get(3).stream().distinct().collect(Collectors.toList());
        Collections.sort(variableNames);

        //plotNo<varName,value>
        HashMap<Integer, HashMap<String, Double>> phenotypeMapping = new HashMap<Integer, HashMap<String, Double>>();

        for (int i = 0; i < phenos.get(0).size(); i++) {
            if (!phenotypeMapping.containsKey(Integer.parseInt(phenos.get(2).get(i)))) {
                phenotypeMapping.put(Integer.parseInt(phenos.get(2).get(i)), new HashMap<String, Double>());
            }
            //phenotypeMapping.get(taxaNames.get(i)).put(phenos.get(3).get(i),Double.parseDouble(phenos.get(4).get(i)));
            phenotypeMapping.get(Integer.parseInt(phenos.get(2).get(i))).put(phenos.get(3).get(i), Double.parseDouble(phenos.get(4).get(i)));

        }

        System.out.println("Size of full phenoMap: " + phenotypeMapping.size());


        HashMap<String, double[]> attributeMap = new HashMap<String, double[]>();
        //Loop through the variables
        for (int i = 0; i < variableNames.size(); i++) {
            double[] currentAttribute = new double[taxaNames.size()];

            for (int j = 0; j < taxaList.size(); j++) {
                currentAttribute[j] = phenotypeMapping.get(plotNoMapping.get(j)).get(variableNames.get(i));
            }

            attributeMap.put(variableNames.get(i), currentAttribute);
        }

        ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(variableNames.size() + 1);
        ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(variableNames.size() + 1);

        attributes.add(new TaxaAttribute(taxaList));
        types.add(ATTRIBUTE_TYPE.taxa);

        for (String variable : variableNames) {
            //NumericAttribute(String name, double[] doubleValues)
            attributes.add(new NumericAttribute(variable, attributeMap.get(variable)));
            types.add(ATTRIBUTE_TYPE.data);
        }

        System.out.println("Number of attributes: " + attributes.size());

        return new CorePhenotype(attributes, types, "B4R_Phenotype");
    }


    //Helper method to make a Phenotype from a list of taxaNames, variableNames and the 2-D dataset
    public static Phenotype createPhenotypeFromTransform(ArrayList<String> taxaListNames, ArrayList<String> variableNames, List<ArrayList<Double>> data, String outputObjectName, ArrayList<ATTRIBUTE_TYPE> types) throws Exception {
        //TODO check to make sure the dimensions of the phenotype match
        boolean taxaCorrectSize = true;
        boolean variableCorrectSize = true;
        if (taxaListNames.size() != data.size()) {
            taxaCorrectSize = false;
        }
        if (variableNames.size() != data.get(0).size()) {
            variableCorrectSize = false;
        }
        //Throw the appropriate Exception to the user if the phenotype object does not have the correct dimensions of the data
        if (!taxaCorrectSize && !variableCorrectSize) {
            throw new Exception("Error Size of Taxa names and Variable names do not match:\nNumber of taxa: " + taxaListNames.size() + " Number of samples: " + data.size() + "\nNumber of variables: " + variableNames.size() + " Number of columns: " + data.get(0).size());
        }
        if (!taxaCorrectSize) {
            throw new Exception("Error: Number of Taxa Names and Number of samples does not match: Number of taxa: " + taxaListNames.size() + " Number of samples: " + data.size());
        }
        if (!variableCorrectSize) {
            throw new Exception("Error: Number of Variable Names and Number of columns does not match: Number of variables: " + variableNames.size() + " Number of columns: " + data.get(0).size());
        }

        //Create a new List of Taxon objects
        ArrayList<Taxon> taxaList = (ArrayList<Taxon>) taxaListNames.stream().map((name) -> Taxon.instance(name)).collect(Collectors.toList());

        ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(variableNames.size() + 1);

        attributes.add(new TaxaAttribute(taxaList));

        for (int col = 0; col < data.get(0).size(); col++) {
            double[] currentColumnValues = new double[data.size()];
            for (int row = 0; row < data.size(); row++) {
                currentColumnValues[row] = data.get(row).get(col);
            }

            attributes.add(new NumericAttribute(variableNames.get(col), currentColumnValues));
        }

        System.out.println("Number of Taxa: " + taxaListNames.size() + " Number of attributes: " + attributes.size());

        return new CorePhenotype(attributes, types, outputObjectName);

    }


}
