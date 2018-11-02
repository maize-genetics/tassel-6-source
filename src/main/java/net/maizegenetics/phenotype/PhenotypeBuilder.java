package net.maizegenetics.phenotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.log4j.Logger;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 * @author Peter Bradbury
 *
 */
public class PhenotypeBuilder {
	private Logger myLogger = Logger.getLogger(PhenotypeBuilder.class);
	
	private enum SOURCE_TYPE{file, phenotype, list, attributes};
	private enum ACTION{importFile, union, intersect, separate, concatenate, keepTaxa, removeTaxa, keepAttributes, changeType, buildFromAttributes, removeMissing};
	
	private SOURCE_TYPE source;
	private List<ACTION> actionList = new ArrayList<PhenotypeBuilder.ACTION>();
	private List<String> filenameList = null;
	private List<Phenotype> phenotypeList = new ArrayList<Phenotype>();
	private String phenotypeName = "Phenotype";
	private List<PhenotypeAttribute> attributeList = null;
	private List<ATTRIBUTE_TYPE> attributeTypeList = null;
	private List<PhenotypeAttribute> separateByList = null;
	
	//filter criteria
	private List<Taxon> taxaToKeep = null;
	private List<Taxon> taxaToRemove = null;
  	private int[] indexOfAttributesToKeep = null;
	private Map<PhenotypeAttribute, ATTRIBUTE_TYPE> attributeChangeMap = null;
	
	private boolean addSourceDataFactor = false;
	
	public PhenotypeBuilder() {
		
	}
	
	/**
	 * @param filename	the name of a file containing phenotype data to be imported
	 * @return	a PhenotypeBuilder that will import a file
	 */
	public PhenotypeBuilder fromFile(String filename) {
		if (filenameList == null) filenameList = new ArrayList<String>();
		filenameList.add(filename);
		source = SOURCE_TYPE.file;
		actionList.add(ACTION.importFile);
		return this;
	}
	
	/**
	 * @param base	the Phenotype from which to create a new Phenotype
	 * @return	a PhenotypeBuilder
	 * Multiple phenotypes can be added to a list using this method.
	 */
	public PhenotypeBuilder fromPhenotype(Phenotype base) {
		phenotypeList.add(base);
		source = SOURCE_TYPE.phenotype;
		return this;
	}
	
	/**
	 * @param phenotypes	a List of Phenotypes
	 * @return	a Phenotype builder that will act on all Phenotypes in the list
	 * A union join returns a Phenotype containing any taxon present in at least one of the Phenotypes to be joined.
	 * An intersect join returns a Phenotype containing only taxa present in all of the Phenotypes to be joined.
	 */
	public PhenotypeBuilder fromPhenotypeList(List<Phenotype> phenotypes) {
		phenotypeList.addAll(phenotypes);
		source = SOURCE_TYPE.list;
		StringBuilder sb = new StringBuilder();
		for (Phenotype pheno : phenotypes) {
			if (sb.length() > 0) sb.append(" + ");
			sb.append(pheno.name());
		}
		phenotypeName = sb.toString();
		return this;
	}
	
	/**
	 * @param attributes	a list of attributes
	 * @param types	a list of types matching the attribute list
	 * @return	a PhenotypeBuilder that will build using these lists
	 * The attribute and type lists must be the same size and the types must be compatible with the attributes
	 */
	public PhenotypeBuilder fromAttributeList(List<PhenotypeAttribute> attributes, List<ATTRIBUTE_TYPE> types) {
		attributeList = attributes;
		attributeTypeList = types;
		source = SOURCE_TYPE.attributes;
		actionList.add(ACTION.buildFromAttributes);
		return this;
	}
	
	/**
	 * @return	a PhenotypeBuilder that will perform a union join of Phenotypes added using one or more of the from methods.
	 *  A union join returns a Phenotype containing Taxa that are in any of the source Phenotypes.
	 */
	public PhenotypeBuilder unionJoin() {
		actionList.add(ACTION.union);
		return this;
	}
	
	/**
	 * @return	a PhenotypeBuilder that will perform an intersect join of Phenotypes added using one or more of the from methods.
	 *  An intersect join returns a Phenotype containing only Taxa that are in all of the source Phenotypes.
	 */
	public PhenotypeBuilder intersectJoin() {
		actionList.add(ACTION.intersect);
		return this;
	}
	
	/**
	 * @return	a builder that will concatenate phenotypes 
	 */
	public PhenotypeBuilder concatenate() {
		actionList.add(ACTION.concatenate);
		return this;
	}
	
	public PhenotypeBuilder assignName(String name) {
		phenotypeName = name;
		return this;
	}
	
	/**
	 * @param taxaToKeep	a list of taxa to be kept from the base Phenotype
	 * @return	a PhenotypeBuilder that will return a FilterPhenotype with taxa in the taxaToKeep list
	 * Only taxa that are in both taxaToKeep and the base Phenotype will be included in the Phenotype that is built.
	 */
	public PhenotypeBuilder keepTaxa(List<Taxon> taxaToKeep) {
		actionList.add(ACTION.keepTaxa);
		this.taxaToKeep = taxaToKeep;
		return this;
	}
	
	/**
	 * @param taxaToRemove	a list of taxa to removed from the base Phenotype
	 * @return	a PhenotypeBuilder that will return a Phenotype with taxa from the supplied list excluded
	 * Any taxa in taxaToRemove but not in the base Phenotype will be ignored.
	 */
	public PhenotypeBuilder removeTaxa(List<Taxon> taxaToRemove)  {
		actionList.add(ACTION.removeTaxa);
		this.taxaToRemove = taxaToRemove;
		return this;
	}
	
	/**
	 * @param attributesToKeep	a list of the attributes to be kept
	 * @return	a PhenotypeBuilder that will return a new Phenotype with only the attributes in the supplied list
	 * Only attributes in both attributesToKeep and the base Phenotype will be included in the Phenotype that is built.
	 * The order of the attributes in the new Phenotype will match that in attributesToKeep.
	 */
	public PhenotypeBuilder keepAttributes(List<PhenotypeAttribute> attributesToKeep) {
		actionList.add(ACTION.keepAttributes);
		attributeList = attributesToKeep;
		return this;
	}
	
	/**
	 * @param indexOfAttributes	the column numbers of the attributes in the base Phenotype to be included in the newly built Phenotype
	 * @return	a PhenotypeBuilder that will build a Phenotype with the specified attributes
	 */
	public PhenotypeBuilder keepAttributes(int[] indexOfAttributes) {
		actionList.add(ACTION.keepAttributes);
		this.indexOfAttributesToKeep = indexOfAttributes;
		return this;
	}
	
	/**
	 * @param attributeIndex	the numeric index (column number) of an attribute in the base Phenotype
	 * @param type	the new type for that attribute
	 * @return	a PhenotypeBuilder that will build a phenotype with the changed attribute type
	 */
	public PhenotypeBuilder changeAttributeType(Map<PhenotypeAttribute, ATTRIBUTE_TYPE> changeMap) {
		if (attributeChangeMap == null) attributeChangeMap = new HashMap<PhenotypeAttribute, Phenotype.ATTRIBUTE_TYPE>();
		attributeChangeMap.putAll(changeMap);
		actionList.add(ACTION.changeType);
		return this;
	}
	
	/**
	 * @param attributeTypes	a list of attribute types for the attributes to be built
	 * @return	a PhenotypeBuilder that will build a Phenotype that will have this list of types
	 * The order of types must be the same as the order of attributes as supplied by the keepAttributes methods if used or in the base Phenotype if the attribute list is not changed.
	 */
	public PhenotypeBuilder typesOfRetainedAttributes(List<ATTRIBUTE_TYPE> attributeTypes) {
		attributeTypeList = attributeTypes;
		return this;
	}
	
	/**
	 * @param factors	a list of attributes on which to pivot the table. These must be CategoricalAttributes.
	 * @return	a PhenotypeBuilder that will create a pivoted Phenotype.
	 *  Pivoting is equivalent to using the separateOn() method followed by a union join.
	 */
	public PhenotypeBuilder pivotOn(List<PhenotypeAttribute> factors) {
		for (PhenotypeAttribute attr : factors) {
			if (!attr.isTypeCompatible(ATTRIBUTE_TYPE.factor)) {
				String msg = "One of the attributes in factors is not a CategoricalAttribute. Only CategoricalAttributes (factors) can be used for pivots.";
				throw new IllegalArgumentException(msg);
			}
		}
		separateByList = factors;
		actionList.add(ACTION.separate);
		actionList.add(ACTION.union);
		return this;
	}
	
	/**
	 * @param factors	a list of attributes used to create separate Phenotypes. These must be CategoricalAttributes.
	 * @return	a PhenotypeBuilder that will create a separate Phenotype for each distinct level of the combined factors.
	 */
	public PhenotypeBuilder separateOn(List<PhenotypeAttribute> factors) {
		for (PhenotypeAttribute attr : factors) {
			if (!attr.isTypeCompatible(ATTRIBUTE_TYPE.factor)) {
				String msg = "One of the attributes in factors is not a CategoricalAttribute. Only CategoricalAttributes (factors) can be used to separate Phenotype observations.";
				throw new IllegalArgumentException(msg);
			}
		}
		separateByList = factors;
		actionList.add(ACTION.separate);
		return this;
	}
	
	/**
	 * @return	a PhenotypeBuilder that will create a Phenotype with no missing values for any attribute by removing observations that have any missing values.
	 */
	public PhenotypeBuilder removeMissingObservations() {
		actionList.add(ACTION.removeMissing);
		return this;
	}
	
	/**
	 * @return	a list of new Phenotypes built with the supplied parameters
	 */
	public List<Phenotype> build() {
		for (ACTION action : actionList) {
			performAction(action);
		}
		return phenotypeList;
	}
	
	private void performAction(ACTION action) {
		switch (action) {
		case importFile:
			importFile();
			break;
		case intersect:
			joinPhenotypes(action);
			break;
		case union:
			joinPhenotypes(action);
			break;
		case concatenate:
			concatenatePhenotypes();
			break;
		case keepTaxa:
			keepTaxaFilter();
			break;
		case keepAttributes:
			keepAttributesFilter();
			break;
		case removeTaxa:
			removeTaxaFilter();
			break;
		case changeType:
			applyAttributeChangeMap();
			break;
		case separate:
			separateByFactors();
			break;
		case buildFromAttributes:
			createPhenotypeFromLists();
			break;
		case removeMissing:
			removeAllObservationsWithMissingValues();
			break;
		}
	}
	
	//private methods  ------------------------------------------------------
	
	private void importFile() {
		while (filenameList.size() > 0) {
			File phenotypeFile = new File(filenameList.remove(0));

			try (BufferedReader br = new BufferedReader(new FileReader(phenotypeFile))) {
				String topline = br.readLine();
				Phenotype myPhenotype;
				if (phenotypeName.equals("Phenotype")) {
					phenotypeName = phenotypeFile.getName();
					if (phenotypeName.endsWith(".txt")) phenotypeName = phenotypeName.substring(0, phenotypeName.length() - 4);
				}

				if (topline.toLowerCase().startsWith("<phenotype")) {
					myPhenotype = importPhenotypeFile(phenotypeFile);
				} else {
					myPhenotype = importTraits(phenotypeFile);
				}
				br.close();
				phenotypeList.add(myPhenotype);
			} catch (IOException e) {
				String errorMsg = "Error reading " + phenotypeFile.getPath() + " in PhenotypeBuilder.importFile()." ;
				myLogger.error(errorMsg);
				e.printStackTrace();
				throw new IllegalArgumentException(errorMsg);
			}
			
		}

	}
	
	private Phenotype importPhenotypeFile(File phenotypeFile) throws IOException {
		Pattern whiteSpace = Pattern.compile("\\s+");
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>();
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>();

		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));
		phenotypeReader.readLine();  //assumes the first line has been read to determine that this is indeed a Phenotype file
		String[] typeString = whiteSpace.split(phenotypeReader.readLine());
		String[] phenoNames = whiteSpace.split(phenotypeReader.readLine());
		
		//find the taxa column
		int taxaCol = 0;
		while (!typeString[taxaCol].toLowerCase().equals("taxa")) taxaCol++;
		
		int nPheno = typeString.length;
		ArrayList<String[]> stringData = new ArrayList<String[]>();
		String inputStr;
		while ((inputStr = phenotypeReader.readLine()) != null) {
			stringData.add(whiteSpace.split(inputStr));
		}

		int nObs = stringData.size();
		for (int pheno = 0; pheno < nPheno; pheno++) {
			if (typeString[pheno].toLowerCase().startsWith("cov") || typeString[pheno].toLowerCase().startsWith("dat")) {
				float[] dataArray = new float[nObs];
				OpenBitSet missing = new OpenBitSet(nObs);
				int obsCount = 0;
				for (String[] inputLine : stringData) {
                                        if (inputLine[pheno].equalsIgnoreCase("NaN")
                                                || inputLine[pheno].equalsIgnoreCase("NA")
                                                || inputLine[pheno].equals(".")) {
                                            dataArray[obsCount] = Float.NaN;
                                            missing.fastSet(obsCount);
                                        } else {
						try {
							dataArray[obsCount] = Float.parseFloat(inputLine[pheno]);
						} catch (NumberFormatException nfe) {
							throw new IllegalArgumentException("PhenotypeBuilder: importPhenotypeFile: at observation " + obsCount + " Taxon: " + inputLine[taxaCol] + " Value: " + inputLine[pheno] + " is not allowed.");
						}
					}
					obsCount++;
				}
				attributes.add(new NumericAttribute(phenoNames[pheno], dataArray, missing));
				if (typeString[pheno].toLowerCase().startsWith("cov")) types.add(ATTRIBUTE_TYPE.covariate);
				else types.add(ATTRIBUTE_TYPE.data);
			} else if (typeString[pheno].toLowerCase().startsWith("tax")) {
				ArrayList<Taxon> taxa = new ArrayList<>();
				for (String[] inputLine : stringData) {
					taxa.add(new Taxon(inputLine[pheno]));
				}
				attributes.add(new TaxaAttribute(taxa, phenoNames[pheno]));
				types.add(ATTRIBUTE_TYPE.taxa);
			} else if (typeString[pheno].toLowerCase().startsWith("fac")) {
				String[] labelArray = new String[nObs];
				int obsCount = 0;
				for (String[] inputLine : stringData) {
					labelArray[obsCount++] = inputLine[pheno];
				}
				attributes.add(new CategoricalAttribute(phenoNames[pheno], labelArray));
				types.add(ATTRIBUTE_TYPE.factor);
			}
		}
		phenotypeReader.close();

		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private Phenotype importTraits(File phenotypeFile) throws IOException {
		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));

		int numberOfDataLines = 0;
		String inputline = phenotypeReader.readLine();
		ArrayList<String> headerLines = new ArrayList<>();
		boolean isFactor = false;
		boolean isCovariate = false;
		boolean hasHeaders = false;
		boolean isTrait = false;
		String[] traitnames = new String[0];
		while (inputline != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				numberOfDataLines++;
			} else if (inputline.toLowerCase().startsWith("<trai")) {
				isTrait = true;
				String[] splitLine = inputline.split("[<>\\s]+");
				traitnames = Arrays.copyOfRange(splitLine, 2, splitLine.length);
			} else if (inputline.toLowerCase().startsWith("<cov")) {
				isCovariate = true;
			} else if (inputline.toLowerCase().startsWith("<fac")) {
				isFactor = true;
			} else if (inputline.toLowerCase().startsWith("<head")) {
				hasHeaders = true;
				headerLines.add(inputline);
			}
			
			inputline = phenotypeReader.readLine();
		}
		phenotypeReader.close();
		
		if (hasHeaders) {
			return processTraitsAndFactors(phenotypeFile, traitnames, numberOfDataLines, isCovariate, headerLines);
		} else if (isFactor) {
			return processFactors(phenotypeFile, traitnames, numberOfDataLines);
		} else if (isTrait) {
			return processTraits(phenotypeFile, traitnames, numberOfDataLines, isCovariate);
		} else throw new IllegalArgumentException("Unrecognized format for a phenotype.");
		
	}
	
	private Phenotype processTraits(File phenotypeFile, String[] traitnames, int numberOfDataLines, boolean isCovariate) throws IOException {
		int ntraits = traitnames.length;
		int nattributes = ntraits + 1;
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(nattributes);
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(nattributes);
		ArrayList<float[]> traitValues = new ArrayList<>(ntraits);
		ArrayList<BitSet> missingList = new ArrayList<>(ntraits);
		ArrayList<Taxon> taxaList = new ArrayList<>();
		
		types.add(ATTRIBUTE_TYPE.taxa);
		ATTRIBUTE_TYPE myAttributeType;
		if (isCovariate) myAttributeType = ATTRIBUTE_TYPE.covariate;
		else myAttributeType = ATTRIBUTE_TYPE.data;
		for (int i = 0; i < ntraits; i++) {
			traitValues.add(new float[numberOfDataLines]);
			missingList.add(new OpenBitSet(numberOfDataLines));
			types.add(myAttributeType);
		}
		
		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));
		int dataCount = 0;
		int lineCount = 1;
		String inputline;
		while ((inputline = phenotypeReader.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				String[] values = inputline.split("\\s+");
				if (values.length != ntraits + 1) {
					String msg = String.format("Incorrect number of values in line %d of %s", lineCount, phenotypeFile.getName());
					phenotypeReader.close();
					throw new IllegalArgumentException(msg);
				}
				taxaList.add(new Taxon(values[0]));
				for (int i = 0; i < ntraits; i++) {
					float val;
					String inval = values[i + 1].trim();
					if (inval.equals("-99") || inval.equals("-999")) {
						val = Float.NaN;
					}
					else {
						try {
							val = Float.parseFloat(inval);
						} catch (NumberFormatException e) {
							val = Float.NaN;
						}
					}
					if (Double.isNaN(val)) missingList.get(i).fastSet(dataCount);
					traitValues.get(i)[dataCount] = val;
				}
				dataCount++;
			}
			
			lineCount++;
		}
		
		phenotypeReader.close();
		
		attributes.add(new TaxaAttribute(taxaList));
		for (int i = 0; i < ntraits; i++) attributes.add(new NumericAttribute(traitnames[i], traitValues.get(i), missingList.get(i)));
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
	private Phenotype processFactors(File phenotypeFile, String[] traitnames, int numberOfDataLines) throws IOException {
		int ntraits = traitnames.length;
		int nattributes = ntraits + 1;
		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(nattributes);
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(nattributes);
		ArrayList<String[]> traitValues = new ArrayList<>(ntraits);
		ArrayList<Taxon> taxaList = new ArrayList<>();
		
		types.add(ATTRIBUTE_TYPE.taxa);
		ATTRIBUTE_TYPE myAttributeType = ATTRIBUTE_TYPE.factor;
		for (int i = 0; i < ntraits; i++) {
			traitValues.add(new String[numberOfDataLines]);
			types.add(myAttributeType);
		}
		
		BufferedReader phenotypeReader = new BufferedReader(new FileReader(phenotypeFile));
		int dataCount = 0;
		int lineCount = 1;
		String inputline;
		while ((inputline = phenotypeReader.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				String[] values = inputline.split("\\s+");
				if (values.length != ntraits + 1) {
					String msg = String.format("Incorrect number of values in line %d of %s", lineCount, phenotypeFile.getName());
					phenotypeReader.close();
					throw new IllegalArgumentException(msg);
				}
				taxaList.add(new Taxon(values[0]));
				for (int i = 0; i < ntraits; i++) {
					traitValues.get(i)[dataCount] = values[i + 1];
				}
				dataCount++;
			}
			
			lineCount++;
		}
		
		phenotypeReader.close();
		
		attributes.add(new TaxaAttribute(taxaList));
		for (int i = 0; i < ntraits; i++) attributes.add(new CategoricalAttribute(traitnames[i], traitValues.get(i)));
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	
//	private void processHeader(String headerLine, Map<String, ArrayList<String>> headerMap) {
//		String[] header = headerLine.split("[<>=\\s]+");
//		if ( header[1].toLowerCase().equals("header") && header[2].toLowerCase().equals("name")) {
//			String name = header[3];
//			ArrayList<String> values = new ArrayList<>();
//			for (int i = 4; i < header.length; i++) {
//				values.add(header[i]);
//			}
//			headerMap.put(name, values);
//		} else throw new IllegalArgumentException("Improperly formatted Header: " + headerLine);
//	}
	
	private Phenotype processTraitsAndFactors(File phenotypeFile, String[] traitnames, int numberOfDataLines, boolean isCovariate, ArrayList<String> headerList)
	throws IOException {
		TreeSet<String> traitSet = new TreeSet<>();
		for (String trait : traitnames) traitSet.add(trait);
		HashMap<String, Integer> traitMap = new HashMap<>();
		int traitCount = 0;
		for (String trait : traitSet) traitMap.put(trait, traitCount++);
		
		int ntraitnames = traitnames.length;
		int ntraits = traitSet.size();
		int nfactors = headerList.size();
		String[] factorNames = new String[nfactors];
		
		//create set of composite headers and get factor names
		String[] splitHeader = headerList.get(0).split("[<>=\\s]+");
		factorNames[0] = splitHeader[3];
		String[] factorValues = Arrays.copyOfRange(splitHeader, 4, splitHeader.length);
		
		for (int i = 1; i < nfactors; i++) {
			splitHeader = headerList.get(i).split("[<>=\\s]+");
			factorNames[i] = splitHeader[3].replace("|", ":");
			for (int j = 0; j < factorValues.length; j++) {
				factorValues[j] += "|" + splitHeader[j + 4].replace("|", ":");
			}
		}
		
		TreeSet<String> factorSet = new TreeSet<>();
		for (String val:factorValues) factorSet.add(val);
		int nCompositeFactors = factorSet.size();
		HashMap<String,Integer> factorMap = new HashMap<>();
		int factorCount = 0;
		for (String factor : factorSet) factorMap.put(factor, factorCount++);
		
		//the length of each attribute array will be numberOfDataLines * number of composite headers
		int ndata = numberOfDataLines * nCompositeFactors;

		//create the factor arrays
		ArrayList<String[]> factorAttributeArrays = new ArrayList<>(nfactors);
		for (int i = 0; i < nfactors; i++) factorAttributeArrays.add(new String[ndata]);
		int fromIndex = 0;
		int toIndex = numberOfDataLines;
		for (String factor : factorSet) {
			String[] subFactor = factor.split("\\|");
			for (int i = 0; i < nfactors; i++) {
				Arrays.fill(factorAttributeArrays.get(i), fromIndex, toIndex, subFactor[i]);
			}
			fromIndex += numberOfDataLines;
			toIndex += numberOfDataLines;
		}
		
		//create the trait arrays and temporary taxa list
		ArrayList<float[]> traitAttributeArrays = new ArrayList<>();
		ArrayList<BitSet> missingList = new ArrayList<>();
		ArrayList<Taxon> tempTaxa = new ArrayList<>();

		for (int i = 0; i < ntraits; i++) {  //initialize the trait arrays to NaN
			float[] traitArray = new float[ndata];
			Arrays.fill(traitArray, Float.NaN);
			traitAttributeArrays.add(traitArray);
			missingList.add(new OpenBitSet(ndata));
		}
		
		BufferedReader br = new BufferedReader(new FileReader(phenotypeFile));
		int dataCount = 0;
		int lineCount = 1;
		String inputline;
		while ((inputline = br.readLine()) != null) {
			inputline = inputline.trim();
			if (inputline.length() > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				String[] values = inputline.split("\\s+");
				if (values.length != ntraitnames + 1) {
					String msg = String.format("Incorrect number of values in line %d of %s", lineCount, phenotypeFile.getName());
					br.close();
					throw new IllegalArgumentException(msg);
				}
				tempTaxa.add(new Taxon(values[0]));
				for (int i = 0; i < ntraitnames; i++) {
					float val;
					int traitnum = traitMap.get(traitnames[i]);
					int factornum = factorMap.get(factorValues[i]);
					int dataIndex = factornum * numberOfDataLines + dataCount;
					try {
						val = Float.parseFloat(values[i + 1]);
					} catch (NumberFormatException e) {
						val = Float.NaN;
						missingList.get(traitnum).fastSet(dataIndex);
					}
					traitAttributeArrays.get(traitnum)[dataIndex] = val;
				}
				dataCount++;
			}
			
			lineCount++;
		}

		br.close();
		
		//create the taxa list
		ArrayList<Taxon> taxaList = new ArrayList<>(ndata);
		for (int i = 0; i < nCompositeFactors; i++) taxaList.addAll(tempTaxa);

		//create the taxa, factor, and trait attributes and types in that order
		ATTRIBUTE_TYPE myAttributeType;
		if (isCovariate) myAttributeType = ATTRIBUTE_TYPE.covariate;
		else myAttributeType = ATTRIBUTE_TYPE.data;

		ArrayList<PhenotypeAttribute> attributes = new ArrayList<>(nfactors + ntraits);
		ArrayList<ATTRIBUTE_TYPE> types = new ArrayList<>(nfactors + ntraits);
		
		attributes.add(new TaxaAttribute(taxaList));
		types.add(ATTRIBUTE_TYPE.taxa);
		
		for (int i = 0; i < nfactors; i++) {
			attributes.add(new CategoricalAttribute(factorNames[i], factorAttributeArrays.get(i)));
			types.add(ATTRIBUTE_TYPE.factor);
		}
		
		traitCount = 0;
		for (String trait : traitSet) {
			attributes.add(new NumericAttribute(trait, traitAttributeArrays.get(traitCount), missingList.get(traitCount)));
			types.add(myAttributeType);
			traitCount++;
		}
		
		return new CorePhenotype(attributes, types, phenotypeName);
	}
	

	private void keepTaxaFilter() {
		List<Phenotype> newList = new ArrayList<Phenotype>();
		for (Phenotype pheno : phenotypeList) {
			newList.add(FilterPhenotype.getInstance(pheno, taxaToKeep, filterPhenotypeName(pheno.name())));
		}
		phenotypeList = newList;
	}
	
	private void keepAttributesFilter() {
		List<Phenotype> newList = new ArrayList<Phenotype>();
		for (Phenotype pheno : phenotypeList) {
			if (attributeList != null) {
				if (attributeTypeList == null) {
					attributeTypeList = new ArrayList<ATTRIBUTE_TYPE>();
					for (PhenotypeAttribute attr:attributeList) {
						attributeTypeList.add(pheno.attributeType(pheno.indexOfAttribute(attr)));
					}
				}
				applyAttributeChangeMap();
				newList.add( new CorePhenotype(attributeList, attributeTypeList, filterPhenotypeName(pheno.name())) );
				
			} else if (indexOfAttributesToKeep != null) {
				if (indexOfAttributesToKeep.length == 1 && indexOfAttributesToKeep[0] == -1) {
					int nattr = pheno.numberOfAttributes();
					indexOfAttributesToKeep = new int[nattr - 1];
					for (int i = 0; i < indexOfAttributesToKeep.length; i++) {
						indexOfAttributesToKeep[i] = i;
					}
				}
				attributeList = new ArrayList<PhenotypeAttribute>();
				attributeTypeList = new ArrayList<ATTRIBUTE_TYPE>();
				for (int attrnum : indexOfAttributesToKeep) {
					attributeList.add(pheno.attribute(attrnum));
				}
				if (attributeTypeList == null || attributeTypeList.size() != attributeList.size()) {
					for (int attrnum : indexOfAttributesToKeep) {
						attributeTypeList.add(pheno.attributeType(attrnum));
					}
				}
				applyAttributeChangeMap();
				newList.add( new CorePhenotype(attributeList, attributeTypeList, filterPhenotypeName(pheno.name())) );
			}
		}
		phenotypeList = newList;

	}
	
	private void removeTaxaFilter() {
		List<Phenotype> newList = new ArrayList<Phenotype>();
		for (Phenotype pheno : phenotypeList) {
			TaxaList myTaxaList = pheno.taxa();
			Iterator<Taxon> taxaIter = myTaxaList.iterator();
			while (taxaIter.hasNext()) {
				if (taxaToRemove.contains(taxaIter.next())) taxaIter.remove();
			}
			newList.add(FilterPhenotype.getInstance(pheno, myTaxaList,  "filtered_" + pheno.name()));
		}
		phenotypeList = newList;
	}
	
	private String filterPhenotypeName(String oldName) {
		if (oldName.toLowerCase().startsWith("filter")) return oldName;
		return "filtered_" + oldName;
	}
	
	private void separateByFactors() {
		List<Phenotype> newList = new ArrayList<Phenotype>();
		for (Phenotype pheno : phenotypeList) {
			newList.addAll(separatePhenotypeByFactors(pheno));
		}
		phenotypeList = newList;
	}
	
	private List<Phenotype> separatePhenotypeByFactors(Phenotype base) {
		List<Phenotype> phenotypeList = new ArrayList<Phenotype>();
		phenotypeList.add(base);
		for (PhenotypeAttribute factor:separateByList) {
			List<Phenotype> subsettedPhenotypeList = new ArrayList<Phenotype>();
			for (Phenotype pheno : phenotypeList) subsettedPhenotypeList.addAll(subsetPhenotypeByOneFactor(pheno, (CategoricalAttribute) factor));
			phenotypeList = subsettedPhenotypeList;
		}
		return phenotypeList;
	}
	
	private List<Phenotype> subsetPhenotypeByOneFactor(Phenotype base, CategoricalAttribute byFactor) {
		List<Phenotype> phenoList = new ArrayList<Phenotype>();
		int nlevels = byFactor.numberOfLevels();
		for (int i = 0; i < nlevels; i++) {
			int[] subset = byFactor.whichObservations(i);
			
			//create attribute and type lists, then create a new Phenotype
			List<PhenotypeAttribute> newAttr = new ArrayList<PhenotypeAttribute>();
			List<ATTRIBUTE_TYPE> newTypes = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
			String factorLevelName = byFactor.name() + "." + byFactor.labelList().get(i);
			String newPhenoName = base.name() + "_" + factorLevelName;
			for (int j = 0; j < base.numberOfAttributes(); j++) {
				if (!base.attribute(j).equals(byFactor)) {
					String newName = null;
					if (base.attributeType(j) == ATTRIBUTE_TYPE.data) newName = base.attribute(j).name() + "_" + factorLevelName;
					newAttr.add(base.attribute(j).subset(subset, newName));
					newTypes.add(base.attributeType(j));
				} 
			}
			phenoList.add(new CorePhenotype(newAttr, newTypes, newPhenoName));
		}
		return phenoList;
	}
	
	private void applyAttributeChangeMap() {
		if (attributeChangeMap != null) {
			List<Phenotype> newList = new ArrayList<Phenotype>();
			for (Phenotype pheno : phenotypeList) {
				List<PhenotypeAttribute> attrList = pheno.attributeListCopy();
				List<ATTRIBUTE_TYPE> typeList = pheno.typeListCopy();
				int nattr = attrList.size();
				for (int i = 0; i < nattr; i++) {
					ATTRIBUTE_TYPE newType = attributeChangeMap.get(attrList.get(i));
					if (newType != null) typeList.set(i, newType);
				}
				newList.add(new CorePhenotype(attrList, typeList, "retyped_" + pheno.name()));
			}
			phenotypeList = newList;
		}
	}
	
	private void createPhenotypeFromLists() {
		if (attributeList.size() != attributeTypeList.size()) throw new IllegalArgumentException("Error building Phenotype: attribute list size not equal to type list size.");
		Iterator<ATTRIBUTE_TYPE> typeIter = attributeTypeList.iterator();
		for (PhenotypeAttribute attr : attributeList) {
			ATTRIBUTE_TYPE type = typeIter.next();
			if (!attr.isTypeCompatible(type)) {
				throw new IllegalArgumentException("Error building Phenotype: types not compatible with attributes.");
			}
		}
		Phenotype newPhenotype = new CorePhenotype(attributeList, attributeTypeList, phenotypeName);
		phenotypeList.add(newPhenotype);
	}
	
	private void joinPhenotypes(ACTION joinAction) {
		if (phenotypeList.size() < 2) throw new IllegalArgumentException("No join will be made because joining phenotypes requires at least two phenotypes.");
		
		if (addSourceDataFactor) {
			mergePhenotypesWithDataSourceFactor();
		} else {
//			phenotypeList = makeAttributeNamesDifferent(phenotypeList);
			Iterator<Phenotype> phenoIter = phenotypeList.iterator();
			Phenotype firstPhenotype = phenoIter.next();
			Phenotype secondPhenotype = phenoIter.next();
			Phenotype mergedPhenotype = mergeTwoPhenotypes(firstPhenotype, secondPhenotype, joinAction);
			while (phenoIter.hasNext()) mergedPhenotype = mergeTwoPhenotypes(mergedPhenotype, phenoIter.next(), joinAction);
			phenotypeList.clear();
			phenotypeList.add(mergedPhenotype);
		}
		
	}
	
	private void concatenatePhenotypes() {
		int nPheno = phenotypeList.size();
		if (nPheno < 2) throw new IllegalArgumentException("No phenotypes to join: must specify at least two phenotypes.");
		attributeList = new ArrayList<PhenotypeAttribute>();
		attributeTypeList = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		
		//create a new TaxaAttribute
		List<Taxon> jointTaxaList = new ArrayList<Taxon>();
		String attrName = phenotypeList.get(0).name();
		for (Phenotype pheno : phenotypeList) {
			TaxaAttribute someTaxa = pheno.taxaAttribute();
			if (someTaxa == null) throw new IllegalArgumentException(String.format("Phenotypes cannot be concatenated because %s has no taxa.", pheno.name()));
			jointTaxaList.addAll(someTaxa.allTaxaAsList());
		}
		TaxaAttribute myTaxaAttribute = new TaxaAttribute(jointTaxaList, attrName);
		attributeList.add(myTaxaAttribute);
		attributeTypeList.add(ATTRIBUTE_TYPE.taxa);
		
		//make a list of attribute names (except for taxa) and find out how many total observations there will be
		HashSet<String> attributeNameSet = new HashSet<>();
		int totalNumberOfObs = 0;
		for (Phenotype pheno : phenotypeList) {
			int nattr = pheno.numberOfAttributes();
			totalNumberOfObs += pheno.numberOfObservations();
			for (int i = 0; i < nattr; i++) {
				if (pheno.attributeType(i) != ATTRIBUTE_TYPE.taxa) attributeNameSet.add(pheno.attributeName(i));
			}
		}
		
		//create a new attribute for each name. Throw an error if attributes of the same name are of incompatible types.
		for (String attributeName : attributeNameSet) {
			ATTRIBUTE_TYPE myType = null;
			//take the type from the first phenotype with this name
			for (Phenotype pheno : phenotypeList) {
				int attrIndex = pheno.attributeIndexForName(attributeName);
				if (attrIndex > -1) {
					myType = pheno.attributeType(attrIndex);
					break;
				}
			}
			
			//create this attribute
			if (myType == ATTRIBUTE_TYPE.factor) {
				String[] values = new String[totalNumberOfObs];
				int nPreviousObs = 0;
				for (Phenotype pheno : phenotypeList) {
					int nObs = pheno.numberOfObservations();
					int attrIndex = pheno.attributeIndexForName(attributeName);
					if (attrIndex > -1) {
						System.arraycopy((String[]) pheno.attribute(attrIndex).allValues(), 0, values, nPreviousObs, nObs);
					} else {
						Arrays.fill(values, nPreviousObs, nPreviousObs + nObs, CategoricalAttribute.missingValue);
					}
					nPreviousObs += nObs;
				}
				attributeList.add(new CategoricalAttribute(attributeName, values));
				attributeTypeList.add(myType);
			} else {
				float[] values = new float[totalNumberOfObs];
				BitSet missing = new OpenBitSet(totalNumberOfObs);
				int nPreviousObs = 0;
				for (Phenotype pheno : phenotypeList) {
					int nObs = pheno.numberOfObservations();
					int attrIndex = pheno.attributeIndexForName(attributeName);
					if (attrIndex > -1) {
						PhenotypeAttribute thisAttribute = pheno.attribute(attrIndex);
						System.arraycopy((float[]) thisAttribute.allValues(), 0, values, nPreviousObs, nObs);
						for (int i = 0; i < nObs; i++) if (thisAttribute.isMissing(i)) missing.fastSet(nPreviousObs + i);
					} else {
						Arrays.fill(values, nPreviousObs, nPreviousObs + nObs, Float.NaN);
						missing.set(nPreviousObs, nPreviousObs + nObs);
					}
					nPreviousObs += nObs;
				}
				attributeList.add(new NumericAttribute(attributeName, values, missing));
				attributeTypeList.add(myType);
			}
		}
		
		sortAttributes();
		phenotypeList.clear();
		phenotypeList.add(new CorePhenotype(attributeList, attributeTypeList, phenotypeName));
	}
	
	private Phenotype mergeTwoPhenotypes(Phenotype pheno1, Phenotype pheno2, ACTION thisAction) {
		
		//build attribute name list for the new phenotype
		//use a TreeSet so that each name is used once and so that the names are sorted
		LinkedHashSet<String> attributeNameSet = new LinkedHashSet<>();
		
		int nAttributes1 = pheno1.numberOfAttributes();
		for (int a = 0; a < nAttributes1; a++) {
			if (pheno1.attributeType(a) != ATTRIBUTE_TYPE.taxa) attributeNameSet.add(pheno1.attributeName(a));
		}

		int nAttributes2 = pheno2.numberOfAttributes();
		for (int a = 0; a < nAttributes2; a++) {
			if (pheno2.attributeType(a) != ATTRIBUTE_TYPE.taxa) attributeNameSet.add(pheno2.attributeName(a));
		}
		
		ArrayList<String> attributeNameList = new ArrayList<>(attributeNameSet);
		ArrayList<ATTRIBUTE_TYPE> typeList = new ArrayList<>();
		
		//create a type list for these attribute names
		//if a name is used in both phenotypes, and the attribute types do not match, throw an error
		for (String name : attributeNameList) {
				 int ndx = pheno1.attributeIndexForName(name);
				 if (ndx > -1) typeList.add(pheno1.attributeType(ndx));
				 else {
					 typeList.add(pheno2.attributeType(pheno2.attributeIndexForName(name)));
			 }
		}
		
		//create a list of taxa to be included in the merged phenotype
		TaxaList outTaxa;
		if (thisAction.equals(ACTION.union)) outTaxa = TaxaListUtils.getAllTaxa(pheno1.taxa(), pheno2.taxa());
		else outTaxa = TaxaListUtils.getCommonTaxa(pheno1.taxa(), pheno2.taxa());
		
		//for each phenotype create a Multimap with Taxon as key, obs number as value
		Multimap<Taxon, Integer> pheno1ObservationMap = ArrayListMultimap.create();
		List<Taxon> pheno1Taxa = pheno1.taxaAttribute().allTaxaAsList();
		int taxonCount = 0;
		for (Taxon taxon:pheno1Taxa) pheno1ObservationMap.put(taxon, taxonCount++);
		
		Multimap<Taxon, Integer> pheno2ObservationMap = ArrayListMultimap.create();
		List<Taxon> pheno2Taxa = pheno2.taxaAttribute().allTaxaAsList();
		taxonCount = 0;
		for (Taxon taxon:pheno2Taxa) pheno2ObservationMap.put(taxon, taxonCount++);
		
		//create a list of new observations
		//first find the factors (if any) in common
		//for each factor in common record index in pheno1 and pheno2
		
		//  if pheno1 and pheno2 have any factors with the same name then these become merge factors
		//  only merge observations that have the same Taxa name and the same values for merge factors
		ArrayList<String> mergeFactors = new ArrayList<String>();
		ArrayList<int[]> mergeFactorIndex = new ArrayList<>();
		int numberOfMergeFactors;
		mergeFactors.addAll(listCommonNames(pheno1.attributeListOfType(ATTRIBUTE_TYPE.factor), pheno2.attributeListOfType(ATTRIBUTE_TYPE.factor)));
		mergeFactors.addAll(listCommonNames(pheno1.attributeListOfType(ATTRIBUTE_TYPE.covariate), pheno2.attributeListOfType(ATTRIBUTE_TYPE.covariate)));
		numberOfMergeFactors = mergeFactors.size();
		
		if (numberOfMergeFactors > 0) {
			Collections.sort(mergeFactors);
			for (String factorName : mergeFactors) {
				mergeFactorIndex.add(new int[]{pheno1.attributeIndexForName(factorName), pheno2.attributeIndexForName(factorName)});
			}
		}
		
		//code for creating a list of new observations
		ArrayList<int[]> mergeObservation = new ArrayList<>();
		ArrayList<Taxon> listOfTaxaForObservations = new ArrayList<>();
		for (Taxon taxon:outTaxa) {
			Collection<Integer> pheno1obs = pheno1ObservationMap.get(taxon);
			Collection<Integer> pheno2obs = pheno2ObservationMap.get(taxon);
			boolean[] wasPheno2ObsUsed = new boolean[pheno2obs.size()];
			Arrays.fill(wasPheno2ObsUsed, false);
			for (Integer obs1 : pheno1obs) {
				boolean wasPheno1ObsUsed = false;
				int obs2Count = 0;
				for (Integer obs2 : pheno2obs) {
					boolean mergeTheseObs = true;
					for (int[] ndx : mergeFactorIndex) {
						if (!pheno1.value(obs1, ndx[0]).equals(pheno2.value(obs2, ndx[1]))) mergeTheseObs = false;
						break;
					}
					if (mergeTheseObs) {
						wasPheno1ObsUsed = true;
						wasPheno2ObsUsed[obs2Count] = true;
						mergeObservation.add(new int[]{obs1, obs2});
						listOfTaxaForObservations.add(taxon);
//						break;
					}
					obs2Count++;
				}
				if (!wasPheno1ObsUsed) {
					mergeObservation.add(new int[]{obs1, -1});
					listOfTaxaForObservations.add(taxon);
				}
					
			}
			int obs2Count = 0;
			for (Integer obs2:pheno2obs) {
				if (!wasPheno2ObsUsed[obs2Count++]) {
					mergeObservation.add(new int[]{-1, obs2});
					listOfTaxaForObservations.add(taxon);
				}
			}
		}
		
		//create the new attribute and type lists. Add the taxa attribute to them.
		TaxaAttribute newTaxaAttribute = new TaxaAttribute(listOfTaxaForObservations, pheno1.taxaAttribute().name());
		ArrayList<PhenotypeAttribute> newAttributes = new ArrayList<>();
		ArrayList<ATTRIBUTE_TYPE> newTypes = new ArrayList<>();
		newAttributes.add(newTaxaAttribute);
		newTypes.add(ATTRIBUTE_TYPE.taxa);
		
		//create and add the other attributes (in the name list) to the new attribute and type lists
		int nAdd = attributeNameList.size();
		int nObs = mergeObservation.size();
		for (int i = 0; i < nAdd; i++) {
			ATTRIBUTE_TYPE myType = typeList.get(i);
			String attrName = attributeNameList.get(i);
			
			//the attrnum array stores the index of this attribute in pheno1 and pheno2
			int[] attrnum = new int[]{pheno1.attributeIndexForName(attrName), pheno2.attributeIndexForName(attrName)};
			
			switch(myType) {
			case data:
				
			case covariate:
				float[] myFloatData = new float[nObs];
				BitSet myMissing = new OpenBitSet(nObs);
				int obsCount = 0;
				for (int[] ndx : mergeObservation) {
					boolean pheno1HasNonmissingValue = attrnum[0] > -1 && ndx[0] > -1 && !pheno1.isMissing(ndx[0], attrnum[0]);
					boolean pheno2HasNonmissingValue = attrnum[1] > -1 && ndx[1] > -1 && !pheno2.isMissing(ndx[1], attrnum[1]);
					if (pheno1HasNonmissingValue && pheno2HasNonmissingValue) {
						throw new IllegalArgumentException("Data sets will not be joined because both phenotypes have values for " + attrName);
					} else if (pheno1HasNonmissingValue) {
						myFloatData[obsCount] = (Float) pheno1.value(ndx[0], attrnum[0]);
					} else if (pheno2HasNonmissingValue) {
						myFloatData[obsCount] = (Float) pheno2.value(ndx[1], attrnum[1]);
					} else {
						myFloatData[obsCount] = Float.NaN;
						myMissing.fastSet(obsCount);
					}
					obsCount++;
				}
				newAttributes.add(new NumericAttribute(attrName, myFloatData, myMissing));
				newTypes.add(myType);
				break;
			case factor:
				boolean isMergeAttribute;
				if (mergeFactors.contains(attrName)) isMergeAttribute = true;
				else isMergeAttribute = false;

				String[] myStringData = new String[nObs];
				obsCount = 0;
				for (int[] ndx : mergeObservation) {
					boolean pheno1HasNonmissingValue = attrnum[0] > -1 && ndx[0] > -1 && !pheno1.isMissing(ndx[0], attrnum[0]);
					boolean pheno2HasNonmissingValue = attrnum[1] > -1 && ndx[1] > -1 && !pheno2.isMissing(ndx[1], attrnum[1]);
					if (pheno1HasNonmissingValue && pheno2HasNonmissingValue) {
						if (isMergeAttribute) myStringData[obsCount] = (String) pheno1.value(ndx[0], attrnum[0]);
						else throw new IllegalArgumentException("Data sets will not be joined because both phenotypes have values for " + attrName);
					} else if (pheno1HasNonmissingValue) {
						myStringData[obsCount] = (String) pheno1.value(ndx[0], attrnum[0]);
					} else if (pheno2HasNonmissingValue) {
						myStringData[obsCount] = (String) pheno2.value(ndx[1], attrnum[1]);
					} else {
						myStringData[obsCount] = CategoricalAttribute.missingValue;
					}
					obsCount++;
				}
				newAttributes.add(new CategoricalAttribute(attrName, myStringData));
				newTypes.add(myType);
				break;
			}
		}
		
		return new CorePhenotype(newAttributes, newTypes, phenotypeName);
	}
	
	private List<Phenotype> makeAttributeNamesDifferent(List<Phenotype> phenoList) {
		int npheno = phenoList.size();
		
		//add all attribute names from phenoList to attrNameSet
		Multiset<String> attrNameSet = HashMultiset.create();
		for (Phenotype pheno : phenoList) {
			for (int i = 0; i < pheno.numberOfAttributes(); i++) {
				if (pheno.attributeType(i) != ATTRIBUTE_TYPE.taxa) attrNameSet.add(pheno.attributeName(i));
			}
		}

		//if an attribute name appears in more than one of phenotypes, append the phenotype number to it for each phenotype that uses it
		ArrayList<HashMap<String, String>> attrNameMapList = new ArrayList<>();
		for ( int i = 0; i < npheno; i++) {
			attrNameMapList.add(new HashMap<>());
		}
		
		for (String attrName : attrNameSet) {
			if (attrNameSet.count(attrName) > 1) {
				for ( int i = 0; i < npheno; i++) {
					int ndx = phenoList.get(i).attributeIndexForName(attrName);
					if (ndx > -1) {
						attrNameMapList.get(i).put(attrName, attrName + "." + i);
					}
				}
			}
		}
		
		//create a new phenotype with changed names as necessary
		List<Phenotype> newPhenotypeList = new ArrayList<Phenotype>();
		for ( int i = 0; i < npheno; i++) {
			if (attrNameMapList.get(i).size() > 0) {
				List<PhenotypeAttribute> attrList = phenoList.get(i).attributeListCopy();
				List<ATTRIBUTE_TYPE> typeList = phenoList.get(i).typeListCopy();
				for (String oldname : attrNameMapList.get(i).keySet()) {
					int ndx = phenoList.get(i).attributeIndexForName(oldname);
					if (ndx > -1) {
						PhenotypeAttribute newAttr = phenoList.get(i).attribute(ndx).changeName(attrNameMapList.get(i).get(oldname));
						attrList.set(ndx, newAttr);
					}
				}
				newPhenotypeList.add(new CorePhenotype(attrList, typeList, phenoList.get(i).name()));
			} else {
				newPhenotypeList.add(phenoList.get(i));
			}
		}
		
		return newPhenotypeList;
	}
	
	private void mergePhenotypesWithDataSourceFactor() {
		//create lists of attribute names and types
		HashMap<String, ATTRIBUTE_TYPE> typeMap = new HashMap<>();
		for (Phenotype pheno : phenotypeList) {
			for (int i = 0; i < pheno.numberOfAttributes(); i++) {
				if (pheno.attributeType(i) != ATTRIBUTE_TYPE.taxa) typeMap.put(pheno.attributeName(i), pheno.attributeType(i));
			}
		}
		
		ArrayList<String> attrNames = new ArrayList<>(typeMap.keySet());
		Collections.sort(attrNames);
		List<ATTRIBUTE_TYPE> newTypeList = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		for (String aname : attrNames) newTypeList.add(typeMap.get(aname));
		
		//number of observations in resulting Phenotype
		int totalObservations = 0;
		for (Phenotype ph :phenotypeList) totalObservations += ph.numberOfObservations();

		ArrayList<Object> attrData = new ArrayList<>();
		ArrayList<OpenBitSet> missingData = new ArrayList<>();
		for (String aname:attrNames) {
			missingData.add(new OpenBitSet(totalObservations)); 
			ATTRIBUTE_TYPE atype = typeMap.get(aname);
			if (atype == ATTRIBUTE_TYPE.covariate || atype == ATTRIBUTE_TYPE.data) {
				attrData.add(new float[0]);
			} else {
				attrData.add(new String[0]);
			}
		}
		
		List<Taxon> taxaList = new ArrayList<Taxon>();
		ArrayList<String> dataSourceNames = new ArrayList<>();
		int startObs = 0;
		for (int p = 0; p < phenotypeList.size(); p++) {
			Phenotype pheno = phenotypeList.get(p);
			int nobs = pheno.numberOfObservations();
			taxaList.addAll(pheno.taxa());
			String sourceName = pheno.name();
			for (int i = 0; i < pheno.numberOfObservations(); i++) dataSourceNames.add(sourceName);
			for (int a = 0; a < attrNames.size(); a++) {
				String aname = attrNames.get(a);
				int ndx = pheno.attributeIndexForName(aname);
				OpenBitSet md = missingData.get(a);
				Object data = attrData.get(a);
				if (ndx > -1) {
					if (data instanceof String[]) {
						String[] stringData = (String[]) data;
						for (int i = startObs; i < nobs; i++) {
							int origObs = i - startObs;
							if (pheno.isMissing(origObs, ndx)) md.fastSet(i);
							stringData[i] = (String) pheno.value(origObs, ndx);
						}
					} else {
						float[] floatData = (float[]) data;
						for (int i = startObs; i < nobs; i++) {
							int origObs = i - startObs;
							if (pheno.isMissing(origObs, ndx)) md.fastSet(i);
							floatData[i] = (Float) pheno.value(origObs, ndx);;
						}
					}
				} else {
					if (data instanceof String[]) {
						String[] stringData = (String[]) data;
						for (int i = startObs; i < nobs; i++) {
							md.fastSet(i);
							stringData[i] = CategoricalAttribute.missingValue;
						}
					} else {
						float[] floatData = (float[]) data;
						for (int i = startObs; i < nobs; i++) {
							md.fastSet(i);
							floatData[i] = Float.NaN;
						}
					}
				}
			}
		}
		
		List<PhenotypeAttribute> newAttrList = new ArrayList<PhenotypeAttribute>();
		for (int a = 0; a < attrNames.size(); a++) {
			Object data = attrData.get(a);
			if (data instanceof String[]) {
				newAttrList.add(new CategoricalAttribute(attrNames.get(a), (String[]) data));
			} else {
				newAttrList.add(new NumericAttribute(attrNames.get(a), (float[]) data, missingData.get(a)));
			}
		}
		phenotypeList.clear();
		StringBuilder newPhenotypeName = new StringBuilder();
		for (Phenotype pheno : phenotypeList) {
			if (newPhenotypeName.length() > 0) newPhenotypeName.append("+");
			newPhenotypeName.append(pheno.name());
		}
		phenotypeList.add(new CorePhenotype(newAttrList, newTypeList, newPhenotypeName.toString()));
	}
	
	private List<String> listCommonNames(List<PhenotypeAttribute> list1, List<PhenotypeAttribute> list2) {
		List<String> outlist = new ArrayList<String>();
		for (PhenotypeAttribute attr1 : list1) {
			String name1 = attr1.name();
			for (PhenotypeAttribute attr2 : list2) {
				if (name1.equalsIgnoreCase(attr2.name())) outlist.add(name1);
			}
		}
		return outlist;
	}
	
	private void sortAttributes() {
		//sort attributes, taxa first, then factors, then covariates, then data. Sorted by name within type
		int nAttr = attributeList.size();
		ArrayList<Integer> index = new ArrayList<>();
		for (int i = 0; i < nAttr; i++) index.add(i);
		
		Collections.sort(index, new Comparator<Integer>() {

			@Override
			public int compare(Integer o1, Integer o2) {
				ATTRIBUTE_TYPE type1 = attributeTypeList.get(o1);
				ATTRIBUTE_TYPE type2 = attributeTypeList.get(o2);
				if (type1 == type2) {
					return attributeList.get(o1).name().compareTo(attributeList.get(o2).name());
				} else {
					if (type1 == ATTRIBUTE_TYPE.taxa) return -1;
					else if (type2 == ATTRIBUTE_TYPE.taxa) return 1;
					else if (type1 == ATTRIBUTE_TYPE.factor) return -1;
					else if (type2 == ATTRIBUTE_TYPE.factor) return 1;
					else if (type1 == ATTRIBUTE_TYPE.covariate) return -1;
					else if (type2 == ATTRIBUTE_TYPE.covariate) return 1;
				}
				return 0;
			}
			
		});
		
		List<PhenotypeAttribute> unsortedAttributeList = attributeList;
		List<ATTRIBUTE_TYPE> unsortedAttributeTypeList = attributeTypeList;
		attributeList = new ArrayList<PhenotypeAttribute>();
		attributeTypeList = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
		for (Integer ndx : index) {
			attributeList.add(unsortedAttributeList.get(ndx));
			attributeTypeList.add(unsortedAttributeTypeList.get(ndx));
		}
		
	}
	
	private void removeAllObservationsWithMissingValues() {
		List<Phenotype> newList = new ArrayList<Phenotype>();
		for (Phenotype pheno : phenotypeList) {
			int nobs = pheno.numberOfObservations();
			BitSet anyMissing = new OpenBitSet(nobs);
			pheno.attributeStream().map(PhenotypeAttribute::missing).reduce(anyMissing, (a,b) -> {a.or(b); return a;});
			int[] nonMissingIndex = IntStream.range(0, nobs).filter(i -> !anyMissing.fastGet(i)).toArray();
			List<PhenotypeAttribute> subsetList = pheno.attributeStream().map(a -> a.subset(nonMissingIndex, a.name())).collect(Collectors.toList());
			Phenotype myNewPhenotype = new CorePhenotype(subsetList, pheno.typeListCopy(), "NonMissing_" + pheno.name());
			newList.add(myNewPhenotype);
		}
		phenotypeList = newList;
	}
	
}
