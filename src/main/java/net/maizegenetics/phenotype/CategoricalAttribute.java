package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;


public class CategoricalAttribute implements PhenotypeAttribute {

	private final String name;
	private final int[] values;
	private final BitSet missing;
	private final String[] categoryNames;
	private final HashMap<String, Integer> nameToIndexMap;
	
	public static final String missingValue = "?";
	private static final List<ATTRIBUTE_TYPE> myAllowedTypes;
	static{
		myAllowedTypes = new ArrayList<ATTRIBUTE_TYPE>();
		myAllowedTypes.add(ATTRIBUTE_TYPE.factor);
	}
	
	public CategoricalAttribute(String name, String[] stringValues) {
		this.name = name;
		int n = stringValues.length;
		missing = new OpenBitSet(n);
		values = new int[n];
		
		TreeSet<String> labelSet = Arrays.stream(stringValues).collect(TreeSet::new, TreeSet::add, TreeSet::addAll);
		
		int nlevels = labelSet.size();
		categoryNames = new String[nlevels];
		labelSet.toArray(categoryNames);
		
		nameToIndexMap = new HashMap<>();
		for (int i = 0; i < nlevels; i++) nameToIndexMap.put(categoryNames[i], i);
		for (int i = 0; i < n; i++) {
			values[i] = nameToIndexMap.get(stringValues[i]);
		}
	}
	
	/**
	 * @param obs	the observation number
	 * @return		zero-based index of category levels for this observation; -1 = missing value
	 */
	public int intValue(int obs) {
		return values[obs];
	}
	
	/**
	 * @return	zero-based category level index for each observation; -1 = missing value
	 */
	public int[] allIntValues() {
		return values;
	}
	
	/**
	 * @param index		a zero-based category level index 
	 * @return			the String label for this level
	 */
	public String attributeLabelForIndex(int index) {
		return categoryNames[index];
	}
	
	/**
	 * @param label		the name of a Category level
	 * @return			the zero-based category level index for this name
	 */
	public int indexForAttrLabel(String label) {
		return nameToIndexMap.get(label);
	}

	/**
	 * @param obs	the observation number
	 * @return		the name of the Category for this observation
	 */
	public String label(int obs) {
		return attributeLabelForIndex(values[obs]);
	}
	
	/**
	 * @return	the Category names for all observations
	 */
	public String[] allLabels() {
		int n = values.length;
		String[] labels = new String[n];
		for (int i = 0; i < n; i++) {
			if (values[i] == -1) labels[i] = missingValue;
			labels[i] = categoryNames[values[i]];
		}
		return labels;
	}
	
	/**
	 * @return	a list each unique Category name ordered by index
	 */
	public List<String> labelList() {
		return new ArrayList<String>(Arrays.asList(categoryNames));
	}
	
	/**
	 * @return	the number of Category levels
	 */
	public int numberOfLevels() {
		return categoryNames.length;
	}
	
	/**
	 * @param 	level
	 * @return	the observation numbers corresponding to this level index
	 */
	public int[] whichObservations(int level) {
		int nvalues = values.length;
		int[] obs = new int[nvalues];
		int obsCount = 0;
		for (int i = 0; i < nvalues; i++) {
			if (values[i] == level) obs[obsCount++] = i;
		}
		return Arrays.copyOf(obs, obsCount);
	}

	/**
	 * @param factorName	the name of a Category level
	 * @return				the observation numbers corresponding to this name
	 */
	public int[] whichObservations(String factorName) {
		int nvalues = values.length;
		int[] obs = new int[nvalues];
		int obsCount = 0;
		for (int i = 0; i < nvalues; i++) {
			if (categoryNames[values[i]].equals(factorName)) obs[obsCount++] = i;
		}
		return Arrays.copyOf(obs, obsCount);
	}

	@Override
	public Object value(int obs) {
		return categoryNames[values[obs]];
	}

	@Override
	public Object allValues() {
		return allLabels();
	}

	@Override
	public PhenotypeAttribute subset(int[] obs, String newName) {
		int n = obs.length;
		String[] labels = new String[n];
		for (int i = 0; i < n; i++) {
			if (values[obs[i]] == -1) labels[i] = missingValue;
			else labels[i] = categoryNames[values[obs[i]]];
		}
		if (newName == null) newName = name;
		return new CategoricalAttribute(newName, labels);
	}

	@Override
	public PhenotypeAttribute changeName(String newName) {
		return new CategoricalAttribute(newName, allLabels());
	}

	@Override
	public boolean isMissing(int obs) {
		return missing.fastGet(obs);
	}

	@Override
	public BitSet missing() {
		return missing;
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public int size() {
		return values.length;
	}

	@Override
	public List<ATTRIBUTE_TYPE> getCompatibleTypes() {
		return myAllowedTypes;
	}

	@Override
	public boolean isTypeCompatible(ATTRIBUTE_TYPE type) {
		return myAllowedTypes.contains(type);
	}

	@Override
	public String toString() {
		return name();
	}

}
