package net.maizegenetics.phenotype;

import java.util.List;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.util.BitSet;

public interface PhenotypeAttribute {

	/**
	 * @param obs	the observation number
	 * @return	the value for this observation
	 */
	Object value(int obs);
	
	/**
	 * The return value will typically be a primitive array whose type depends on the sub class
	 * @return	the values of this Attribute for all observations, in order by observation number
	 */
	Object allValues();
	
	/**
	 * @param obs		an array of observation numbers
	 * @param newName	The name for the new PhenotypeAttribute
	 * @return			a new PhenotypeAttribute equivalent to this one but with a subset of observations specified by obs
	 */
	PhenotypeAttribute subset(int[] obs, String newName);
	
	/**
	 * @param newName	a new name for this PhenotypeAttribute
	 * @return			a copy of this PhenotypeAttribute with newName as its name
	 */
	PhenotypeAttribute changeName(String newName);
	
	/**
	 * @param obs	the observation number
	 * @return	if the value of the observation is missing, true, otherwise false.
	 */
	boolean isMissing(int obs);
	
	/**
	 * @return an array whose elements are true for each missing observation, false for observations with valid values.
	 */
	BitSet missing();
	
	/**
	 * @return	the name of this Attribute
	 */
	String name();
	
	/**
	 * @return	the number of observations in this Attribute
	 */
	int size();
	
	/**
	 * @return	a list of attribute types that are compatible with this attribute
	 * For instance, a CategoricalAttribute can only have a type of factor while a NumericPhenotype can have a type of data or covariate.
	 */
	List<ATTRIBUTE_TYPE> getCompatibleTypes();
	
	/**
	 * @param type	an attribute type
	 * @return	true if type is compatible with this attribute
	 */
	boolean isTypeCompatible(ATTRIBUTE_TYPE type);
	
}
