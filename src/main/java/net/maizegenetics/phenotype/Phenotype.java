package net.maizegenetics.phenotype;

import java.util.List;
import java.util.stream.Stream;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.TableReport;

/**
 * Phenotype represents phenotype data as a two dimensional matrix. Rows are observations. Columns are attributes.
 * 
 * @author pbradbury
 *
 */
public interface Phenotype  extends TableReport {
	public enum ATTRIBUTE_TYPE {data, covariate, factor, taxa};

	/**
	 * @param obs	an observation number
	 * @param attrnum	the index or column number of the attribute
	 * @return	the attribute value for the observation
	 */
	public Object value(int obs, int attrnum);
	
	/**
	 * @param obs	an observation number
	 * @param attrnum	the index or column number of the attribute
	 * @return	true if the phenotype value of attrnum is missing for obs
	 */
	public boolean isMissing(int obs, int attrnum);
	
	/**
	 * @param attrnum	the index or column number of the attribute
	 * @return	the PhenotypeAttribute 
	 */
	public PhenotypeAttribute attribute(int attrnum);

	/**
	 * @param attribute	a PhenotypeAttribute
	 * @return	the index of attribute in this Phenotype, -1 if the attribute is not contained in this Phenotype
	 */
	public int indexOfAttribute(PhenotypeAttribute attribute);
	
	/**
	 * @return	the number of attributes or columns in the Phenotype
	 */
	public int numberOfAttributes();
	
	/**
	 * @return	the number of observations or rows in this phenotype
	 */
	public int numberOfObservations();

	/**
	 * @return	the unique set of taxa in this Phenotype. 
	 * Each taxon will be in the TaxaList only once though there may be more than one observation of the taxon in the Phenotype.
	 */
	public TaxaList taxa();

	/**
	 * @param type	an attribute type
	 * @return	the number of attributes of this type
	 */
	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type);
	
	/**
	 * @param type	an attribute type
	 * @return	an array of the indices of all the attributes of this type
	 */
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type);
	
	/**
	 * @param attrnum	the index or column number of an attribute
	 * @return	the type of the attribute
	 */
	public ATTRIBUTE_TYPE attributeType(int attrnum);
	
	/**
	 * @param attrnum	the index or column number of the attribute
	 * @return	the name of the attribute
	 */
	public String attributeName(int attrnum);
	
	/**
	 * @return	the name of this Phenotype
	 */
	public String name();
	
	/**
	 * @return	true, if this Phenotype type has one and only one TaxaAttribute
	 */
	public boolean hasTaxaAttribute();
	
	/**
	 * @return	this Phenotype's TaxaAttribute, which is the Taxa represented by the observations in this Phenotype
	 */
	public TaxaAttribute taxaAttribute();
	
	/**
	 * @return	a shallow copy of the attribute list for this Phenotype
	 */
	public List<PhenotypeAttribute> attributeListCopy();
	
	/**
	 * @param type	an attribute type
	 * @return	a List of all attributes of the requested type
	 */
	public List<PhenotypeAttribute> attributeListOfType(ATTRIBUTE_TYPE type);
	
	/**
	 * @return	a sequential Stream of all PhenotypeAttributes
	 */
	public Stream<PhenotypeAttribute> attributeStream();
	
	/**
	 * @return	a sequential Stream of all data type PhenotypeAttributes as NumericAttributes
	 */
	public Stream<NumericAttribute> dataAttributeStream();
	
	/**
	 * @return	a sequential Stream of all covariate type PhenotypeAttributes as NumericAttributes
	 */
	public Stream<NumericAttribute> covariateAttributeStream();
	
	/**
	 * @return	a sequential Stream of all factor type PhenotypeAttributes as CategoricalAttributes
	 */
	public Stream<CategoricalAttribute> factorAttributeStream();
	
	/**
	 * @return	a shallow copy of the attribute type list for this Phenotype
	 */
	public List<ATTRIBUTE_TYPE> typeListCopy();

	/**
	 * @param name	the name of an attribute
	 * @return	the index of the attribute with this name, -1 if there is no attribute with this name
	 */
	public int attributeIndexForName(String name);
	
	/**
	 * @return true, if the number of observations is greater than the number of unique taxa, false otherwise
	 */
	public boolean areTaxaReplicated();
	
	/**
	 * @return	this Phenotype as a CorePhenotype.
	 * <p>
	 * CorePhenotype may be more efficient in some cases than FilterPhenotype, though the difference will probably be small.
	 * For very large filtered data sets consider converting to CorePhenotype. The operation performed on a CorePhenotype will return
	 * the same Phenotype.
	 */
	public Phenotype asCorePhenotype();
}