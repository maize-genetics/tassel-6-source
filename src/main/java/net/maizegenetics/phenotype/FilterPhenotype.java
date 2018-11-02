package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

public class FilterPhenotype implements Phenotype {
	private int[] myRowRedirect;
	private int numberOfObservations;
	private CorePhenotype basePhenotype;
	private String name;
	
	FilterPhenotype(CorePhenotype basePheno, List<Taxon> taxaToKeep, String name) {
		this.basePhenotype = basePheno;
		TaxaAttribute myTaxaAttribute = basePheno.taxaAttribute();
		List<Taxon> baseTaxa = myTaxaAttribute.allTaxaAsList();
		int numberOfBaseObservations = basePheno.numberOfObservations();
		myRowRedirect = new int[numberOfBaseObservations];
		int obsCount = 0;
		int baseTaxonCount = 0;
		for (Taxon baseTaxon:baseTaxa) {
			if (taxaToKeep.contains(baseTaxon)) myRowRedirect[obsCount++] = baseTaxonCount;
			baseTaxonCount++;		
		}
		myRowRedirect = Arrays.copyOf(myRowRedirect, obsCount);
		numberOfObservations = obsCount;
		this.name = name;
	}
	
	FilterPhenotype(FilterPhenotype basePheno, List<Taxon> taxaToKeep, String name) {
		int numberOfBaseObservations = basePheno.numberOfObservations();
		TaxaAttribute myTaxaAttribute = basePheno.taxaAttribute();
		myRowRedirect = new int[numberOfBaseObservations];
		List<Taxon> baseTaxa = myTaxaAttribute.allTaxaAsList();
		int obsCount = 0;
		int baseTaxonCount = 0;
		for (Taxon baseTaxon:baseTaxa) {
			if (taxaToKeep.contains(baseTaxon)) myRowRedirect[obsCount++] = basePheno.myRowRedirect[baseTaxonCount];
			baseTaxonCount++;		
		}
		myRowRedirect = Arrays.copyOf(myRowRedirect, obsCount);
		this.basePhenotype = basePheno.basePhenotype;
		numberOfObservations = obsCount;
		this.name = name;
	}
	
	static FilterPhenotype getInstance(Phenotype basePheno, List<Taxon> taxaToKeep, String name) {
		if (basePheno instanceof CorePhenotype) return new FilterPhenotype((CorePhenotype) basePheno, taxaToKeep, name);
		if (basePheno instanceof FilterPhenotype) return new FilterPhenotype((FilterPhenotype) basePheno, taxaToKeep, name);
		else return null; 
	}
	
	//TableReport methods
	@Override
	public Object[] getTableColumnNames() {
		return basePhenotype.getTableColumnNames();
	}

	@Override
	public String getTableTitle() {
		return name;
	}

	@Override
	public int getColumnCount() {
		return basePhenotype.getColumnCount();
	}

	@Override
	public long getRowCount() {
		return numberOfObservations;
	}

	@Override
	public long getElementCount() {
		return getRowCount() * getColumnCount();
	}

	@Override
	public Object[] getRow(long row) {
		return basePhenotype.getRow(myRowRedirect[(int) row]);
	}

	@Override
	public Object getValueAt(long row, int col) {
		return basePhenotype.getValueAt(myRowRedirect[(int) row], col);
	}

	//Phenotype methods
	
	@Override
	public Object value(int obs, int attrnum) {
		return basePhenotype.value(myRowRedirect[obs], attrnum);
	}

	@Override
	public boolean isMissing(int obs, int attrnum) {
		return basePhenotype.isMissing(myRowRedirect[obs], attrnum);
	}

	@Override
	public PhenotypeAttribute attribute(int attrnum) {
		return basePhenotype.attribute(attrnum).subset(myRowRedirect, null);
	}

	@Override
	public int indexOfAttribute(PhenotypeAttribute attribute) {
		return basePhenotype.indexOfAttribute(attribute);
	}

	@Override
	public List<PhenotypeAttribute> attributeListCopy() {
		List<PhenotypeAttribute> attrList = new ArrayList<PhenotypeAttribute>();
		for (int i = 0; i < numberOfAttributes(); i++) attrList.add(attribute(i));
		return attrList;
	}

	@Override
	public List<PhenotypeAttribute> attributeListOfType(ATTRIBUTE_TYPE type) {
		int[] indices = attributeIndicesOfType(type);
		ArrayList<PhenotypeAttribute> attrList = new ArrayList<>();
		for (int ndx:indices) attrList.add(attribute(ndx));
		return attrList;
	}

	@Override
	public Stream<PhenotypeAttribute> attributeStream() {
		return IntStream.range(0, numberOfAttributes()).mapToObj(i -> attribute(i));
	}

	@Override
	public Stream<NumericAttribute> dataAttributeStream() {
		return IntStream.iterate(0, i -> i + 1).limit(numberOfAttributes())
				.filter(i -> attributeType(i) == ATTRIBUTE_TYPE.data)
				.mapToObj(i -> (NumericAttribute) attribute(i));
	}

	@Override
	public Stream<NumericAttribute> covariateAttributeStream() {
		return IntStream.iterate(0, i -> i + 1).limit(numberOfAttributes())
				.filter(i -> attributeType(i) == ATTRIBUTE_TYPE.covariate)
				.mapToObj(i -> (NumericAttribute) attribute(i));
	}

	@Override
	public Stream<CategoricalAttribute> factorAttributeStream() {
		return IntStream.iterate(0, i -> i + 1).limit(numberOfAttributes())
				.filter(i -> attributeType(i) == ATTRIBUTE_TYPE.factor)
				.mapToObj(i -> (CategoricalAttribute) attribute(i));
	}

	@Override
	public List<ATTRIBUTE_TYPE> typeListCopy() {
		return basePhenotype.typeListCopy();
	}

	@Override
	public int numberOfAttributes() {
		return basePhenotype.numberOfAttributes;
	}

	@Override
	public int numberOfObservations() {
		return numberOfObservations;
	}

	@Override
	public TaxaList taxa() {
		if (!basePhenotype.hasTaxaAttribute()) return null;
		Taxon[] taxaArray = basePhenotype.taxaAttribute().allTaxa();
		TreeSet<Taxon> taxaSet = new TreeSet<>();
		for (int ndx : myRowRedirect) taxaSet.add(taxaArray[ndx]);
		
		return new TaxaListBuilder().addAll(taxaSet).build();
	}

	@Override
	public int numberOfAttributesOfType(ATTRIBUTE_TYPE type) {
		return basePhenotype.numberOfAttributesOfType(type);
	}

	@Override
	public int[] attributeIndicesOfType(ATTRIBUTE_TYPE type) {
		return basePhenotype.attributeIndicesOfType(type);
	}

	@Override
	public ATTRIBUTE_TYPE attributeType(int attrnum) {
		return basePhenotype.attributeType(attrnum);
		
	}

	@Override
	public String attributeName(int attrnum) {
		return basePhenotype.attributeName(attrnum);
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public boolean hasTaxaAttribute() {
		return basePhenotype.hasTaxaAttribute();
	}

	@Override
	public TaxaAttribute taxaAttribute() {
		return (TaxaAttribute) basePhenotype.taxaAttribute().subset(myRowRedirect, null);
	}

	@Override
	public int attributeIndexForName(String name) {
		return basePhenotype.attributeIndexForName(name);
	}

	@Override
	public boolean areTaxaReplicated() {
		int numberOfUniqueTaxa = taxa().size();
		return (numberOfUniqueTaxa > numberOfObservations);
	}

	@Override
	public Phenotype asCorePhenotype() {
		List<PhenotypeAttribute> allAttributes = IntStream.iterate(0, i-> i + 1).limit(numberOfAttributes())
				.mapToObj(i -> attribute(i)).collect(Collectors.toList());
		return new PhenotypeBuilder().fromAttributeList(allAttributes, basePhenotype.typeListCopy()).build().get(0);
	}

}
