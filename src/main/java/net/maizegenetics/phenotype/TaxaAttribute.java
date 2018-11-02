package net.maizegenetics.phenotype;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 * @author pbradbury
 *
 */
public class TaxaAttribute implements PhenotypeAttribute {
	private final ArrayList<Taxon> taxaList;
	private final int numberOfTaxa;
	private final String name;
	public static final String DEFAULT_NAME = "Taxa";
	
	private static final List<ATTRIBUTE_TYPE> myAllowedTypes;
	static{
		myAllowedTypes = new ArrayList<ATTRIBUTE_TYPE>();
		myAllowedTypes.add(ATTRIBUTE_TYPE.taxa);
	}

	public TaxaAttribute(List<Taxon> taxa, String name) {
		this.name = DEFAULT_NAME;
		if (taxa instanceof ArrayList) taxaList = (ArrayList<Taxon>) taxa;
		else taxaList = new ArrayList<>(taxa);
		numberOfTaxa = taxa.size();
	}
	
	public TaxaAttribute(List<Taxon> taxa) {
		this(taxa, DEFAULT_NAME);
	}
	
	public Taxon[] allTaxa() {
		return taxaList.toArray(new Taxon[numberOfTaxa]);
	}
	
	public List<Taxon> allTaxaAsList() {
		return new ArrayList<Taxon>(taxaList);
	}
	
	public Taxon taxon(int obs) {
		return taxaList.get(obs);
	}
	
	@Override
	public Object value(int obs) {
		return taxaList.get(obs);
	}

	@Override
	public Object allValues() {
		Taxon[] taxaArray = new Taxon[numberOfTaxa];
		return taxaList.toArray(taxaArray);
	}

	@Override
	public PhenotypeAttribute subset(int[] obs, String newName) {
		int n = obs.length;
		ArrayList<Taxon> subset = new ArrayList<Taxon>();
		for (int i = 0; i < n; i++) subset.add(taxaList.get(obs[i]));
		if (newName == null) newName = name;
		return new TaxaAttribute(subset, newName);
	}

	@Override
	public PhenotypeAttribute changeName(String newName) {
		return new TaxaAttribute(taxaList, newName);
	}

	@Override
	public boolean isMissing(int obs) {
		return false;
	}

	@Override
	public BitSet missing() {
		return new OpenBitSet(numberOfTaxa);
	}

	@Override
	public String name() {
		return name;
	}

	@Override
	public int size() {
		return numberOfTaxa;
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
