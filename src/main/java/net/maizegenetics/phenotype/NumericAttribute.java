package net.maizegenetics.phenotype;



import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class NumericAttribute implements PhenotypeAttribute {
	private final String name;
	private final float[] values;
	private final BitSet missing;
	
	private static final List<ATTRIBUTE_TYPE> myAllowedTypes;
	static{
		myAllowedTypes = new ArrayList<ATTRIBUTE_TYPE>();
		myAllowedTypes.add(ATTRIBUTE_TYPE.data);
		myAllowedTypes.add(ATTRIBUTE_TYPE.covariate);
	}

	public NumericAttribute(String name, float[] values, BitSet missing) {
		this.name = name;
		this.values = values;
		this.missing = missing;
	}
	
	public NumericAttribute(String name, float[] values) {
            this.name = name;
            this.values = values;
            int n = values.length;
            missing = new OpenBitSet(n);
            for (int i = 0; i < n; i++) {
                if (!Float.isFinite(values[i])) {
                    missing.set(i);
                }
            }
	}
	
        public NumericAttribute(String name, double[] doubleValues) {
            this.name = name;
            int n = doubleValues.length;
            missing = new OpenBitSet(n);
            values = new float[n];
            for (int i = 0; i < n; i++) {
                if (!Double.isFinite(doubleValues[i])) {
                    missing.set(i);
                    values[i] = Float.NaN;
                } else values[i] = (float) doubleValues[i];
            }
        }
        
	
	public float floatValue(int obs) {
		return values[obs];
	}
	
	public float[] floatValues() {
		return values;
	}
	
	public double doubleValue(int obs) {
		return (double) values[obs];
	}
	
	public double[] doubleValues() {
		double[] dblVal = new double[values.length];
		for (int i = 0; i < values.length; i++) {
			dblVal[i] = (double) values[i];
		}
		return dblVal;
	}
	
	@Override
	public Object value(int obs) {
		return new Float(values[obs]);
	}

	@Override
	public Object allValues() {
		return values;
	}

	@Override
	public PhenotypeAttribute subset(int[] obs, String newName) {
		int n = obs.length;
		float[] valueSubset = new float[n];
		OpenBitSet missingSubset = new OpenBitSet(n);
		for (int i = 0; i < n; i++) {
			valueSubset[i] = values[obs[i]];
			if (missing.fastGet(obs[i])) missingSubset.fastSet(i);
		}
		if (newName == null) newName = name;
		return new NumericAttribute(newName, valueSubset, missingSubset);
	}

	@Override
	public PhenotypeAttribute changeName(String newName) {
		return new NumericAttribute(newName, values, missing);
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
