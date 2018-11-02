/*
 *  FilterTaxa
 *
 *  Created on Nov 29, 2016
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.taxa.TaxaList;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * @author Terry Casstevens
 */
final public class FilterTaxa implements Filter {

    public static enum FILTER_TAXA_ATTRIBUTES {

        filterName,
        minNotMissing, minHeterozygous, maxHeterozygous,
        includeTaxa, taxaList;
    }

    private final Map<FILTER_TAXA_ATTRIBUTES, Object> myAttributes = new LinkedHashMap<>();

    public FilterTaxa(Map<String, Object> values) {
        for (Map.Entry<String, Object> current : values.entrySet()) {
            if (current.getValue() == null) {
                continue;
            }
            if ((current.getValue() instanceof String) && (((String) current.getValue()).isEmpty())) {
                continue;
            }
            FILTER_TAXA_ATTRIBUTES attribute = FILTER_TAXA_ATTRIBUTES.valueOf(current.getKey());
            myAttributes.put(attribute, current.getValue());
        }

        if ((filterName() == null) || (filterName().length() == 0)) {
            myAttributes.put(FILTER_TAXA_ATTRIBUTES.filterName, "Filter");
        }

    }

    public Map<FILTER_TAXA_ATTRIBUTES, Object> attributes() {
        return Collections.unmodifiableMap(myAttributes);
    }

    @Override
    public int numAttributes() {
        return myAttributes.size();
    }

    public String filterName() {
        return (String) myAttributes.get(FILTER_TAXA_ATTRIBUTES.filterName);
    }

    public double minNotMissing() {
        Double value = (Double) myAttributes.get(FILTER_TAXA_ATTRIBUTES.minNotMissing);
        if (value == null) {
            return 0.0;
        } else {
            return value;
        }
    }

    public double minHeterozygous() {
        Double value = (Double) myAttributes.get(FILTER_TAXA_ATTRIBUTES.minHeterozygous);
        if (value == null) {
            return 0.0;
        } else {
            return value;
        }
    }

    public double maxHeterozygous() {
        Double value = (Double) myAttributes.get(FILTER_TAXA_ATTRIBUTES.maxHeterozygous);
        if (value == null) {
            return 1.0;
        } else {
            return value;
        }
    }

    public boolean includeTaxa() {
        Boolean value = (Boolean) myAttributes.get(FILTER_TAXA_ATTRIBUTES.includeTaxa);
        if (value == null) {
            return true;
        } else {
            return value;
        }
    }

    public TaxaList taxaList() {
        return (TaxaList) myAttributes.get(FILTER_TAXA_ATTRIBUTES.taxaList);
    }

}
