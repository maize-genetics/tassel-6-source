/*
 *  FilterSite
 *
 *  Created on Jun 30, 2014
 */
package net.maizegenetics.dna.snp;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;

/**
 * @author Terry Casstevens
 */
final public class FilterSite implements Filter {

    public enum FILTER_SITES_ATTRIBUTES {
        filterName,
        siteMinCount, siteMinAlleleFreq, siteMaxAlleleFreq, siteRangeFilterType,
        startSite, endSite, startChr, startPos, endChr,
        endPos, includeSites, siteNames, chrPosFile,
        positionList, removeMinorSNPStates, removeSitesWithIndels, bedFile,
        minHeterozygous, maxHeterozygous;
    }

    public enum SITE_RANGE_FILTER_TYPES {
        NONE, SITES, POSITIONS
    }

    private final Map<FILTER_SITES_ATTRIBUTES, Object> myAttributes = new LinkedHashMap<>();

    public FilterSite(Map<String, Object> values) {
        for (Map.Entry<String, Object> current : values.entrySet()) {
            if (current.getValue() == null) {
                continue;
            }
            if ((current.getValue() instanceof String) && (((String) current.getValue()).isEmpty())) {
                continue;
            }
            FILTER_SITES_ATTRIBUTES attribute = FILTER_SITES_ATTRIBUTES.valueOf(current.getKey());
            myAttributes.put(attribute, current.getValue());
        }

        if ((filterName() == null) || (filterName().length() == 0)) {
            myAttributes.put(FILTER_SITES_ATTRIBUTES.filterName, "Filter");
        }

        if (siteFilterType() == SITE_RANGE_FILTER_TYPES.SITES) {
            System.out.println("startSite: " + startSite() + "  endSite: " + endSite());
            if ((startSite() == -1) || (endSite() == -1)) {
                throw new IllegalArgumentException("Filter: init: both start site and end site must be specified.");
            } else if (startSite() > endSite()) {
                throw new IllegalArgumentException("Filter: init: start site can't be larger than end site.");
            }
        }

        if (siteFilterType() == SITE_RANGE_FILTER_TYPES.POSITIONS) {
            if (startChr() == null || endChr() == null) {
                throw new IllegalArgumentException("Filter: init: start chr and end chr must be specified.");
            }
        }

        if (siteMinCount() < 0) {
            throw new IllegalArgumentException("Filter: init: site min count can't be less than zero.");
        }

        if ((siteMinAlleleFreq() < 0.0) || (siteMinAlleleFreq() > 1.0)) {
            throw new IllegalArgumentException("Filter: init: site min allele freq must be between 0.0 and 1.0");
        }

        if ((siteMaxAlleleFreq() < 0.0) || (siteMaxAlleleFreq() > 1.0)) {
            throw new IllegalArgumentException("Filter: init: site max allele freq must be between 0.0 and 1.0");
        }

    }

    public Map<FILTER_SITES_ATTRIBUTES, Object> attributes() {
        return Collections.unmodifiableMap(myAttributes);
    }

    @Override
    public int numAttributes() {
        return myAttributes.size();
    }

    public String filterName() {
        return (String) myAttributes.get(FILTER_SITES_ATTRIBUTES.filterName);
    }

    public int siteMinCount() {
        Integer value = (Integer) myAttributes.get(FILTER_SITES_ATTRIBUTES.siteMinCount);
        if (value == null) {
            return 0;
        } else {
            return value;
        }
    }

    public double siteMinAlleleFreq() {
        Double value = (Double) myAttributes.get(FILTER_SITES_ATTRIBUTES.siteMinAlleleFreq);
        if (value == null) {
            return 0.0;
        } else {
            return value;
        }
    }

    public double siteMaxAlleleFreq() {
        Double value = (Double) myAttributes.get(FILTER_SITES_ATTRIBUTES.siteMaxAlleleFreq);
        if (value == null) {
            return 1.0;
        } else {
            return value;
        }
    }

    public double minHeterozygous() {
        Double value = (Double) myAttributes.get(FILTER_SITES_ATTRIBUTES.minHeterozygous);
        if (value == null) {
            return 0.0;
        } else {
            return value;
        }
    }

    public double maxHeterozygous() {
        Double value = (Double) myAttributes.get(FILTER_SITES_ATTRIBUTES.maxHeterozygous);
        if (value == null) {
            return 1.0;
        } else {
            return value;
        }
    }

    public boolean removeMinorSNPStates() {
        Boolean value = (Boolean) myAttributes.get(FILTER_SITES_ATTRIBUTES.removeMinorSNPStates);
        if (value == null) {
            return false;
        } else {
            return value;
        }
    }

    public boolean removeSitesWithIndels() {
        Boolean value = (Boolean) myAttributes.get(FILTER_SITES_ATTRIBUTES.removeSitesWithIndels);
        if (value == null) {
            return false;
        } else {
            return value;
        }
    }

    public SITE_RANGE_FILTER_TYPES siteFilterType() {
        SITE_RANGE_FILTER_TYPES value = (SITE_RANGE_FILTER_TYPES) myAttributes.get(FILTER_SITES_ATTRIBUTES.siteRangeFilterType);
        if (value == null) {
            return SITE_RANGE_FILTER_TYPES.NONE;
        } else {
            return value;
        }
    }

    public int startSite() {
        Integer value = (Integer) myAttributes.get(FILTER_SITES_ATTRIBUTES.startSite);
        if (value == null) {
            System.out.println("start site is null");
            return -1;
        } else {
            return value;
        }
    }

    public int endSite() {
        Integer value = (Integer) myAttributes.get(FILTER_SITES_ATTRIBUTES.endSite);
        if (value == null) {
            return -1;
        } else {
            return value;
        }
    }

    public Chromosome startChr() {
        Chromosome value = (Chromosome) myAttributes.get(FILTER_SITES_ATTRIBUTES.startChr);
        if (value == null) {
            return null;
        } else {
            return value;
        }
    }

    public int startPos() {
        Integer value = (Integer) myAttributes.get(FILTER_SITES_ATTRIBUTES.startPos);
        if (value == null) {
            return -1;
        } else {
            return value;
        }
    }

    public Chromosome endChr() {
        Chromosome value = (Chromosome) myAttributes.get(FILTER_SITES_ATTRIBUTES.endChr);
        if (value == null) {
            return null;
        } else {
            return value;
        }
    }

    public int endPos() {
        Integer value = (Integer) myAttributes.get(FILTER_SITES_ATTRIBUTES.endPos);
        if (value == null) {
            return -1;
        } else {
            return value;
        }
    }

    public boolean includeSites() {
        Boolean value = (Boolean) myAttributes.get(FILTER_SITES_ATTRIBUTES.includeSites);
        if (value == null) {
            return true;
        } else {
            return value;
        }
    }

    public List<String> siteNames() {
        List<String> value = (List<String>) myAttributes.get(FILTER_SITES_ATTRIBUTES.siteNames);
        if (value == null) {
            return null;
        } else {
            return Collections.unmodifiableList(value);
        }
    }

    public String chrPosFile() {
        return (String) myAttributes.get(FILTER_SITES_ATTRIBUTES.chrPosFile);
    }

    public PositionList positionList() {
        return (PositionList) myAttributes.get(FILTER_SITES_ATTRIBUTES.positionList);
    }

    public String bedFile() {
        return (String) myAttributes.get(FILTER_SITES_ATTRIBUTES.bedFile);
    }

}
