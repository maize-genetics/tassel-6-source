/*
 *  FilterJSONUtils
 * 
 *  Created on Apr 2, 2015
 */
package net.maizegenetics.dna.snp.io;

import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.Map;
import javax.json.Json;
import javax.json.stream.JsonGenerator;
import javax.json.stream.JsonGeneratorFactory;
import net.maizegenetics.dna.snp.Filter;
import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.FilterSite;
import net.maizegenetics.dna.snp.FilterTaxa;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class FilterJSONUtils {

    private static final Logger myLogger = Logger.getLogger(FilterJSONUtils.class);

    private FilterJSONUtils() {
        // utility
    }

    public static String exportFilterToJSON(FilterList filters, String filename) {
        filename = Utils.addGzSuffixIfNeeded(filename, ".json");
        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            Map<String, Object> properties = new HashMap<>(1);
            properties.put(JsonGenerator.PRETTY_PRINTING, true);
            JsonGeneratorFactory factory = Json.createGeneratorFactory(properties);
            try (JsonGenerator generator = factory.createGenerator(writer)) {
                generator.writeStartObject();
                generator.writeStartArray("Filter");
                for (Filter filter : filters) {
                    if (filter instanceof FilterSite) {
                        filterToJSON((FilterSite) filter, generator);
                    } else if (filter instanceof FilterTaxa) {
                        filterToJSON((FilterTaxa) filter, generator);
                    }
                }
                generator.writeEnd();
                generator.writeEnd();
            }
            return filename;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("FilterJSONUtils: exportFilterToJSON: problem saving file: " + filename + "\n" + e.getMessage());
        }
    }

    private static void filterToJSON(FilterSite filter, JsonGenerator generator) {
        generator.writeStartObject();
        generator.writeStartObject("FilterSite");
        Map<FilterSite.FILTER_SITES_ATTRIBUTES, Object> attributes = filter.attributes();
        for (FilterSite.FILTER_SITES_ATTRIBUTES current : attributes.keySet()) {
            generator.write(current.name(), attributes.get(current).toString());
        }
        generator.writeEnd();
        generator.writeEnd();
    }

    private static void filterToJSON(FilterTaxa filter, JsonGenerator generator) {
        generator.writeStartObject();
        generator.writeStartObject("FilterTaxa");
        Map<FilterTaxa.FILTER_TAXA_ATTRIBUTES, Object> attributes = filter.attributes();
        for (FilterTaxa.FILTER_TAXA_ATTRIBUTES current : attributes.keySet()) {
            generator.write(current.name(), attributes.get(current).toString());
        }
        generator.writeEnd();
        generator.writeEnd();
    }

    public static FilterList importJSONToFilter(String filename) {
        return null;
    }

}
