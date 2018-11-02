/*
 *  JSONUtils
 *
 *  Created on Mar 6, 2015
 */
package net.maizegenetics.dna.snp.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import javax.json.Json;
import javax.json.JsonArray;
import javax.json.JsonObject;
import javax.json.JsonReader;
import javax.json.JsonString;
import javax.json.stream.JsonGenerator;
import javax.json.stream.JsonGeneratorFactory;
import javax.json.stream.JsonParser;
import javax.json.stream.JsonParser.Event;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 */
public class JSONUtils {

    private static final Logger myLogger = Logger.getLogger(JSONUtils.class);

    private JSONUtils() {
        // utility
    }

    public static void taxaListToJSON(TaxaList taxa, JsonGenerator generator) {
        generator.writeStartArray("TaxaList");
        taxa.stream().forEach((current) -> {
            taxaToJSON(current, generator);
        });
        generator.writeEnd();
    }

    private static void taxaToJSON(Taxon taxon, JsonGenerator generator) {
        generator.writeStartObject();
        generator.write("name", taxon.getName());
        generalAnnotationToJSON(taxon.getAnnotation(), generator);
        generator.writeEnd();
    }

    private static void generalAnnotationToJSON(GeneralAnnotation annotation, JsonGenerator generator) {
        if (annotation == null) {
            return;
        }
        Set<String> keys = annotation.getAnnotationKeys();
        if (keys.isEmpty()) {
            return;
        }
        generator.writeStartObject("anno");
        keys.stream().forEach((key) -> {
            String[] values = annotation.getTextAnnotation(key);
            if (values.length == 1) {
                generator.write(key, values[0]);
            } else {
                generator.writeStartArray(key);
                for (String current : values) {
                    generator.write(current);
                }
                generator.writeEnd();
            }
        });
        generator.writeEnd();
    }

    /**
     * Exports given taxa list to JSON file.
     *
     * @param taxa taxa list
     * @param filename file name (adds .json if needed)
     *
     * @return final filename
     */
    public static String exportTaxaListToJSON(TaxaList taxa, String filename) {
        filename = Utils.addGzSuffixIfNeeded(filename, ".json");
        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            Map<String, Object> properties = new HashMap<>(1);
            properties.put(JsonGenerator.PRETTY_PRINTING, true);
            JsonGeneratorFactory factory = Json.createGeneratorFactory(properties);
            try (JsonGenerator generator = factory.createGenerator(writer)) {
                generator.writeStartObject();
                taxaListToJSON(taxa, generator);
                generator.writeEnd();
            }
            return filename;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("JSONUtils: exportTaxaListToJSON: problem saving file: " + filename + "\n" + e.getMessage());
        }
    }

    /**
     * Imports taxa list from JSON file.
     *
     * @param filename filename
     *
     * @return taxa list
     */
    public static TaxaList importTaxaListFromJSON(String filename) {
        try (BufferedReader reader = Utils.getBufferedReader(filename)) {
            try (JsonReader jsonReader = Json.createReader(reader)) {
                JsonObject obj = jsonReader.readObject();
                JsonArray taxaList = obj.getJsonArray("TaxaList");
                if (taxaList == null) {
                    throw new IllegalArgumentException("JSONUtils: importTaxaListFromJSON: There is no TaxaList in this file: " + filename);
                }
                TaxaList result = taxaListFromJSON(taxaList);
                myLogger.info("JSONUtils: importTaxaListFromJSON: Imported: " + result.numberOfTaxa() + " taxa from: " + filename);
                return result;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("JSONUtils: importTaxaListFromJSON: problem reading file: " + filename + "\n" + e.getMessage());
        }
    }

    private static TaxaList taxaListFromJSON(JsonArray taxaList) {
        if (taxaList.isEmpty()) {
            return null;
        }
        TaxaListBuilder builder = new TaxaListBuilder();
        taxaList.forEach((current) -> {
            builder.add(taxonFromJSON((JsonObject) current));
        });
        return builder.build();
    }

    private static Taxon taxonFromJSON(JsonObject json) {
        String name = json.getString("name");
        if (name == null) {
            throw new IllegalStateException("JSONUtils: taxonFromJSON: All Taxa must have a name.");
        }

        JsonObject anno = json.getJsonObject("anno");
        if (anno != null) {
            return new Taxon(name, generalAnnotationFromJSON(anno));
        } else {
            return new Taxon(name);
        }
    }

    private static GeneralAnnotationStorage.Builder generalAnnotationBuilderFromJSON(JsonObject json) {
        GeneralAnnotationStorage.Builder builder = GeneralAnnotationStorage.getBuilder();
        if (json != null) {
            json.forEach((key, value) -> {
                if (value instanceof JsonArray) {
                    JsonArray jsonArray = (JsonArray) value;
                    jsonArray.forEach((str) -> {
                        builder.addAnnotation(key, ((JsonString) str).getString());
                    });
                } else if (value instanceof JsonString) {
                    builder.addAnnotation(key, ((JsonString) value).getString());
                } else {
                    throw new IllegalArgumentException("JSONUtils: generalAnnotationBuilderFromJSON: unknown value type: " + value.getClass().getName());
                }
            });
        }
        return builder;
    }

    private static GeneralAnnotation generalAnnotationFromJSON(JsonObject json) {
        if (json == null) {
            return GeneralAnnotationStorage.EMPTY_ANNOTATION_STORAGE;
        }
        GeneralAnnotationStorage.Builder builder = generalAnnotationBuilderFromJSON(json);
        return builder.build();
    }

    public static String exportPositionListToJSON(PositionList positions, String filename) {
        filename = Utils.addGzSuffixIfNeeded(filename, ".json");
        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {
            Map<String, Object> properties = new HashMap<>(1);
            properties.put(JsonGenerator.PRETTY_PRINTING, true);
            JsonGeneratorFactory factory = Json.createGeneratorFactory(properties);
            try (JsonGenerator generator = factory.createGenerator(writer)) {
                generator.writeStartObject();
                positionListToJSON(positions, generator);
                generator.writeEnd();
            }
            return filename;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("JSONUtils: exportTaxaListToJSON: problem saving file: " + filename + "\n" + e.getMessage());
        }
    }

    private static void positionListToJSON(PositionList positions, JsonGenerator generator) {
        generator.writeStartArray("PositionList");
        positions.stream().forEach((current) -> {
            positionToJSON(current, generator);
        });
        generator.writeEnd();
    }

    private static void positionToJSON(Position position, JsonGenerator generator) {
        generator.writeStartObject();

        String snpID = position.getActualSNPID();
        if (snpID != null) {
            generator.write("SNPID", snpID);
        }

        chromosomeToJSON(position.getChromosome(), generator);
        generator.write("position", position.getPosition());

        short subPos = position.getInsertionPosition();
        if (subPos != 0) {
            generator.write("subPos", subPos);
        }

        float globalMAF = position.getGlobalMAF();
        if (!Float.isNaN(globalMAF)) {
            generator.write("globalMAF", globalMAF);
        }

        float globalSiteCoverage = position.getGlobalSiteCoverage();
        if (!Float.isNaN(globalSiteCoverage)) {
            generator.write("globalSiteCoverage", globalSiteCoverage);
        }

        generator.write("strand", Position.getStrand(position.getStrand()));

        generalAnnotationToJSON(position.getAnnotation(), generator);

        generator.writeEnd();
    }

    public static void chromosomeToJSON(Chromosome chromosome, JsonGenerator generator) {
        if (chromosome == null) {
            return;
        }
        generator.writeStartObject("chr");
        generator.write("name", chromosome.getName());
        int length = chromosome.getLength();
        if (length != -1) {
            generator.write("length", length);
        }
        generalAnnotationToJSON(chromosome.getAnnotation(), generator);
        generator.writeEnd();
    }

    public static PositionList importPositionListFromJSON(String filename) {
        try (InputStream reader = Utils.getInputStream(filename)) {
            try (JsonParser parser = Json.createParser(reader)) {
                Event temp = parser.next();
                if (temp != Event.START_OBJECT) {
                    throw new IllegalStateException("JSONUtils: importPositionListFromJSON: Position List must start with JSON Object: " + filename);
                }
                PositionList result = positionListFromJSON(parser);
                myLogger.info("JSONUtils: importPositionListFromJSON: Imported: " + result.numberOfSites() + " positions from: " + filename);
                return result;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("JSONUtils: importPositionListFromJSON: problem reading file: " + filename + "\n" + e.getMessage());
        }
    }

    private static PositionList positionListFromJSON(JsonParser parser) {
        Event current = parser.next();
        if ((current != Event.KEY_NAME) || (!parser.getString().equals("PositionList"))) {
            throw new IllegalStateException("JSOnUtils: positionListFromJSON: Expecting KEY_NAME PositionList");
        }
        current = parser.next();
        if (current != Event.START_ARRAY) {
            throw new IllegalStateException("JSONUtils: positionListFromJSON: Position List must start with JSON Array");
        }
        PositionListBuilder builder = new PositionListBuilder();
        current = parser.next();
        while (current == Event.START_OBJECT) {
            builder.add(positionFromJSON(parser));
            current = parser.next();
        }
        return builder.build();
    }

    private static Position positionFromJSON(JsonParser parser) {

        Map<String, Object> values = new HashMap<>();
        Event current = parser.next();
        while (current != Event.END_OBJECT) {
            if (current == Event.KEY_NAME) {
                String key = parser.getString();
                current = parser.next();
                if (current == Event.START_OBJECT) {
                    if (key.equals("chr")) {
                        values.put(key, chromosomeFromJSON(parser));
                    } else if (key.equals("anno")) {
                        values.put(key, generalAnnotationBuilderFromJSON(parser));
                    } else {
                        throw new IllegalStateException("JSONUtils: positionFromJSON: Unknown Object value key: " + key);
                    }
                } else {
                    values.put(key, parser.getString());
                }
            } else {
                throw new IllegalStateException("JSONUtils: positionFromJSON: Don't know how to handle Event type: " + current.name());
            }
            current = parser.next();
        }

        Chromosome chr = (Chromosome) values.remove("chr");
        if (chr == null) {
            throw new IllegalStateException("JSONUtils: positionFromJSON: No chromosome defined.");
        }
        Object positionStr = values.remove("position");
        if (positionStr == null) {
            throw new IllegalStateException("JSONUtils: positionFromJSON: No position defined.");
        }
        int position = Integer.parseInt((String) positionStr);
        GeneralPosition.Builder builder;
        GeneralAnnotationStorage.Builder anno = (GeneralAnnotationStorage.Builder) values.remove("anno");
        if (anno == null) {
            builder = new GeneralPosition.Builder(chr, position);
        } else {
            builder = new GeneralPosition.Builder(chr, position, anno);
        }

        for (Map.Entry<String, Object> entry : values.entrySet()) {
            String key = entry.getKey();
            Object value = entry.getValue();
            if (key.equals("SNPID")) {
                builder.snpName((String) value);
            } else if (key.equals("subPos")) {
                builder.insertionPosition(Short.parseShort((String) value));
            } else if (key.equals("globalMAF")) {
                builder.maf(Float.parseFloat((String) value));
            } else if (key.equals("globalSiteCoverage")) {
                builder.siteCoverage(Float.parseFloat((String) value));
            } else if (key.equals("strand")) {
                builder.strand((String) value);
            } else {
                throw new IllegalStateException("JSONUtils: positionFromJSON: Unknown key: " + key);
            }
        }

        return builder.build();

    }

    private static GeneralAnnotationStorage.Builder generalAnnotationBuilderFromJSON(JsonParser parser) {

        GeneralAnnotationStorage.Builder builder = GeneralAnnotationStorage.getBuilder();
        Event current = parser.next();
        while (current != Event.END_OBJECT) {
            if (current == Event.START_ARRAY) {
                String key = parser.getString();
                current = parser.next();
                while (current != Event.END_ARRAY) {
                    builder.addAnnotation(key, parser.getString());
                    current = parser.next();
                }
            } else if (current == Event.KEY_NAME) {
                String key = parser.getString();
                parser.next();
                builder.addAnnotation(key, parser.getString());
            } else {

            }
            current = parser.next();
        }

        return builder;

    }

    private static GeneralAnnotation generalAnnotationFromJSON(JsonParser parser) {
        if (!parser.getString().equals("anno")) {
            return GeneralAnnotationStorage.EMPTY_ANNOTATION_STORAGE;
        }
        GeneralAnnotationStorage.Builder builder = generalAnnotationBuilderFromJSON(parser);
        return builder.build();
    }

    public static Chromosome chromosomeFromJSON(JsonParser parser) {
        String name = null;
        int length = -1;
        GeneralAnnotation anno = null;
        Event current = parser.next();
        while (current != Event.END_OBJECT) {
            if (current == Event.KEY_NAME) {
                String key = parser.getString();
                parser.next();
                if (key.equals("name")) {
                    name = parser.getString();
                } else if (key.equals("length")) {
                    length = parser.getInt();
                } else {
                    throw new IllegalStateException("JSONUtils: chromosomeFromJSON: Unknown key: " + key);
                }
            } else if (current == Event.START_OBJECT) {
                if (parser.getString().equals("anno")) {
                    anno = generalAnnotationFromJSON(parser);
                }
            }
            current = parser.next();
        }
        if (length == -1 && anno == null) {
            return Chromosome.instance(name);
        }
        return new Chromosome(name, length, anno);
    }

}
