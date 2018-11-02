/*
 *  SiteScoresIO
 * 
 *  Created on Feb 27, 2015
 */
package net.maizegenetics.dna.snp.io;

import java.io.BufferedWriter;

import java.text.NumberFormat;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class SiteScoresIO {

    private static final Logger myLogger = Logger.getLogger(SiteScoresIO.class);

    private static final String DELIMITER = "\t";

    private static final NumberFormat DECIMAL_FORMAT = NumberFormat.getNumberInstance();

    static {
        DECIMAL_FORMAT.setMaximumFractionDigits(3);
    }

    private SiteScoresIO() {
        // utility
    }

    public static String writeReferenceProbability(GenotypeTable genotypeTable, String filename) {

        filename = Utils.addSuffixIfNeeded(filename, ".txt");

        ReferenceProbability probability = genotypeTable.referenceProbability();
        if (probability == null) {
            throw new IllegalStateException("SiteScoresIO: writeReferenceProbability: this genotype table has no reference probability.");
        }

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {

            writer.write("<Numeric>\n");

            writer.write("<Marker>");
            PositionList positions = genotypeTable.positions();
            for (Position current : positions) {
                writer.write(DELIMITER);
                writer.write(current.getSNPID());
            }
            writer.write("\n");

            TaxaList taxa = genotypeTable.taxa();
            for (int r = 0, n = probability.numTaxa(); r < n; r++) {
                writer.write(taxa.get(r).getName());
                for (int i = 0; i < probability.numSites(); i++) {
                    writer.write(DELIMITER);
                    float value = probability.value(r, i);
                    if (Float.isNaN(value)) {
                        writer.write("NA");
                    } else {
                        writer.write(DECIMAL_FORMAT.format(probability.value(r, i)));
                    }
                }
                writer.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("SiteScoresIO: writeReferenceProbability: problem saving file: " + filename);
        }

        return filename;

    }

    public static String writeDepth(GenotypeTable genotypeTable, String filename) {

        filename = Utils.addSuffixIfNeeded(filename, ".txt");

        AlleleDepth depth = genotypeTable.depth();
        if (depth == null) {
            throw new IllegalStateException("SiteScoresIO: writeDepth: this genotype table has no depth.");
        }

        try (BufferedWriter writer = Utils.getBufferedWriter(filename)) {

            writer.write("Taxa");
            PositionList positions = genotypeTable.positions();
            for (Position current : positions) {
                writer.write(DELIMITER);
                writer.write(current.getSNPID());
            }
            writer.write("\n");

            TaxaList taxa = genotypeTable.taxa();
            for (int t = 0, n = depth.numTaxa(); t < n; t++) {
                writer.write(taxa.get(t).getName());
                for (int s = 0, m = depth.numSites(); s < m; s++) {
                    byte genotype = genotypeTable.genotype(t, s);
                    writer.write(DELIMITER);
                    int value = 0;
                    byte allele1 = (byte) ((genotype >>> 4) & 0xf);
                    if (allele1 < 6) {
                        value = depth.depthForAllele(t, s, allele1);
                    }
                    byte allele2 = (byte) (genotype & 0xf);
                    if ((allele2 < 6) && (allele2 != allele1)) {
                        value += depth.depthForAllele(t, s, allele2);
                    }
                    writer.write(String.valueOf(value));
                }
                writer.write("\n");
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("SiteScoresIO: writeDepth: problem saving file: " + filename);
        }

        return filename;

    }
}
