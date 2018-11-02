/**
 *
 */
package net.maizegenetics.dna.map;

import java.io.BufferedReader;
import net.maizegenetics.util.Utils;

/**
 * Utilities for reading and writing Position Lists
 *
 * @author lcj34
 *
 */
public class PositionListIOUtils {

    private PositionListIOUtils() {
    }

    /**
     * Returns a PositionList from a tab-delimited text SNP Conserve file. The
     * input file has 2 tab-delimited fields indicating Chromosome Number and
     * Position A header row is the first line and looks like this: #CHROM	POS
     * The remaining rows contains integer values as below: 9	18234
     *
     * @param fileName with complete path
     * @return PositionList
     */
    public static PositionList readSNPConserveFile(String fileName) {
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            PositionListBuilder plb = new PositionListBuilder();
            //parse SNP position rows
            // First value is Chromosome number, second is position
            String line = fileIn.readLine(); // read/skip header
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t");
                if (tokens.length != 2) {
                    System.err.println("Error in SNP Conserve File format:" + fileName);
                    System.err.println("Expecting tab-delimited file with 2 integer values per row");
                }
                Chromosome chrom = new Chromosome(tokens[0]);
                int pos = Integer.parseInt(tokens[1]);
                Position position = new GeneralPosition.Builder(chrom, pos).build();
                plb.add(position);
            }
            return plb.build();
        } catch (Exception e) {
            System.err.println("Error in Reading SNP Conserve File:" + fileName);
            e.printStackTrace();
        }
        return null;
    }
    
    /**
     * Returns a PositionList from a tab-delimited text SNP Quality Score file. The
     * input file has 3 tab-delimited fields indicating Chromosome Number and Quality Score
     * Position A header row is the first line and looks like this: CHROM	POS	QUALITYSCORE
     * The remaining rows contains integer values as below: 9	18234	15.5
     * NOTE: the CHROM field is a string.
     *
     * @param fileName with complete path
     * @return PositionList
     */
    public static PositionList readQualityScoreFile(String fileName) {
    	if (fileName == null) return null;
        try {
            BufferedReader fileIn = Utils.getBufferedReader(fileName, 1000000);
            PositionListBuilder plb = new PositionListBuilder();
            //parse SNP position rows
            // First value is Chromosome string, second is position, third is qualityScore
            String line = fileIn.readLine(); // read/skip header
            while ((line = fileIn.readLine()) != null) {
                String[] tokens = line.split("\\t");
                if (tokens.length < 3) {
                    System.err.println("Error in SNP Position QualityScore file format:" + fileName);
                    System.err.println("Expecting tab-delimited file where first 3 values in each row are type: integer, integer, float  "
                    		+ " with header values CHROM POS QUALITYSCORE");
                }
                Chromosome chrom = new Chromosome(tokens[0]);
                int pos = Integer.parseInt(tokens[1]);
                double qscore = Double.parseDouble(tokens[2]);
                Position position = new GeneralPosition.Builder(chrom, pos)
                		.addAnno("QualityScore",qscore).build();
                plb.add(position);
            }
            return plb.build();
        } catch (Exception e) {
            System.err.println("Error in Reading Quality Score File:" + fileName);
            e.printStackTrace();
        }
        return null;   	
    }

}
