package net.maizegenetics.dna.snp.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Optional;
import java.util.regex.Pattern;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.Utils;

public class ReadNumericMarkerUtils {

    //prevents instantiation
    private ReadNumericMarkerUtils() {
    }

    /**
     * @param inputFile	the input file with TASSEL v3 annotations or with no
     * input directives
     * @return
     * @throws IOException
     */
    public static GenotypeTable readNumericMarkerFile(String inputFile) throws IOException {

        BufferedReader br = Utils.getBufferedReader(inputFile);
        String inputline = br.readLine();
        Pattern sep = Pattern.compile("\\s+");

        String[] markerName = null;
        int numberOfColumns = 0;

        //process header rows and count the non-blank rows
        int numberOfDataLines = 0;
        while (inputline != null) {
            inputline = inputline.trim();
            String[] parsedline = sep.split(inputline);
            if (parsedline.length > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
                numberOfDataLines++;
            } else if (parsedline[0].toUpperCase().equals("<MARKER>")) {
                markerName = processHeader(numberOfColumns, parsedline, inputFile);
                numberOfColumns = markerName.length;
            }
            inputline = br.readLine();
        }
        br.close();

        if (numberOfDataLines == 0) {
            StringBuilder msg = new StringBuilder("Error in ");
            msg.append(inputFile);
            msg.append(": Missing taxa values.");
            throw new IllegalArgumentException(msg.toString());
        }
        if (numberOfColumns == 0) {
            StringBuilder msg = new StringBuilder("Error in ");
            msg.append(inputFile);
            msg.append(": Missing taxa data values.");
            throw new IllegalArgumentException(msg.toString());
        }

        //process body of data:  we needed numberOfColumns and numberOfDataLines from above
        String[][] textdata = new String[numberOfColumns][numberOfDataLines];
        String[] taxanames = new String[numberOfDataLines];
        br = Utils.getBufferedReader(inputFile);
        inputline = br.readLine();
        int totallines = 0;
        int linecount = 0;
        while (inputline != null) {
            totallines++;
            inputline = inputline.trim();
            String[] parsedline = sep.split(inputline);
            if (parsedline.length > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
                if (parsedline.length != numberOfColumns + 1) {
                    StringBuilder msg = new StringBuilder("Error in ");
                    msg.append(inputFile);
                    msg.append(" line ").append(totallines);
                    msg.append(": Incorrect number of data values.");
                    throw new IllegalArgumentException(msg.toString());
                }
                taxanames[linecount] = parsedline[0];
                for (int c = 0; c < numberOfColumns; c++) {
                    textdata[c][linecount] = parsedline[c + 1];
                }
                linecount++;
            }
            inputline = br.readLine();
        }
        br.close();

        TaxaList tL = new TaxaListBuilder().addAll(taxanames).build();
        ReferenceProbabilityBuilder rpb = ReferenceProbabilityBuilder.getInstance(numberOfDataLines,
                numberOfColumns, tL);

        //Create a list of values for each taxon.
        for (int indexR = 0; indexR < numberOfDataLines; indexR++) {
            float[] fvalues = new float[numberOfColumns];
            for (int indexC = 0; indexC < numberOfColumns; indexC++) {
                // Create array of floats - these are values for each taxon
                // Note when we read a line at a time when processing the file, the 
                // text array was defined as textData[columns][rows], sa access it thus
                if (textdata[indexC][indexR].equalsIgnoreCase("NaN")
                        || textdata[indexC][indexR].equalsIgnoreCase("NA")
                        || textdata[indexC][indexR].equals(".")) {
                    fvalues[indexC] = Float.NaN;
                } else {
                    try {
                        fvalues[indexC] = Float.parseFloat(textdata[indexC][indexR]);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("ReadNumericMarkerUtils: readNumericMarkerFile: Can't convert: " + textdata[indexC][indexR] + " to a number on data line: " + indexR);
                    }
                }
            }
            rpb.addTaxon(indexR, fvalues); //  taxon is the row.
        }

        ReferenceProbability rp = rpb.build();  // build does the "new"

        // Build PositionList for GenotypeTable
        PositionListBuilder posBuilder = new PositionListBuilder();
        for (int mNum = 0; mNum < numberOfColumns; mNum++) {
            String snpname = markerName[mNum];
            Optional<Object[]> info = decodeSnpName(snpname);
            if (info.isPresent()) {
                posBuilder.add(new GeneralPosition.Builder((Chromosome) info.get()[0], (Integer) info.get()[1]).snpName(snpname).build());
            } else {
                posBuilder.add(new GeneralPosition.Builder(Chromosome.UNKNOWN, mNum).snpName(snpname).build());
            }
        }

        PositionList pl = posBuilder.build();

        return GenotypeTableBuilder.getInstance(null, pl, tL, null, null, rp, null, null);
    }

    private static String[] processHeader(int numberOfColumns, String[] parsedline, String filename) {
        if (parsedline[0].equalsIgnoreCase("<Header")) {
            String headername = parsedline[1].split("[=>\\s]")[1];
            if (!parsedline[1].contains("name=") || headername.length() == 0) {
                StringBuilder msg = new StringBuilder("Error in ");
                msg.append(filename);
                msg.append(": Improperly formatted <Header name=> line.");
                throw new IllegalArgumentException(msg.toString());
            }
            int finalBracketPosition = 0;
            while (!parsedline[finalBracketPosition].contains(">")) {
                finalBracketPosition++;
            }
            if (numberOfColumns == 0) {
                numberOfColumns = parsedline.length - finalBracketPosition - 1;
            } else if (numberOfColumns != parsedline.length - finalBracketPosition - 1) {
                StringBuilder msg = new StringBuilder("Error in ");
                msg.append(filename);
                msg.append(": The number of ");
                msg.append(parsedline[0]).append(" ").append(parsedline[1]);
                msg.append(" columns does not match the number of columns in previous header rows");
                throw new IllegalArgumentException(msg.toString());
            }

            String[] contents = new String[numberOfColumns + 1];
            contents[0] = headername;
            for (int i = 0; i < numberOfColumns; i++) {
                contents[i + 1] = parsedline[i + finalBracketPosition + 1];
            }
            return contents;
        } else {
            if (numberOfColumns == 0) {
                numberOfColumns = parsedline.length - 1;
            } else if (numberOfColumns != parsedline.length - 1) {
                StringBuilder msg = new StringBuilder("Error in ");
                msg.append(filename);
                msg.append(": The number of ");
                msg.append(parsedline[0]);
                msg.append(" columns does not match the number of columns in previous header rows");
                throw new IllegalArgumentException(msg.toString());
            }
            String[] contents = new String[numberOfColumns];
            for (int i = 0; i < numberOfColumns; i++) {
                contents[i] = parsedline[i + 1];
            }
            return contents;
        }
    }
    
    public static Optional<Object[]> decodeSnpName(String name) {
        if (!name.startsWith("S")) return Optional.empty();
        try {
            name = name.substring(1);
            int ndx = name.indexOf('_');
            return Optional.of(new Object[]{new Chromosome(name.substring(0, ndx)) ,Integer.parseInt(name.substring(ndx + 1))});
        } catch (Exception e) {
            return Optional.empty();
        }
    }
}
