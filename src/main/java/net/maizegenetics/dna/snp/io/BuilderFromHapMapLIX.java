/*
 *  BuilderFromHapMapLIX
 * 
 *  Created on Aug 29, 2015
 */
package net.maizegenetics.dna.snp.io;

import com.google.common.collect.SetMultimap;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.tribble.util.ParsingUtils;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.LineIndexHapmapGenotypeCallTable;
import static net.maizegenetics.dna.snp.io.BuilderFromHapMap.processTaxa;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class BuilderFromHapMapLIX {

    private static final Logger myLogger = Logger.getLogger(BuilderFromHapMapLIX.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private static final int SNPID_INDEX = 0;
    private static final int VARIANT_INDEX = 1;
    private static final int CHROMOSOME_INDEX = 2;
    private static final int POSITION_INDEX = 3;

    private BuilderFromHapMapLIX() {
    }

    public static GenotypeTable build(String hapmapFileBGZip) {
        return build(hapmapFileBGZip, ParsingUtils.appendToPath(hapmapFileBGZip, LineIndexBuilder.LINE_INDEX_FILE_EXTENSION));
    }

    public static GenotypeTable build(String hapmapFileBGZip, String indexFilename) {

        Tuple<LineIndex, String[]> indexPositionInfo = LineIndexBuilder.readIndex(indexFilename);

        Map<String, Chromosome> chromosomeLookup = new HashMap<>();
        PositionListBuilder positions = new PositionListBuilder();
        for (String current : indexPositionInfo.y) {

            String[] tokens = current.split("\t");

            String chrName = tokens[CHROMOSOME_INDEX];
            Chromosome currChr = chromosomeLookup.get(chrName);
            if (currChr == null) {
                currChr = new Chromosome(new String(chrName));
                chromosomeLookup.put(chrName, currChr);
            }

            String variants = tokens[VARIANT_INDEX];

            int physicalPos;
            try {
                physicalPos = Integer.parseInt(tokens[POSITION_INDEX]);
            } catch (Exception ex) {
                throw new IllegalArgumentException("BuilderFromHapMapLIX: Position must be an integer: " + tokens[POSITION_INDEX]);
            }

            GeneralPosition.Builder positionBuilder = new GeneralPosition.Builder(currChr, physicalPos)
                    .snpName(tokens[SNPID_INDEX])
                    .knownVariants(variants);

            byte glbMajor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(0));
            positionBuilder.allele(WHICH_ALLELE.GlobalMajor, glbMajor);
            if (variants.length() == 3) {
                byte glbMinor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(2));
                positionBuilder.allele(WHICH_ALLELE.GlobalMinor, glbMinor);
            }

            positions.add(positionBuilder.build());
        }
        PositionList positionList = positions.build();

        TaxaList taxaList = null;
        boolean isOneLetter = false;
        try (BlockCompressedInputStream reader = new BlockCompressedInputStream(new File(hapmapFileBGZip))) {

            Map<String, SetMultimap<String, String>> sampAnnoBuild = new TreeMap<>();

            String currLine = reader.readLine();
            while ((currLine != null) && currLine.startsWith("##")) {

                String[] cat = currLine.split("=", 2);
                if (cat.length < 2) {
                    continue;
                }
                if (cat[0].startsWith("##SAMPLE")) {

                    SetMultimap<String, String> mapOfAnno = TaxaListIOUtils.parseVCFHeadersIntoMap(cat[1]);
                    String taxaID = mapOfAnno.get("ID").iterator().next();
                    if (taxaID != null) {
                        sampAnnoBuild.put(taxaID, mapOfAnno);
                    }

                }

                currLine = reader.readLine();
            }

            taxaList = processTaxa(currLine, sampAnnoBuild).build();
            int numTaxa = taxaList.numberOfTaxa();

            currLine = reader.readLine();
            String[] tokens = WHITESPACE_PATTERN.split(currLine, NUM_HAPMAP_NON_TAXA_HEADERS + 1);
            if (tokens.length <= NUM_HAPMAP_NON_TAXA_HEADERS) {
                throw new IllegalStateException("BuilderFromHapMapLIX: Header Incorrectly Formatted: See:\nhttps://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-hapmap");
            }
            double avg = (double) (tokens[NUM_HAPMAP_NON_TAXA_HEADERS].length() + 1) / (double) numTaxa;
            if ((avg > 1.99) && (avg < 2.01)) {
                isOneLetter = true;
            } else if ((avg > 2.99) && (avg < 3.01)) {
                isOneLetter = false;
            } else {
                throw new IllegalStateException("BuilderFromHapMapLIX: Genotype coded wrong use 1 or 2 letters per genotype.  Or first site has incorrect number of values.");
            }

        } catch (Exception e) {
            throw new IllegalStateException("BuilderFromHapMapLIX: Problem opening file: " + hapmapFileBGZip + "\n" + e.getMessage());
        }

        return GenotypeTableBuilder.getInstance(LineIndexHapmapGenotypeCallTable.getInstance(taxaList.numberOfTaxa(), positionList.numberOfSites(), false, isOneLetter, indexPositionInfo.x, hapmapFileBGZip), positionList, taxaList);

    }

}
