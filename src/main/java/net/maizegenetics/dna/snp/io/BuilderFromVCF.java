package net.maizegenetics.dna.snp.io;

import com.google.common.base.Splitter;
import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.score.AlleleDepthBuilder;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

/**
 * Create an alignment based on VCF format file (either .txt or compressed).  Alleles are set as global reference.
 * e.g. code <p></p>
 * {@code
 * Alignment a=BuilderFromVCF.getBuilder(infileName).build();
 * }
 * <p></p>
 * TODO:  Add filtering while reading, provide an option to define the alleles as reference and alternate
 *
 * @author Ed Buckler
 */
public class BuilderFromVCF {

    private static final Logger myLogger = LogManager.getLogger(BuilderFromVCF.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("[\\s]+");
    private static final Pattern TAB_PATTERN = Pattern.compile("[\\t]+");
    private HeaderPositions hp = null;
    private final String infile;
    private boolean includeDepth = false;
    private final ProgressListener myProgressListener;

    private BuilderFromVCF(String infile, ProgressListener listener) {
        this.infile = infile;
        this.myProgressListener = listener;
    }


    /**
     * Create a builder for loading a VCF file into memory
     *
     * @param infile name of the VCF file to be load
     *
     * @return a builder
     */
    public static BuilderFromVCF getBuilder(String infile) {
        return new BuilderFromVCF(infile, null);
    }

    /**
     * Create a builder for loading a VCF file into memory
     *
     * @param infile name of the VCF file to be load
     * @param listener Progress listener for GUI
     *
     * @return a builder
     */
    public static BuilderFromVCF getBuilder(String infile, ProgressListener listener) {
        return new BuilderFromVCF(infile, listener);
    }

    public BuilderFromVCF keepDepth() {
        includeDepth = true;
        return this;
    }

    public GenotypeTable buildAndSortInMemory() {
        return buildEngine(true);
    }

    public GenotypeTable build() {
        return buildEngine(false);
    }

    //TODO provide options on caching to use, read only some sites, etc.
    private GenotypeTable buildEngine(boolean fullSort) {
        long time = System.nanoTime();
        GenotypeTable result = null;
        int totalSites = -1;//unknown
        GenotypeTableBuilder gtbDiskBuild = null;
        ExecutorService pool = null;
        try {

            int numThreads = Runtime.getRuntime().availableProcessors();
            pool = Executors.newFixedThreadPool(numThreads);

            BufferedReader r = Utils.getBufferedReader(infile, -1);
            //Read the ## annotation rows
            String currLine;
            Map<String, String> infoMap = new HashMap<>();
            Map<String, String> formatMap = new HashMap<>();
            Map<String, SetMultimap<String, String>> sampAnnoBuild = new TreeMap<>();
            currLine = parseVCFHeadersIntoMaps(infoMap, formatMap, sampAnnoBuild, r);

            TaxaList taxaList = processTaxa(currLine, sampAnnoBuild);
            int linesAtTime = 1 << 12;
            //  int linesAtTime=1<<8;  //better for with lots of taxa.
            ArrayList<String> txtLines = new ArrayList<>(linesAtTime);
            ArrayList<ProcessVCFBlock> pbs = new ArrayList<>();
            List<Future<ProcessVCFBlock>> futures = new ArrayList<>();
            int sitesRead = 0;
            while ((currLine = r.readLine()) != null) {
                if (currLine.startsWith("#")) continue;
                txtLines.add(currLine);
                sitesRead++;
                if (sitesRead % linesAtTime == 0) {
                    ProcessVCFBlock pb;
                    pb = ProcessVCFBlock.getInstance(taxaList.numberOfTaxa(), hp, txtLines, includeDepth);
                    //pbs.add(pb);
                    //     pb.run(); //used for testing
                    try {
                        futures.add(pool.submit(pb));
                    } catch (Exception e) {
                        myLogger.debug(e.getMessage(), e);
                        throw new IllegalStateException(e.getMessage());
                    }
                    txtLines = new ArrayList<>(linesAtTime);
                }
            }
            r.close();
            //Handle whatever is left over in the file
            if (txtLines.size() > 0) {
                ProcessVCFBlock pb;//=ProcessVCFBlock.getInstance(taxaList.numberOfTaxa(), hp, txtLines);
                pb = ProcessVCFBlock.getInstance(taxaList.numberOfTaxa(), hp, txtLines, includeDepth);
                //pbs.add(pb);
                //  pb.run(); //used for testing

                try {
                    futures.add(pool.submit(pb));
                } catch (Exception e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException(e.getMessage());
                }
            }
            int numFutures = futures.size();
            int count = 0;
            for (Future<ProcessVCFBlock> future : futures) {
                try {
                    pbs.add(future.get());
                } catch (Exception e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException(e.getMessage());
                }
                if (myProgressListener != null) {
                    count++;
                    myProgressListener.progress(count * 100 / numFutures, null);
                }
            }
            pool.shutdown();

            //result=completeInMemoryBuilding(futures, taxaList, sitesRead, includeDepth, fullSort);
            result = completeInMemoryBuilding(pbs, taxaList, sitesRead, includeDepth, fullSort);

//            int currentSite=0;
//            PositionListBuilder posBuild=new PositionListBuilder();
//            GenotypeCallTableBuilder gb=GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.numberOfTaxa(), lines);
//            AlleleDepthBuilder db=null;
//            if(includeDepth) db=AlleleDepthBuilder.getInstance(taxaList.numberOfTaxa(),lines,6);
//            for (ProcessVCFBlock pb : pbs) {
//                posBuild.addAll(pb.getBlkPosList());
//                byte[][] bgTS=pb.getGenoTS();
//                for (int t=0; t<bgTS.length; t++) {
//                    gb.setBaseRangeForTaxon(t, currentSite, bgTS[t]);
//                }
//                if(includeDepth) {
//                    byte[][][] bdTS=pb.getDepthTS();
//                    for (int t=0; t<bgTS.length; t++) {
//                        db.setDepthRangeForTaxon(t, currentSite, bdTS[t]);
//                    }
//                }
//                currentSite+=pb.getSiteNumber();
//            }
//            if (posBuild.validateOrdering()==false) {
//                throw new IllegalStateException("BuilderFromVCF: Ordering incorrect HapMap must be ordered by position");
//            }
//            GenotypeCallTable g=gb.build();
//            if(includeDepth) {result=GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList, null, db.build());}
//            else {result=GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList);}
        } catch (IOException e) {
            e.printStackTrace();
        } catch (IllegalStateException e) {
            if (pool != null) {
                pool.shutdown();
            }
            //e.printStackTrace();
            throw e;
        }
        long totalTime = System.nanoTime() - time;
        System.out.printf("BuilderFromVCF data timing %gs %n", totalTime / 1e9);
        return result;
    }

    private static GenotypeTable completeInMemoryBuilding(List<ProcessVCFBlock> pbs, TaxaList taxaList, int numberOfSites, boolean includeDepth, boolean fullSort) {
        int currentSite = 0;
        PositionListBuilder posBuild = new PositionListBuilder();
        GenotypeCallTableBuilder gb = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(taxaList.numberOfTaxa(), numberOfSites);
        AlleleDepthBuilder db = null;

        if (includeDepth) db = AlleleDepthBuilder.getInstance(taxaList.numberOfTaxa(), numberOfSites, taxaList);


        for (ProcessVCFBlock pb : pbs) {
            posBuild.addAll(pb.getBlkPosList());
            byte[][] bgTS = pb.getGenoTS();
            for (int t = 0; t < bgTS.length; t++) {
                gb.setBaseRangeForTaxon(t, currentSite, bgTS[t]);
            }
            if (includeDepth) {
                byte[][][] bdTS = pb.getDepthTS();
                for (int t = 0; t < bgTS.length; t++) {
                    db.setDepthRangeForTaxon(t, currentSite, bdTS[t]);
                }

            }
            currentSite += pb.getSiteNumber();
        }


        //Check that result is in correct order. If not, either try to sort or just throw an error (determined by what was passed to fullSort)
        if (posBuild.validateOrdering() == false) {
            if (fullSort) {
                int[] siteRedirect = posBuild.sort();
                gb.reorderPositions(siteRedirect);
                if (includeDepth) {
                    db.reorderPositions(siteRedirect);
                }
                if (posBuild.validateOrdering() == false) {   //Double-check post-sort ordering. Should never happen, but just to be safe
                    throw new IllegalStateException("BuilderFromVCF: Ordering of VCF file held in memory failed.");
                }
            } else {
                throw new IllegalStateException("BuilderFromVCF: Ordering incorrect. VCF file must be ordered by position. Please first use SortGenotypeFilePlugin to correctly order the file.");
            }
        }
        GenotypeCallTable g = gb.build();
        if (includeDepth) {
            return GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList, db.build());
        } else {
            return GenotypeTableBuilder.getInstance(g, posBuild.build(), taxaList);
        }

    }

    private static String parseVCFHeadersIntoMaps(Map<String, String> infoMap, Map<String, String> formatMap,
                                                  Map<String, SetMultimap<String, String>> sampAnnoBuild, BufferedReader r) throws IOException {
        String currLine;
        while (((currLine = r.readLine()) != null) && (currLine.startsWith("##"))) {
            String[] cat = currLine.split("=", 2);
            if (cat.length < 2) continue;
            switch (cat[0]) {
//                    case "##INFO":
//                        infoMap.put(mapOfAnno.get("ID"), mapOfAnno.get("Description"));
//                        break;
//                    case "##FILTER":break;
//                    case "##FORMAT":
//                        formatMap.put(mapOfAnno.get("ID"),mapOfAnno.get("Description"));
//                        break;
                case "##SAMPLE":
                    SetMultimap<String, String> mapOfAnno = TaxaListIOUtils.parseVCFHeadersIntoMap(cat[1]);
                    String taxaID = mapOfAnno.get("ID").iterator().next();
                    if (taxaID == null) break;
                    sampAnnoBuild.put(taxaID, mapOfAnno);
                    break;
                case "##PEDIGREE":
                    break;
                default:
                    break;
            }

//                System.out.println(currLine);
        }
        return currLine;
    }

    private static String getReplaceCommaWithinQuote(String s) {
        StringBuilder sb = new StringBuilder(s);
        boolean inQuote = false;
        for (int i = 0; i < sb.length(); i++) {
            if (sb.charAt(i) == '\"') inQuote = (!inQuote);
            if (inQuote && sb.charAt(i) == ',') sb.setCharAt(i, (char) ((int) ',' + 256));//(char)167);
            if (inQuote && sb.charAt(i) == '=') sb.setCharAt(i, (char) ((int) '=' + 256));//(char)167);
        }
        return sb.toString();
    }

    private TaxaList processTaxa(String readLn, Map<String, SetMultimap<String, String>> taxaAnnotation) {
        String[] header = TAB_PATTERN.split(readLn);
        hp = new HeaderPositions(header);
        int numTaxa = header.length - hp.NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (int i = 0; i < numTaxa; i++) {
            String taxonID = header[i + hp.NUM_HAPMAP_NON_TAXA_HEADERS];
            Taxon.Builder at = new Taxon.Builder(taxonID);
            SetMultimap<String, String> taMap = taxaAnnotation.get(taxonID);
            if (taMap != null) {
                for (Map.Entry<String, String> en : taMap.entries()) {
                    if (en.getKey().equals("ID")) continue; //skip the IDs as these became the name
                    String s = en.getValue().replace("\"", "");
                    at.addAnno(en.getKey(), s);
                }
            }
            tlb.add(at.build());
        }
        return tlb.build();
    }


}

class HeaderPositions {
    final int NUM_HAPMAP_NON_TAXA_HEADERS;
    final int GENOIDX;
    final int SNPID_INDEX;
    //  final int VARIANT_INDEX;
    final int FILTER_INDEX;
    final int QUAL_INDEX;
    final int CHROMOSOME_INDEX;
    final int POSITION_INDEX;
    final int REF_INDEX;
    final int ALT_INDEX;
    final int INFO_INDEX;
    final int FORMAT_INDEX;

    public HeaderPositions(String[] header) {
        int chrIdx = firstEqualIndex(header, "#CHROM");
        if (chrIdx < 0) chrIdx = firstEqualIndex(header, "#CHR");
        CHROMOSOME_INDEX = chrIdx;
        POSITION_INDEX = firstEqualIndex(header, "POS");
        SNPID_INDEX = firstEqualIndex(header, "ID");
        REF_INDEX = firstEqualIndex(header, "REF");
        ALT_INDEX = firstEqualIndex(header, "ALT");
        QUAL_INDEX = firstEqualIndex(header, "QUAL");
        FILTER_INDEX = firstEqualIndex(header, "FILTER");
        INFO_INDEX = firstEqualIndex(header, "INFO");
        FORMAT_INDEX = firstEqualIndex(header, "FORMAT");

        NUM_HAPMAP_NON_TAXA_HEADERS = Math.max(INFO_INDEX, FORMAT_INDEX) + 1;
        GENOIDX = NUM_HAPMAP_NON_TAXA_HEADERS;
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i = 0; i < sa.length; i++) {
            if (sa[i].equals(match)) return i;
        }
        return -1;
    }

}

//class ProcessVCFBlock implements Runnable {
class ProcessVCFBlock implements Callable<ProcessVCFBlock> {

    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    private static final Pattern SLASH_PATTERN = Pattern.compile("/");
    private final HeaderPositions hp;
    private final int taxaN;
    private final int siteN;
    private final int startSite; //if unknown Int.Mini
    private ArrayList<String> txtL;
    private byte[][] gTS;  //genotypes
    private byte[][][] dTS; //depth
    private final ArrayList<Position> blkPosList;
    private final boolean keepDepth;


    private ProcessVCFBlock(int taxaN, HeaderPositions hp, ArrayList<String> txtL, int startSite,
                            boolean keepDepth) {
        this.taxaN = taxaN;
        this.siteN = txtL.size();
        this.txtL = txtL;
        this.hp = hp;
        blkPosList = new ArrayList<>(siteN);
        this.startSite = startSite;
        this.keepDepth = keepDepth;
    }

    /*Used to process VCF blocks and return the result for a in memory GenotypeTable*/
    static ProcessVCFBlock getInstance(int taxaN, HeaderPositions hp, ArrayList<String> txtL, boolean keepDepth) {
        return new ProcessVCFBlock(taxaN, hp, txtL, Integer.MIN_VALUE, keepDepth);
    }

    @Override
    public ProcessVCFBlock call() throws Exception {
        Map<String, Chromosome> chromosomeLookup = new HashMap<>();
        gTS = new byte[taxaN][siteN];
        if (keepDepth == true) dTS = new byte[taxaN][6][siteN];
        for (int s = 0; s < siteN; s++) {
            //really needs to use a Splitter iterator to make this cleaner if it is performant
            String input = txtL.get(s);
            try {
                int[] tabPos = new int[hp.NUM_HAPMAP_NON_TAXA_HEADERS + taxaN];
                int tabIndex = 0;
                int len = input.length();
                for (int i = 0; (tabIndex < hp.NUM_HAPMAP_NON_TAXA_HEADERS + taxaN) && (i < len); i++) {
                    if (input.charAt(i) == '\t') {
                        tabPos[tabIndex++] = i;
                    }
                }
                String chrName = input.substring(0, tabPos[hp.CHROMOSOME_INDEX]);
                Chromosome currChr = chromosomeLookup.get(chrName);
                if (currChr == null) {
                    currChr = new Chromosome(new String(chrName));
                    chromosomeLookup.put(chrName, currChr);
                }
                String snpID = null;
                if (hp.SNPID_INDEX > 0) snpID = input.substring(tabPos[hp.SNPID_INDEX - 1] + 1, tabPos[hp.SNPID_INDEX]);
                String refS = input.substring(tabPos[hp.REF_INDEX - 1] + 1, tabPos[hp.REF_INDEX]);
                String alt = input.substring(tabPos[hp.ALT_INDEX - 1] + 1, tabPos[hp.ALT_INDEX]);
                String variants;
                if (alt.equals(".")) {
                    variants = refS;
                } else {
                    variants = (refS + "/" + alt).replace(',', '/')
                            .replace("<INS>", "+").replace('I', '+')
                            .replace("<DEL>", "-").replace('D', '-')
                            .replace("*", "N");
                }

                //GeneralPosition.Builder apb=new GeneralPosition.Builder(currChr, currentPosition)
                //                                               .knownVariants(variants); //TODO strand, variants,
                //ZRM 8_26
                GeneralPosition.Builder apb = new GeneralPosition.Builder(currChr, Integer.parseInt(input.substring(tabPos[hp.POSITION_INDEX - 1] + 1, tabPos[hp.POSITION_INDEX])))
                        .knownVariants(variants) //TODO strand, variants,
                        ;
                if (snpID != null && !snpID.equals(".")) {
                    apb.snpName(snpID);
                }
                //byte[] alleles=new byte[(variants.length()+1)/2];
                byte[] alleles = new byte[variants.split("/").length];
                for (int i = 0, varInd = 0; i < alleles.length; i++, varInd += 2) {
                    alleles[i] = NucleotideAlignmentConstants.getNucleotideAlleleByte(variants.charAt(varInd));
                }
                /***ZRM 8_27 New code ***/
                String[] variantList = variants.split("/");
                if (variantList[0].length() > 1) {
                    String[] parsedVariantList = new String[variantList.length];
                    //alt deletion
                    for (int i = 0; i < variantList.length; i++) {
                        //Pull off the first character if it exists
                        if (variantList[i].length() > 1) {
                            parsedVariantList[i] = variantList[i].substring(1);
                            if (parsedVariantList[i].length() == 0) {
                                parsedVariantList[i] = "-";
                            }
                        } else {
                            //Mark as deletion
                            parsedVariantList[i] = "-";
                        }
                    }
                    for (int i = 0; i < parsedVariantList.length; i++) {
                        alleles[i] = NucleotideAlignmentConstants.getNucleotideAlleleByte(parsedVariantList[i].charAt(0));
                    }
                } else {
                    //Check for reference deletion(insertion)
                    //Loop through all variants to see if one alt is longer than the ref
                    boolean isIndel = false;
                    for (int i = 1; i < variantList.length; i++) {
                        if (variantList[i].length() > variantList[0].length()) {
                            isIndel = true;
                            break;
                        }
                    }
                    if (isIndel) {
                        String[] parsedVariantList = new String[variantList.length];
                        //ref+alt deletion
                        for (int i = 0; i < variantList.length; i++) {
                            //Pull off the first character if it exists
                            if (variantList[i].length() > 1) {
                                parsedVariantList[i] = variantList[i].substring(1);
                                if (parsedVariantList[i].length() == 0) {
                                    parsedVariantList[i] = "-";
                                }
                            } else {
                                //Mark as deletion
                                parsedVariantList[i] = "-";
                            }
                        }
                        for (int i = 0; i < parsedVariantList.length; i++) {
                            alleles[i] = NucleotideAlignmentConstants.getNucleotideAlleleByte(parsedVariantList[i].charAt(0));
                        }
                    } else {
                        //if not just put it in the allele array
                        for (int i = 0; i < variantList.length; i++) {
                            alleles[i] = NucleotideAlignmentConstants.getNucleotideAlleleByte(variantList[i].charAt(0));
                        }
                    }
                }
                /***ZRM 8_27 New code end ***/
                apb.allele(WHICH_ALLELE.Reference, alleles[0]);
                if (alleles.length > 1) {
                    apb.allele(WHICH_ALLELE.Alternate, alleles[1]);
                }
                for (String annoS : Splitter.on(";").split(input.substring(tabPos[hp.INFO_INDEX - 1] + 1, tabPos[hp.INFO_INDEX]))) {
                    apb.addAnno(annoS);
                }
                blkPosList.add(apb.build());
                final int iGT = 0; //genotype index
                int iAD = -1, iDP = -1, iGQ = -1, iPL = -1;  //alleleDepth, overall depth, genotypeQuality, phredGenotypeLikelihoods
                if (hp.FORMAT_INDEX >= 0) {
                    //Check to see if FORMAT tag is missing. Only applicable for single taxa files
                    if (tabPos[hp.FORMAT_INDEX] == 0) {
                        throw new IllegalStateException("Error Processing VCF: Missing FORMAT tag.");
                    }
                    String unsplitInput = input.substring(tabPos[hp.FORMAT_INDEX - 1] + 1, tabPos[hp.FORMAT_INDEX]);
                    if (unsplitInput.length() == 0 || !unsplitInput.startsWith("GT")) {
                        //Check to see it has the GT field
                        if (unsplitInput.contains("GT")) {
                            throw new IllegalStateException("Error Processing VCF Block: GT field is not in first position of FORMAT.");
                        }
                        //If GT isnt in, we assume that it is missing FORMAT
                        else {
                            throw new IllegalStateException("Error Processing VCF Block: Missing FORMAT tag.");
                        }
                    }
                    String[] formatS = unsplitInput.split(":");

                    iAD = firstEqualIndex(formatS, "AD");
                }
                int t = 0;
                for (String taxaAllG : Splitter.on("\t").split(input.substring(tabPos[hp.NUM_HAPMAP_NON_TAXA_HEADERS - 1] + 1))) {
                    int f = 0;
                    //if(taxaAllG.equals(".")) {
                    if (taxaAllG.startsWith(".")) {
                        gTS[t][s] = GenotypeTable.UNKNOWN_GENOTYPE;
                        //need to still move up the taxa counter by one as we have covered this now.
                        t++;
                        continue;
                    }
                    for (String fieldS : Splitter.on(":").split(taxaAllG)) {
                        if (f == iGT) {
                            //String "[.0-9]\\/[.0-9]||[.0-9]\\|[.0-9]" will match a valid diploid
                            if (!fieldS.equals(".")) { //[TAS-509] Check to make sure we are using diploids in the form 0/1 or 0|0
                                int a1 = 0;
                                int a2 = 0;
                                //Should be this, but for speed we just check the character
                                //if(fieldS.split("|/").length ==1) {
                                if (fieldS.length() == 1) {
                                    a1 = fieldS.charAt(0) - '0';
                                    a2 = fieldS.charAt(0) - '0';
                                } else {
                                    //zrm22 Jul 31, 2017
                                    //int a1 = fieldS.charAt(0) - '0';
                                    //int a2 = fieldS.charAt(2) - '0';

                                    a1 = fieldS.charAt(0) - '0';
                                    a2 = fieldS.charAt(2) - '0';
                                }
                                if (a1 > alleles.length - 1 || a2 > alleles.length - 1) {
                                    Position pos = blkPosList.get(blkPosList.size() - 1);
                                    throw new IllegalStateException("\nError Processing VCF block: Mismatch of alleles.\n  At Chromosome " + pos.getChromosome().getName() + ", Position " + pos.getPosition() + ".\nAllele ID larger than number of alleles");
                                }
                                if (a1 < 0 || a2 < 0) {
                                    gTS[t][s] = GenotypeTable.UNKNOWN_GENOTYPE;
                                } else {
                                    gTS[t][s] =
                                            GenotypeTableUtils.getDiploidValue(alleles[a1], alleles[a2]);
                                }
                            } else {    //[TAS-509] if it isnt a diploid error out early
                                throw new IllegalStateException("Error Processing VCF block: Found haploid information for the element: "
                                        + taxaAllG + ".\nExpected a diploid entry.");
                            }
                        } else if ((f == iAD) && keepDepth) {
                            // System.out.println("GTS: "+gTS[t][s] + " "+fieldS);
                            if (gTS[t][s] != GenotypeTable.UNKNOWN_GENOTYPE) {
                                int i = 0;
                                for (String ad : Splitter.on(",").split(fieldS)) {
                                    if (i >= alleles.length) continue;
                                    if (alleles[i] == GenotypeTable.UNKNOWN_ALLELE || ad.equals(".") || alleles[i] == NucleotideAlignmentConstants.UNDEFINED_ALLELE ||
                                            alleles[i] == NucleotideAlignmentConstants.UNDEFINED_HOMOZYGOUS) {  //no position for depth of unknown alleles or depth is set to missing, so skip
                                        //Uncomment when converted
                                        //dTS[t][alleles[i++]][s] = AlleleDepthUtil.depthIntToByte(AlleleDepthUtil.DEPTH_MISSING);
                                        //Comment next two lines when converted
                                        i++;
                                        continue;
                                    }

                                    int adInt = Integer.parseInt(ad);
                                    dTS[t][alleles[i++]][s] = AlleleDepthUtil.depthIntToByte(adInt);
                                }
                            }
                        }
                        f++;
                    }
                    t++;
                }
            } catch (IllegalStateException e) {
                e.printStackTrace();
                throw e;
            } catch (Exception e) {
                System.err.println("Err Site Number:" + s);
                //System.err.println("Err Site Number:"+blkPosList.get(blkPosList.size()-1).toString());
                //System.err.println("Err:"+input);
                e.printStackTrace();
                throw e;
            }

        }
        txtL = null;
        return this;
    }

    int getSiteNumber() {
        return siteN;
    }

    byte[][] getGenoTS() {
        return gTS;
    }

    byte[][][] getDepthTS() {
        return dTS;
    }

    ArrayList<Position> getBlkPosList() {
        return blkPosList;
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i = 0; i < sa.length; i++) {
            if (sa[i].equals(match)) return i;
        }
        return -1;
    }
}

