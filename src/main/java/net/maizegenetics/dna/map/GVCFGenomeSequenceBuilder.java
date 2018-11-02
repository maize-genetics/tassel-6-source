package net.maizegenetics.dna.map;

import com.google.common.base.Splitter;
import com.google.common.collect.*;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;


/**
 * Created by zrm22 on 3/27/17.
 */
public class GVCFGenomeSequenceBuilder extends GenomeSequenceBuilder {
    private static final Pattern TAB_PATTERN = Pattern.compile("[\\t]+");

    /**
     * Builds GenomeSequence from a fasta file and a GVCF file.
     *
     * @param fastaFileName full path to fasta file
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(String fastaFileName, String gvcfFileName) throws Exception{
        Function<Character, Character> charConversion = (c) -> c;
        return instance(fastaFileName, charConversion, gvcfFileName);
    }

    /**
     * Builds GenomeSequence from a fasta file.  The char conversion provide a mechanism to convert upper and lower case
     * or convert one case to N.  This is useful if a case if used to define a certain class of bases
     *
     * @param fastaFileName  full path to fasta file
     * @param charConversion lambda Function to convert characters
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(String fastaFileName, Function<Character, Character> charConversion, String gvcfFileName) throws Exception{
        long chrPosMapStartTime = System.currentTimeMillis();
        Map<Chromosome, byte[]> chromPositionMap = readReferenceGenomeChr(fastaFileName, charConversion);
        System.out.println("Done Reading in Reference: Total Time Taken: "+(System.currentTimeMillis()-chrPosMapStartTime)+"ms");
        //need to create a Map<chr,Map<range,Map<AnnotationName,Value>>>
//        Map<Chromosome, RangeMap<Integer,GeneralAnnotationStorage>> gvcfAnnotationsAndCalls = readGVCFFile(gvcfFileName);
        System.out.println("Generating Position List");
        long startTime = System.currentTimeMillis();
        PositionList gvcfPositionsAndAnnotations = readGVCFFilePositionList(gvcfFileName);
        System.out.println("Done Creating Pos list. Total Time Taken: "+(System.currentTimeMillis() - startTime)+"ms");
        return new HalfByteGenomeSequenceGVCF(chromPositionMap,gvcfPositionsAndAnnotations);
    }

    public static GenomeSequence instance(GVCFGenomeSequence base, BitSet maskedBitSet, BitSet filteredBitSet) throws Exception{
        return new HalfByteGenomeSequenceGVCF(base.getChrPosMap(), base.getGVCFPositions(), maskedBitSet, filteredBitSet);

//        Map<Chromosome, byte[]> chromPositionMap = readReferenceGenomeChr(fastaFileName, charConversion);
//        //need to create a Map<chr,Map<range,Map<AnnotationName,Value>>>
//        Map<Chromosome, RangeMap<Integer,GeneralAnnotationStorage>> gvcfAnnotationsAndCalls = readGVCFFile(gvcfFileName);
//        PositionList gvcfPositionsAndAnnotations = readGVCFFilePositionList(gvcfFileName);
//        return new HalfByteGenomeSequenceGVCF(chromPositionMap,gvcfPositionsAndAnnotations);
    }


    private static Map<Chromosome,RangeMap<Integer,GeneralAnnotationStorage>> readGVCFFile(String gvcfFileName) {
        HashMap<Chromosome, RangeMap<Integer,GeneralAnnotationStorage>> chromosomeRangeMapHashMap = new HashMap<>();
        //TODO multithread using similar code to BuilderFromVCF

        return chromosomeRangeMapHashMap;
    }


    private static PositionList readGVCFFilePositionList(String gvcfFileName) throws Exception{
        ArrayList<Position> positionArrayList = new ArrayList<>();

//        BufferedReader gvcfFileReader = new BufferedReader(new FileReader(gvcfFileName));
        BufferedReader gvcfFileReader = Utils.getBufferedReader(gvcfFileName,-1);
        //Loop through the headers
        String currentLine = "";
        while(!(currentLine = gvcfFileReader.readLine()).startsWith("#CHROM")) {

        }
        //parse the header
        String[] header = TAB_PATTERN.split(currentLine);
        HeaderPositions hp= new HeaderPositions(header);
        GVCFPositionRecord gvcfPositionRecord = new GVCFPositionRecord(hp);

        while((currentLine = gvcfFileReader.readLine())!=null) {
            positionArrayList.add(gvcfPositionRecord.parseGVCFRecords(currentLine));
        }
        PositionList instance = new PositionArrayList(positionArrayList, "AGPv4");
        return instance;
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

    public HeaderPositions(String[] header){
        int chrIdx=firstEqualIndex(header,"#CHROM");
        if(chrIdx<0) chrIdx=firstEqualIndex(header,"#CHR");
        CHROMOSOME_INDEX=chrIdx;
        POSITION_INDEX=firstEqualIndex(header,"POS");
        SNPID_INDEX=firstEqualIndex(header,"ID");
        REF_INDEX=firstEqualIndex(header,"REF");
        ALT_INDEX=firstEqualIndex(header,"ALT");
        QUAL_INDEX=firstEqualIndex(header,"QUAL");
        FILTER_INDEX=firstEqualIndex(header,"FILTER");
        INFO_INDEX=firstEqualIndex(header,"INFO");
        FORMAT_INDEX=firstEqualIndex(header,"FORMAT");

        NUM_HAPMAP_NON_TAXA_HEADERS=Math.max(INFO_INDEX,FORMAT_INDEX)+1;
        GENOIDX=NUM_HAPMAP_NON_TAXA_HEADERS;
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i=0; i<sa.length; i++) {
            if(sa[i].equals(match)) return i;
        }
        return -1;
    }

}




class GVCFRecord {
    //Simple class to hold a gvcf record in a storage efficient manner
    //Need the following:
    //the call
    int[] call = new int[2];
    //phasing
    boolean phased = false;
    //the known variants
    ArrayList<String> knownVariants = new ArrayList<>();
    //the depth
    int depth = 0;
    //GQ score
    int gqScore = 0;
    //GeneralAnnotations for any other in INFO tag
    GeneralAnnotationStorage generalAnnotationStorage;


}

class GVCFPositionRecord {
    GeneralPosition ap= new GeneralPosition.Builder(new Chromosome("1"),1232)
            .maf(0.05f)
            .build();
    HeaderPositions hp;

    public GVCFPositionRecord(HeaderPositions hp) {
        this.hp = hp;

    }
    public Position parseGVCFRecords(String input) {
        //Figure out the tab positioning for the header columns
        int[] tabPos=new int[hp.NUM_HAPMAP_NON_TAXA_HEADERS+1];
        int tabIndex=0;
        int len=input.length();
        for (int i=0; (tabIndex<hp.NUM_HAPMAP_NON_TAXA_HEADERS+1)&&(i<len); i++) {
            if (input.charAt(i)=='\t') {
                tabPos[tabIndex++]=i;
            }
        }
        String chrName=input.substring(0, tabPos[hp.CHROMOSOME_INDEX]);
        Chromosome currChr=new Chromosome(new String(chrName));

        String snpID=null;
        if(hp.SNPID_INDEX>0) snpID=input.substring(tabPos[hp.SNPID_INDEX-1]+1, tabPos[hp.SNPID_INDEX]);

        //create the General position
        GeneralPosition.Builder apb=new GeneralPosition.Builder(currChr, Integer.parseInt(input.substring(tabPos[hp.POSITION_INDEX-1]+1, tabPos[hp.POSITION_INDEX])));
        if(snpID!=null && !snpID.equals(".")) {
            apb.snpName(snpID);
        }



        String refS=input.substring(tabPos[hp.REF_INDEX-1]+1, tabPos[hp.REF_INDEX]);
        String alt=input.substring(tabPos[hp.ALT_INDEX-1]+1, tabPos[hp.ALT_INDEX]);
        //create an String to hold the list of variants, we will have to parse this on export
        String allAlleles = refS+"/"+alt.replace(",","/");

        apb = apb.knownVariants(allAlleles);

        //loop through the INFO tag and save off any of those values into positions
        for(String annoS: Splitter.on(";").split(input.substring(tabPos[hp.INFO_INDEX-1]+1, tabPos[hp.INFO_INDEX]))) {
            apb.addAnno(annoS);
        }


        final int iGT=0; //genotype index
        int iAD=-1,iDP=-1,iGQ=-1, iPL=-1;  //alleleDepth, overall depth, genotypeQuality, phredGenotypeLikelihoods
        if(hp.FORMAT_INDEX>=0) {
            //Check to see if FORMAT tag is missing. Only applicable for single taxa files
            if(tabPos[hp.FORMAT_INDEX]==0) {
                throw new IllegalStateException("Error Processing VCF: Missing FORMAT tag.");
            }
            String unsplitInput = input.substring(tabPos[hp.FORMAT_INDEX-1]+1, tabPos[hp.FORMAT_INDEX]);
            if(unsplitInput.length()==0|| !unsplitInput.startsWith("GT")) {
                //Check to see it has the GT field
                if(unsplitInput.contains("GT")) {
                    throw new IllegalStateException("Error Processing VCF Block: GT field is not in first position of FORMAT.");
                }
                //If GT isnt in, we assume that it is missing FORMAT
                else {
                    throw new IllegalStateException("Error Processing VCF Block: Missing FORMAT tag.");
                }
            }
            String[] formatS = unsplitInput.split(":");

            iAD=firstEqualIndex(formatS,"AD");
            iDP=firstEqualIndex(formatS,"DP");
            iGQ=firstEqualIndex(formatS,"GQ");
        }



        //after info is recorded we need to record the call section
        String taxaAllG = input.substring(tabPos[hp.NUM_HAPMAP_NON_TAXA_HEADERS-1]+1);
        int f = 0;
        for(String fieldS: Splitter.on(":").split(taxaAllG)) {
            if (f == iGT) {
                apb.addAnno("GT",fieldS);
            }
            if(f == iAD) {
                apb.addAnno("AD",fieldS);
            }
            if(f == iDP) {
                apb.addAnno("DP",fieldS);
            }
            if(f == iGQ) {
                apb.addAnno("GQ",fieldS);
            }

            f++;
        }



        return apb.build();
    }

    private static int firstEqualIndex(String[] sa, String match) {
        for (int i=0; i<sa.length; i++) {
            if(sa[i].equals(match)) return i;
        }
        return -1;
    }
}
/**
 * ReferenceGenomeSequence class.  This class is used to read chromosome sequences
 * from fasta files.  Data is stored as half-bytes packed into a byte array.
 * This byte array comprises the "value" for a hash map whose key is a
 * Chromosome object.
 *
 * The class also contains methods to obtain a full or partial genome sequence for a
 * specified stored chromosome.
 *
 * @author Lynn Johnson, Zack Miller
 *
 */
//Try to use an index to quickly grab only the useful positions from the position list.
//Trying to multithread this.
class HalfByteGenomeSequenceGVCF implements GVCFGenomeSequence{
    private Map<Chromosome, byte[]> chromPositionMap;
    private Map<Chromosome, Integer> chromLengthLookup=new HashMap<>();
    private RangeMap<Long,Chromosome> wholeGenomeIndexMap= TreeRangeMap.create();
    private  PositionList gvcfAnnotationsAndCalls;
    private final long genomeSize;
    private BitSet maskBitSet;
    private BitSet filterBitSet;

    //Timing stuff
    long totalQueryTime = 0;


    HashMap<Chromosome,int[]> positionIndex = new HashMap<>();
    int indexScaleFactor = 10000;
    protected HalfByteGenomeSequenceGVCF(Map<Chromosome, byte[]>chromPositionMap, PositionList gvcfAnnotationsAndCalls) {
        this.chromPositionMap = chromPositionMap;
        this.gvcfAnnotationsAndCalls = gvcfAnnotationsAndCalls;
        long initialSequenceSetupStart = System.currentTimeMillis();
        chromPositionMap.entrySet().stream()
                .forEach(e -> chromLengthLookup.put(e.getKey(),e.getKey().getLength()));
        LongAdder genomeIndex=new LongAdder();
        chromosomes().stream().sorted()
                .forEach(chrom -> {
                            int length=chromLengthLookup.get(chrom);
                            wholeGenomeIndexMap.put(Range.closed(genomeIndex.longValue(),
                                    genomeIndex.longValue()+length-1),chrom);
                            genomeIndex.add(length);
                            int[] currentIndexArray = new int[length/indexScaleFactor+2];
                            Arrays.fill(currentIndexArray,-1);
                            positionIndex.put(chrom,currentIndexArray);
                        }
                );
        genomeSize=genomeIndex.longValue();
        maskBitSet = new OpenBitSet(gvcfAnnotationsAndCalls.size());
        filterBitSet = new OpenBitSet(gvcfAnnotationsAndCalls.size());
        System.out.println("Initial Startup time for GVCF Sequence: "+(System.currentTimeMillis() - initialSequenceSetupStart)+"ms");
        long indexInitStart = System.currentTimeMillis();
        for(int i = 0; i < gvcfAnnotationsAndCalls.size(); i++) {
            int currentBin = gvcfAnnotationsAndCalls.get(i).getPosition()/indexScaleFactor;
            //new record.  add it to the index
            if(positionIndex.get(gvcfAnnotationsAndCalls.get(i).getChromosome())[currentBin] == -1) {
                positionIndex.get(gvcfAnnotationsAndCalls.get(i).getChromosome())[currentBin] = i;
            }
        }
        System.out.println("Done Creating index: "+(System.currentTimeMillis() - indexInitStart)+"ms");
    }

    protected HalfByteGenomeSequenceGVCF(Map<Chromosome, byte[]>chromPositionMap, PositionList gvcfAnnotationsAndCalls,BitSet maskBitSet, BitSet filterBitSet) {
        this.chromPositionMap = chromPositionMap;
        this.gvcfAnnotationsAndCalls = gvcfAnnotationsAndCalls;
        chromPositionMap.entrySet().stream()
                .forEach(e -> chromLengthLookup.put(e.getKey(),e.getKey().getLength()));
        LongAdder genomeIndex=new LongAdder();
        chromosomes().stream().sorted()
                .forEach(chrom -> {
                            int length=chromLengthLookup.get(chrom);
                            wholeGenomeIndexMap.put(Range.closed(genomeIndex.longValue(),
                                    genomeIndex.longValue()+length-1),chrom);
                            genomeIndex.add(length);
                            int[] currentIndexArray = new int[length/indexScaleFactor+2];
                            Arrays.fill(currentIndexArray,-1);
                            positionIndex.put(chrom,currentIndexArray);
                        }
                );
        genomeSize=genomeIndex.longValue();
        this.maskBitSet = maskBitSet;
        this.filterBitSet = filterBitSet;

        for(int i = 0; i < gvcfAnnotationsAndCalls.size(); i++) {
            int currentBin = gvcfAnnotationsAndCalls.get(i).getPosition()/indexScaleFactor;
            //new record.  add it to the index
            if(positionIndex.get(gvcfAnnotationsAndCalls.get(i).getChromosome())[currentBin] == -1) {
                positionIndex.get(gvcfAnnotationsAndCalls.get(i).getChromosome())[currentBin] = i;
            }
        }
    }
    @Override
    public Set<Chromosome> chromosomes() {
        return chromPositionMap.keySet();
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom) {
        return chromosomeSequence(chrom,1,chromLengthLookup.get(chrom));
    }

    @Override
    // Code assumes 1-based coordinates have been passed.  It will catch and return
    // null if the startSite is 0.  Otherwise, the user is on their own to ensure
    // input is 1-based.
    //TODO add in functionality using GVCF annotations
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int lastSite) {
        final int startSiteFinal = startSite;
        final int lastSiteFinal = lastSite;

        //add one as we grab the last site as well.
//        regionTotalReferenceBPCount = lastSite - startSite+1;

        //zrm22 new for getting gvcf records.
        //Basically the loop now changes to loop through the positionList to check if the allele should be the ref or not
        gvcfAnnotationsAndCalls.chromosomeSiteCount(chrom);

        //use the index to get the currentList of positions
        int[] currentIndex = positionIndex.get(chrom);
        int startPositionListIndex = currentIndex[startSite/indexScaleFactor];
        int endPositionListIndex = currentIndex[lastSite/indexScaleFactor+1];

        int startPositionListCounter = 1;
        while(startPositionListIndex==-1) {
            if (startSite / indexScaleFactor - startPositionListCounter >= 0) {
                startPositionListIndex = currentIndex[startSite / indexScaleFactor - startPositionListCounter];
                startPositionListCounter++;
            }
            else {
                startPositionListIndex = 0;
                break;
            }
        }
        ArrayList<Position> listOfChrPositions = new ArrayList<>();
        if(endPositionListIndex==-1) {
            for(int i = startPositionListIndex; i < startPositionListIndex+indexScaleFactor && i < gvcfAnnotationsAndCalls.size(); i++) {

                if(gvcfAnnotationsAndCalls.get(i).getPosition()>=startSiteFinal && gvcfAnnotationsAndCalls.get(i).getPosition()<=lastSiteFinal) {
                    listOfChrPositions.add(gvcfAnnotationsAndCalls.get(i));
                }
            }
        }
        else {
            for (int i = startPositionListIndex; i <= endPositionListIndex; i++) {
                if (gvcfAnnotationsAndCalls.get(i).getPosition() >= startSiteFinal && gvcfAnnotationsAndCalls.get(i).getPosition() <= lastSiteFinal) {
                    listOfChrPositions.add(gvcfAnnotationsAndCalls.get(i));
                }
            }
        }

        //sort things to make sure we can iterate properly
        Collections.sort(listOfChrPositions);


        int actualStartPos = 0;
        int startIndex = 0;
        int actualEndPos = 0;
        int endIndex = 0;

        boolean addingAtStart = false;
        if(listOfChrPositions.size()>0) {
            for (int i = 1; i < gvcfAnnotationsAndCalls.size(); i++) {
                if (gvcfAnnotationsAndCalls.get(i).getPosition() == listOfChrPositions.get(0).getPosition()) {
                    if (checkPositionCoverage(gvcfAnnotationsAndCalls.get(i - 1), startSiteFinal)) {
                        startIndex = i-1;
                        listOfChrPositions.add(0,gvcfAnnotationsAndCalls.get(startIndex));
                        addingAtStart = true;
                    }
                    break;
                }
            }
        }

        //loop through from startSite till the first position in the GVCF export these alleles directly from the ref
        int listOfChrPositionsCounter = 0;
        ArrayList<Byte> byteList = new ArrayList<>();

        if(listOfChrPositions.size()==0) {
            //Export a single N for a missing sequence
            byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
        }
        else {
            for (int siteCounter = startSite; siteCounter <= lastSite && listOfChrPositionsCounter < listOfChrPositions.size(); siteCounter++) {
                GeneralAnnotation ga = listOfChrPositions.get(listOfChrPositionsCounter).getAnnotation();
                Set<String> annotationKeys = ga.getAnnotationKeys();
                if(addingAtStart) {
                    if (listOfChrPositionsCounter != 0) {
                        if (filterBitSet.fastGet(listOfChrPositionsCounter-1)) {
                            //get the size of the position
                            if (annotationKeys.contains("END")) {
                                int lengthOfPosBlock = Integer.parseInt((String) ga.getTextAnnotation("END")[0]) - listOfChrPositions.get(listOfChrPositionsCounter).getPosition();
                                siteCounter+=lengthOfPosBlock;
                            }
                            else{
                                //check insertion or deletion
                                String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                                siteCounter+=variants[0].length()-1;
                            }
                            listOfChrPositionsCounter++;
                            continue;
                        }
                    }
                }
                else {
                    if (filterBitSet.fastGet(listOfChrPositionsCounter)) {
                        if (annotationKeys.contains("END")) {
                            int lengthOfPosBlock = Integer.parseInt((String) ga.getTextAnnotation("END")[0]) - listOfChrPositions.get(listOfChrPositionsCounter).getPosition();
                            siteCounter+=lengthOfPosBlock;
                        }
                        else{
                            //check insertion or deletion
                            String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                            siteCounter+=variants[0].length()-1;
                        }
                        listOfChrPositionsCounter++;
                        continue;
                    }
                }
                if (maskBitSet.fastGet(listOfChrPositionsCounter)) {
                    //if true it means we have to mask all of these positions
                    //Check to see if we have an END anno
                }
                if (listOfChrPositions.get(listOfChrPositionsCounter).getPosition() <= siteCounter) {
                    if (annotationKeys.contains("END")) {
                        //this means we have a block
                        //Check the call and grab the corresponding allele values
                        String call = ga.getTextAnnotation("GT")[0];
                        boolean phased = true;
                        boolean haploid = false;
                        if (call.contains("/")) {
                            phased = false;
                        } else if (call.contains("|")) {
                            phased = true;
                        } else {
                            haploid = true;
                        }
                        String[] callSplit = phased ? call.split("|") : call.split("/");
                        int leftAllele = Integer.parseInt(callSplit[0]);
                        int rightAllele = leftAllele;
                        if (!haploid) {
                            rightAllele = Integer.parseInt(callSplit[1]);
                        }
                        int endPoint = Integer.parseInt(ga.getTextAnnotation("END")[0]) - 1;
                        //check to see if our requested end point is smaller than the region
                        if(lastSite-1<endPoint) {
                            endPoint = lastSite-1;
                        }
                        //Get the known Variants
                        int startSiteShifted = siteCounter - 1;  //shift over to zero bas

                        if (startSite < 0)
                            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: starting parameter is less than 1 for 1-based method");
                        ; // method needs 1-based coordinates.
                        byte[] packedBytes = chromPositionMap.get(chrom);
                        if (packedBytes == null)
                            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); // chromosome not found
                        if (startSiteShifted > packedBytes.length * 2 || endPoint > packedBytes.length * 2) {
                            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested sequence is out of range"); // requested sequence is out of range
                        }
                        String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                        String leftAlleleString = variants[leftAllele];

                        //check depth
//                        String dpFullString = (String) annos.get("DP").toArray()[0];
                        String dpFullString = ga.getTextAnnotation("DP")[0];

                        // String minDPString = (String) annos.get("MIN_DP").toArray()[0];
                        //TODO check to see if we should assume only Ref or<NON_REF> call for blocks
                        if (maskBitSet.fastGet(listOfChrPositionsCounter) || dpFullString.equals("0")) {

                            //if true we have a masking
                            for (int i = startSiteShifted; i <= endPoint; i++) {
//                                regionDepthCount+=0;
//                                regionZeroCoverageCount++;

                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                        } else if (leftAllele == 0) {
                            //fill in with reference sequence
                            //pull the sequence from siteCounter till you get to endPoint
                            //Zrm Apr17 fix for default rules
                            //check if homozygous or het based on GT calls
                            int depth = 0;
                            int minDepth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            for (int i = startSiteShifted; i <= endPoint; i++) {
//                                regionDepthCount+=depth;
//                                regionMinDepthCount+=minDepth;
                                byteList.add((byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F)));
                            }

                        } else {
                            int depth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            for (int i = startSiteShifted; i <= endPoint; i++) {
//                                regionDepthCount+=depth;
//                                regionZeroCoverageCount++;
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                        }
                        //shift up the siteCounter to match the end
                        siteCounter = endPoint + 1;
                    } else {
                        //it is likely a snp
                        //check the call and grab the correct allele
                        //siteCounter will increment correctly
                        String call = ga.getTextAnnotation("GT")[0];
                        boolean phased = true;
                        boolean haploid = false;
                        if (call.contains("/")) {
                            phased = false;
                        } else if (call.contains("|")) {
                            phased = true;
                        } else {
                            //haploid
                            haploid = true;
                        }
                        String[] callSplit = phased ? call.split("|") : call.split("/");
                        int leftAllele = Integer.parseInt(callSplit[0]);
                        int rightAllele = leftAllele;
                        if (!haploid) {
                            rightAllele = Integer.parseInt(callSplit[1]);
                        }

                        //Get the AD field and split it on commas
                        String adsFullString = ga.getTextAnnotation("AD")[0];
                        String dpFullString = ga.getTextAnnotation("DP")[0];
                        String[] adSplit = adsFullString.split(",");
                        //Get the known Variants
                        String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                        boolean isHet = calcHet(adSplit);
                        if (isHet) {
                            //TODO handle hets correctly for indels and such
                            String leftAlleleString = variants[0];
                            for(int alleleCounter = 0; alleleCounter < leftAlleleString.length(); alleleCounter++) {
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                            siteCounter+=leftAlleleString.length()-1;
//                            for(int i = 0; i < leftAlleleString.length();i++) {
//                                regionHetCount++;
//                            }
                        }
                        //ZRM April17 do calls based on depth only
                        //This block is not actually used as it has already been handled if it is a het or a homozygous ref(in an allele block)
                        else if (adSplit[0].equals("1") || adSplit[0].equals("2")) {
                            if (adSplit[1].equals("0")) {
                                //export the ref
                                String leftAlleleString = variants[0];

                                int depth = 0;
                                if(!dpFullString.equals("")) {
                                    depth = Integer.parseInt(dpFullString);
                                }

                                if (maskBitSet.fastGet(listOfChrPositionsCounter)) {
                                    //if true we have to mask it
                                    for (int i = 0; i < leftAlleleString.length() && siteCounter+i < lastSite-1; i++) {
//                                        regionDepthCount+=depth;
//                                        regionNumberHomoRef1Or2Count++;
                                        byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                                    }
                                } else {
                                    for (int i = 0; i < leftAlleleString.length() && siteCounter+i < lastSite-1; i++) {
//                                        regionDepthCount+=depth;
//                                        regionNumberHomoRef1Or2Count++;
                                        byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                                    }
                                }
                            }

                        } else if (adSplit[1].equals("1") || adSplit[1].equals("2")) { //ZRM22 Jun19 ChangeforChristy
//                        } else if (adSplit[1].equals("0")) {
                            //call Ns
                            //export the ref
                            String leftAlleleString = variants[1];
                            int depth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            for (int i = 0; i < leftAlleleString.length(); i++) {
//                                regionDepthCount+=depth;
//                                regionNumberHomoAlt1Or2Count++;
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                        } else if (Integer.parseInt(adSplit[1]) >= 3) { //ZRM22 Jun19 change for christy
//                        } else if (Integer.parseInt(adSplit[1]) >= 1) {
                            int depth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            String leftAlleleString = variants[1];
                            if (maskBitSet.fastGet(listOfChrPositionsCounter)) {
                                //if true we have to mask it
                                for (int i = 0; i < leftAlleleString.length(); i++) {
//                                    regionAltCount++;
//                                    regionNumberHomoAlt3PlusCount++;
//                                    regionDepthCount+=depth;
                                    byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                                }
                            } else {
                                for (int i = 0; i < leftAlleleString.length(); i++) {
//                                    regionAltCount++;
//                                    regionNumberHomoAlt3PlusCount++;
//                                    regionDepthCount+=depth;
                                    byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                                }
                                siteCounter+=variants[0].length()-1;
                            }
                        }
                        else {
                            //in the case where we have 0 depth for each allele it will fall through all the cases
                            //lets stick in Ns for the number of reference basepairs  fix the old method header as well.
                            String leftAlleleString = variants[0];

                            for (int i = 0; i < leftAlleleString.length(); i++) {
//                                    regionAltCount++;
//                                    regionNumberHomoAlt3PlusCount++;
//                                    regionDepthCount+=depth;
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                            }
                            siteCounter+=variants[0].length()-1;
                        }
                    }
                    listOfChrPositionsCounter++;
                } else {
                    //grab the reference as we dont have a gvcf for the requested site
                    //Should this be the breaking point of the sequence??
                    //TODO should we mark these with Ns?

                    //With 0 coverage we want to export an N
                    byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
//                    regionDepthCount+=0;
//                    regionZeroCoverageCount++;

                }
            }
        }
//        regionTotalExportedBPCount = byteList.size();
        byte[] fullBytes = new byte[byteList.size()];
        for(int i = 0; i < fullBytes.length; i++) {
            fullBytes[i] = byteList.get(i);
        }
        return fullBytes;
    }


    @Override
    //Code to get the Sequence in a string as well as the various stats we have collected.
    // This should allow for safe multithreading for the stats.
    public HashMap<String, String> chromosomeSequenceAndStats(Chromosome chrom, int startSite, int lastSite) {
        long queryTimeStart = System.currentTimeMillis();

        int regionTotalReferenceBPCount = 0;
        int regionTotalExportedBPCount = 0;
        int regionHetCount = 0;
        int regionAltCount = 0;
        int regionDepthCount = 0;
        int regionGQCount = 0;
        int regionMinDepthCount = 0;
        int regionZeroCoverageCount = 0;
        int regionNumberHomoRef1Or2Count = 0;
        int regionNumberHomoRef3PlusCount = 0;
        int regionNumberHomoAlt1Or2Count = 0;
        int regionNumberHomoAlt3PlusCount = 0;

        final int startSiteFinal = startSite;
        final int lastSiteFinal = lastSite;

        //add one as we grab the last site as well.
        regionTotalReferenceBPCount = lastSite - startSite+1;

        //zrm22 new for getting gvcf records.
        //Basically the loop now changes to loop through the positionList to check if the allele should be the ref or not
//        gvcfAnnotationsAndCalls.chromosomeSiteCount(chrom);

        long getCurrentSubsetPositionListTimeStart = System.currentTimeMillis();
        //use the index to get the currentList of positions
        int[] currentIndex = positionIndex.get(chrom);
        int startPositionListIndex = currentIndex[startSite/indexScaleFactor];
        int endPositionListIndex = currentIndex[lastSite/indexScaleFactor+1];

        int startPositionListCounter = 1;
        while(startPositionListIndex==-1) {
            startPositionListIndex = currentIndex[startSite/indexScaleFactor-startPositionListCounter];
            startPositionListCounter++;
        }
        ArrayList<Position> listOfChrPositions = new ArrayList<>();
        int actualStartIndex = Integer.MAX_VALUE;
        if(endPositionListIndex==-1) {
            for(int i = startPositionListIndex; i < startPositionListIndex+indexScaleFactor && i < gvcfAnnotationsAndCalls.size(); i++) {

                if(gvcfAnnotationsAndCalls.get(i).getPosition()>=startSiteFinal && gvcfAnnotationsAndCalls.get(i).getPosition()<=lastSiteFinal) {
                    if(i < actualStartIndex) {
                        actualStartIndex = i;
                    }
                    listOfChrPositions.add(gvcfAnnotationsAndCalls.get(i));
                }
            }
        }
        else {
            for (int i = startPositionListIndex; i <= endPositionListIndex; i++) {
                if (gvcfAnnotationsAndCalls.get(i).getPosition() >= startSiteFinal && gvcfAnnotationsAndCalls.get(i).getPosition() <= lastSiteFinal) {
                    if(i < actualStartIndex) {
                        actualStartIndex = i;
                    }
                    listOfChrPositions.add(gvcfAnnotationsAndCalls.get(i));
                }
            }
        }
        long totalTimeToSubsetPosList = System.currentTimeMillis() - getCurrentSubsetPositionListTimeStart;
        //sort things to make sure we can iterate properly
        long sortPosTimeStart = System.currentTimeMillis();
        Collections.sort(listOfChrPositions);
        long totalTimeToSortPosList = System.currentTimeMillis() - sortPosTimeStart;

        int actualStartPos = 0;
        int startIndex = 0;
        int actualEndPos = 0;
        int endIndex = 0;

        long addAPositionToStartTime = System.currentTimeMillis();
        boolean addingAtStart = false;
        if(listOfChrPositions.size()>0 && startPositionListIndex!=0) {
            if(checkPositionCoverage(gvcfAnnotationsAndCalls.get(actualStartIndex-1),startSiteFinal)) {
                startIndex = actualStartIndex -1;
                listOfChrPositions.add(0,gvcfAnnotationsAndCalls.get(startIndex));
                addingAtStart = true;
            }


//            for (int i = 1; i < gvcfAnnotationsAndCalls.size(); i++) {
//            for (int i = startPositionListIndex; i < gvcfAnnotationsAndCalls.size(); i++) {
//                if (gvcfAnnotationsAndCalls.get(i).getPosition() == listOfChrPositions.get(0).getPosition()) {
//                    if (checkPositionCoverage(gvcfAnnotationsAndCalls.get(i - 1), startSiteFinal)) {
//                        startIndex = i-1;
//                        listOfChrPositions.add(0,gvcfAnnotationsAndCalls.get(startIndex));
//                        addingAtStart = true;
//                    }
//                    break;
//                }
//            }
        }
        long totalTimeToAddStartPosition = System.currentTimeMillis() - addAPositionToStartTime;

        //loop through from startSite till the first position in the GVCF export these alleles directly from the ref
        int listOfChrPositionsCounter = 0;
        ArrayList<Byte> byteList = new ArrayList<>();

        if(listOfChrPositions.size()==0) {
            //Export a single N for a missing sequence
            byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
        }
        else {
            for (int siteCounter = startSite; siteCounter <= lastSite && listOfChrPositionsCounter < listOfChrPositions.size(); siteCounter++) {
                GeneralAnnotation ga = listOfChrPositions.get(listOfChrPositionsCounter).getAnnotation();
                Set<String> annotationKeys = ga.getAnnotationKeys();
                if(addingAtStart) {
                    if (listOfChrPositionsCounter != 0) {
                        if (filterBitSet.fastGet(listOfChrPositionsCounter-1)) {
                            //get the size of the position
                            if (annotationKeys.contains("END")) {
                                int lengthOfPosBlock = Integer.parseInt((String) ga.getTextAnnotation("END")[0]) - listOfChrPositions.get(listOfChrPositionsCounter).getPosition();
                                siteCounter+=lengthOfPosBlock;
                            }
                            else{
                                //check insertion or deletion
                                String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                                siteCounter+=variants[0].length()-1;
                            }
                            listOfChrPositionsCounter++;
                            continue;
                        }
                    }
                }
                else {
                    if (filterBitSet.fastGet(listOfChrPositionsCounter)) {
                        if (annotationKeys.contains("END")) {
                            int lengthOfPosBlock = Integer.parseInt((String) ga.getTextAnnotation("END")[0]) - listOfChrPositions.get(listOfChrPositionsCounter).getPosition();
                            siteCounter+=lengthOfPosBlock;
                        }
                        else{
                            //check insertion or deletion
                            String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                            siteCounter+=variants[0].length()-1;
                        }
                        listOfChrPositionsCounter++;
                        continue;
                    }
                }
                if (maskBitSet.fastGet(listOfChrPositionsCounter)) {
                    //if true it means we have to mask all of these positions
                    //Check to see if we have an END anno
                }
                if (listOfChrPositions.get(listOfChrPositionsCounter).getPosition() <= siteCounter) {
                    if (annotationKeys.contains("END")) {
                        //this means we have a block
                        //Check the call and grab the corresponding allele values
                        String call = ga.getTextAnnotation("GT")[0];
                        boolean phased = true;
                        boolean haploid = false;
                        if (call.contains("/")) {
                            phased = false;
                        } else if (call.contains("|")) {
                            phased = true;
                        } else {
                            haploid = true;
                        }
                        String[] callSplit = phased ? call.split("|") : call.split("/");
                        int leftAllele = Integer.parseInt(callSplit[0]);
                        int rightAllele = leftAllele;
                        if (!haploid) {
                            rightAllele = Integer.parseInt(callSplit[1]);
                        }
                        int endPoint = Integer.parseInt(ga.getTextAnnotation("END")[0]) - 1;
                        //check to see if our requested end point is smaller than the region
                        if(lastSite-1<endPoint) {
                            endPoint = lastSite-1;
                        }
                        //Get the known Variants
                        int startSiteShifted = siteCounter - 1;  //shift over to zero bas

                        if (startSite < 0)
                            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: starting parameter is less than 1 for 1-based method");
                        ; // method needs 1-based coordinates.
                        byte[] packedBytes = chromPositionMap.get(chrom);
                        if (packedBytes == null)
                            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); // chromosome not found
                        if (startSiteShifted > packedBytes.length * 2 || endPoint > packedBytes.length * 2) {
                            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested sequence is out of range"); // requested sequence is out of range
                        }
                        String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                        String leftAlleleString = variants[leftAllele];

                        //check depth
//                        String dpFullString = (String) annos.get("DP").toArray()[0];
                        String dpFullString = ga.getTextAnnotation("DP")[0];

                        // String minDPString = (String) annos.get("MIN_DP").toArray()[0];
                        //TODO check to see if we should assume only Ref or<NON_REF> call for blocks
                        if (maskBitSet.fastGet(listOfChrPositionsCounter) || dpFullString.equals("0")) {

                            //if true we have a masking
                            for (int i = startSiteShifted; i <= endPoint; i++) {
                                regionDepthCount+=0;
                                regionZeroCoverageCount++;

                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                        } else if (leftAllele == 0) {
                            //fill in with reference sequence
                            //pull the sequence from siteCounter till you get to endPoint
                            //Zrm Apr17 fix for default rules
                            //check if homozygous or het based on GT calls
                            int depth = 0;
                            int minDepth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            for (int i = startSiteShifted; i <= endPoint; i++) {
                                regionDepthCount+=depth;
                                regionMinDepthCount+=minDepth;
                                byteList.add((byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F)));
                            }

                        } else {
                            int depth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            for (int i = startSiteShifted; i <= endPoint; i++) {
                                regionDepthCount+=depth;
                                regionZeroCoverageCount++;
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                        }
                        //shift up the siteCounter to match the end
                        siteCounter = endPoint + 1;
                    } else {
                        //it is likely a snp
                        //check the call and grab the correct allele
                        //siteCounter will increment correctly
                        String call = ga.getTextAnnotation("GT")[0];
                        boolean phased = true;
                        boolean haploid = false;
                        if (call.contains("/")) {
                            phased = false;
                        } else if (call.contains("|")) {
                            phased = true;
                        } else {
                            //haploid
                            haploid = true;
                        }
                        String[] callSplit = phased ? call.split("|") : call.split("/");
                        int leftAllele = Integer.parseInt(callSplit[0]);
                        int rightAllele = leftAllele;
                        if (!haploid) {
                            rightAllele = Integer.parseInt(callSplit[1]);
                        }

                        //Get the AD field and split it on commas
                        String adsFullString = ga.getTextAnnotation("AD")[0];
                        String dpFullString = ga.getTextAnnotation("DP")[0];
                        String[] adSplit = adsFullString.split(",");
                        //Get the known Variants
                        String[] variants = listOfChrPositions.get(listOfChrPositionsCounter).getKnownVariants();
                        boolean isHet = calcHet(adSplit);
                        if (isHet) {
                            //TODO handle hets correctly for indels and such
                            String leftAlleleString = variants[0];
                            siteCounter+=leftAlleleString.length()-1;
                            for(int i = 0; i < leftAlleleString.length();i++) {
                                regionHetCount++;
                            }
                        }
                        //ZRM April17 do calls based on depth only
                        else if (adSplit[0].equals("1") || adSplit[0].equals("2")) {
                            if (adSplit[1].equals("0")) {
                                //export the ref
                                String leftAlleleString = variants[0];

                                int depth = 0;
                                if(!dpFullString.equals("")) {
                                    depth = Integer.parseInt(dpFullString);
                                }

                                if (maskBitSet.fastGet(listOfChrPositionsCounter)) {
                                    //if true we have to mask it
                                    for (int i = 0; i < leftAlleleString.length() && siteCounter+i < lastSite-1; i++) {
                                        regionDepthCount+=depth;
                                        regionNumberHomoRef1Or2Count++;
                                        byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                                    }
                                } else {
                                    for (int i = 0; i < leftAlleleString.length() && siteCounter+i < lastSite-1; i++) {
                                        regionDepthCount+=depth;
                                        regionNumberHomoRef1Or2Count++;
                                        byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                                    }
                                }
                            }

                        } else if (adSplit[1].equals("1") || adSplit[1].equals("2")) {
                            //call Ns
                            //export the ref
                            String leftAlleleString = variants[1];
                            int depth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            for (int i = 0; i < leftAlleleString.length(); i++) {
                                regionDepthCount+=depth;
                                regionNumberHomoAlt1Or2Count++;
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                            }
                        } else if (Integer.parseInt(adSplit[1]) >= 3) {
                            int depth = 0;
                            if(!dpFullString.equals("")) {
                                depth = Integer.parseInt(dpFullString);
                            }
                            String leftAlleleString = variants[1];
                            if (maskBitSet.fastGet(listOfChrPositionsCounter)) {
                                //if true we have to mask it
                                for (int i = 0; i < leftAlleleString.length(); i++) {
                                    regionAltCount++;
                                    regionNumberHomoAlt3PlusCount++;
                                    regionDepthCount+=depth;
                                    byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                                }
                            } else {
                                for (int i = 0; i < leftAlleleString.length(); i++) {
                                    regionAltCount++;
                                    regionNumberHomoAlt3PlusCount++;
                                    regionDepthCount+=depth;
                                    byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                                }
                                siteCounter+=variants[0].length()-1;
                            }
                        }
                        else {
                            //in the case where we have 0 depth for each allele it will fall through all the cases
                            //lets stick in Ns for the number of reference basepairs  fix the old method header as well.
                            String leftAlleleString = variants[0];

                            for (int i = 0; i < leftAlleleString.length(); i++) {
//                                    regionAltCount++;
//                                    regionNumberHomoAlt3PlusCount++;
//                                    regionDepthCount+=depth;
                                byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte(leftAlleleString.charAt(i)));
                            }
                            siteCounter+=variants[0].length()-1;
                        }
                    }
                    listOfChrPositionsCounter++;
                } else {
                    //grab the reference as we dont have a gvcf for the requested site
                    //Should this be the breaking point of the sequence??
                    //TODO should we mark these with Ns?

                    //With 0 coverage we want to export an N
                    byteList.add(NucleotideAlignmentConstants.getNucleotideAlleleByte("N"));
                    regionDepthCount+=0;
                    regionZeroCoverageCount++;

                }
            }
        }
        regionTotalExportedBPCount = byteList.size();
        byte[] fullBytes = new byte[byteList.size()];
        for(int i = 0; i < fullBytes.length; i++) {
            fullBytes[i] = byteList.get(i);
        }

        //Compile the HashMap<String,String>();
        HashMap<String,String> sequenceAndStats = new HashMap<>();
        sequenceAndStats.put("Sequence",NucleotideAlignmentConstants.nucleotideBytetoString(fullBytes));

        sequenceAndStats.put("RefSize",""+regionTotalReferenceBPCount);
        sequenceAndStats.put("Size", ""+regionTotalExportedBPCount);
        sequenceAndStats.put("HetCount",""+regionHetCount);
        sequenceAndStats.put("AltCount",""+regionAltCount);
        sequenceAndStats.put("Depth",""+regionDepthCount);
        sequenceAndStats.put("GQ",""+regionGQCount);
        sequenceAndStats.put("Min_Depth",""+regionMinDepthCount);
        sequenceAndStats.put("ZeroCoverageCount",""+regionZeroCoverageCount);
        sequenceAndStats.put("HomoRefCount",""+regionNumberHomoRef1Or2Count);
        sequenceAndStats.put("HomoAltLowDepthCount",""+regionNumberHomoAlt1Or2Count);
        sequenceAndStats.put("HomoAltHighDepthCount",""+regionNumberHomoAlt3PlusCount);

        long queryTimeEnd = System.currentTimeMillis();
//        System.out.println("Current Query Time:"+(queryTimeEnd - queryTimeStart)+"ms\n" +
//                "Time Taken to subset Pos List"+totalTimeToSubsetPosList+"ms\n" +
//                "Time Taken to sort Pos List"+totalTimeToSortPosList+"ms\n" +
//                "Time Taken to add a Position to start of List:"+totalTimeToAddStartPosition+"ms\n");

//        System.out.println("Current Query Time:"+(queryTimeEnd - queryTimeStart)+"ms");


        return sequenceAndStats;


    }

    @Override
    //TODO add in functionality to handle GVCF annotations
    public byte[] genomeSequence(long startSite, long lastSite) {
        if(lastSite-startSite>Integer.MAX_VALUE) throw
                new IllegalArgumentException("Less than "+Integer.MAX_VALUE+" sites must be requested at a time");
        byte[] fullBytes=new byte[(int)(lastSite-startSite+1)];
        long currentSiteToGet=startSite;
        while(currentSiteToGet<lastSite) {
            Map.Entry<Range<Long>,Chromosome> rangeChromEntry=wholeGenomeIndexMap.getEntry(currentSiteToGet);
            int chrStart=(int)(currentSiteToGet-rangeChromEntry.getKey().lowerEndpoint());
            int chrLast=(int)Math.min(rangeChromEntry.getKey().upperEndpoint()-rangeChromEntry.getKey().lowerEndpoint(),lastSite-rangeChromEntry.getKey().lowerEndpoint());
            byte[] chromoSeq=chromosomeSequence(rangeChromEntry.getValue(), chrStart+1,chrLast+1);  //+1 for 0 based genome, 1 based chromosomes
            System.arraycopy(chromoSeq,0,fullBytes,(int)(currentSiteToGet-startSite),chromoSeq.length);
            currentSiteToGet+=chromoSeq.length;
        }
        return fullBytes;
    }

    @Override
    public int chromosomeSize(Chromosome chromosome) {
        return chromLengthLookup.get(chromosome);
    }

    @Override
    public long genomeSize() {
        return genomeSize;
    }

    @Override
    public int numberOfChromosomes() {
        return chromPositionMap.size();
    }

    @Override
    //TODO see if we need to fix this for GVCF annotations
    public Map<Long, Tuple<Chromosome, Integer>> fullRefCoordinateToChromCoordinate(ArrayList<Long> coordinates) {
        // Returns 0-based value from Chromosome array (values are stored as 0-based)
        Map<Long, Tuple<Chromosome, Integer>> mappedCoordinates = new ConcurrentHashMap<Long, Tuple<Chromosome, Integer>>();
        coordinates.stream().parallel().forEach(coordinate -> {
            Map.Entry<Range<Long>,Chromosome> rangeChromEntry=wholeGenomeIndexMap.getEntry(coordinate);
            Chromosome curChrom = rangeChromEntry.getValue();
            long chromCoordinate = coordinate - rangeChromEntry.getKey().lowerEndpoint();
            Tuple<Chromosome, Integer> chromWithCoordinate = new Tuple<>(curChrom, (int)chromCoordinate);
            mappedCoordinates.put(coordinate, chromWithCoordinate);
        });
        return mappedCoordinates;
    }

    private boolean checkPositionCoverage(Position pos, int site) {
//        SetMultimap<String,String> annos = pos.getAnnotation().getAnnotationAsMap();
        GeneralAnnotation ga = pos.getAnnotation();
//        if(annos.containsKey("END")) {
//            int endPoint = Integer.parseInt((String)annos.get("END").toArray()[0]);
//            if(pos.getPosition()<=site && site <= endPoint) {
//                return true;
//            }
//            else {
//                return false;
//            }
//
//        }

        String[] endArray = ga.getTextAnnotation("END");

//        if(ga.getAnnotationKeys().contains("END")) {
        if(endArray.length!=0) {
//            int endPoint = Integer.parseInt(ga.getTextAnnotation("END")[0]);
            int endPoint = Integer.parseInt(endArray[0]);
            if(pos.getPosition()<=site && site <= endPoint) {
                return true;
            }
            else {
                return false;
            }

        }
        else if(pos.getPosition()==site) {
            return true;
        }
        else {
            String[] alleles = pos.getKnownVariants();
            //check to see if we have a deletion and the reference covers the whole
            if(alleles[0].length()+pos.getPosition()>=site) {
                return true;
            }
            else {
                return false;
            }
        }
    }

    public Map<Chromosome, byte[]> getChrPosMap() {
        return chromPositionMap;
    }
    public HashMap<Chromosome,ArrayList<ArrayList<Integer>>> getConsecutiveRegions() {
        HashMap<Chromosome,ArrayList<ArrayList<Integer>>> consecRegions = new HashMap<>();
        Set<Chromosome> chromosomeSet = chromosomes();
        //Loop through each chromosome add to a Map<RangeMap>
        HashMap<Chromosome,RangeSet<Integer>> rangeMaps = new HashMap<>();
        for(Chromosome chr : chromosomeSet) {
            rangeMaps.put(chr,TreeRangeSet.create());
            consecRegions.put(chr,new ArrayList<>());
        }

        for(int i = 0 ; i < gvcfAnnotationsAndCalls.size(); i++) {
            if(filterBitSet.fastGet(i)) {
                //if we have a position to filter out, we should probably break the fasta sequence
                //to do this, just do not add in the range... RangeMap will handle it for you
                continue;
            }

            Position currentPosition = gvcfAnnotationsAndCalls.get(i);
            SetMultimap<String,String> annos = currentPosition.getAnnotation().getAnnotationAsMap();
            if(annos.containsKey("END")) {
                int endPoint = Integer.parseInt((String)annos.get("END").toArray()[0]);
                rangeMaps.get(currentPosition.getChromosome()).add(Range.closed(currentPosition.getPosition(),endPoint+1));
            }
            else {
                //check to make sure that the reference allele is more than one allele(A deletion)
                String[] variants = currentPosition.getKnownVariants();
                //TODO handle hets correctly for indels and such
                String refAlleleString = variants[0];
                if(refAlleleString.length()>1) {
                    //we have a deletion or an insertion+deletion as such our END point will need to take into account the size of the deletion
                    rangeMaps.get(currentPosition.getChromosome()).add(Range.closed(currentPosition.getPosition(),currentPosition.getPosition()+refAlleleString.length()));
                }
                else {
                    //This will cause the index of the next line in the GVCF file to be non consecutive with the previous
                    rangeMaps.get(currentPosition.getChromosome()).add(Range.closed(currentPosition.getPosition(), currentPosition.getPosition() + 1));

                }
            }

        }

        for(Chromosome chr : chromosomeSet) {
            Set<Range<Integer>> rangeSet = rangeMaps.get(chr).asRanges();
            for(Range<Integer> currentRange : rangeSet) {
                ArrayList<Integer> currentList = new ArrayList<Integer>();
                currentList.add(currentRange.lowerEndpoint());
                //We need to subtract 1 point so it will be [inclusive,inclusive] instead of [inclusive,exclusive)
                currentList.add(currentRange.upperEndpoint()-1);
                consecRegions.get(chr).add(currentList);
            }
        }

        return consecRegions;
    }
    //TODO move this to a different class
    public void writeFASTA(String fileName){
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
            HashMap<Chromosome,ArrayList<ArrayList<Integer>>> consecutiveRegions = getConsecutiveRegions();
            Set<Chromosome> chromosomes = chromosomes();
            for(Chromosome chr : chromosomes) {
                for(ArrayList<Integer> bounds : consecutiveRegions.get(chr)) {
                    writer.write(">Chr_"+chr.getChromosomeNumber()+"_StartSite_"+bounds.get(0)+"_EndSite_"+bounds.get(1));
                    writer.newLine();
                    writer.write(""+NucleotideAlignmentConstants.nucleotideBytetoString(chromosomeSequence(chr,bounds.get(0),bounds.get(1))));
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }

    }

    public BitSet getMaskBitSet() {
        return maskBitSet;
    }
    public void setMaskBitSet(BitSet newMaskBitSet) {
        maskBitSet = newMaskBitSet;
    }

    public BitSet getFilterBitSet() {
        return filterBitSet;
    }
    public void setFilterBitSet(BitSet newFilterBitSet) {
        filterBitSet = newFilterBitSet;
    }

    public void flipMaskBit(int index) {
        maskBitSet.fastFlip(index);
    }
    public void flipFilterBit(int index) {
        filterBitSet.fastFlip(index);
    }
    public PositionList getGVCFPositions() {
        return gvcfAnnotationsAndCalls;
    }

    private boolean calcHet(String[] adSplit) {
        boolean isHet = false;
        int refDepth = Integer.parseInt(adSplit[0]);
        int altDepth = Integer.parseInt(adSplit[1]);
        if(refDepth>0 && altDepth>0) {
            isHet = true;
        }
        return isHet;
    }

    @Override
    public void resetCounters() {
//        regionTotalReferenceBPCount = 0;
//        regionTotalExportedBPCount = 0;
//        regionHetCount = 0;
//        regionAltCount = 0;
//        regionDepthCount = 0;
//        regionGQCount = 0;
//        regionMinDepthCount = 0;
//        regionZeroCoverageCount = 0;
//        regionNumberHomoRef1Or2Count = 0;
//        regionNumberHomoAlt1Or2Count = 0;
//        regionNumberHomoAlt3PlusCount = 0;
    }

    @Override
    public HashMap<String,Integer> getPreviousRegionStats() {
        HashMap<String,Integer> stats = new HashMap<>();
//        stats.put("RefSize",regionTotalReferenceBPCount);
//        stats.put("Size", regionTotalExportedBPCount);
//        stats.put("HetCount",regionHetCount);
//        stats.put("AltCount",regionAltCount);
//        stats.put("Depth",regionDepthCount);
//        stats.put("GQ",regionGQCount);
//        stats.put("Min_Depth",regionMinDepthCount);
//        stats.put("ZeroCoverageCount",regionZeroCoverageCount);
//        stats.put("HomoRefCount",regionNumberHomoRef1Or2Count);
//        stats.put("HomoAltLowDepthCount",regionNumberHomoAlt1Or2Count);
//        stats.put("HomoAltHighDepthCount",regionNumberHomoAlt3PlusCount);

        return stats;
    }

    @Override
    public byte genotype(Chromosome chrom, int position) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public byte genotype(Chromosome chrom, Position positionObject) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public String genotypeAsString(Chromosome chrom, int position) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String genotypeAsString(Chromosome chrom, Position positionObject) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String genotypeAsString(Chromosome chrom, int startSite, int endSite) {
        // TODO Auto-generated method stub
        return null;
    }
}