package net.maizegenetics.dna.pd;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.HapMapHDF5Constants;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class PDAnnotation {

    private static final String HAS_DATA = "HasData"; // summary index where if any trait has a value at that location, value is set to 1

    private String hapMapFile_prefix = "/maizeHapMapV2_B73RefGenV2_201203028_";
    private String hapMapFile_suffix = "h.hmp.txt.gz";
    private static final int PHYSICAL_POSITION_COLUMN = 0;
    private static final int MINOR_ALLELE_FREQUENCY_COLUMN = 1;
    private static final int COLUMN_OFFSET = 1; // one for physical position column and another for minor allele frequency column

    private int traitIndex = 0;
    private  int chrIndex = 1;
    private int physPosIndex = 2;
    private int resultIndex = 5;

    private String[] allTraits;
    private int[][] allPositions;
    private float[][] allResults;
    private int[][] featurePositions;
    private String[] allFeatures;



    public PDAnnotation(String hapMapPath, String pathGwas, String annoPath, String outputFile,
            int startChr, int endChr) {
        File aHapMapDir = new File(hapMapPath);
        File aGwasDir = new File(pathGwas);

        // Ed
//        loadGWAS(aGwasDir, outputFile, startChr);        // previous version - for comparison/testing

        // Dallas
        createAnnoHDF5WithHapMap(aHapMapDir, outputFile, startChr);
        loadGWAS(aGwasDir, outputFile, startChr, "\t");  // new version
//        File annoFile = new File(annoPath);
        //instatiate annotations in the HDF5 file - GENE (String), DistToGene (Integer), 
        //downstream_gene_variant, 3_prime_UTR_variant, missense_variant, synonymous_variant
//        loadAnnotationsByFile(annoFile, "\t");
    }



    public void createAnnoHDF5WithHapMap(File hapMapDir, String outputFile, int currChr) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(outputFile);
        config.overwrite();
        IHDF5Writer writer = config.writer();
        int[] hasData = null;  //recorder if there is any GWAS data for a site; TODO perhaps increment
        String chromosomeFile = hapMapDir + hapMapFile_prefix + "chr" + currChr + hapMapFile_suffix;
        System.out.println("Loading:" + chromosomeFile);
        GenotypeTable bna = ImportUtils.readFromHapmap(chromosomeFile, null /*progressListener*/);
        //System.out.printf("Sites:%d StartPosition:%d EndPosition:%d %n", bna.numberOfSites(), bna.chromosomalPosition(0), bna.chromosomalPosition(bna.numberOfSites() - 1));
        int siteCnt = bna.numberOfSites();
        int[] alignmentPhysPos = bna.physicalPositions();
        hasData = new int[siteCnt];
        float[] maf = new float[siteCnt];
        byte[] mjAllele = new byte[siteCnt];
        byte[] mnAllele = new byte[siteCnt];
        for (int j = 0; j < siteCnt; j++) {
            mjAllele[j] = bna.majorAlleleAsString(j).getBytes()[0];
            mnAllele[j] = bna.minorAlleleAsString(j).getBytes()[0];
            maf[j] = (float) bna.minorAlleleFrequency(j);
        }
        //write positions to hdf "pos"+chromosome
        String chrGroup = "chr" + currChr + "/";
        writer.object().createGroup(chrGroup);
        writer.int32().setAttr(chrGroup, HapMapHDF5Constants.NUM_SITES, siteCnt);
        writer.int32().createArray(chrGroup + HapMapHDF5Constants.POSITIONS, alignmentPhysPos.length);
        writer.writeIntArray(chrGroup + HapMapHDF5Constants.POSITIONS, alignmentPhysPos);

        //write alleles to hdf "allele"+chromosome
        // which version? String[][] ?
        writer.int8().createArray(chrGroup + HapMapHDF5Constants.MAJOR_ALLELE, mjAllele.length);
        writer.writeByteArray(chrGroup + HapMapHDF5Constants.MAJOR_ALLELE, mjAllele);
        writer.int8().createArray(chrGroup + HapMapHDF5Constants.MINOR_ALLELE, mnAllele.length);
        writer.writeByteArray(chrGroup + HapMapHDF5Constants.MINOR_ALLELE, mnAllele);
        // write minor allele frequencies
        writer.float32().createArray(chrGroup + HapMapHDF5Constants.MAF, maf.length);
        writer.float32().writeArray(chrGroup + HapMapHDF5Constants.MAF, maf);

        writer.object().createGroup(chrGroup + HapMapHDF5Constants.GWAS);
        writer.object().createGroup(chrGroup + HapMapHDF5Constants.GENOMIC);
        writer.object().createGroup(chrGroup + HapMapHDF5Constants.POP_GEN);

        writer.close();
    }

    /**
     * For the provided chromosome, writes out all positions and gwas results on a trait-by-trait basis
     *
     * @param gwasFileIn
     * @param outputFile
     * @param currChr
     * @param delimiter
     */
    private void loadGWAS(File gwasFileIn, String outputFile, int currChr, String delimiter){
        IHDF5Writer writer = HDF5Factory.open(outputFile);
        System.out.println("Loading GWAS by chromosome");
        //1. ArrayList<String> traitsInFile=getTraitListFromGWASInputFile();
        //2. Evaluate whether the traits already exist in the HDF5 file, if not create
        //3. Add GWAS results to the HDF5 file
        
        loadDataByChromosome(gwasFileIn, currChr,  delimiter);
        String chrGroup = "chr" + currChr + "/";

        int[] positions = writer.readIntArray(chrGroup + HapMapHDF5Constants.POSITIONS);

        for(int i =0; i < allTraits.length; i++){
            int posMatch = 0, posMisMatch = 0;
            float[] rmip = new float[positions.length];
            for(int j = 0; j < allPositions[i].length; j++) {
                if(allPositions[i][j]>3600000) continue;  // TODO: remove after testing

                // for the current traits, transfer result values to array
                int[] aInt = allPositions[i];
                int site = Arrays.binarySearch(positions, allPositions[i][j]);
                if (site < 0) {
                    System.out.println("Error Position not found:" + allPositions[i][j]);
                    posMisMatch++;
                } else {
                    posMatch++;
                    rmip[site] = allResults[i][j];
                    System.out.printf("Hit Chr:%d Trait:%s Position:%d site:%d rmip:%f %n ", currChr, allTraits[i], allPositions[i][j], site, allResults[i][j]);
                }
            }
            System.out.printf("Position matches:%d errors:%d %n", posMatch, posMisMatch);
            String dataSetName = chrGroup + HapMapHDF5Constants.GWAS + "/" + allTraits[i];
            writer.float32().createArray(dataSetName, rmip.length);
            writer.writeFloatArray(dataSetName, rmip);
        }
    }

    // Original - deprecated
    private void loadGWAS(File gwasFileIn, String outputFile, int currChr) {
        IHDF5Writer writer = HDF5Factory.open(outputFile);

        String[] traits =getGWASTraits(gwasFileIn, traitIndex, "\t");

        String chrGroup = "chr" + currChr + "/";
        //read in all chromosome position
        //create a method to hold this memory


        int[] positions = writer.readIntArray(chrGroup + HapMapHDF5Constants.POSITIONS);

        for (int j = 0; j < traits.length; j++) {
           
            System.out.println(gwasFileIn.toString());
            // pull out the physical location and p-value
            int posMatch = 0, posMisMatch = 0;
            float[] rmip = new float[positions.length];
            Arrays.fill(rmip, Float.NaN);
            try {
                BufferedReader fileIn = Utils.getBufferedReader(gwasFileIn, 1000000);
                String s;
                while ((s = fileIn.readLine()) != null) {
                    String[] fields = s.split("\t");
                    try {
                        int theChr = Integer.parseInt(fields[chrIndex]);
                        int position = Integer.parseInt(fields[physPosIndex]);
                        float rmipValue = Float.parseFloat(fields[resultIndex]);
                        if(theChr!=9) continue;
                        if(position>7600000) continue;
                        //int site = Arrays.binarySearch(allPositions[theChr-1], position);
                        int site = Arrays.binarySearch(positions, position); 
                        if (site < 0) {
                            System.out.println("Error Position not found:"+position);
                            posMisMatch++;
                        } else {
                            posMatch++;
                            rmip[site] = rmipValue;
                            System.out.printf("Hit Chr:%d Trait:%s Position:%d site:%d rmip:%d %n ",theChr, traits[j], position, site, rmipValue);
                        }

                    } catch (Exception e) {
                        //                     System.out.println("Header");
                    }
                }
            } catch (IOException e) {
                System.out.println("IOError");
                e.printStackTrace();
            }
            System.out.printf("Position matches:%d errors:%d %n", posMatch, posMisMatch);
            String dataSetName = chrGroup + HapMapHDF5Constants.GWAS + "/" + traits[j];
            writer.float32().createArray(dataSetName, rmip.length);
            writer.writeFloatArray(dataSetName, rmip);
        } // end of traits loop
    }

    // Only used in deprecated version of loadGWAS
    private String[] getGWASTraits(File gwasResults, int traitIndex, String delimiter){
        BufferedReader br = Utils.getBufferedReader(gwasResults, 1000000);
        String line = null;
        Set<String> aSet = new HashSet();
        try{
            while((line =  br.readLine()) != null){
                String[] fields = line.split(delimiter);

                aSet.add(fields[traitIndex]);
            }
        }catch(IOException ioe){
            ioe.printStackTrace();
        }
    
        String[] result = new String[aSet.size()];
        aSet.toArray(result);
        return result;
    }

    /**
     * For a given chromosome, loads all positions and results for all traits
     *
     * @param gwasFileIn
     * @param currChr
     * @param delimiter
     */
    private void loadDataByChromosome(File gwasFileIn, int currChr, String delimiter){
        BufferedReader br = Utils.getBufferedReader(gwasFileIn, 1000000);
        String line = null;

        Map<String, List> traitPosition = new HashMap<String, List>();
        Map<String, List> traitResult = new HashMap<String, List>();
        try{
            while(( line = br.readLine()) != null) {
                String[] fields = line.split(delimiter);

                // for the current chromosome and trait, hold the positions
                try {
                    int chromosome = Integer.parseInt(fields[chrIndex]);
                    if (currChr != chromosome) continue;

                    String aTrait = fields[traitIndex].trim();
                    if (traitPosition.containsKey(aTrait)) {
                        traitPosition.get(aTrait).add(fields[physPosIndex]);
                        List l = traitResult.get(aTrait);
                        l.add(fields[resultIndex]);
                        traitPosition.put(aTrait,l);
                    } else {
                        List<String> position = new ArrayList();
                        position.add(fields[physPosIndex]);
                        
                        Object o=traitPosition.put(aTrait, position);
                        System.out.printf("%s %s %d %n", o, aTrait, position);
                        List<String> result = new ArrayList();
                        result.add(fields[resultIndex]);
                        traitResult.put(fields[traitIndex], result);
                    }
                } catch (Exception e) {
                    System.out.println("Header");
                }
            }
        }catch(IOException ioe){
            ioe.printStackTrace();
        }

        allTraits = new String[traitPosition.size()];

        // create a two-dimensional array of positions for each trait
        // first dimension index is shared with allTraits
        traitPosition.keySet().toArray(allTraits);

        int traitCount = allTraits.length;
        allPositions = new int[traitCount][];
        allResults = new float[traitCount][];
        for(int i = 0; i < traitCount; i++){
            List posList = traitPosition.get(allTraits[i]);
            int posCount = posList.size();
            allPositions[i] = new int[posCount];
            Iterator posIt = posList.iterator();
            int count = 0;
            while(posIt.hasNext()){
                allPositions[i][count++] = Integer.parseInt((String)posIt.next());
            }

            List resList = traitResult.get(allTraits[i]);
            int resCount = resList.size();
            allResults[i] = new float[resCount];
            Iterator resIt = resList.iterator();
            count = 0;
            while(resIt.hasNext()){
                allResults[i][count++] = Float.parseFloat((String)resIt.next());
            }
        }
        System.out.println("Breakpoint");
    }

    private void loadAnnotations(){

    }

    // annotations files have been organized by chromosome
    private void loadAnnotationsByFile(File annoFile, String delimiter){

        int snpIdIndex = 0;
        int locationIndex = 1;      // location is  specified as <chr>:<position>, e.g., 9:30
                                    // TODO: how to handle range locations? e.g., 9:513883-513884
        String locationDelimiter = ":";
        int featureIndex = 6;     // Feature_typeConsequence

        Map<String, List> featurePosition = new HashMap<String, List>();

        BufferedReader br = Utils.getBufferedReader(annoFile, 1000000);
        String line = null;
        try {
            while ((line = br.readLine()) != null) {
                String[] fields = line.split(delimiter);

                try {
                    String aLoc =  fields[locationIndex];
                    int aPosition = getPosition(aLoc, locationDelimiter);
                    String feature = fields[featureIndex].trim();
                    if (featurePosition.containsKey(feature)) {
                        List l = featurePosition.get(feature);
                        l.add(aPosition);
                        featurePosition.put(feature, l);
                    } else {
                        List<Integer> l = new ArrayList<Integer>();
                        l.add(aPosition);
                        featurePosition.put(feature, l);
                    }

                } catch (Exception e) {
                    System.out.println("Header");
                }
            }
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        // convert to a two dimensional array
        int featureCount = featurePosition.size();
        featurePositions = new int[featureCount][];
        allFeatures = new String[featureCount];
        featurePosition.keySet().toArray(allFeatures);
        for(int i = 0; i < featureCount; i++){
            List posList = featurePosition.get(allFeatures[i]);
            int posCount = posList.size();
            featurePositions[i] = new int[posCount];
            Iterator<Integer> iterator = posList.iterator();
            for(int j = 0; j < posCount; j++){
                featurePositions[i][j] = iterator.next().intValue();
            }
        }
    }


    private int getPosition(String location, String delimiter){
        String[] fields = location.split(delimiter);
        int position = Integer.parseInt(fields[1]);
        return position;
    }


    public static void main(String[] args) {
        // Dallas
//        String hapMapPath = "C:\\Documents and Settings\\dkroon\\My Documents\\PD\\HapMap\\compressed";
//        String pathGwas = "C:\\Documents and Settings\\dkroon\\My Documents\\PD\\gwas\\gwas_hits_all.txt";
//        String PDfile = "C:\\Documents and Settings\\dkroon\\My Documents\\PD\\out\\testPD.h5";
//        String annoPath = "C:\\Documents and Settings\\dkroon\\My Documents\\PD\\Annotations\\20130522_SnpAnnotations_FromJason\\maizeHapMapV2_B73RefGenV2_201203028_chr9h.WorstPerSnp.vcf";

        // Ed
         String hapMapPath = "/Volumes/LaCie/HapMapV2/compressed/";
        String pathGwas = "/Volumes/LaCie/PolymorphismDescriptors/gwas_hits_all.txt";
        String PDfile = "/Volumes/LaCie/PolymorphismDescriptors/XtestPD.h5";
        String annoPath = "/Volumes/LaCie/PolymorphismDescriptors/maizeHapMapV2_B73RefGenV2_201203028_chr9h.WorstPerSnp.vcf";

        PDAnnotation p = new PDAnnotation(hapMapPath, pathGwas, annoPath, PDfile, 9, 9);

    }
}
