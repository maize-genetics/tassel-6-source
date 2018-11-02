/*
 * ExportUtils
 */
package net.maizegenetics.dna.snp;

import com.google.common.base.Joiner;
import com.google.common.collect.Multimap;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.io.VCFUtil;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.*;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;

/**
 * Exports Genotype Tables to various file formats.
 *
 * @author Jon Zhang
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class ExportUtils {

    private static final Logger myLogger = Logger.getLogger(ExportUtils.class);
    private static FormattedOutput format = FormattedOutput.getInstance();

    private ExportUtils() {
        // Utility Class - do not instantiate.
    }

    public static String writeGenotypeHDF5(GenotypeTable a, String newHDF5file) {
        return writeGenotypeHDF5(a, newHDF5file, null, true);
    }

    public static String writeGenotypeHDF5(GenotypeTable a, String newHDF5file, boolean keepDepth) {
        return writeGenotypeHDF5(a, newHDF5file, null, keepDepth);
    }

    /**
     * Exports a alignment into the Byte HDF5 format.
     *
     * @param a alignment to be exported
     * @param newHDF5file filename for the new file (should end with "hmp.h5")
     * @param exportTaxa subset of taxa (if null exports ALL taxa)
     * @return
     */
    public static String writeGenotypeHDF5(GenotypeTable a, String newHDF5file, TaxaList exportTaxa, boolean keepDepth) {
        GenotypeTableBuilder aB = GenotypeTableBuilder.getTaxaIncremental(a.positions(), newHDF5file);
        if ((exportTaxa != null) && (exportTaxa.numberOfTaxa() == 0)) {
            aB.build();
            return newHDF5file;
        }
        for (int t = 0; t < a.numberOfTaxa(); t++) {
            if ((exportTaxa != null) && (!exportTaxa.contains(a.taxa().get(t)))) {
                continue;  //taxon not in export list
            }
            byte[] bases = a.genotypeAllSites(t);
            if (a.hasDepth() == false || keepDepth == false) {
                aB.addTaxon(a.taxa().get(t), bases, null);
            } else {
                aB.addTaxon(a.taxa().get(t), bases, a.depth().valuesForTaxonByte(t));
            }
        }
        aB.build();
        return newHDF5file;
    }

    /**
     * Write a GenotypeTable to HapMap format with standard settings - unphased
     * single character, tab delimiter, and no progress tracking.
     *
     * @param alignment genotype table
     * @param filename outfile name (will add ".hmp.txt" if needed)
     * @return name of the outfile with the appropriate suffix
     */
    public static String writeToHapmap(GenotypeTable alignment, String filename) {
        return writeToHapmap(alignment, false, filename, '\t', null);
    }

    /**
     * Write a GenotypeTable to HapMap format.
     *
     * @param alignment genotype table
     * @param diploid true uses phased two letter encoding, false one letter
     * unphased
     * @param filename outfile name (will add ".hmp.txt" if needed)
     * @param delimChar delimiter character normally tab
     * @param listener progress listener, (null if unneeded)
     * @return name of the outfile with the appropriate suffix
     */
    public static String writeToHapmap(GenotypeTable alignment, boolean diploid, String filename, char delimChar, ProgressListener listener) {
        return writeToHapmap(alignment, diploid, filename, delimChar, true, listener);
    }

    public static String writeToHapmap(GenotypeTable alignment, boolean diploid, String filename, char delimChar, boolean includeTaxaAnnotations, ProgressListener listener) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(filename, ".hmp.txt", new String[]{".hmp.txt", ".hmp.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
            if (includeTaxaAnnotations) {
                for (Taxon taxon : alignment.taxa()) {
                    GeneralAnnotation annotation = taxon.getAnnotation();
                    if ((annotation == null) || (annotation.numAnnotations() == 0)) {
                        continue;
                    }
                    bw.write("##SAMPLE=" + taxon.toStringWithVCFAnnotation() + "\n");
                }
            }
            bw.write(Joiner.on(delimChar).join("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID",
                    "assayLSID", "panelLSID", "QCcode"));
            bw.write(delimChar);
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String sequenceID = alignment.taxaName(taxa).trim();
                bw.write(sequenceID);
                if (taxa != numTaxa - 1) {
                    bw.write(delimChar);
                }
            }
            bw.write("\n");
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                bw.write(alignment.siteName(site));
                bw.write(delimChar);
                byte[] genotypes = alignment.genotypeAllTaxa(site);
                // which alleles are present among the genotypes
                int[][] sortedAlleles = AlleleFreqCache.allelesSortedByFrequencyNucleotide(genotypes);
                int numAlleles = sortedAlleles[0].length;
                if (numAlleles == 0) {
                    bw.write("NA"); //if data does not exist
                } else if (numAlleles == 1) {
                    bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                } else {
                    bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][0]));
                    for (int allele = 1; allele < sortedAlleles[0].length; allele++) {
                        if (sortedAlleles[0][allele] != GenotypeTable.UNKNOWN_ALLELE) {
                            bw.write('/');
                            bw.write(alignment.genotypeAsString(site, (byte) sortedAlleles[0][allele]));  // will write out a third allele if it exists
                        }
                    }
                }
                bw.write(delimChar);
                bw.write(Joiner.on(delimChar).join(alignment.chromosomeName(site), String.valueOf(alignment.chromosomalPosition(site)),
                        "+", "NA", "NA", "NA", "NA", "NA", "NA"));
                bw.write(delimChar);
                for (int taxa = 0; taxa < numTaxa; taxa++) {
                    if (diploid == false) {
                        String baseIUPAC = null;
                        try {
                            baseIUPAC = alignment.diploidAsString(site, genotypes[taxa]);
                        } catch (Exception e) {
                            String[] b = alignment.genotypeAsStringArray(taxa, site);
                            bw.close();
                            myLogger.debug(e.getMessage(), e);
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                        }
                        if ((baseIUPAC == null) || baseIUPAC.equals("?")) {
                            String[] b = alignment.genotypeAsStringArray(taxa, site);
                            bw.close();
                            throw new IllegalArgumentException("There is no String representation for diploid values: " + b[0] + ":" + b[1] + " getBase(): 0x" + Integer.toHexString(alignment.genotype(taxa, site)) + "\nTry Exporting as Diploid Values.");
                        }
                        bw.write(baseIUPAC);
                    } else {
                        byte[] temp = GenotypeTableUtils.getDiploidValues(genotypes[taxa]);
                        bw.write(alignment.genotypeAsString(site, temp[0]));
                        bw.write(alignment.genotypeAsString(site, temp[1]));
                    }
                    if (taxa != (numTaxa - 1)) {
                        bw.write(delimChar);
                    }
                }
                bw.write("\n");

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }
            return fullFileName;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalArgumentException("Error writing Hapmap file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Writes given alignment to a VCF file
     *
     * @param gt
     * @param filename
     * @return
     */
    public static String writeToVCF(GenotypeTable gt, String filename, boolean keepDepth) {
        return writeToVCF(gt,filename,keepDepth,null);
    }
    public static String writeToVCF(GenotypeTable gt, String filename, boolean keepDepth, ProgressListener listener) {
        final char delimChar = '\t';
        boolean hasDepth = gt.hasDepth() && keepDepth;
        try {

            filename = Utils.addSuffixIfNeeded(filename, ".vcf", new String[]{".vcf", ".vcf.gz"});
            BufferedWriter bw = Utils.getBufferedWriter(filename);
            bw.write("##fileformat=VCFv4.0");
            bw.newLine();
            if (!gt.hasReference()) {
                bw.write("##Tassel=<ID=GenotypeTable,Version=5,Description=\"Reference allele is not known. The major allele was used as reference allele\">");
                bw.newLine();
            }
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (only filtered reads used for calling)\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
            bw.newLine();
            bw.write("##FORMAT=<ID=PL,Number=.,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic\">");
            bw.newLine();
            bw.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
            bw.newLine();
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
            bw.newLine();
            bw.write("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">");
            bw.newLine();
            writeVCFSampleAnnotationToWriter(gt, bw);
            bw.write("#CHROM" + delimChar + "POS" + delimChar + "ID" + delimChar + "REF" + delimChar + "ALT" + delimChar + "QUAL" + delimChar + "FILTER" + delimChar + "INFO" + delimChar + "FORMAT");
            for (int taxa = 0; taxa < gt.numberOfTaxa(); taxa++) {
                String taxonName = gt.taxaName(taxa).trim();
                bw.write(delimChar + taxonName);
            }
            bw.newLine();

            int noAlleles = 0;
            for (int site = 0; site < gt.numberOfSites(); site++) {
                int numSites = gt.numberOfSites();
                Position p = gt.positions().get(site);
                String[] knownVariants = p.getKnownVariants();
                byte refAllele = p.getAllele(WHICH_ALLELE.Reference);
                int[] sortedAlleles = gt.allelesSortedByFrequency(site)[0]; // which alleles are actually present among the genotypes



                //ZRM22 March 18 2016 move to add reference into sortedAlleles array if its missing
                int[] sortedAllelesTemp = VCFUtil.resolveRefSorted(sortedAlleles, refAllele);

                //ZRM22 June 6 2016 fix variants with ref allele


                sortedAlleles = sortedAllelesTemp;
                //ZRM22 Jan7 Remake
                //If knownVariants.length is greater than 0 its either from a VCF file or Hapmap
                if(knownVariants.length>0) {

                    //ReOrder based on variant alleles
                    //Store a tempSortedAlleles so we can appropriately handle hapmap to vcf
                    //int[] tempSortedAlleles = new int[knownVariants.length];

                    //ArrayList to hold the Sorted Alleles Indices Temporarily as the ordering will change
                    ArrayList<Integer> tempSortedAlleles = new ArrayList<Integer>();

                    //Loop through all the knownVariants and check to see if we have an indel
                    boolean knownVariantIndel = VCFUtil.indelInKnownVariant(knownVariants);


                    //If we do have an indel, we can add the variants after picking off the first character to the tempSortedAlleles
                    if(knownVariantIndel) {
                        //Loop through the variants
                        for(int i = 0; i < knownVariants.length; i++) {
                            //Pull off the first character if it exists
                            if(knownVariants[i].length()>1) {
                                String parsedVariant = knownVariants[i].substring(1);
                                tempSortedAlleles.add((int)NucleotideAlignmentConstants.getNucleotideAlleleByte(parsedVariant.charAt(0)));
                            }
                            else {
                                //Mark as deletion
                                tempSortedAlleles.add((int)NucleotideAlignmentConstants.getNucleotideAlleleByte('-'));
                            }
                        }
                    } else {
                        //If we dont have an indel, we can add it to the allele array
                        if(sortedAlleles.length<knownVariants.length){
                            //Clear it out, we probably dont need to do this
                            tempSortedAlleles = new ArrayList<Integer>();
                        }
                        int nIndex = -1;
                        for(int i = 0; i<knownVariants.length; i++) {
                            //ZRM22 Mar 22
                            if(knownVariants[i].charAt(0)!='N') {
                                tempSortedAlleles.add((int)NucleotideAlignmentConstants.getNucleotideAlleleByte(knownVariants[i].charAt(0)));
                            }
                            else {
                                //If N is in our known Variants list but we do not have an indel, we need to remove it
                                nIndex = i;
                            }
                        }
                        if(nIndex != -1) {
                            //if we have an N we need to resize KnownVariants
                            String[] knownVariantsSmall = new String[knownVariants.length-1];
                            for(int i = 0; i<knownVariants.length; i++) {
                                if(i < nIndex) {
                                    knownVariantsSmall[i] = knownVariants[i];
                                }
                                else if(i > nIndex) {
                                    knownVariantsSmall[i-1] = knownVariants[i];
                                }
                            }
                            knownVariants = knownVariantsSmall;
                        }
                    }
                    //END ZRM22 Jan7

                    //Make a copy of KnownVaraints in case we need to add some
                    ArrayList<String> knownVariantsList = new ArrayList<String>();
                    boolean indelsExist = false;
                    boolean indelsInKnownVariants = VCFUtil.indelInKnownVariant(knownVariants);
                    if(indelsInKnownVariants) {
                        indelsExist = true;
                    }

                    //Go through sorted alleles and also check for indels
                    for(int i = 0 ;i<sortedAlleles.length; i++) {
                        if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]).equals("-")) {
                            indelsExist = true;
                        }
                    }
                    //Move To Function/

                    for(String variant:knownVariants) {
                        if(indelsExist && !indelsInKnownVariants) {
                            knownVariantsList.add("N"+variant);
                        }
                        else {
                            knownVariantsList.add(variant);
                        }
                    }
                    //ZRM Jun6 fix to force Ref annotated alleles to stay in REF for export
                    //Need to reorder the variants based on the original sorting
                    ArrayList<Integer> sortedAllelesList = new ArrayList<Integer>();
                    HashMap<Integer,String> sortedAlleleKnownVariantMap = new HashMap<Integer, String>();
                    for(int i = 0; i < sortedAlleles.length; i++) {
                        //Add it to the new sorted list
                        sortedAllelesList.add(sortedAlleles[i]);
                        if(!tempSortedAlleles.contains(sortedAlleles[i])) {
                            //Check for an indel
                            if(indelsExist) {
                                if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]).equals("-")) {
                                    //Add an Entry to the sortedAllele, knownVariant mapping
                                    sortedAlleleKnownVariantMap.put(sortedAlleles[i],"N");
                                }
                                else {
                                    sortedAlleleKnownVariantMap.put(sortedAlleles[i],NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
                                }
                            }
                            else {
                                sortedAlleleKnownVariantMap.put(sortedAlleles[i],NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
                            }
                        }
                        else {
                            //Find the index in tempSortedAlleles
                            int variantIndex = tempSortedAlleles.indexOf(sortedAlleles[i]);
                            //Use it to get the correct KnownVariant
                            sortedAlleleKnownVariantMap.put(sortedAlleles[i],knownVariants[variantIndex]);
                        }
                    }
                    //loop through tempSortedAlleles and make sure we have them all
                    //Else add to the end
                    for(int i = 0; i < tempSortedAlleles.size(); i++) {
                        if(!sortedAllelesList.contains(tempSortedAlleles.get(i))) {
                            sortedAllelesList.add(tempSortedAlleles.get(i));
                            sortedAlleleKnownVariantMap.put(tempSortedAlleles.get(i),knownVariantsList.get(i));
                        }
                    }
                    int[] sortedAllelesExtended = new int[sortedAllelesList.size()];
                    for(int i = 0; i < sortedAllelesExtended.length; i++) {
                        sortedAllelesExtended[i] = sortedAllelesList.get(i);
                    }
                    sortedAlleles = sortedAllelesExtended;

                    String[] knownVariantsExtended = new String[sortedAllelesList.size()];
                    for(int i = 0; i < knownVariantsExtended.length; i++) {
                        knownVariantsExtended[i] = sortedAlleleKnownVariantMap.get(sortedAllelesList.get(i));
                    }
                    knownVariants = knownVariantsExtended;
                    //TODO Cleanup
//                    //Go through sorted alleles
//                    for(int i = 0 ;i<sortedAlleles.length; i++) {
//                    //If a sorted allele is not in tempSortedAlleles,
//                        if(!tempSortedAlleles.contains(sortedAlleles[i])) {
//                            //if its not add it to sorted alleles and knownVariants
//                            tempSortedAlleles.add(sortedAlleles[i]);
//                            //Check for an indel
//                            if(indelsExist) {
//                                if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]).equals("-")) {
//                                    knownVariantsList.add("N");
//                                }
//                                else {
//                                    knownVariantsList.add("N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
//                                }
////                                knownVariantsList.add("N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
//                            }
//                            else {
//                                knownVariantsList.add(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
//                            }
//                        }
//                    }
//                    //reset knownVariants and sortedAlleles to reflect the changes
//                    String[] knownVariantsExtended = new String[knownVariantsList.size()];
//                    for(int i = 0; i < knownVariantsExtended.length; i++) {
//                        knownVariantsExtended[i] = knownVariantsList.get(i);
//                    }
//                    knownVariants = knownVariantsExtended;
//
//                    int[] sortedAllelesExtended = new int[tempSortedAlleles.size()];
//                    for(int i = 0; i < sortedAllelesExtended.length; i++) {
//                        sortedAllelesExtended[i] = tempSortedAlleles.get(i);
//                    }
//                    sortedAlleles = sortedAllelesExtended;
//                    //sortedAlleles = tempSortedAlleles.toArray(new int[tempSortedAlleles.size()]);
                }
                else {
                    //No known variants, but we need to handle indels
                    int indelIndex = -1;
                    //loop through sorted alleles
                    for(int i = 0; i<sortedAlleles.length; i++) {
                        //if we find an indel mark the index and set a boolean
                        if(sortedAlleles[i] == (int)NucleotideAlignmentConstants.getNucleotideAlleleByte('-')) {
                            indelIndex = i;
                            break;
                        }
                    }

                    knownVariants = new String[sortedAlleles.length];
                    for(int i = 0; i<knownVariants.length; i++) {
                        if(indelIndex==-1) {
                            knownVariants[i] = ""+ NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]);
                        }
                        else {
                            if(indelIndex == i) {
                                knownVariants[i] = "N";
                            }
                            else {
                                knownVariants[i] = "N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]);
                            }
                        }
                    }

                }

                int nAlleles = sortedAlleles.length;
                HashMap<String,Integer> alleleRedirectMap = new HashMap<String,Integer>();
                String[] alleleRedirect = new String[16];
                Arrays.fill(alleleRedirect, ".");
                for (int i = 0; i < sortedAlleles.length; i++) {
                    alleleRedirect[sortedAlleles[i]] = "" + i;
                    alleleRedirectMap.put(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]), i);
                }

                bw.write(gt.chromosomeName(site)); // chromosome
                bw.write(delimChar);
                bw.write(gt.chromosomalPosition(site) + ""); // position
                bw.write(delimChar);
                bw.write(gt.siteName(site)); // site name
                bw.write(delimChar);
                if (nAlleles == 0) {                                                  //used to be ==0
                    //System.out.println("A0:"+gt.chromosomeName(site)+":"+gt.chromosomalPosition(site));
                    noAlleles++;
                    bw.write(".\t.\t.\tPASS\t.\tGT");
                    for (int taxa = 0; taxa < gt.numberOfTaxa(); taxa++) {
                        bw.write("\t./.");
                    }
                    bw.newLine();
                    continue;
                }
                //bw.write(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[0])); // ref allele
                //Fix for indels 8_27
                if(knownVariants.length==0) {
                    bw.write(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[0])); // ref allele
                }
                else {
                    bw.write(knownVariants[0]);
                }
                bw.write(delimChar);

                StringBuilder altAllelesBuilder = new StringBuilder("");

                //ZRM 8_27
                String altString = "";
                int indelIndex = -1;

                if(knownVariants.length==0 || knownVariants.length<sortedAlleles.length) {
                    ArrayList<String> altAlleles = new ArrayList<String>();
                    for(int aa = 1; aa<sortedAlleles.length; aa++) {
                        //Ramu Fix
                        //altAlleles.add(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[aa]));
                        //UNCOMMENT BEFORE COMMIT
                        if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[aa]) != "-") {
                            altAlleles.add(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[aa]));
                        }
                        else {
                            indelIndex = aa;
                        }
                    }
                    altString = altAlleles.stream().collect(Collectors.joining(","));
                }
                else {
                    altString = Arrays.stream(knownVariants, 1, knownVariants.length).collect(Collectors.joining(","));
                }

                if(altString.length()==0) {
                    altString = ".";
                }

                ////bw.write(altAllelesBuilder.toString()); // alt alleles
                bw.write(altString);
                bw.write(delimChar);

                bw.write("."); // qual score
                bw.write(delimChar);

                bw.write("PASS"); // filter
                bw.write(delimChar);

                //INFO
                GeneralAnnotation ga = p.getAnnotation();
                String annotationHolder=ga.getAnnotationKeys().stream().sorted()
                        .filter(k->!k.equals("VARIANT"))
                        .filter(k -> !(k.equals("DP") && hasDepth)) //Get rid of the DP tag if we have depth in the genotype table already.  
                        .map(key->{
                            String[] annos=ga.getTextAnnotation(key);
                            if(annos[0].equals("TRUE")) return key;
                            return key+Arrays.stream(annos).collect(Collectors.joining(",","=",""));
                        })
                        .collect(Collectors.joining(";"));
                if (hasDepth) {
                    //bw.write("DP=" + gt.depth().depthForSite(site)); // DP
                    //To Fix bug where ";DP=100" string would occur
                    if(annotationHolder.equals("")) {
                        annotationHolder += "DP=" + gt.depth().depthForSite(site);
                    }
                    else if(annotationHolder.equals(".")) {
                        //Fix bug where we have .;DP=100 showing up.
                        annotationHolder = "DP="+ gt.depth().depthForSite(site);
                    }
                    else {
                        annotationHolder += ";DP=" + gt.depth().depthForSite(site);
                    }
                }
                if(!annotationHolder.equals("")) {
                    bw.write(annotationHolder);
                }
                else {
                    bw.write("."); // DP
                }
                bw.write(delimChar);

                if (hasDepth) {
                    bw.write("GT:AD:DP:GQ:PL");
                } else {
                    bw.write("GT");
                }
                for (int taxa = 0; taxa < gt.numberOfTaxa(); taxa++) {
                    bw.write(delimChar);
                    // GT = genotype
                    byte[] values = gt.genotypeArray(taxa, site);
                    if(knownVariants.length>0) {
                        bw.write(alleleRedirect[values[0]]+"/"+alleleRedirect[values[1]]);
                    }
                    else {
                        //Ramu Fix
//                        if(alleleRedirect[values[0]].equals(".")) {
//                            bw.write(alleleRedirect[values[0]]);
//                        }
//                        else {
//                            if(indelIndex != -1 && Integer.parseInt(alleleRedirect[values[0]]) > indelIndex) {
//                                bw.write(""+(Integer.parseInt(alleleRedirect[values[0]]) -1 ));
//                            }
//                            else {
//                                bw.write(alleleRedirect[values[0]]);
//                            }
//                        }

                        //handle if no Known Variants(from a different file type)
                        if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)values[0]).equals("-")) {
                            //TODO handle Missing better
                            bw.write(".");
                        }
                        else {
                            if(alleleRedirect[values[0]].equals(".")) {
                                bw.write(alleleRedirect[values[0]]);
                            }
                            else {
                                if(indelIndex != -1 && Integer.parseInt(alleleRedirect[values[0]]) > indelIndex) {
                                    bw.write(""+(Integer.parseInt(alleleRedirect[values[0]]) -1 ));
                                }
                                else {
                                    bw.write(alleleRedirect[values[0]]);
                                }
                            }
                        }

                        bw.write("/");

                        //Ramu Fix
//                        if(alleleRedirect[values[1]].equals(".")) {
//                            bw.write(alleleRedirect[values[1]]);
//                        }
//                        else {
//                            if(indelIndex != -1 && Integer.parseInt(alleleRedirect[values[1]]) > indelIndex) {
//                                bw.write(""+(Integer.parseInt(alleleRedirect[values[1]]) -1));
//                            }
//                            else {
//                                bw.write(alleleRedirect[values[1]]);
//                            }
//                        }

                        if(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)values[1]).equals("-")) {
                            bw.write(".");
                        }
                        else {
                            if(alleleRedirect[values[1]].equals(".")) {
                                bw.write(alleleRedirect[values[1]]);
                            }
                            else {
                                if(indelIndex != -1 && Integer.parseInt(alleleRedirect[values[1]]) > indelIndex) {
                                    bw.write(""+(Integer.parseInt(alleleRedirect[values[1]]) -1));
                                }
                                else {
                                    bw.write(alleleRedirect[values[1]]);
                                }
                            }
                        }
                    }
                    if (!(hasDepth)) {
                        continue;
                    }
                    bw.write(":");

                    // AD
                    int[] siteAlleleDepths = gt.depthForAlleles(taxa, site);
                    int siteTotalDepth = 0;

                    ArrayList<Integer> depthsList = new ArrayList<Integer>();
                    //Fix missing commas in depth information
                    for(int ss = 0; ss < sortedAlleles.length; ss++) {
                        if(ss!=indelIndex && sortedAlleles[ss]<siteAlleleDepths.length) {
                            try {
                                depthsList.add(siteAlleleDepths[sortedAlleles[ss]]);
                                siteTotalDepth += siteAlleleDepths[sortedAlleles[ss]];
                            }
                            catch(Exception e) {
                                System.out.println(Arrays.toString(alleleRedirect));
                                System.out.println(altString);
                                System.out.println(Arrays.toString(siteAlleleDepths));
                                System.out.println(Arrays.toString(sortedAlleles));
                                System.out.println(ss);
                                throw e;
                            }
                            //TODO Cleanup
//                            depthsList.add(AlleleDepthUtil.depthByteToInt((byte)siteAlleleDepths[sortedAlleles[ss]]));
//                            siteTotalDepth += AlleleDepthUtil.depthByteToInt((byte)siteAlleleDepths[sortedAlleles[ss]]);
                        }
                    }

                    bw.write(depthsList.stream().map((depth)->""+depth).collect(Collectors.joining(",")));
//
//                    for (int ss = 0; ss < sortedAlleles.length; ss++) {
//                        //bw.write("" + AlleleDepthUtil.decode(siteAlleleDepths[sortedAlleles[ss]]));
//                        if(ss!=indelIndex) {
//                            bw.write("" + siteAlleleDepths[sortedAlleles[ss]]);
//                            if (ss < sortedAlleles.length - 1 && ss+1!=indelIndex) {
//                                bw.write(',');
//                            }
//                            siteTotalDepth += siteAlleleDepths[sortedAlleles[ss]];
//                        }
//
//                    }
                    bw.write(":");
                    // DP
                    bw.write(siteTotalDepth + "");

                    int[] scores = new int[]{-1, -1, -1, -1};
                    if (values[0] != GenotypeTable.UNKNOWN_ALLELE) {
                        int altDepth = (sortedAlleles.length < 2 || sortedAlleles[1]>=siteAlleleDepths.length) ? 0 : siteAlleleDepths[sortedAlleles[1]];
                        altDepth = (altDepth<0) ? 0 : altDepth;
                        //int refDepth = (siteAlleleDepths[sortedAlleles[0]]==-1) ? 0 : siteAlleleDepths[sortedAlleles[0]];

                        //Check to see if either the major or alt allele has depth
                        if(siteAlleleDepths[sortedAlleles[0]] >= 0 && altDepth >= 0) {
                            scores = VCFUtil.getScore(siteAlleleDepths[sortedAlleles[0]], altDepth);
                            bw.write(":");
                            // GQ
                            bw.write(scores[3] + "");
                            bw.write(":");
                            // PL
                            int k = sortedAlleles.length - 1;
                            int[] fullPL = new int[(k * (k+1)/2)+k+1];


                            //Set all the values to 255 as Higher PL means its less likely to be correct
                            //Zero PL means the probability of error is 0
                            Arrays.fill(fullPL,255);

                            //Leaving these indicies in expanded form so we know its correct
                            //it should really just be in positions 0,1, and 2 regardless of number of sites
                            //(k*(k+1)/2)+j
                            //If we only have 1 allele we should only have 1 likelihood
                            if(fullPL.length==1) {
                                fullPL[0] = scores[0];
                            }
                            else {
                                fullPL[(0 * (0 + 1)/2) + 0] = scores[0];
                                fullPL[(1 * (1 + 1)/2) + 0] = scores[1];
                                fullPL[(1 * (1 + 1)/2) + 1] = scores[2];
                            }
                            for(int i = 0; i < fullPL.length-1; i++) {
                                bw.write(fullPL[i] + ",");
                            }
                            bw.write(""+fullPL[fullPL.length-1]);

//                            //Leaving these indicies in expanded form so we know its correct
//                            //it should really just be in positions 0,1, and 2 regardless of number of sites
//                            //(k*(k+1)/2)+j
//                            fullPL[(0 * (0 + 1)/2) + 0] = scores[0];
//                            fullPL[(1 * (1 + 1)/2) + 0] = scores[1];
//                            fullPL[(1 * (1 + 1)/2) + 1] = scores[2];
//                            for(int i = 0; i < fullPL.length-1; i++) {
//                                bw.write(fullPL[i] + ",");
//                            }
//                            bw.write(""+fullPL[fullPL.length-1]);
//                            //
//                            //bw.write(scores[0] + "," + scores[1] + "," + scores[2]);
                        }
                    }
//                    else {
//                        //If unknown just write out :0:0:0,0,0
//                        bw.write(":0:0,0,0");
//                    }
                }
                bw.newLine();
                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }
            if (noAlleles > 0) {
                myLogger.warn("Warning: " + noAlleles + " sites have no alleles.");
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing VCF file: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        }
        return filename;
    }

    private static void writeVCFSampleAnnotationToWriter(GenotypeTable gt, BufferedWriter bw) throws IOException {
        for (Taxon taxon : gt.taxa()) {
            GeneralAnnotation annotation = taxon.getAnnotation();
            if ((annotation == null) || (annotation.numAnnotations() == 0)) {
                continue;
            }
            Multimap annoMap = taxon.getAnnotation().getAnnotationAsMap();
            String annoString = Joiner.on(',').withKeyValueSeparator("=").join(annoMap.entries());
            bw.write("##SAMPLE=<ID=" + taxon.getName() + "," + annoString + ">");
            bw.newLine();
        }
    }

    /**
     * Writes given set of alignments to a set of Plink files
     *
     * @param alignment
     * @param filename
     * @param delimChar
     */
    public static String writeToPlink(GenotypeTable alignment, String filename, char delimChar) {
        if (delimChar != ' ' && delimChar != '\t') {
            throw new IllegalArgumentException("Delimiter charater must be either a blank space or a tab.");
        }

        BufferedWriter MAPbw = null;
        BufferedWriter PEDbw = null;
        String mapFileName = Utils.addSuffixIfNeeded(filename, ".plk.map");
        String pedFileName = Utils.addSuffixIfNeeded(filename, ".plk.ped");
        try {
            MAPbw = new BufferedWriter(new FileWriter(mapFileName), 1000000);
            int numSites = alignment.numberOfSites();
            for (int site = 0; site < numSites; site++) {
                MAPbw.write(alignment.chromosomeName(site)); // chromosome name
                MAPbw.write(delimChar);
                MAPbw.write(alignment.siteName(site)); // rs#
                MAPbw.write(delimChar);
                MAPbw.write("-9"); // genetic distance unavailable
                MAPbw.write(delimChar);
                MAPbw.write(Integer.toString(alignment.chromosomalPosition(site))); // position
                MAPbw.write("\n");
            }
            MAPbw.close();

            PEDbw = new BufferedWriter(new FileWriter(pedFileName), 1000000);
            // Compiled : Pattern
            Pattern splitter = Pattern.compile(":");
            int numTaxa = alignment.numberOfTaxa();
            for (int taxa = 0; taxa < numTaxa; taxa++) {
                String[] name = splitter.split(alignment.taxaName(taxa).trim());
                if (name.length != 1) {
                    PEDbw.write(name[1]); // namelvl 1 if is available
                } else {
                    PEDbw.write("-9");
                }
                PEDbw.write(delimChar);
                PEDbw.write(alignment.taxaName(taxa).trim()); // namelvl 0
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // paternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // maternal ID unavailable
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // gender is both
                PEDbw.write(delimChar);
                PEDbw.write("-9"); // phenotype unavailable, might have to change to "-9" for missing affection status
                PEDbw.write(delimChar);
                for (int site = 0; site < numSites; site++) {
                    String[] b = getSNPValueForPlink(alignment.genotypeAsStringArray(taxa, site));
                    PEDbw.write(b[0]);
                    PEDbw.write(delimChar);
                    PEDbw.write(b[b.length - 1]);
                    if (site != numSites - 1) {
                        PEDbw.write(delimChar);
                    }
                }
                PEDbw.write("\n");
            }
            PEDbw.close();
            return mapFileName + " and " + pedFileName;
        } catch (Exception e) {
            myLogger.error("Error writing Plink files: " + mapFileName + " and " + pedFileName + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Plink files: " + mapFileName + " and " + pedFileName + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                PEDbw.close();
            } catch (Exception e) {
                // do nothing
            }
            try {
                MAPbw.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    private static String[] getSNPValueForPlink(String[] base) {
        for (int i = 0; i < base.length; i++) {
            if (base[i].equals("N")) {
                base[i] = "0";
            } else if (base[i].equals("0")) {
                base[i] = "D";
            }
        }
        return base;
    }

    public static String saveDelimitedAlignment(GenotypeTable theAlignment, String delimit, String saveFile) {

        if ((saveFile == null) || (saveFile.length() == 0)) {
            return null;
        }
        saveFile = Utils.addSuffixIfNeeded(saveFile, ".txt");
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(new File(saveFile));
            bw = new BufferedWriter(fw);

            bw.write("Taxa");
            int numSites = theAlignment.numberOfSites();
            for (int j = 0; j < numSites; j++) {
                bw.write(delimit);
                bw.write(String.valueOf(theAlignment.chromosomalPosition(j)));
            }
            bw.write("\n");

            for (int r = 0, n = theAlignment.numberOfTaxa(); r < n; r++) {
                bw.write(theAlignment.taxaName(r));
                for (int i = 0; i < numSites; i++) {
                    bw.write(delimit);
                    bw.write(theAlignment.genotypeAsString(r, i));
                }
                bw.write("\n");
            }

            return saveFile;

        } catch (Exception e) {
            myLogger.error("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    /**
     * print alignment (in PHYLIP SEQUENTIAL format)
     */
    public static void printSequential(GenotypeTable a, Writer out) throws IOException {
        // PHYLIP header line
        out.write("  " + a.numberOfTaxa() + " " + a.numberOfSites() + "  S" + "\n");

        // Print sequences
        for (int s = 0; s < a.numberOfTaxa(); s++) {
            int n = 0;
            while (n < a.numberOfSites()) {
                if (n == 0) {
                    format.displayLabel(out, a.taxaName(s), 10);
                    out.write("     ");
                } else {
                    out.write("               ");
                }
                printNextSites(a, out, false, s, n, 50);
                out.write("\n");
                n += 50;
            }
        }
    }

    /**
     * print alignment (in PHYLIP 3.4 INTERLEAVED format)
     */
    public static void printInterleaved(GenotypeTable a, Writer out) throws IOException {
        int n = 0;

        // PHYLIP header line
        out.write("  " + a.numberOfTaxa() + " " + a.numberOfSites() + "\n");

        // Print sequences
        while (n < a.numberOfSites()) {
            for (int s = 0; s < a.numberOfTaxa(); s++) {
                if (n == 0) {
                    format.displayLabel(out, a.taxaName(s), 10);
                    out.write("     ");
                } else {
                    out.write("               ");
                }
                printNextSites(a, out, true, s, n, 50);
                out.write("\n");
            }
            out.write("\n");
            n += 50;
        }
    }

    /**
     * Print alignment (in CLUSTAL W format)
     */
    public static void printCLUSTALW(GenotypeTable a, Writer out) throws IOException {
        int n = 0;

        // CLUSTAL W header line
        out.write("CLUSTAL W multiple sequence alignment\n\n");

        // Print sequences
        while (n < a.numberOfSites()) {
            out.write("\n");
            for (int s = 0; s < a.numberOfTaxa(); s++) {
                format.displayLabel(out, a.taxaName(s), 10);
                out.write("     ");

                printNextSites(a, out, false, s, n, 50);
                out.write("\n");
            }
            // Blanks in status line are necessary for some parsers)
            out.write("               \n");
            n += 50;
        }
    }

    private static void printNextSites(GenotypeTable a, Writer out, boolean chunked, int seq, int start, int num) throws IOException {
        // Print next num characters
        for (int i = 0; (i < num) && (start + i < a.numberOfSites()); i++) {
            // Chunks of 10 characters
            if (i % 10 == 0 && i != 0 && chunked) {
                out.write(' ');
            }
            out.write(a.genotypeAsString(seq, start + i));
        }
    }

    public static String writeAlignmentToSerialGZ(GenotypeTable sba, String outFile) {

        long time = System.currentTimeMillis();

        File theFile = null;
        FileOutputStream fos = null;
        GZIPOutputStream gz = null;
        ObjectOutputStream oos = null;
        try {
            theFile = new File(Utils.addSuffixIfNeeded(outFile, ".serial.gz"));
            fos = new FileOutputStream(theFile);
            gz = new GZIPOutputStream(fos);
            oos = new ObjectOutputStream(gz);
            oos.writeObject(sba);
            return theFile.getName();
        } catch (Exception e) {
            e.printStackTrace();
            myLogger.error("Error writing Serial GZ: " + theFile.getName() + ": " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing Serial GZ: " + theFile.getName() + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                oos.flush();
                oos.close();
                gz.close();
                fos.close();
            } catch (Exception e) {
                // do nothing
            }
            myLogger.info("writeAlignmentToSerialGZ: " + theFile.toString() + "  Time: " + (System.currentTimeMillis() - time));
        }

    }

    /**
     * Simple method to export a GenotypeTable to a fasta file.
     * @param gt Import GenotypeTable which we want to export.
     * @param out FileWriter which will be used to export the FASTA file.
     */
    public static void writeFasta(GenotypeTable gt, Writer out) {
        try {
            TaxaList tl = gt.taxa();
            for(int i = 0; i < tl.size(); i++) {
                out.write(">");
                out.write(tl.get(i).getName());
                out.write("\n");
                for(int j = 0; j < gt.positions().size(); j++) {
                    byte call = gt.genotype(i,j);
                    out.write(NucleotideAlignmentConstants.getNucleotideIUPAC(call));
                }
                out.write("\n");
            }
        }
        catch(Exception e) {
            e.printStackTrace();
            myLogger.error("Error writing FASTA file: " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing FastaAlignment: "+ ExceptionUtils.getExceptionCauses(e));
        }
    }

    /**
     * Simple method to export a GenotypeTable to a fasta file while removing Gap characters.
     * @param gt GenotypeTable which we want to export.
     * @param out FileWriter which will be used to export the FASTA file.
     */
    public static void writeFastaNoGaps(GenotypeTable gt, Writer out) {
        try {
            TaxaList tl = gt.taxa();
            for(int i = 0; i < tl.size(); i++) {
                out.write(">");
                out.write(tl.get(i).getName());
                out.write("\n");
                for(int j = 0; j < gt.positions().size(); j++) {
                    byte call = gt.genotype(i,j);
                    if(!NucleotideAlignmentConstants.getNucleotideIUPAC(call).equals("-")) {
                        out.write(NucleotideAlignmentConstants.getNucleotideIUPAC(call));
                    }
                }
                out.write("\n");
            }
        }
        catch(Exception e) {
            e.printStackTrace();
            myLogger.error("Error writing FASTA file: " + ExceptionUtils.getExceptionCauses(e));
            throw new IllegalArgumentException("Error writing FastaAlignment: "+ ExceptionUtils.getExceptionCauses(e));
        }
    }
}