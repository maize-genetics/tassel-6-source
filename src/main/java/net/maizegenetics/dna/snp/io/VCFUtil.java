/*
 * VCFUtil
 */
package net.maizegenetics.dna.snp.io;

import com.google.common.primitives.Ints;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Qi Sun
 * @author Zack Miller
 */
public class VCFUtil {

    private static final Logger myLogger = Logger.getLogger(VCFUtil.class);

    // variables for calculating OS and PL for VCF, might not be in the correct class
    private static double error;
    private static double v1;
    private static double v2;
    private static double v3;
    private static int[][][] myGenoScoreMap;
    public static final int VCF_DEFAULT_MAX_NUM_ALLELES = 3;

    static {
        error = 0.001; //TODO this seems low, is this the standard
        v1 = Math.log10(1.0 - error * 3.0 / 4.0);
        v2 = Math.log10(error / 4);
        v3 = Math.log10(0.5 - (error / 4.0));
        myGenoScoreMap = new int[128][128][];
        for (int i = 0; i < 128; i++) {
            for (int j = 0; j < 128; j++) {
                myGenoScoreMap[i][j] = calcScore(i, j);
            }
        }
    }

    private VCFUtil() {
        // utility
    }


    public static int[] getScore(int i, int j) {
        if (i > 127 || j > 127) return calcScore(i, j);
        return myGenoScoreMap[i][j];
    }

    // Calculate QS and PL for VCF might not be in the correct class
    private static int[] calcScore(int a, int b) {
        int[] results = new int[4];
        int n = a + b;
        int m = a;
        if (b > m) {
            m = b;
        }

        double fact = 0;
        if (n > m) {
            for (int i = n; i > m; i--) {
                fact += Math.log10(i);
            }
            for (int i = 1; i <= (n - m); i++) {
                fact -= Math.log10(i);
            }
        }
        double aad = Math.pow(10, fact + (double) a * v1 + (double) b * v2);
        double abd = Math.pow(10, fact + (double) n * v3);
        double bbd = Math.pow(10, fact + (double) b * v1 + (double) a * v2);
        double md = aad;
        if (md < abd) {
            md = abd;
        }
        if (md < bbd) {
            md = bbd;
        }
        int gq = 0;
        if ((aad + abd + bbd) > 0) {
            gq = (int) (md / (aad + abd + bbd) * 100);
        }

        int aa = (int) (-10 * (fact + (double) a * v1 + (double) b * v2));
        int ab = (int) (-10 * (fact + (double) n * v3));
        int bb = (int) (-10 * (fact + (double) b * v1 + (double) a * v2));

        m = aa;
        if (m > ab) {
            m = ab;
        }
        if (m > bb) {
            m = bb;
        }
        aa -= m;
        ab -= m;
        bb -= m;
        results[0] = aa > 255 ? 255 : aa;
        results[1] = ab > 255 ? 255 : ab;
        results[2] = bb > 255 ? 255 : bb;
        results[3] = gq;

        return results;
    }

    public static byte resolveVCFGeno(byte[] alleles, int[][] allelesInTaxa, int tx) {
        int[] alleleDepth = new int[allelesInTaxa.length];
        for (int i = 0; i < allelesInTaxa.length; i++) {
            alleleDepth[i] = allelesInTaxa[i][tx];
        }
        return resolveVCFGeno(alleles, alleleDepth);
    }

    public static byte resolveVCFGeno(byte[] alleles, int[] alleleDepth) {
        int depth = 0;
        for (int i = 0; i < alleleDepth.length; i++) {
            depth += alleleDepth[i];
        }
        if (depth == 0) {
            return (byte) ((GenotypeTable.UNKNOWN_ALLELE << 4) | GenotypeTable.UNKNOWN_ALLELE);
        }
        int max = 0;
        byte maxAllele = GenotypeTable.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = GenotypeTable.UNKNOWN_ALLELE;
        for (int i = 0; i < alleles.length; i++) {
            if (alleleDepth[i] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = alleleDepth[i];
                maxAllele = alleles[i];
            } else if (alleleDepth[i] > nextMax) {
                nextMax = alleleDepth[i];
                nextMaxAllele = alleles[i];
            }
        }
        if (alleles.length == 1) {
            return (byte) ((alleles[0] << 4) | alleles[0]);
        } else {
            max = (max > 127) ? 127 : max;
            nextMax = (nextMax > 127) ? 127 : nextMax;
            int[] scores = getScore(max, nextMax);
            if ((scores[1] <= scores[0]) && (scores[1] <= scores[2])) {
                return (byte) ((maxAllele << 4) | nextMaxAllele);
            } else if ((scores[0] <= scores[1]) && (scores[0] <= scores[2])) {
                return (byte) ((maxAllele << 4) | maxAllele);
            } else {
                return (byte) ((nextMaxAllele << 4) | nextMaxAllele);
            }
        }
    }


    public static int[] resolveRefSorted(int[] sortedAlleles, byte refAllele) {
        int[] sortedAllelesResolved = new int[sortedAlleles.length];
        int indexOfRefAllele = Ints.indexOf(sortedAlleles, refAllele);

        //If indexOfRefAllele is -1 the refAllele is not contained in sortedAlleles
        //We need to add it
        if (indexOfRefAllele < 0) {
            //Set index to 0 as we want the Ref in the first position
            indexOfRefAllele = 0;

            if (refAllele != GenotypeTable.UNKNOWN_ALLELE) {
                int[] sortedAllelesExpanded = new int[sortedAlleles.length + 1];
                sortedAllelesExpanded[0] = refAllele;
                for (int i = 0; i < sortedAlleles.length; i++) {
                    sortedAllelesExpanded[i + 1] = sortedAlleles[i];
                }
                sortedAllelesResolved = sortedAllelesExpanded;
            } else {
                for (int i = 0; i < sortedAllelesResolved.length; i++) {
                    sortedAllelesResolved[i] = sortedAlleles[i];
                }
            }
        }
        //Resort sorted Alleles if ref is not first in array
        else if (indexOfRefAllele != 0) {
            sortedAllelesResolved[0] = refAllele;

            for (int i = indexOfRefAllele; i > 0; i--) {
                sortedAllelesResolved[i] = sortedAlleles[i - 1];
            }
            for (int i = indexOfRefAllele + 1; i < sortedAllelesResolved.length; i++) {
                sortedAllelesResolved[i] = sortedAlleles[i];
            }
        } else {
            for (int i = 0; i < sortedAllelesResolved.length; i++) {
                sortedAllelesResolved[i] = sortedAlleles[i];
            }
        }

        return sortedAllelesResolved;
    }

    public static boolean indelInKnownVariant(String[] knownVariants) {
        for (String variant : knownVariants) {
            if (variant.length() > 1) {
                return true;
            }
        }
        return false;
    }

    public static boolean indelMinusInKnownVariant(String[] knownVariants) {
        for (String variant : knownVariants) {
            if (variant.equals("-") || variant.equals("+")) {
                return true;
            }
        }
        return false;
    }

    public static Tuple<int[], String[]> resolveSortedAndKnownVariantsExport(int[] sortedAllelesInput, String[] knownVariantsInput) {
        int[] sortedAlleles = Arrays.copyOf(sortedAllelesInput, sortedAllelesInput.length);
        String[] knownVariants = Arrays.copyOf(knownVariantsInput, knownVariantsInput.length);
        if (knownVariants.length > 0) {
            //ReOrder based on variant alleles
            //Store a tempSortedAlleles so we can appropriately handle hapmap to vcf
            //int[] tempSortedAlleles = new int[knownVariants.length];

            //ArrayList to hold the Sorted Alleles Indices Temporarily as the ordering will change
            ArrayList<Integer> tempSortedAlleles = new ArrayList<Integer>();

            //Loop through all the knownVariants and check to see if we have an indel
            boolean knownVariantIndel = VCFUtil.indelInKnownVariant(knownVariants);

            //If we do have an indel, we can add the variants after picking off the first character to the tempSortedAlleles
            if (knownVariantIndel) {
                //Loop through the variants
                for (int i = 0; i < knownVariants.length; i++) {
                    //Pull off the first character if it exists
                    if (knownVariants[i].length() > 1) {
                        String parsedVariant = knownVariants[i].substring(1);
                        tempSortedAlleles.add((int) NucleotideAlignmentConstants.getNucleotideAlleleByte(parsedVariant.charAt(0)));
                    } else {
                        //Mark as deletion
                        tempSortedAlleles.add((int) NucleotideAlignmentConstants.getNucleotideAlleleByte('-'));
                    }
                }
            } else {
                //If we dont have an indel, we can add it to the allele array
                if (sortedAlleles.length < knownVariants.length) {
                    //Clear it out, we probably dont need to do this
                    tempSortedAlleles = new ArrayList<Integer>();
                }
                int nIndex = -1;
                for (int i = 0; i < knownVariants.length; i++) {
                    //ZRM22 Mar 22
                    if (knownVariants[i].charAt(0) != 'N') {
                        tempSortedAlleles.add((int) NucleotideAlignmentConstants.getNucleotideAlleleByte(knownVariants[i].charAt(0)));
                    } else {
                        //If N is in our known Variants list but we do not have an indel, we need to remove it
                        nIndex = i;
                    }
                }
                if (nIndex != -1) {
                    //if we have an N we need to resize KnownVariants
                    String[] knownVariantsSmall = new String[knownVariants.length - 1];
                    for (int i = 0; i < knownVariants.length; i++) {
                        if (i < nIndex) {
                            knownVariantsSmall[i] = knownVariants[i];
                        } else if (i > nIndex) {
                            knownVariantsSmall[i - 1] = knownVariants[i];
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
            if (indelsInKnownVariants) {
                indelsExist = true;
            }

            //Go through sorted alleles and also check for indels
            for (int i = 0; i < sortedAlleles.length; i++) {
                if (NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]).equals("-")) {
                    indelsExist = true;
                }
            }
            //Move To Function/

            for (String variant : knownVariants) {
                if (indelsExist && !indelsInKnownVariants) {
                    knownVariantsList.add("N" + variant);
                } else {
                    knownVariantsList.add(variant);
                }
            }

            //Go through sorted alleles
            for (int i = 0; i < sortedAlleles.length; i++) {
                //If a sorted allele is not in tempSortedAlleles,
                if (!tempSortedAlleles.contains(sortedAlleles[i])) {
                    //if its not add it to sorted alleles and knownVariants
                    tempSortedAlleles.add(sortedAlleles[i]);
                    //Check for an indel
                    if (indelsExist) {
                        if (NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]).equals("-")) {
                            knownVariantsList.add("N");
                        } else {
                            knownVariantsList.add("N" + NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]));
                        }
                        //                     knownVariantsList.add("N"+NucleotideAlignmentConstants.getHaplotypeNucleotide((byte)sortedAlleles[i]));
                    } else {
                        knownVariantsList.add(NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]));
                    }
                }
            }
            //reset knownVariants and sortedAlleles to reflect the changes
            String[] knownVariantsExtended = new String[knownVariantsList.size()];
            for (int i = 0; i < knownVariantsExtended.length; i++) {
                knownVariantsExtended[i] = knownVariantsList.get(i);
            }
            knownVariants = knownVariantsExtended;

            int[] sortedAllelesExtended = new int[tempSortedAlleles.size()];
            for (int i = 0; i < sortedAllelesExtended.length; i++) {
                sortedAllelesExtended[i] = tempSortedAlleles.get(i);
            }
            sortedAlleles = sortedAllelesExtended;
            //sortedAlleles = tempSortedAlleles.toArray(new int[tempSortedAlleles.size()]);
        } else {
            //No known variants, but we need to handle indels
            int indelIndex = -1;
            //loop through sorted alleles
            for (int i = 0; i < sortedAlleles.length; i++) {
                //if we find an indel mark the index and set a boolean
                if (sortedAlleles[i] == (int) NucleotideAlignmentConstants.getNucleotideAlleleByte('-')) {
                    indelIndex = i;
                    break;
                }
            }

            knownVariants = new String[sortedAlleles.length];
            for (int i = 0; i < knownVariants.length; i++) {
                if (indelIndex == -1) {
                    knownVariants[i] = "" + NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]);
                } else {
                    if (indelIndex == i) {
                        knownVariants[i] = "N";
                    } else {
                        knownVariants[i] = "N" + NucleotideAlignmentConstants.getHaplotypeNucleotide((byte) sortedAlleles[i]);
                    }
                }
            }
        }
        return new Tuple<int[], String[]>(sortedAlleles, knownVariants);
    }

    //Utility method to find the top 2 depths from the depth array
    public static int[][] calcTop2Depths(int[] siteAlleleDepths) {
        //loop through the allele depths and figure out which one
        int highestDepth = -1;
        int highestIndex = -1;
        int secondHighestDepth = -1;
        int secondHighestIndex = -1;

        for (int depthCounter = 0; depthCounter < siteAlleleDepths.length; depthCounter++) {
            if (siteAlleleDepths[depthCounter] > highestDepth) {
                //Take the Current Highest Depth and push it to the second highest
                secondHighestDepth = highestDepth;
                secondHighestIndex = highestIndex;
                //Then set the highestDepth to the current value
                highestDepth = siteAlleleDepths[depthCounter];
                highestIndex = depthCounter;
            } else {
                //If not the best check for the second highest
                if (siteAlleleDepths[depthCounter] > secondHighestDepth) {
                    //Then reset the second highest value with the current
                    secondHighestDepth = siteAlleleDepths[depthCounter];
                    secondHighestIndex = depthCounter;
                }
                //Otherwise let it fall through
            }
        }
        int[][] top2Depths = {{highestIndex, highestDepth}, {secondHighestIndex, secondHighestDepth}};
        return top2Depths;
    }


    /**
     * Method to convert a genotype table to a List of VariantContexts
     *
     * This is done by walking through the genotypeTable's positions in order.
     * To Handle indels, we need to process a set of positions at a time
     *
     * It will start with the first position and hold its index temporarily
     * Combine consecutive indel positions into this index list.
     * A consecutive position is one where the chromosomes match, and the current position is less than 1 physical position away from the previous
     * We need to check less than one as an insertion has the same physical reference position
     * If the consecutive position also contains an indel character "+" or "-" we need to add it to the block
     * Otherwise(either non indel or non-consecutive) process the previous block and setup a new block with the current site
     *
     * @param genotypeTable
     *
     * @return
     */
    public static List<VariantContext> convertGenotypeTableToVariantContextList(GenotypeTable genotypeTable) {
        List<VariantContext> variantContextList = new ArrayList<>();

        List<Integer> indelSiteBlockList = new ArrayList<>();
        indelSiteBlockList.add(0);
        for (int i = 1; i < genotypeTable.numberOfSites(); i++) {
            //Check to see if the current position is consecutive to the previous
            //Note we need to check that the difference in physical position is less than 1 in case of insertions
            if (genotypeTable.positions().get(i - 1).getChromosome().equals(genotypeTable.positions().get(i).getChromosome()) &&
                    (genotypeTable.positions().get(i - 1).getPosition() + 1) >= genotypeTable.positions().get(i).getPosition()) {

                //Figure out if we have an insertion or deletion allele
                boolean isIndel = false;

                if(genotypeTable.referenceAllele(i)== NucleotideAlignmentConstants.INSERT_ALLELE || genotypeTable.referenceAllele(i)== NucleotideAlignmentConstants.GAP_ALLELE ) {
                    indelSiteBlockList.add(i);
                    continue;
                }

                for (byte allele : genotypeTable.alleles(i)) {
                    if (allele == NucleotideAlignmentConstants.INSERT_ALLELE || allele == NucleotideAlignmentConstants.GAP_ALLELE) {
                        isIndel = true;
                        break;
                    }
                }
                if (isIndel) {
                    //If we know it is an indel we want to add the site to the block
                    indelSiteBlockList.add(i);
                    continue;
                }
            }

            //Process the current block and add it to the variantContextList
            variantContextList.add(convertGenotypeTableSiteToVariantContext(genotypeTable, indelSiteBlockList));
            //Reset the block
            indelSiteBlockList = new ArrayList<>();
            indelSiteBlockList.add(i);

        }

        //Process the remaining block
        variantContextList.add(convertGenotypeTableSiteToVariantContext(genotypeTable, indelSiteBlockList));


        return variantContextList;
    }

    /**
     * Method to convert a list of positions into a single variantContext record
     *
     * If sites only has one element, it will just create genotypes for that position(Fixing indels with missing information)
     * If sites has multiple consecutive elements, it is likely an indel
     * This needs to walk through the sites and combine the alleles together ignoring - and + nucleotides
     *
     * @param genotypeTable
     * @param sites
     *
     * @return
     */
    public static VariantContext convertGenotypeTableSiteToVariantContext(GenotypeTable genotypeTable, List<Integer> sites) {

        //Get the start position for use later
        int startPos = genotypeTable.positions().get(sites.get(0)).getPosition();


        //Create a StringBuilder object for the left and the right alleles.
        //The StringBuilders are building the alleleString in order of sites
        List<Tuple<StringBuilder, StringBuilder>> alleleStrings = IntStream.range(0, genotypeTable.numberOfTaxa())
                .boxed()
                .map(index -> new Tuple<>(new StringBuilder(), new StringBuilder()))
                .collect(Collectors.toList());

        //Set up the reference Allele String builder
        StringBuilder referenceAlleleStringBuilder = new StringBuilder();

        Set<String> knownVariantSet = new LinkedHashSet<>();

        //Loop through each site and figure out the allele strings
        //If it is a - or + we ignore it
        for (Integer currentSite : sites) {
            byte refAllele = genotypeTable.referenceAllele(currentSite);
            if (refAllele != NucleotideAlignmentConstants.GAP_ALLELE && refAllele != NucleotideAlignmentConstants.INSERT_ALLELE) {
                //We add it because it is an actual allele value
                referenceAlleleStringBuilder.append(NucleotideAlignmentConstants.getHaplotypeNucleotide(refAllele));
            }

            String[] knownVariantsAtSite = genotypeTable.positions().get(currentSite).getKnownVariants();
            if(knownVariantsAtSite!=null) {
                knownVariantSet.addAll(Arrays.asList(knownVariantsAtSite));
            }

            //loop through each taxon
            for (int taxonIndex = 0; taxonIndex < genotypeTable.numberOfTaxa(); taxonIndex++) {
                String[] genotypeCalls = genotypeTable.genotypeAsStringArray(taxonIndex,currentSite);

                if(genotypeCalls[0].equals(NucleotideAlignmentConstants.UNDEFINED_ALLELE_STR)) {
                    genotypeCalls[0] = GenotypeTable.UNKNOWN_ALLELE_STR;
                }
                if(genotypeCalls[1].equals(NucleotideAlignmentConstants.UNDEFINED_ALLELE_STR)) {
                    genotypeCalls[1] = GenotypeTable.UNKNOWN_ALLELE_STR;
                }

                //Add the calls to the string builders for this taxon
                //We only want to add in the alleles if they are not + or minus as those are not valid VCF characters
                if (!genotypeCalls[0].equals(NucleotideAlignmentConstants.GAP_ALLELE_STR) && !genotypeCalls[0].equals(NucleotideAlignmentConstants.INSERT_ALLELE_STR)) {
                    alleleStrings.get(taxonIndex).getX().append(genotypeCalls[0]);
                }
                if (!genotypeCalls[1].equals(NucleotideAlignmentConstants.GAP_ALLELE_STR) && !genotypeCalls[1].equals(NucleotideAlignmentConstants.INSERT_ALLELE_STR)) {
                    alleleStrings.get(taxonIndex).getY().append(genotypeCalls[1]);
                }

            }
        }

        //Because we are not going to change the String Builders any more, convert them to String objects
        List<Tuple<String, String>> alleleStringVals = alleleStrings.stream()
                .map(singleAlleleTuple -> new Tuple<>(singleAlleleTuple.getX().toString(), singleAlleleTuple.getY().toString()))
                .collect(Collectors.toList());

        //We need to check to see if we have to add in Ns for indels where the previous position is not known
        boolean fixIndels = sites.size() == 1 && alleleStringVals.stream().flatMap(alleleTuple -> Arrays.asList(alleleTuple.getX(), alleleTuple.getY()).stream()).anyMatch(alleleString -> alleleString.equals(""));


        //Fix the allele values if we have an indel and no additional information
        if (fixIndels) {
            //Add an N to the reference allele
            referenceAlleleStringBuilder.insert(0, "N");
            //Subtract the startPosition by 1 bp as we are adding an allele to the start of the string
            startPos--;
            //Fix the indel positions by adding in the Ns
            alleleStringVals = fixIndelPositions(alleleStringVals);
        }

        //We now have all the possible allele strings figured out we can get the Allele objects
        Map<String, Allele> alleleStringToObjMap = createAlleleStringToObjMap(referenceAlleleStringBuilder.toString(), alleleStringVals,knownVariantSet, genotypeTable.positions().get(sites.get(0)));


        //Convert the allele objects into a list so we can add them to the VariantContextBuilder
        List<Allele> alleleObjectList = alleleStringToObjMap.keySet().stream()
                .map(alleleString -> alleleStringToObjMap.get(alleleString))
                .filter(alleleObject -> !alleleObject.equals(Allele.NO_CALL)) //Remove No_Call(missing) alleles as they will break the variantContext
                .distinct()//need this as X, N and * all map to Allele("N",false)
                .collect(Collectors.toList());

        //Get the genotype Calls for each taxon
        List<Genotype> genotypeList = getGenotypes(genotypeTable, alleleStringVals, alleleStringToObjMap);

        //Create the variant Context for this record
        VariantContextBuilder vcb = new VariantContextBuilder(".", //Is there a better source to put here?
                genotypeTable.positions().chromosomeName(sites.get(0)),
                startPos,
                startPos + referenceAlleleStringBuilder.toString().length() - 1,
                alleleObjectList)
                .genotypes(genotypeList);

        //Make the VaraintContext
        return vcb.make();
    }

    /**
     * Fix the indel positions for the allele strings
     * Basically just add in an N at the start of each String
     *
     * TODO handle indels at the start of the chromosome(Add N to the end)
     *
     * @param alleleStringVals
     *
     * @return
     */
    private static List<Tuple<String, String>> fixIndelPositions(List<Tuple<String, String>> alleleStringVals) {
        return alleleStringVals.stream()
                .map(alleleTuple -> new Tuple<>("N" + alleleTuple.getX(), "N" + alleleTuple.getY()))
                .collect(Collectors.toList());
    }

    /**
     * Convert the alleleStrings into HTSJDK Allele objects and create a Mapping of String to Allele Object
     *
     * @param referenceString
     * @param alleleStrings
     *
     * @return
     */
    private static Map<String, Allele> createAlleleStringToObjMap(String referenceString, List<Tuple<String, String>> alleleStrings, Set<String> knownVariantSet,Position position) {
        //We need to loop through each possible allele string and make a new HTSJDK Allele object for each one
        try {

            //TODO setup a cache for the common ACGT ref and non Ref alleles for speed
            Map<String, Allele> stringValueToAlleleMap = new HashMap<>();

            if (referenceString == null || referenceString.equals("")) {
                myLogger.warn("NULL reference allele found: " + position.toString() + " Putting N in for ref.");
                stringValueToAlleleMap.put(GenotypeTable.UNKNOWN_ALLELE_STR, Allele.create("N", true));
            } else {
                //Setup the reference allele
                stringValueToAlleleMap.put(referenceString, Allele.create(referenceString, true));
            }
            Map<String, Allele> uniqueAlternateAlleles = alleleStrings.stream()
                    .flatMap(sbTuple -> Arrays.asList(sbTuple.getX(), sbTuple.getY()).stream()) //FlatMap the tuples
                    .filter(alleleString -> !stringValueToAlleleMap.containsKey(alleleString))//Make sure we dont include the reference allele as we do not want to overwrite
                    .distinct() //Remove duplicates
                    .filter(alleleString -> alleleString != null)
                    .collect(Collectors.toMap(alleleString -> alleleString, alleleString -> {
                        if(alleleString.equals("*") || alleleString.equals(NucleotideAlignmentConstants.UNDEFINED_ALLELE_STR)) {
                            return Allele.NO_CALL;
                        }
                        else {
                            return Allele.create(alleleString, false);
                        }
                    })); //Convert each string to an Allele Object and add it to a map

            //Add in all the unique alternate alleles to the map
            stringValueToAlleleMap.putAll(uniqueAlternateAlleles);

            Map<String, Allele> uniqueKnownVariantAlleles = knownVariantSet.stream()
                    .filter(alleleString -> !stringValueToAlleleMap.containsKey(alleleString))//Make sure we dont include the reference allele as we do not want to overwrite
                    .distinct() //Remove duplicates
                    .filter(alleleString -> alleleString != null)
                    .collect(Collectors.toMap(alleleString -> alleleString, alleleString -> {
                        if(alleleString.equals("*") || alleleString.equals(NucleotideAlignmentConstants.UNDEFINED_ALLELE_STR)) {
                            return Allele.NO_CALL;
                        }
                        else {
                            return Allele.create(alleleString, false);
                        }
                    })); //Convert each string to an Allele Object and add it to a map

            stringValueToAlleleMap.putAll(uniqueKnownVariantAlleles);


            //we need to check to see if we have an indel in our allele list.  If we do not, we need to remove the N key object as these genotypes will not be in the GenotypeContext
            boolean hasIndel = false;
            int refLength = referenceString.length();
            for (String alleleValue : stringValueToAlleleMap.keySet()) {
                if (alleleValue.length() != refLength) {
                    hasIndel = true;
                    break;
                }
            }

            if (!hasIndel && !referenceString.equals(GenotypeTable.UNKNOWN_ALLELE_STR) && stringValueToAlleleMap.containsKey(GenotypeTable.UNKNOWN_ALLELE_STR)) {
                //remove the N allele from the map as it is not valid
                stringValueToAlleleMap.put(GenotypeTable.UNKNOWN_ALLELE_STR, Allele.NO_CALL);
            }

            return stringValueToAlleleMap;
        }
        catch(Exception e) {
            throw new IllegalStateException(e);
        }
    }

    /**
     * Method to create a list of genotype objects for the sites
     *
     * @param genotypeTable
     * @param listOfAlleleCalls
     * @param callToAlleleObjMap
     *
     * @return
     */
    private static List<Genotype> getGenotypes(GenotypeTable genotypeTable, List<Tuple<String, String>> listOfAlleleCalls, Map<String, Allele> callToAlleleObjMap) {
        List<Genotype> genotypeList = new ArrayList<>();

        for (int taxonIndex = 0; taxonIndex < genotypeTable.numberOfTaxa(); taxonIndex++) {
            List<Allele> currentTaxonsAlleles = Arrays.asList(callToAlleleObjMap.get(listOfAlleleCalls.get(taxonIndex).getX()),
                    callToAlleleObjMap.get(listOfAlleleCalls.get(taxonIndex).getY()));


            GenotypeBuilder currentGenotypeBuilder = new GenotypeBuilder(genotypeTable.taxaName(taxonIndex), currentTaxonsAlleles);
            // .AD(new int[0]) // For now we do not care about depth.  TODO Figure out a better way to handle allele depth
            // .DP(0)
            // .GQ(0)
            // .PL(new int[0]);

            genotypeList.add(currentGenotypeBuilder.make());
        }

        return genotypeList;
    }

}
