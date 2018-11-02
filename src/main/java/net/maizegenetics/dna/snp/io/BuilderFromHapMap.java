/*
 *  BuilderFromHapMap
 * 
 *  Created on Aug 1, 2014
 */
package net.maizegenetics.dna.snp.io;

import com.google.common.collect.SetMultimap;

import java.io.BufferedReader;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListIOUtils;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class BuilderFromHapMap {

    private static final Logger myLogger = Logger.getLogger(BuilderFromHapMap.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\t");
    private static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    private static final int SNPID_INDEX = 0;
    private static final int VARIANT_INDEX = 1;
    private static final int CHROMOSOME_INDEX = 2;
    private static final int POSITION_INDEX = 3;
    private static final int NUM_VALUES_PROCESSED_TOGETHER = 7 << 20;

    private final String myHapmapFile;
    private boolean mySortTaxaAlphabetically = false;
    private boolean mySortPositions = false;
    private final ProgressListener myProgressListener;

    private BuilderFromHapMap(String hapmapFile, ProgressListener listener) {
        myHapmapFile = hapmapFile;
        myProgressListener = listener;
    }

    public static BuilderFromHapMap getBuilder(String hapmapFile) {
        return new BuilderFromHapMap(hapmapFile, null);
    }

    public static BuilderFromHapMap getBuilder(String hapmapFile, ProgressListener listener) {
        return new BuilderFromHapMap(hapmapFile, listener);
    }

    public GenotypeTable build() {

        ExecutorService pool = null;
        try (BufferedReader reader = Utils.getBufferedReader(myHapmapFile, 1 << 20)) {

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

            TaxaListBuilder taxaList = processTaxa(currLine, sampAnnoBuild);
            int numTaxa = taxaList.numberOfTaxa();

            Map<String, Chromosome> chromosomeLookup = new ConcurrentHashMap<>();

            currLine = reader.readLine();

            boolean isOneLetter = false;
            String[] tokens = WHITESPACE_PATTERN.split(currLine, NUM_HAPMAP_NON_TAXA_HEADERS + 1);
            if (tokens.length <= NUM_HAPMAP_NON_TAXA_HEADERS) {
                throw new IllegalStateException("BuilderFromHapMap: Header Incorrectly Formatted: See:\nhttps://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-hapmap");
            }
            double avg = (double) (tokens[NUM_HAPMAP_NON_TAXA_HEADERS].length() + 1) / (double) numTaxa;
            if ((avg > 1.99) && (avg < 2.01)) {
                isOneLetter = true;
            } else if ((avg > 2.99) && (avg < 3.01)) {
                isOneLetter = false;
            } else {
                throw new IllegalStateException("BuilderFromHapMap: Genotype coded wrong use 1 or 2 letters per genotype. Average chars including tab: " + avg + "  Or first site has incorrect number of values. Number of taxa: " + numTaxa);
            }

            int numThreads = Runtime.getRuntime().availableProcessors();
            pool = Executors.newFixedThreadPool(numThreads);
            List<Future<ProcessHapmapBlock>> futures = new ArrayList<>();

            int numSitesToProcessTogether = NUM_VALUES_PROCESSED_TOGETHER / numTaxa;
            numSitesToProcessTogether = Math.min(1 << 16, numSitesToProcessTogether);
            numSitesToProcessTogether = Math.max(512, numSitesToProcessTogether);

            ArrayList<String> textLines = new ArrayList<>(numSitesToProcessTogether);
            int numLines = 0;
            while (currLine != null) {
                textLines.add(currLine);
                numLines++;
                if (numLines % numSitesToProcessTogether == 0) {
                    ProcessHapmapBlock processBlock = new ProcessHapmapBlock(textLines, numTaxa, chromosomeLookup, isOneLetter);
                    futures.add(pool.submit(processBlock));
                    textLines = new ArrayList<>(numSitesToProcessTogether);
                }
                currLine = reader.readLine();
            }

            if (textLines.size() > 0) {
                ProcessHapmapBlock processBlock = new ProcessHapmapBlock(textLines, numTaxa, chromosomeLookup, isOneLetter);
                futures.add(pool.submit(processBlock));
            }

            int currentSite = 0;
            PositionListBuilder positions = new PositionListBuilder();
            GenotypeCallTableBuilder genotypes = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numTaxa, numLines);

            int numFutures = futures.size();
            int count = 0;
            for (Future<ProcessHapmapBlock> future : futures) {
                ProcessHapmapBlock pb = future.get();
                positions.addAll(pb.getPositions());
                SuperByteMatrix bgTS = pb.getGenotypes();
                for (int t = 0; t < bgTS.getNumRows(); t++) {
                    for (int s = 0; s < bgTS.getNumColumns(); s++) {
                        genotypes.setBase(t, currentSite + s, bgTS.get(t, s));
                    }
                }
                currentSite += pb.getNumberSitesProcessed();
                if (myProgressListener != null) {
                    count++;
                    myProgressListener.progress(count * 100 / numFutures, null);
                }
            }
            pool.shutdown();

            if (mySortTaxaAlphabetically) {
                taxaList.sortTaxaAlphabetically(genotypes);
            }

            if (mySortPositions) {
                positions.sortPositions(genotypes);
            }

            if (positions.validateOrdering() == false) {
                throw new IllegalStateException("BuilderFromHapMap: Ordering incorrect. HapMap must be ordered by position. Please first use SortGenotypeFilePlugin to correctly order the file.");
            }

            return GenotypeTableBuilder.getInstance(genotypes.build(), positions.build(), taxaList.build());

        } catch (Exception e) {
            if (pool != null) {
                pool.shutdown();
            }
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException(e.getMessage());
        }

    }

    private class ProcessHapmapBlock implements Callable<ProcessHapmapBlock> {

        private List<String> myInputLines;
        private final Map<String, Chromosome> myChromosomeLookup;
        private final boolean myIsOneLetter;
        private final List<Position> myPositionList;
        private final int myNumSitesToProcess;
        private final int myNumTaxa;
        private SuperByteMatrix myGenotypes;

        public ProcessHapmapBlock(List<String> inputLines, int numTaxa, Map<String, Chromosome> chromosomeLookup, boolean isOneLetter) {
            myInputLines = inputLines;
            myChromosomeLookup = chromosomeLookup;
            myNumTaxa = numTaxa;
            myIsOneLetter = isOneLetter;
            myNumSitesToProcess = inputLines.size();
            myPositionList = new ArrayList<>(myNumSitesToProcess);
        }

        @Override
        public ProcessHapmapBlock call() throws Exception {

            myGenotypes = SuperByteMatrixBuilder.getInstance(myNumTaxa, myNumSitesToProcess);
            for (int site = 0; site < myNumSitesToProcess; site++) {
                String input = myInputLines.get(site);
                try {
                    int[] tabPos = new int[NUM_HAPMAP_NON_TAXA_HEADERS];
                    int tabIndex = 0;
                    int len = input.length();
                    for (int i = 0; (tabIndex < NUM_HAPMAP_NON_TAXA_HEADERS) && (i < len); i++) {
                        if (input.charAt(i) == '\t') {
                            tabPos[tabIndex++] = i;
                        }
                    }
                    String chrName = input.substring(tabPos[CHROMOSOME_INDEX - 1] + 1, tabPos[CHROMOSOME_INDEX]);
                    Chromosome currChr = myChromosomeLookup.get(chrName);
                    if (currChr == null) {
                        currChr = new Chromosome(new String(chrName));
                        myChromosomeLookup.put(chrName, currChr);
                    }
                    String variants = input.substring(tabPos[VARIANT_INDEX - 1] + 1, tabPos[VARIANT_INDEX]);
                    int physicalPos;
                    try {
                        physicalPos = Integer.parseInt(input.substring(tabPos[POSITION_INDEX - 1] + 1, tabPos[POSITION_INDEX]));
                    } catch (Exception ex) {
                        throw new IllegalArgumentException("BuilderFromHapMap: Position must be an integer: " + input.substring(tabPos[POSITION_INDEX - 1] + 1, tabPos[POSITION_INDEX]).trim());
                    }
                    GeneralPosition.Builder apb = new GeneralPosition.Builder(currChr, physicalPos)
                            .snpName(input.substring(0, tabPos[SNPID_INDEX]))
                            .knownVariants(variants) //TODO   strand, variants,
                            ;

                    byte glbMajor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(0));
                    apb.allele(WHICH_ALLELE.GlobalMajor, glbMajor);
                    if (variants.length() == 3) {
                        byte glbMinor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants.charAt(2));
                        apb.allele(WHICH_ALLELE.GlobalMinor, glbMinor);
                    }

                    myPositionList.add(apb.build());
                    int offset = tabPos[NUM_HAPMAP_NON_TAXA_HEADERS - 1] + 1;

                    int taxon = 0;
                    if (myIsOneLetter) {
                        for (int i = offset; i < len; i += 2) {
                            if (taxon >= myNumTaxa) {
                                throw new IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList.get(myPositionList.size() - 1).getSNPID() + " has too many values.");
                            }
                            byte value = NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i));
                            if (value == NucleotideAlignmentConstants.UNDEFINED_DIPLOID_ALLELE) {
                                throw new IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList.get(myPositionList.size() - 1).getSNPID() + " has illegal value: " + input.charAt(i));
                            }
                            myGenotypes.set(taxon++, site, value);
                        }
                    } else {
                        for (int i = offset; i < len; i += 3) {
                            if (taxon >= myNumTaxa) {
                                throw new IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList.get(myPositionList.size() - 1).getSNPID() + " has too many values.");
                            }
                            // there is a phasing conflict with the existing import approach
                            byte value = GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i + 1)),
                                    NucleotideAlignmentConstants.getNucleotideDiploidByte(input.charAt(i)));
                            if (value == NucleotideAlignmentConstants.UNDEFINED_DIPLOID_ALLELE) {
                                throw new IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList.get(myPositionList.size() - 1).getSNPID() + " has illegal value: " + input.charAt(i) + input.charAt(i + 1));
                            }
                            myGenotypes.set(taxon++, site, value);
                        }
                    }
                    if (taxon != myNumTaxa) {
                        throw new IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList.get(myPositionList.size() - 1).getSNPID() + " has too few values.");
                    }

                    swapSitesIfOutOfOrder(site);

                } catch (Exception e) {
                    myLogger.error("Error parsing this row " + input);
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException("BuilderFromHapMap: Error Parsing Line: " + input.substring(0, Math.min(25, input.length())) + "...\n" + e.getMessage());
                }
            }
            myInputLines = null;

            return this;

        }

        // Swap adjacent misordered sites, often caused by two sites at the same positions with a different name order
        private void swapSitesIfOutOfOrder(int site) {
            if (site < 1) {
                return;
            }
            if (myPositionList.get(site - 1).compareTo(myPositionList.get(site)) > 0) {
                //swap
                Position tempP = myPositionList.get(site - 1);
                myLogger.warn("Swapping:" + tempP.toString() + " <-> " + myPositionList.get(site).toString());
                myPositionList.set(site - 1, myPositionList.get(site));
                myPositionList.set(site, tempP);
                for (int t = 0; t < myGenotypes.getNumRows(); t++) {
                    byte tempG = myGenotypes.get(t, site - 1);
                    myGenotypes.set(t, site - 1, myGenotypes.get(t, site));
                    myGenotypes.set(t, site, tempG);
                }
            }

        }

        public int getNumberSitesProcessed() {
            return myNumSitesToProcess;
        }

        public SuperByteMatrix getGenotypes() {
            return myGenotypes;
        }

        public List<Position> getPositions() {
            return myPositionList;
        }

    }

    static TaxaListBuilder processTaxa(String readLn, Map<String, SetMultimap<String, String>> taxaAnnotation) {

        String[] header = WHITESPACE_PATTERN.split(readLn);
        int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
        TaxaListBuilder tlb = new TaxaListBuilder();
        for (int i = 0; i < numTaxa; i++) {
            String taxonID = header[i + NUM_HAPMAP_NON_TAXA_HEADERS];
            if (taxonID == null || taxonID.isEmpty()) {
                throw new IllegalStateException("BuilderFromHapMap: processTaxa: Taxa names should be separated by a single tab and contain no spaces.");
            }
            Taxon.Builder at = new Taxon.Builder(taxonID);
            SetMultimap<String, String> taMap = taxaAnnotation.get(taxonID);
            if (taMap != null) {
                for (Map.Entry<String, String> en : taMap.entries()) {
                    if (en.getKey().equals("ID")) {
                        continue; //skip the IDs as these became the name
                    }
                    String s = en.getValue().replace("\"", "");
                    at.addAnno(en.getKey(), s);
                }
            }
            tlb.add(at.build());
        }
        return tlb;

    }

    /**
     * Set the builder so that when built it will sort the taxa alphabetically.
     */
    public BuilderFromHapMap sortTaxa() {
        mySortTaxaAlphabetically = true;
        return this;
    }

    /**
     * Set the builder so that when built it will sort positions.
     */
    public BuilderFromHapMap sortPositions() {
        mySortPositions = true;
        return this;
    }
}
