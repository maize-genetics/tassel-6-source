package net.maizegenetics.dna.snp.io;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.score.AlleleDepthBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

/**
 * This builder reads a VCF formatted file and creates a {@link GenotypeTable}
 *
 * @author Terry Casstevens
 * Created February 05, 2018
 */
public class BuilderFromVCFUsingHTSJDK {

    private static final Logger myLogger = Logger.getLogger(BuilderFromVCFUsingHTSJDK.class);

    private static final int NUM_VARIANT_CONTEXT_PER_THREAD = 100;

    private final VCFFileReader myReader;
    private final VCFHeader myHeader;
    private final Iterator<VariantContext> myVariants;
    private final TaxaList myTaxa;
    private boolean myKeepDepth = true;
    private ProgressListener myProgressListener = null;

    private BuilderFromVCFUsingHTSJDK(String vcf) {
        File vcfFile = new File(vcf);
        if (!vcfFile.isFile()) {
            throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: init: file doesn't exist: " + vcf);
        }

        myReader = new VCFFileReader(vcfFile, false);
        myHeader = myReader.getFileHeader();
        myVariants = myReader.iterator();
        myTaxa = taxa();
    }

    private BuilderFromVCFUsingHTSJDK(VCFHeader header, Iterator<VariantContext> variants) {
        myReader = null;
        myHeader = header;
        myVariants = variants;
        myTaxa = taxa();
    }

    private TaxaList taxa() {
        return myHeader.getGenotypeSamples().stream()
                .map(name -> new Taxon(name))
                .collect(TaxaList.collect());
    }

    private void close() {

        try {
            if (myVariants instanceof CloseableIterator) {
                ((CloseableIterator) myVariants).close();
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            // do nothing, trying to close variants iterator if possible
        }

        try {
            if (myReader != null) {
                myReader.close();
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            // do nothing, trying to close reader if possible
        }

    }

    public static BuilderFromVCFUsingHTSJDK instance(String vcf) {
        return new BuilderFromVCFUsingHTSJDK(vcf);
    }

    public static BuilderFromVCFUsingHTSJDK instance(VCFHeader header, Iterator<VariantContext> variants) {
        return new BuilderFromVCFUsingHTSJDK(header, variants);
    }

    public static BuilderFromVCFUsingHTSJDK instance(VCFHeader header, List<VariantContext> variants) {
        return new BuilderFromVCFUsingHTSJDK(header, variants.iterator());
    }

    public static GenotypeTable read(String vcf) {
        return new BuilderFromVCFUsingHTSJDK(vcf).build();
    }

    public static GenotypeTable read(VCFHeader header, List<VariantContext> variants) {
        return new BuilderFromVCFUsingHTSJDK(header, variants.iterator()).build();
    }

    public BuilderFromVCFUsingHTSJDK progressListener(ProgressListener listener) {
        myProgressListener = listener;
        return this;
    }

    public BuilderFromVCFUsingHTSJDK keepDepth(boolean keep) {
        myKeepDepth = keep;
        return this;
    }

    public GenotypeTable build() {

        ForkJoinPool threadPool = ForkJoinPool.commonPool();

        int numTaxa = myTaxa.numberOfTaxa();

        try {

            List<Future<ProcessVariantContext>> futures = new ArrayList<>();

            PositionListBuilder positionListBuilder = new PositionListBuilder();
            List<VariantContext> contextsToProcess = new ArrayList<>();
            int numContextsToProcess = 0;
            List<Tuple<Integer, Short>> positionsToProcess = new ArrayList<>();
            int previousEnd = -1;
            Position previousPosition = null;
            Chromosome previousChromosome = null;
            while (myVariants.hasNext()) {

                numContextsToProcess++;

                VariantContext context = myVariants.next();
                // This explicitly loads the variant context in the main thread, so that
                // it will work correctly in the processing threads
                if (context.getGenotypes() instanceof LazyGenotypesContext) {
                    ((LazyGenotypesContext) context.getGenotypes()).decode();
                }

                Chromosome chr = Chromosome.instance(context.getContig());
                if (previousChromosome == null) {
                    previousChromosome = chr;
                }
                if (!chr.equals(previousChromosome)) {
                    previousChromosome = chr;
                    previousEnd = -1;
                    previousPosition = null;

                    ProcessVariantContext process = new ProcessVariantContext(contextsToProcess, numTaxa, positionsToProcess);
                    futures.add(threadPool.submit(process));
                    numContextsToProcess = 0;
                    positionsToProcess = new ArrayList<>();
                    contextsToProcess = new ArrayList<>();
                }
                int start = context.getStart();
                int end = context.getEnd();

                int refLength = context.getLengthOnReference();
                List<Allele> possibleAlleles = context.getAlleles();
                List<String> knownVariants = possibleAlleles.stream()
                        .map(Allele::getBaseString)
                        .collect(Collectors.toList());
                short maxLength = knownVariants.stream()
                        .map(String::length)
                        .max(Integer::compare)
                        .get().shortValue();

                if (refLength >= maxLength) {

                    for (int position = start; position <= end; position++) {

                        GeneralPosition.Builder posBuilder = new GeneralPosition.Builder(chr, position);
                        if (!context.getID().equals(".")) {
                            posBuilder.snpName(context.getID());
                        }
                        if (position == start) {
                            posBuilder.knownVariants(knownVariants.toArray(new String[knownVariants.size()]));
                        }
                        Position pos = posBuilder.build();

                        if (previousPosition == null || pos.compareTo(previousPosition) > 0) {
                            positionListBuilder.add(pos);
                            previousPosition = pos;
                            positionsToProcess.add(new Tuple<>(position, (short) 0));
                        }

                    }

                } else {

                    if (refLength != 1) {
                        throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: build: expected reference to be length 1 for insertion record: chr: " + chr.getName() + " start position: " + context.getStart());
                    }

                    int position = context.getStart();
                    GeneralPosition.Builder posBuilder = new GeneralPosition.Builder(chr, position);
                    if (!context.getID().equals(".")) {
                        posBuilder.snpName(context.getID());
                    }
                    posBuilder.knownVariants(knownVariants.toArray(new String[knownVariants.size()]));
                    Position pos = posBuilder.build();

                    if (previousPosition == null || pos.compareTo(previousPosition) > 0) {
                        positionListBuilder.add(pos);
                        previousPosition = pos;
                        positionsToProcess.add(new Tuple<>(position, (short) 0));
                    }

                    for (short i = 1; i < maxLength; i++) {
                        Position insertionPos = new GeneralPosition.Builder(pos).insertionPosition(i).build();
                        if (previousPosition == null || insertionPos.compareTo(previousPosition) > 0) {
                            positionListBuilder.add(insertionPos);
                            previousPosition = insertionPos;
                            positionsToProcess.add(new Tuple<>(position, i));
                        }
                    }

                }

                contextsToProcess.add(context);

                if (numContextsToProcess >= NUM_VARIANT_CONTEXT_PER_THREAD && start > previousEnd) {
                    ProcessVariantContext process = new ProcessVariantContext(contextsToProcess, numTaxa, positionsToProcess);
                    futures.add(threadPool.submit(process));
                    numContextsToProcess = 0;
                    positionsToProcess = new ArrayList<>();
                    contextsToProcess = new ArrayList<>();
                }

                previousEnd = end;

            }

            if (!contextsToProcess.isEmpty()) {
                ProcessVariantContext process = new ProcessVariantContext(contextsToProcess, numTaxa, positionsToProcess);
                futures.add(threadPool.submit(process));
            }

            PositionList positions = positionListBuilder.build();
            int numSites = positions.numberOfSites();

            GenotypeCallTableBuilder genotypes = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numTaxa, numSites);

            AlleleDepthBuilder depthBuilder = null;
            if (myKeepDepth) {
                depthBuilder = AlleleDepthBuilder.getInstance(numTaxa, numSites, myTaxa);
            }

            int numFutures = futures.size();
            int count = 0;
            int currentSite = 0;
            for (Future<ProcessVariantContext> future : futures) {

                ProcessVariantContext process = future.get();

                SuperByteMatrix genotypesBlock = process.genotypes();
                for (int t = 0; t < genotypesBlock.getNumRows(); t++) {
                    for (int s = 0; s < genotypesBlock.getNumColumns(); s++) {
                        genotypes.setBase(t, currentSite + s, genotypesBlock.get(t, s));
                    }
                }

                if (myKeepDepth) {
                    byte[][][] depth = process.depths();
                    for (int t = 0; t < depth.length; t++) {
                        depthBuilder.setDepthRangeForTaxon(t, currentSite, depth[t]);
                    }
                }

                currentSite += genotypesBlock.getNumColumns();

                if (myProgressListener != null) {
                    count++;
                    myProgressListener.progress(count * 100 / numFutures, null);
                }
            }

            if (myKeepDepth) {
                return GenotypeTableBuilder.getInstance(genotypes.build(), positions, myTaxa, depthBuilder.build());
            } else {
                return GenotypeTableBuilder.getInstance(genotypes.build(), positions, myTaxa);
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: build: problem building GenotypeTable\n" + e.getMessage());
        } finally {
            close();
        }

    }

    /**
     * This class processes a list of {@link VariantContext}
     * producing the genotype calls and depths formats for {@link GenotypeTable}
     */
    private class ProcessVariantContext implements Callable<ProcessVariantContext> {

        private final List<VariantContext> myContexts;
        private final int myNumTaxa;
        private final int myNumSitesToProcess;
        private final List<Tuple<Integer, Short>> myPositionsToProcess;
        private SuperByteMatrix myGenotypes = null;
        private byte[][][] myDepths = null;

        public ProcessVariantContext(List<VariantContext> contexts, int numTaxa, List<Tuple<Integer, Short>> positionsToProcess) {
            myContexts = contexts;
            myNumTaxa = numTaxa;
            myPositionsToProcess = positionsToProcess;
            myNumSitesToProcess = positionsToProcess.size();
        }

        @Override
        public ProcessVariantContext call() {

            myGenotypes = SuperByteMatrixBuilder.getInstance(myNumTaxa, myNumSitesToProcess);
            myGenotypes.setAll(GenotypeTable.UNKNOWN_DIPLOID_ALLELE);

            if (myKeepDepth) {
                myDepths = new byte[myNumTaxa][6][myNumSitesToProcess];
            }

            for (VariantContext context : myContexts) {

                try {

                    if (myNumTaxa != context.getNSamples()) {
                        throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: call: number of taxa: " + myNumTaxa + " should equal samples in context: " + context.getNSamples());
                    }

                    int currentPosition = context.getStart();

                    List<Allele> possibleAlleles = context.getAlleles();

                    // Get maximum length of possible allele for this variant context
                    short maxLength = possibleAlleles.stream()
                            .map(Allele::getBaseString)
                            .map(String::length)
                            .max(Integer::compare)
                            .get().shortValue();

                    // This holds translation from HTSJDK Allele to
                    // TASSEL nucleotide byte codes
                    Map<Allele, byte[]> alleleToBytes = new HashMap<>();

                    for (Allele allele : possibleAlleles) {
                        byte[] result = new byte[maxLength];
                        String value = allele.getBaseString();
                        if (value.equals("*")) {
                            Arrays.fill(result, NucleotideAlignmentConstants.GAP_ALLELE);
                        } else {
                            for (int i = 0; i < value.length(); i++) {
                                result[i] = NucleotideAlignmentConstants.getNucleotideAlleleByte(value.charAt(i));
                            }
                            for (int i = value.length(); i < maxLength; i++) {
                                result[i] = NucleotideAlignmentConstants.GAP_ALLELE;
                            }
                        }
                        alleleToBytes.put(allele, result);
                    }

                    byte[] unknownAllele = new byte[maxLength];
                    Arrays.fill(unknownAllele, GenotypeTable.UNKNOWN_ALLELE);

                    for (int t = 0; t < myNumTaxa; t++) {

                        List<Allele> alleles = context.getGenotype(t).getAlleles();
                        int numAlleles = alleles.size();
                        if (numAlleles != 1 && numAlleles != 2) {
                            throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: call: not haploid or diploid: id: " + context.getID() + " start: " + context.getStart() + " taxon: " + t + ": " + myTaxa.get(t).getName() + " allele size: " + numAlleles);
                        }

                        byte[] first = alleleToBytes.get(alleles.get(0));
                        if (first == null) {
                            first = unknownAllele;
                        }
                        byte[] second = null;
                        if (numAlleles == 2) {
                            second = alleleToBytes.get(alleles.get(1));
                            if (second == null) {
                                second = unknownAllele;
                            }
                        } else {
                            second = first;
                        }

                        int index = myPositionsToProcess.indexOf(new Tuple<>(currentPosition, (short) 0));
                        for (short i = 0; i < maxLength; i++) {
                            myGenotypes.set(t, index, (byte) ((first[i] << 4) | second[i]));
                            if (myKeepDepth) {
                                int[] alleleDepths = context.getGenotype(t).getAD();
                                if (alleleDepths != null) {
                                    if (alleleDepths.length != possibleAlleles.size()) {
                                        throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: call: number allele depths (AD): " + alleleDepths.length + " doesn't equal number alleles: " + possibleAlleles.size() + " position: " + context.getStart() + " taxa: " + t + " depths: " + Arrays.toString(alleleDepths));
                                    }
                                    for (int d = 0; d < alleleDepths.length; d++) {
                                        myDepths[t][alleleToBytes.get(possibleAlleles.get(d))[i]][index] += alleleDepths[d];
                                    }
                                }
                            }
                            index++;
                        }

                    }

                } catch (Exception e) {
                    myLogger.debug(e.getMessage(), e);
                    throw new IllegalStateException("BuilderFromVCFUsingHTSJDK: call: problem with id: " + context.getID() + "  start position: " + context.getStart() + "\n" + e.getMessage());
                }

            }

            return this;

        }

        public SuperByteMatrix genotypes() {
            return myGenotypes;
        }

        public byte[][][] depths() {
            return myDepths;
        }
    }
}
