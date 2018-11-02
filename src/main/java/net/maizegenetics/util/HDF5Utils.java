package net.maizegenetics.util;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import com.google.common.base.Splitter;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.taxa.Taxon;

import java.util.*;

import static net.maizegenetics.util.Tassel5HDF5Constants.*;

/**
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public final class HDF5Utils {

    private HDF5Utils() {
    }

    public static boolean isTASSEL4HDF5Format(IHDF5Reader reader) {
        return reader.exists(net.maizegenetics.dna.snp.HapMapHDF5Constants.LOCI);
    }

    // TAXA Module
    public static void createHDF5TaxaModule(IHDF5Writer h5w) {
        h5w.object().createGroup(Tassel5HDF5Constants.TAXA_MODULE);
        h5w.bool().setAttr(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_LOCKED, false);
        h5w.string().createArray(Tassel5HDF5Constants.TAXA_ORDER, 256, 0, 1);
    }

    public static void lockHDF5TaxaModule(IHDF5Writer h5w) {
        h5w.bool().setAttr(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_LOCKED, true);
    }

    public static void unlockHDF5TaxaModule(IHDF5Writer h5w) {
        h5w.bool().setAttr(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_LOCKED, false);
    }

    public static boolean doesTaxaModuleExist(IHDF5Reader reader) {
        return reader.exists(Tassel5HDF5Constants.TAXA_MODULE);
    }

    public static boolean doTaxonCallsExist(IHDF5Reader reader, String taxonName) {
        return reader.exists(Tassel5HDF5Constants.getGenotypesCallsPath(taxonName));
    }

    public static boolean doTaxonCallsExist(IHDF5Reader reader, Taxon taxon) {
        return reader.exists(Tassel5HDF5Constants.getGenotypesCallsPath(taxon.getName()));
    }

    public static boolean isTaxaLocked(IHDF5Reader reader) {
        return reader.bool().getAttr(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_LOCKED);
    }

    /**
     * Adds a taxon to the taxon module
     *
     * @param h5w
     * @param taxon
     * @return true if new add, false if already exists
     */
    public static boolean addTaxon(IHDF5Writer h5w, Taxon taxon) {
        if (isTaxaLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        if (!h5w.exists(Tassel5HDF5Constants.TAXA_ORDER)) {
            createTaxaOrder(h5w);
        }
        String path = Tassel5HDF5Constants.getTaxonPath(taxon.getName());
        if (h5w.exists(path)) {
            return false;
        }
        h5w.object().createGroup(path);
        writeHDF5Annotation(h5w, path, taxon.getAnnotation());
        long size = h5w.getDataSetInformation(Tassel5HDF5Constants.TAXA_ORDER).getNumberOfElements();
        h5w.string().writeArrayBlockWithOffset(Tassel5HDF5Constants.TAXA_ORDER, new String[]{taxon.getName()}, 1, size);
        return true;
    }

    private static void createTaxaOrder(IHDF5Writer h5w) {
        List<String> taxaNames = getAllTaxaNames(h5w);
        h5w.string().createArray(Tassel5HDF5Constants.TAXA_ORDER, 256, 0, 1);
        for (int i = 0; i < taxaNames.size(); i++) {
            h5w.string().writeArrayBlockWithOffset(Tassel5HDF5Constants.TAXA_ORDER, new String[]{taxaNames.get(i)}, 1, i);
        }
    }

    public static void writeHDF5Annotation(IHDF5Writer writer, String path, GeneralAnnotation annotations) {
        if (annotations == null) {
            return;
        }
        Iterator<Map.Entry<String, String>> itr = annotations.getConcatenatedTextAnnotations().entrySet().iterator();
        while (itr.hasNext()) {
            Map.Entry<String, String> current = itr.next();
            writer.string().setAttr(path, current.getKey(), current.getValue());
        }
    }

    public static void replaceTaxonAnnotations(IHDF5Writer h5w, Taxon modifiedTaxon) {
        if (isTaxaLocked(h5w)) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        String path = Tassel5HDF5Constants.getTaxonPath(modifiedTaxon.getName());
        if (!h5w.exists(path)) {
            throw new IllegalStateException("HDF5Utils: replaceTaxonAnnotations: Taxon does not already exist: " + modifiedTaxon.getName());
        }
        writeHDF5Annotation(h5w, path, modifiedTaxon.getAnnotation());
    }

    public static GeneralAnnotationStorage readHDF5Annotation(IHDF5Reader reader, String path) {
        return readHDF5Annotation(reader, path, null);
    }

    public static GeneralAnnotationStorage readHDF5Annotation(IHDF5Reader reader, String path, String[] annotationKeys) {
        GeneralAnnotationStorage.Builder builder = GeneralAnnotationStorage.getBuilder();
        if (annotationKeys == null) {
            reader.object().getAllAttributeNames(path).stream().forEach((key) -> {
                for (String value : Splitter.on(",").split(reader.string().getAttr(path, key))) {
                    builder.addAnnotation(key, value);
                }
            });
        } else {
            for (String key : annotationKeys) {
                if (reader.hasAttribute(path, key)) {
                    for (String value : Splitter.on(",").split(reader.string().getAttr(path, key))) {
                        builder.addAnnotation(key, value);
                    }
                }
            }
        }
        return builder.build();
    }

    public static Taxon getTaxon(IHDF5Reader reader, String taxonName) {
        String taxonPath = Tassel5HDF5Constants.getTaxonPath(taxonName);
        if (!reader.exists(taxonPath)) {
            return null;
        }
        Taxon.Builder tb = new Taxon.Builder(taxonName);
        for (String a : reader.object().getAllAttributeNames(taxonPath)) {
            for (String s : Splitter.on(",").split(reader.string().getAttr(taxonPath, a))) {
                tb.addAnno(a, s);
            }
        }
        return tb.build();
    }

    public static List<String> getAllTaxaNames(IHDF5Reader reader) {
        List<String> taxaNames = new ArrayList<>();
        if (reader.exists(Tassel5HDF5Constants.TAXA_ORDER)) {
            for (String s : reader.readStringArray(Tassel5HDF5Constants.TAXA_ORDER)) {
                taxaNames.add(s);
            }
        } else {
            List<HDF5LinkInformation> fields = reader.object().getAllGroupMemberInformation(Tassel5HDF5Constants.TAXA_MODULE, true);
            for (HDF5LinkInformation is : fields) {
                if (!is.isGroup()) {
                    continue;
                }
                taxaNames.add(is.getName());
            }
        }
        return taxaNames;
    }

    public static void writeHDF5TaxaNumTaxa(IHDF5Writer h5w, int numTaxa) {
        h5w.int32().setAttr(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_NUM_TAXA, numTaxa);
    }

    public static int getHDF5TaxaNumTaxa(IHDF5Reader reader) {
        return reader.int32().getAttr(Tassel5HDF5Constants.TAXA_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAXA_NUM_TAXA);
    }

    // GENOTYPE Module
    public static int getHDF5GenotypeTaxaNumber(IHDF5Reader reader) {
        return reader.int32().getAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA);
    }

    public static void createHDF5GenotypeModule(IHDF5Writer h5w) {
        if (h5w.exists(Tassel5HDF5Constants.GENOTYPES_MODULE)) {
            throw new UnsupportedOperationException("Genotypes module already exists in HDF5 file");
        }
        h5w.object().createGroup(Tassel5HDF5Constants.GENOTYPES_MODULE);
        h5w.bool().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_LOCKED, false);
    }

    public static void lockHDF5GenotypeModule(IHDF5Writer h5w) {
        h5w.bool().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_LOCKED, true);
    }

    public static void unlockHDF5GenotypeModule(IHDF5Writer h5w) {
        h5w.bool().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_LOCKED, false);
    }

    public static boolean doesGenotypeModuleExist(IHDF5Reader reader) {
        return reader.exists(Tassel5HDF5Constants.GENOTYPES_MODULE);
    }

    public static boolean isHDF5GenotypeLocked(IHDF5Reader reader) {
        if (reader.exists(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH) == false) {
            return false;
        }
        return reader.bool().getAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_LOCKED);
    }

    public static void writeHDF5GenotypesMaxNumAlleles(IHDF5Writer h5w, int maxNumAlleles) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        h5w.int32().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, maxNumAlleles);
    }

    public static void writeHDF5GenotypesRetainRareAlleles(IHDF5Writer h5w, boolean retainRareAlleles) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        h5w.bool().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_RETAIN_RARE_ALLELES, retainRareAlleles);
    }

    public static void writeHDF5GenotypesNumTaxa(IHDF5Writer h5w, int numTaxa) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        h5w.int32().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA, numTaxa);
    }

    public static void writeHDF5GenotypesScoreType(IHDF5Writer h5w, String scoreType) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        h5w.string().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, scoreType);
    }

    public static void writeHDF5GenotypesAlleleStates(IHDF5Writer h5w, String[][] aEncodings) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        int numEncodings = aEncodings.length;
        int numStates = aEncodings[0].length;
        MDArray<String> alleleEncodings = new MDArray<>(String.class, new int[]{numEncodings, numStates});
        for (int s = 0; s < numEncodings; s++) {
            for (int x = 0; x < numStates; x++) {
                alleleEncodings.set(aEncodings[s][x], s, x);
            }
        }
        h5w.string().createMDArray(Tassel5HDF5Constants.GENOTYPES_ALLELE_STATES, 100, new int[]{numEncodings, numStates});
        h5w.string().writeMDArray(Tassel5HDF5Constants.GENOTYPES_ALLELE_STATES, alleleEncodings);
    }

    public static byte[] getHDF5GenotypesCalls(IHDF5Reader reader, String taxon) {
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        return reader.readAsByteArray(callsPath);
    }

    public static void writeHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if (h5w.exists(callsPath)) {
            throw new IllegalStateException("Taxa Calls Already Exists:" + taxon);
        }
        h5w.int8().createArray(callsPath, calls.length, Math.min(Tassel5HDF5Constants.BLOCK_SIZE, calls.length), Tassel5HDF5Constants.intDeflation);
        writeHDF5EntireArray(callsPath, h5w, calls.length, Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static void replaceHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if (!h5w.exists(callsPath)) {
            throw new IllegalStateException("Taxa Calls Do Not Already Exists to replace");
        }
        writeHDF5EntireArray(callsPath, h5w, calls.length, Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static void replaceHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, int startSite, byte[] calls) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        String callsPath = Tassel5HDF5Constants.getGenotypesCallsPath(taxon);
        if (!h5w.exists(callsPath)) {
            throw new IllegalStateException("Taxa Calls Do Not Already Exists to replace");
        }
        if (startSite % Tassel5HDF5Constants.BLOCK_SIZE != 0) {
            throw new IllegalStateException("Taxa Calls Start Site not a multiple of the block size");
        }
        writeHDF5Block(callsPath, h5w, Tassel5HDF5Constants.BLOCK_SIZE, startSite / Tassel5HDF5Constants.BLOCK_SIZE, calls);
    }

    public static byte[][] getHDF5GenotypesDepth(IHDF5Reader reader, String taxon) {
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if (reader.exists(callsPath)) {
            return reader.int8().readMatrix(callsPath);
        } else {
            return null;
        }
    }

    public static boolean doesGenotypeDepthExist(IHDF5Reader reader) {
        List<String> taxaNames = getAllTaxaNames(reader);
        boolean depthExist = false;
        for (String taxon : taxaNames) {
            String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
            if (reader.exists(callsPath)) {
                depthExist = true;
                break;
            }
        }
        return depthExist;
    }

    public static void writeHDF5GenotypesDepth(IHDF5Writer h5w, String taxon, byte[][] depth) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if (h5w.exists(callsPath)) {
            throw new IllegalStateException("Taxa Depth Already Exists:" + taxon);
        }
        h5w.int8().createMatrix(callsPath, depth.length, depth[0].length, 6, Math.min(Tassel5HDF5Constants.BLOCK_SIZE, depth[0].length), Tassel5HDF5Constants.intDeflation);
        writeHDF5EntireArray(callsPath, h5w, depth[0].length, Tassel5HDF5Constants.BLOCK_SIZE, depth);
    }

    public static void replaceHDF5GenotypesDepth(IHDF5Writer h5w, String taxon, byte[][] depth) {
        if (isHDF5GenotypeLocked(h5w) == true) {
            throw new UnsupportedOperationException("Trying to write to a locked HDF5 file");
        }
        String callsPath = Tassel5HDF5Constants.getGenotypesDepthPath(taxon);
        if (!h5w.exists(callsPath)) {
            throw new IllegalStateException("Taxa Depth Does Not Already Exists to Replace");
        }
        writeHDF5EntireArray(callsPath, h5w, depth[0].length, Tassel5HDF5Constants.BLOCK_SIZE, depth);
    }

    public static byte[] getHDF5Alleles(IHDF5Reader reader, WHICH_ALLELE allele) {
        return reader.int8().readMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 1, getHDF5PositionNumber(reader),
                (long) allele.index(), 0)[0];
    }

    public static byte[] getHDF5GenotypeSiteScores(IHDF5Reader reader, String taxon, String siteScoreType) {
        String path = Tassel5HDF5Constants.getGenotypesSiteScorePath(taxon, siteScoreType);
        if (reader.exists(path)) {
            return reader.int8().readArray(path);
        } else {
            return null;
        }
    }

    public static void writeHDF5GenotypeSiteScores(IHDF5Writer writer, String taxon, String siteScoreType, byte[] values) {
        String path = Tassel5HDF5Constants.getGenotypesSiteScorePath(taxon, siteScoreType);
        if (writer.exists(path)) {
            throw new IllegalStateException("HDF5Utils: writeHDF5GenotypeSiteScores: path already exists: " + path);
        } else {
            writer.int8().createArray(path, values.length, Math.min(Tassel5HDF5Constants.BLOCK_SIZE, values.length), Tassel5HDF5Constants.intDeflation);
            writeHDF5EntireArray(path, writer, values.length, Tassel5HDF5Constants.BLOCK_SIZE, values);
        }
    }

    // Positions/numSites
    public static int getHDF5PositionNumber(IHDF5Reader reader) {
        return reader.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
    }

    public static void createHDF5PositionModule(IHDF5Writer h5w) {
        h5w.object().createGroup(Tassel5HDF5Constants.POSITION_MODULE);
    }

    public static void writeHDF5PositionNumSite(IHDF5Writer h5w, int numSites) {
        h5w.int32().setAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES, numSites);
    }

    public static byte[] getHDF5ReferenceAlleles(IHDF5Reader reader) {
        return getHDF5Alleles(reader, Tassel5HDF5Constants.REF_ALLELES, 0, getHDF5PositionNumber(reader));
    }

    public static byte[] getHDF5ReferenceAlleles(IHDF5Reader reader, int startSite, int numSites) {
        return getHDF5Alleles(reader, Tassel5HDF5Constants.REF_ALLELES, startSite, numSites);
    }

    public static byte[] getHDF5AncestralAlleles(IHDF5Reader reader) {
        return getHDF5Alleles(reader, Tassel5HDF5Constants.ANC_ALLELES, 0, getHDF5PositionNumber(reader));
    }

    public static byte[] getHDF5AncestralAlleles(IHDF5Reader reader, int startSite, int numSites) {
        return getHDF5Alleles(reader, Tassel5HDF5Constants.ANC_ALLELES, startSite, numSites);
    }

    private static byte[] getHDF5Alleles(IHDF5Reader reader, String allelePath, int startSite, int numSites) {
        if (reader.exists(allelePath)) {
            return reader.int8().readArrayBlockWithOffset(allelePath, numSites, startSite);
        }
        byte[] unknown = new byte[numSites];
        Arrays.fill(unknown, GenotypeTable.UNKNOWN_ALLELE);
        return unknown;
    }

    // Tags Module
    public static void createHDF5TagModule(IHDF5Writer h5w, int tagLengthInLong) {
        if (h5w.exists(Tassel5HDF5Constants.TAG_MODULE)) {
            throw new UnsupportedOperationException("Tag module already exists in HDF5 file");
        }
        h5w.object().createGroup(Tassel5HDF5Constants.TAG_MODULE);
        h5w.bool().setAttr(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LOCKED, false);
        h5w.int32().setAttr(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LENGTH_LONG, tagLengthInLong);
        h5w.int32().setAttr(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_COUNT, 0);
    }

    public static boolean isHDF5TagLocked(IHDF5Reader reader) {
        if (reader.exists(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH + "/" + Tassel5HDF5Constants.TAG_LOCKED) == false) {
            return false;
        }
        return reader.bool().getAttr(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LOCKED);
    }

    public static int getHDF5TagCount(IHDF5Reader reader) {
        return reader.int32().getAttr(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_COUNT);
    }

    public static int getHDF5TagLengthInLong(IHDF5Reader reader) {
        return reader.int32().getAttr(Tassel5HDF5Constants.TAG_ATTRIBUTES_PATH, Tassel5HDF5Constants.TAG_LENGTH_LONG);
    }

    public static boolean doTagsExist(IHDF5Reader reader) {
        return reader.exists(Tassel5HDF5Constants.TAGS);
    }

    public static String getTagPath(Tag tag) {
        return null;
    }

    public static long[][] getTags(IHDF5Reader reader) {
        return reader.readLongMatrix(Tassel5HDF5Constants.TAGS);
    }

    //TODO invoke all is kicking this off but it is stopping the process before complete.  Only with synchronized is
    //it completing
    public static synchronized void writeTagDistributionBucket(IHDF5Writer h5w, int bucket, long[][] tags, short[] length,
            int[] encodedTaxaDist, int maxTaxa, int[] tagDistOffset) {
        String path = Tassel5HDF5Constants.TAG_MODULE + "/" + bucket + "/";
        h5w.object().createGroup(path);
        h5w.int64().writeMatrix(path + TAG_SEQ, tags, intDeflation);
        h5w.int16().writeArray(path + TAG_LENGTHS, length, intDeflation);
        h5w.int32().createArray(path + TAG_DIST, encodedTaxaDist.length, Math.min(BLOCK_SIZE, encodedTaxaDist.length), intDeflation);
        h5w.int32().writeArray(path + TAG_DIST, encodedTaxaDist, intDeflation);
        h5w.int32().setAttr(path + TAG_DIST, "MaxTaxa", maxTaxa);
        h5w.int32().createArray(path + TAG_DIST_OFFSETS, tagDistOffset.length, intDeflation);
        h5w.int32().writeArray(path + TAG_DIST_OFFSETS, tagDistOffset, intDeflation);
    }

    public static boolean doTagsByTaxaExist(IHDF5Reader reader) {
        throw new UnsupportedOperationException("Not implemented yet");
        //return reader.exists(Tassel5HDF5Constants.TAG_DIST);
    }

    /**
     *
     * @param objectPath
     * @param myWriter
     * @param objMaxLength
     * @param blockSize
     * @param val
     */
    public static void writeHDF5EntireArray(String objectPath, IHDF5Writer myWriter, int objMaxLength, int blockSize, Object val) {
        int blocks = ((objMaxLength - 1) / blockSize) + 1;
        for (int block = 0; block < blocks; block++) {
            int startPos = block * blockSize;
            int length = Math.min(objMaxLength - startPos, blockSize);
            if (val instanceof byte[][]) {
                byte[][] oval = (byte[][]) val;
                byte[][] sval = new byte[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof int[][]) {
                int[][] oval = (int[][]) val;
                int[][] sval = new int[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof long[][]) {
                long[][] oval = (long[][]) val;
                long[][] sval = new long[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof byte[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((byte[]) val, startPos, startPos + length));
            } else if (val instanceof float[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((float[]) val, startPos, startPos + length));
            } else if (val instanceof int[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((int[]) val, startPos, startPos + length));
            } else if (val instanceof String[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((String[]) val, startPos, startPos + length));
            }
        }
    }

    /**
     *
     * @param objectPath
     * @param myWriter
     * @param blockSize
     * @param block
     * @param val
     */
    public static void writeHDF5Block(String objectPath, IHDF5Writer myWriter, int blockSize, int block, Object val) {
        int startPos = block * blockSize;
        if (val instanceof byte[][]) {
            byte[][] bval = (byte[][]) val;
            myWriter.int8().writeMatrixBlockWithOffset(objectPath, bval, bval.length, bval[0].length, 0l, (long) startPos);
        } else if (val instanceof byte[]) {
            byte[] fval = (byte[]) val;
            myWriter.int8().writeArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof float[]) {
            float[] fval = (float[]) val;
            myWriter.float32().writeArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[]) {
            int[] fval = (int[]) val;
            myWriter.int32().writeArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[][]) {
            int[][] ival = (int[][]) val;
            myWriter.int32().writeMatrixBlockWithOffset(objectPath, ival, ival.length, ival[0].length, 0l, (long) startPos);
        } else if (val instanceof long[][]) {
            long[][] lval = (long[][]) val;
            myWriter.int64().writeMatrixBlockWithOffset(objectPath, lval, lval.length, lval[0].length, 0l, (long) startPos);
        } else if (val instanceof String[]) {
            String[] sval = (String[]) val;
            myWriter.string().writeArrayBlockWithOffset(objectPath, sval, sval.length, (long) startPos);
        }
    }

}
