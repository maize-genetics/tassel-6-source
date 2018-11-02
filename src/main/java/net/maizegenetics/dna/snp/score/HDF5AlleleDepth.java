/*
 *  HDF5AlleleDepth
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import java.util.Collection;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.LinkedHashMap;
import java.util.Map;
import net.maizegenetics.dna.snp.byte2d.Byte2D;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5AlleleDepth extends AlleleDepth {

    private static int MAX_CACHE_SIZE = 1 << 16;
    private static final int HDF5_BLOCK = 1 << 16;
    private final Map<Long, byte[][]> myDepthCache = new LinkedHashMap<Long, byte[][]>((3 * MAX_CACHE_SIZE) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };

    private final IHDF5Reader myReader;
    private final int myNumSites;
    private final TaxaList myTaxa;

    HDF5AlleleDepth(IHDF5Reader reader) {
        super(reader.int32().getAttr(Tassel5HDF5Constants.GENOTYPES_MODULE, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA),
                reader.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES)
        );
        myReader = reader;
        myNumSites = reader.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        myTaxa = new TaxaListBuilder().buildFromHDF5Genotypes(reader);
    }

    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_BLOCK);
    }

    /**
     * Returns the depth values (byte representation) of all nucleotides at
     * given taxon and site. Depth values are stored in bytes and translated to
     * integer using AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     *
     * @return depths
     */
    @Override
    public byte[] valuesByte(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        byte[][] data = myDepthCache.get(key);
        if (data == null) {
            data = cacheDepthBlock(taxon, site, key);
        }
        byte[] result = new byte[6];
        for (int i = 0; i < 6; i++) {
            result[i] = data[i][site % MAX_CACHE_SIZE];
        }
        return result;
    }

    /**
     * Returns depth values (byte representation) of all nucleotides and sites
     * for given taxon. The first dimension of returned array is nucleotides
     * (ALLELE_DEPTH_TYPES) and second dimension is sites.
     *
     * @param taxon taxon
     *
     * @return depths
     */
    @Override
    public byte[][] valuesForTaxonByte(int taxon) {
        byte[][] result = new byte[6][myNumSites];
        for (int site = 0; site < myNumSites; site++) {
            long key = getCacheKey(taxon, site);
            byte[][] data = myDepthCache.get(key);
            if (data == null) {
                data = cacheDepthBlock(taxon, site, key);
            }
            for (int i = 0; i < 6; i++) {
                result[i][site] = data[i][site % MAX_CACHE_SIZE];
            }
        }
        return result;
    }

    private byte[][] cacheDepthBlock(int taxon, int site, long key) {
        int start = (site / MAX_CACHE_SIZE) * MAX_CACHE_SIZE;
        int realSiteCache = (myNumSites - start < MAX_CACHE_SIZE) ? myNumSites - start : MAX_CACHE_SIZE;
        byte[][] data = myReader.int8().readMatrixBlockWithOffset(Tassel5HDF5Constants.getGenotypesDepthPath(myTaxa.taxaName(taxon)), 6, realSiteCache, 0, start);
        if (data == null) {
            return null;
        }
        myDepthCache.put(key, data);
        return data;
    }

    /**
     * Returns the depth of nucleotide (scoreType) at given taxon and site.
     * Depth values are stored in bytes and translated to integer using
     * AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     * @param scoreType nucleotide
     *
     * @return depth
     */
    @Override
    public int value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        return AlleleDepthUtil.depthByteToInt(valueByte(taxon, site, scoreType));
    }

    /**
     * Returns the depth values of all nucleotides at given taxon and site.
     * Depth values are stored in bytes and translated to integer using
     * AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     *
     * @return depths
     */
    @Override
    public int[] values(int taxon, int site) {
        return AlleleDepthUtil.depthByteToInt(valuesByte(taxon, site));
    }

    /**
     * Returns the depth (byte representation) of nucleotide (scoreType) at
     * given taxon and site. Depth values are stored in bytes and translated to
     * integer using AlleleDepthUtil.depthByteToInt().
     *
     * @param taxon taxon
     * @param site site
     * @param scoreType nucleotide
     *
     * @return depth
     */
    @Override
    public byte valueByte(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        long key = getCacheKey(taxon, site);
        byte[][] data = myDepthCache.get(key);
        if (data == null) {
            data = cacheDepthBlock(taxon, site, key);
        }
        return data[scoreType.getIndex()][site % MAX_CACHE_SIZE];
    }

    @Override
    Collection<Byte2D> byteStorage() {
        return null;
    }

}
