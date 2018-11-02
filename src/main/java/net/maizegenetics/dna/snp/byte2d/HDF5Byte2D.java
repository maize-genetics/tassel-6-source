/*
 *  HDF5Byte2D
 */
package net.maizegenetics.dna.snp.byte2d;

import ch.systemsx.cisd.hdf5.IHDF5Reader;

import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.Tassel5HDF5Constants;
import net.maizegenetics.util.HDF5Utils;

import java.util.LinkedHashMap;
import java.util.Map;
import net.maizegenetics.dna.snp.score.SiteScore;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5Byte2D extends AbstractByte2D {

    private static final int MAX_CACHE_SIZE = 1 << 16;
    private static final int HDF5_BLOCK = 1 << 16;
    private final Map<Long, byte[]> myCache = new LinkedHashMap<Long, byte[]>((3 * MAX_CACHE_SIZE) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry<Long, byte[]> eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };

    private final IHDF5Reader myReader;
    private final int myNumSites;
    private final TaxaList myTaxa;

    HDF5Byte2D(IHDF5Reader reader, SiteScore.SITE_SCORE_TYPE siteScoreType) {
        super(siteScoreType, reader.int32().getAttr(Tassel5HDF5Constants.GENOTYPES_MODULE, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA),
                reader.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES)
        );
        myReader = reader;
        myNumSites = reader.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        //TODO - Maybe pass in taxa list?
        myTaxa = new TaxaListBuilder().buildFromHDF5(reader);
    }

    private static long getCacheKey(int taxon, int site) {
        return (long) taxon << 33;
    }

    private byte[] cacheValues(int taxon, long key) {
        byte[] data = HDF5Utils.getHDF5GenotypeSiteScores(myReader, myTaxa.taxaName(taxon), siteScoreType().name());
        if (data == null) {
            return null;
        } else {
            myCache.put(key, data);
            return data;
        }
    }

    @Override
    public byte valueForAllele(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        byte[] data = myCache.get(key);
        if (data == null) {
            data = cacheValues(taxon, key);
        }
        return data[site];
    }

}
