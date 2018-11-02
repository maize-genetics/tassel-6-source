package net.maizegenetics.dna.snp.bit;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.cache.Weigher;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;

/**
 * Provides rapid conversion routines and caching from byte encoding of
 * nucleotides to bit encoding. Only two alleles are supported for each scope
 * (e.g. Major and minor, or Reference and Alternate).
 * <p>
 * </p>
 * The cache is designed to support multiple scopes, but currently scope must be
 * passed in at construction.
 * <p>
 * </p>
 * It is not clear that site or taxa optimization is needed. The code should be
 * highly parallelizable as long as the gets are not for adjacent sites.
 *
 * @author Ed Buckler
 */
public class DynamicBitStorage implements BitStorage {

    private GenotypeCallTable myGenotype;
    private final WHICH_ALLELE myWhichAllele;
    private final byte[] myPrefAllele;
    private final int myTaxaCount;
    private final int mySiteCount;
    private static final int SBoff = 58;

    private enum SB {

        TAXA(0), SITE(1);
        public final int index;

        SB(int index) {
            this.index = index;
        }
    };
    Weigher<Long, BitSet> weighByLength = new Weigher<Long, BitSet>() {
        @Override
        public int weigh(Long key, BitSet value) {
            return value.getNumWords() * 8;
        } //words are longs so *8 converts to bytes
    };
    private LoadingCache<Long, BitSet> bitCache;
    private CacheLoader<Long, BitSet> bitLoader = new CacheLoader<Long, BitSet>() {
        @Override
        public BitSet load(Long key) {
            if (getDirectionFromKey(key) == SB.TAXA) {
                int taxon = getSiteOrTaxonFromKey(key);
                if ((myPrefAllele.length == 1) && (myPrefAllele[0] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE)) {
                    return UnmodifiableBitSet.getInstance(GenotypeTableUtils.calcBitUnknownPresenceFromGenotype(myGenotype.genotypeAllSites(taxon)));
                } else {
                    return UnmodifiableBitSet.getInstance(GenotypeTableUtils.calcBitPresenceFromGenotype(myGenotype.genotypeAllSites(taxon), myPrefAllele)); //allele comp
                }
            } else {
                ArrayList<Long> toFill = new ArrayList<>();
                toFill.add(key);
                try {
                    bitCache.putAll(loadAll(toFill));
                    return bitCache.get(key);
                } catch (Exception e) {
                    e.printStackTrace();
                    return null;
                }
            }
        }

        @Override
        public Map<Long, BitSet> loadAll(Iterable<? extends Long> keys) throws Exception {
            long key = keys.iterator().next();
            //This pivoting code is needed if myGenotype is store in taxa direction
            //It runs about 7 times faster than getting base sequentially across taxa.
            HashMap<Long, BitSet> result = new HashMap<Long, BitSet>(64);
            int site = getSiteOrTaxonFromKey(key);
            int length = (mySiteCount - site < 64) ? mySiteCount - site : 64;
            byte[][] genotypeTBlock = new byte[length][myTaxaCount];
            for (int t = 0; t < myTaxaCount; t++) {
                for (int s = 0; s < genotypeTBlock.length; s++) {
                    genotypeTBlock[s][t] = myGenotype.genotype(t, site + s);
                }
            }

            if ((myPrefAllele.length == 1) && (myPrefAllele[0] == GenotypeTable.UNKNOWN_DIPLOID_ALLELE)) {
                for (int i = 0; i < length; i++) {
                    BitSet bs = UnmodifiableBitSet.getInstance(GenotypeTableUtils.calcBitUnknownPresenceFromGenotype(genotypeTBlock[i]));
                    result.put(getKey(SB.SITE, site + i), bs);
                }
            } else {
                for (int i = 0; i < length; i++) {
                    byte a1 = myPrefAllele[site + i];
                    BitSet bs = UnmodifiableBitSet.getInstance(GenotypeTableUtils.calcBitPresenceFromGenotype(genotypeTBlock[i], a1));
                    result.put(getKey(SB.SITE, site + i), bs);
                }
            }
            return result;
        }
    };

    private long getKey(SB direction, int siteOrTaxon) {
        return ((long) direction.index << SBoff) | (long) siteOrTaxon;
    }

    private int getSiteOrTaxonFromKey(long key) {
        return (int) ((key << 32) >>> 32);
    }

    private SB getDirectionFromKey(long key) {
        if (key >>> SBoff == SB.TAXA.index) {
            return SB.TAXA;
        }
        return SB.SITE;
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon) {
        try {
            return bitCache.get(getKey(SB.TAXA, taxon));
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site) {
        try {
            return bitCache.get(getKey(SB.SITE, site));
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, int startBlock, int endBlock) {
        BitSet result = allelePresenceForAllSites(taxon);
        if (result == null) {
            return new long[0];
        }
        return result.getBits(startBlock, endBlock - 1);   //BitSet is inclusive, while this method is exclusive.
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public DynamicBitStorage(GenotypeCallTable genotype, WHICH_ALLELE allele, byte[] prefAllele) {
        myGenotype = genotype;
        myWhichAllele = allele;
        mySiteCount = myGenotype.numberOfSites();
        myTaxaCount = myGenotype.numberOfTaxa();
        myPrefAllele = Arrays.copyOf(prefAllele, prefAllele.length);
        bitCache = CacheBuilder.newBuilder()
                .maximumWeight(2_000_000_000L)
                .weigher(weighByLength)
                .build(bitLoader);
    }

    public static DynamicBitStorage getUnknownInstance(GenotypeCallTable genotype) {
        return new DynamicBitStorage(genotype, WHICH_ALLELE.Unknown, new byte[]{GenotypeTable.UNKNOWN_DIPLOID_ALLELE});
    }

}
