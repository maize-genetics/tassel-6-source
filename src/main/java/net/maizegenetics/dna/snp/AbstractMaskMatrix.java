/*
 *  AbstractMaskMatrix
 *
 *  Created on Jan 6, 2017
 */
package net.maizegenetics.dna.snp;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import net.maizegenetics.util.BitSet;

import java.util.concurrent.ForkJoinPool;

/**
 * @author Terry Casstevens
 */
public abstract class AbstractMaskMatrix implements MaskMatrix {

    private static final int NUM_SITES_PER_BLOCK = 10;

    protected final int myNumTaxa;
    protected final int myNumSites;
    private final Cache<Integer, BitSet> myCache;
    private final Cache<Integer, BitSet> mySmallCache;
    private final ForkJoinPool myThreadPool = ForkJoinPool.commonPool();

    AbstractMaskMatrix(int numTaxa, int numSites) {
        myNumTaxa = numTaxa;
        myNumSites = numSites;
        myCache = CacheBuilder.newBuilder()
                .initialCapacity(1000)
                .maximumSize(1000)
                .build();
        mySmallCache = CacheBuilder.newBuilder()
                .initialCapacity(100)
                .maximumSize(100)
                .build();
    }

    protected abstract BitSet siteMask(int site);

    protected abstract BitSet taxonMask(int taxon);

    protected abstract boolean isMasked(int taxon, int site);

    private BitSet getFromCache(int site) {
        BitSet result = myCache.getIfPresent(site);
        if (result != null) {
            myThreadPool.submit(new ProcessSite(siteBlockKey(site + NUM_SITES_PER_BLOCK)));
            return result;
        } else {
            myThreadPool.submit(new ProcessSite(siteBlockKey(site)));
            myThreadPool.submit(new ProcessSite(siteBlockKey(site + NUM_SITES_PER_BLOCK)));
            result = siteMask(site);
            myCache.put(site, result);
            return result;
        }
    }

    private int siteBlockKey(int site) {
        return site / NUM_SITES_PER_BLOCK * NUM_SITES_PER_BLOCK;
    }

    @Override
    public boolean get(int taxon, int site) {
        BitSet temp = mySmallCache.getIfPresent(site);
        if (temp == null) {
            temp = getFromCache(site);
            mySmallCache.put(site, temp);
        }
        return temp.fastGet(taxon);
    }

    @Override
    public boolean isTaxonMaskedHint(int taxon) {
        return true;
    }

    @Override
    public BitSet maskForTaxon(int taxon) {
        return taxonMask(taxon);
    }

    @Override
    public boolean isSiteMaskedHint(int site) {
        return getFromCache(site).cardinality() != 0;
    }

    @Override
    public BitSet maskForSite(int site) {
        return siteMask(site);
    }

    @Override
    public int numTaxa() {
        return myNumTaxa;
    }

    @Override
    public int numSites() {
        return myNumSites;
    }

    @Override
    public boolean isSiteOptimized() {
        return true;
    }

    private class ProcessSite implements Runnable {

        private final int mySite;

        public ProcessSite(int site) {
            mySite = site;
        }

        @Override
        public void run() {
            int end = Math.min(mySite + NUM_SITES_PER_BLOCK, myNumSites);
            for (int s = mySite; s < end; s++) {
                if (myCache.getIfPresent(s) == null) {
                    myCache.put(s, siteMask(s));
                } else {
                    return;
                }
            }
        }

    }

}
