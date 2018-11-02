/*
 *  FilterAlleleDepth
 * 
 *  Created on Dec 16, 2016
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.Translate;

/**
 * @author Terry Casstevens
 */
public class FilterAlleleDepth extends AlleleDepth {

    final AlleleDepth myBase;
    final Translate myTranslate;

    FilterAlleleDepth(AlleleDepth alleleDepth, Translate translate) {
        super(translate.numTaxa(), translate.numSites());
        myBase = alleleDepth;
        myTranslate = translate;
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
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return new byte[NUM_ALLELE_DEPTH_TYPES];
        }
        return myBase.valuesByte((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF));
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
        if (!myTranslate.hasSiteTranslations()) {
            int translatedTaxon = myTranslate.taxon(taxon);
            if (translatedTaxon == -1) {
                return new byte[NUM_ALLELE_DEPTH_TYPES][numSites()];
            }
            return myBase.valuesForTaxonByte(translatedTaxon);
        } else {
            return super.valuesForTaxonByte(taxon);
        }
    }

    /**
     * Returns depth values (byte representation) of all nucleotides and taxa
     * for given site. The first dimension of returned array is nucleotides
     * (ALLELE_DEPTH_TYPES) and second dimension is taxa.
     *
     * @param site site
     *
     * @return depths
     */
    @Override
    public byte[][] valuesForSiteByte(int site) {
        if (!myTranslate.hasTaxaTranslations()) {
            int translatedSite = myTranslate.site(site);
            if (translatedSite == -1) {
                return new byte[NUM_ALLELE_DEPTH_TYPES][numTaxa()];
            }
            return myBase.valuesForSiteByte(myTranslate.site(site));
        } else {
            return super.valuesForSiteByte(site);
        }
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
    public byte valueByte(int taxon, int site, SiteScore.SITE_SCORE_TYPE scoreType) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return 0;
        }
        return myBase.valueByte((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF), scoreType);
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
    public int value(int taxon, int site, SITE_SCORE_TYPE scoreType) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return 0;
        }
        return myBase.value((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF), scoreType);
    }

}

