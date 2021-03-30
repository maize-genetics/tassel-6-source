/*
 *  DosageBuilder
 */
package net.maizegenetics.dna.snp.score;

import net.maizegenetics.dna.snp.MaskMatrix;
import net.maizegenetics.dna.snp.Translate;
import net.maizegenetics.dna.snp.TranslateBuilder;
import net.maizegenetics.dna.snp.byte2d.Byte2D;
import net.maizegenetics.dna.snp.byte2d.Byte2DBuilder;
import net.maizegenetics.taxa.TaxaList;

/**
 * @author Terry Casstevens
 */
public class DosageBuilder {

    private Byte2DBuilder myBuilder;
    private final int myNumSites;

    private DosageBuilder(int numTaxa, int numSites, TaxaList taxaList) {
        myBuilder = Byte2DBuilder.getInstance(numTaxa, numSites, SiteScore.SITE_SCORE_TYPE.Dosage, taxaList);
        myNumSites = numSites;
    }

    public static DosageBuilder getInstance(int numTaxa, int numSites, TaxaList taxaList) {
        return new DosageBuilder(numTaxa, numSites, taxaList);
    }

    public static Dosage getFilteredInstance(Dosage base, Translate translate) {
        if (base instanceof FilterDosage) {
            FilterDosage filter = (FilterDosage) base;
            Translate merged = TranslateBuilder.getInstance(filter.myTranslate, translate);
            return new FilterDosage(filter.myBase, merged);
        }
        return new FilterDosage(base, translate);
    }

    public static Dosage getMaskInstance(Dosage base, MaskMatrix mask) {
        return new MaskDosage(base, mask);
    }

    public DosageBuilder addTaxon(int taxon, byte[] values) {
        if (myNumSites != values.length) {
            throw new IllegalArgumentException("DosageBuilder: addTaxon: number of values: " + values.length + " doesn't equal number of sites: " + myNumSites);
        }
        myBuilder.addTaxon(taxon, values);
        return this;
    }

    public Dosage build() {
        Byte2D input = myBuilder.build();
        myBuilder = null;
        return new Dosage(input);
    }

}
