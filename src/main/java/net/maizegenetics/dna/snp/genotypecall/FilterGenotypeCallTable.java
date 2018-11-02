/*
 *  FilterGenotypeCallTable
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.util.Arrays;
import java.util.Spliterator;
import static java.util.Spliterator.IMMUTABLE;
import static java.util.Spliterator.ORDERED;
import static java.util.Spliterator.SIZED;
import static java.util.Spliterator.SUBSIZED;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.Translate;

/**
 * Filtering class for GenotypeCallTable. This class is generally never used
 * directly, but rather used through FilterGenotypeTable.
 *
 * @see net.maizegenetics.dna.snp.FilterGenotypeTable
 * @author Terry Casstevens
 */
class FilterGenotypeCallTable extends AbstractGenotypeCallTable {

    final GenotypeCallTable myBaseGenotype;
    final Translate myTranslate;

    FilterGenotypeCallTable(GenotypeCallTable genotype, Translate translate) {
        super(translate.numTaxa(), translate.numSites(), genotype.isPhased(), null, genotype.maxNumAlleles());
        myBaseGenotype = genotype;
        myTranslate = translate;
    }

    @Override
    public byte genotype(int taxon, int site) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        } else {
            return myBaseGenotype.genotype((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF));
        }
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        long taxonSite = myTranslate.taxonSite(taxon, site);
        if (taxonSite == -1) {
            return GenotypeTable.UNKNOWN_ALLELE_STR;
        } else {
            return myBaseGenotype.genotypeAsString((int) (taxonSite >>> 32), (int) (taxonSite & 0xFFFFFFFF));
        }
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myBaseGenotype.diploidAsString(site, value);
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        if (myTranslate.site(site) == -1) {
            byte[] result = new byte[numberOfTaxa()];
            Arrays.fill(result, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
            return result;
        } else if (!myTranslate.hasTaxaTranslations()) {
            return myBaseGenotype.genotypeForAllTaxa(myTranslate.site(site));
        } else {
            byte[] orig = myBaseGenotype.genotypeForAllTaxa(myTranslate.site(site));
            int[] translations = myTranslate.taxaTranslations();
            int numTaxa = myTranslate.numTaxa();
            byte[] result = new byte[numTaxa];
            for (int i = 0; i < numTaxa; i++) {
                if (translations[i] == -1) {
                    result[i] = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
                } else {
                    result[i] = orig[translations[i]];
                }
            }
            return result;
        }
    }

    @Override
    public byte[] genotypeForAllSites(int taxon) {
        if (myTranslate.taxon(taxon) == -1) {
            byte[] result = new byte[numberOfSites()];
            Arrays.fill(result, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
            return result;
        } else if (!myTranslate.hasSiteTranslations()) {
            return myBaseGenotype.genotypeForAllSites(myTranslate.taxon(taxon));
        } else {
            return super.genotypeForAllSites(taxon);
        }
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        if (myTranslate.site(site) == -1) {
            return new int[0][0];
        } else if (!myTranslate.hasTaxaTranslations()) {
            return myBaseGenotype.allelesSortedByFrequency(myTranslate.site(site));
        } else {
            return super.allelesSortedByFrequency(site);
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isSiteOptimized() {
        return myBaseGenotype.isSiteOptimized();
    }

    @Override
    public String[][] alleleDefinitions() {
        String[][] encodings = myBaseGenotype.alleleDefinitions();
        if (encodings.length == 1) {
            return encodings;
        } else if (myTranslate.hasSiteTranslations()) {
            int numSites = numberOfSites();
            String[][] result = new String[numSites][];
            for (int i = 0; i < numSites; i++) {
                result[i] = alleleDefinitions(i);
            }
            return result;
        } else {
            return encodings;
        }
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return myBaseGenotype.alleleDefinitions(myTranslate.site(site));
    }

    @Override
    public int maxNumAlleles() {
        return myBaseGenotype.maxNumAlleles();
    }

    @Override
    public Stream<Byte> stream() {
        return StreamSupport.stream(spliterator(), true);
    }

    @Override
    public Stream<Byte> stream(int taxon) {
        return StreamSupport.stream(new FilterGenotypeCallTableSpliterator<>(taxon, 0, numberOfSites(), taxon, numberOfSites()), true);
    }

    @Override
    public Spliterator<Byte> spliterator() {
        return new FilterGenotypeCallTableSpliterator<>(0, 0, numberOfSites(), numberOfTaxa() - 1, numberOfSites());
    }

    class FilterGenotypeCallTableSpliterator<T extends Byte> implements Spliterator<Byte> {

        protected int myTaxaOrigin;
        protected int mySiteOrigin;
        protected final int myNumSites;
        protected final int myTaxaFence;
        protected final int mySiteFence;

        FilterGenotypeCallTableSpliterator(int taxaOrigin, int siteOrigin, int numSites, int taxaFence, int siteFence) {
            myTaxaOrigin = taxaOrigin;
            mySiteOrigin = siteOrigin;
            myNumSites = numSites;
            myTaxaFence = taxaFence;
            mySiteFence = siteFence;
        }

        @Override
        public void forEachRemaining(Consumer<? super Byte> action) {
            if (myTranslate.hasSiteTranslations()) {
                for (; myTaxaOrigin < myTaxaFence; myTaxaOrigin++) {
                    int taxaIndex = myTranslate.taxon(myTaxaOrigin);
                    if (taxaIndex == -1) {
                        for (; mySiteOrigin < myNumSites; mySiteOrigin++) {
                            action.accept(GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                        }
                    } else {
                        for (; mySiteOrigin < myNumSites; mySiteOrigin++) {
                            action.accept(myBaseGenotype.genotype(taxaIndex, myTranslate.site(mySiteOrigin)));
                        }
                    }
                    mySiteOrigin = 0;
                }
                int taxaIndex = myTranslate.taxon(myTaxaOrigin);
                if (taxaIndex == -1) {
                    for (; mySiteOrigin < mySiteFence; mySiteOrigin++) {
                        action.accept(GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                    }
                } else {
                    for (; mySiteOrigin < mySiteFence; mySiteOrigin++) {
                        action.accept(myBaseGenotype.genotype(taxaIndex, myTranslate.site(mySiteOrigin)));
                    }
                }
            } else {
                for (; myTaxaOrigin < myTaxaFence; myTaxaOrigin++) {
                    int taxaIndex = myTranslate.taxon(myTaxaOrigin);
                    if (taxaIndex == -1) {
                        for (; mySiteOrigin < myNumSites; mySiteOrigin++) {
                            action.accept(GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                        }
                    } else {
                        for (; mySiteOrigin < myNumSites; mySiteOrigin++) {
                            action.accept(myBaseGenotype.genotype(taxaIndex, mySiteOrigin));
                        }
                    }
                    mySiteOrigin = 0;
                }
                int taxaIndex = myTranslate.taxon(myTaxaOrigin);
                if (taxaIndex == -1) {
                    for (; mySiteOrigin < mySiteFence; mySiteOrigin++) {
                        action.accept(GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                    }
                } else {
                    for (; mySiteOrigin < mySiteFence; mySiteOrigin++) {
                        action.accept(myBaseGenotype.genotype(taxaIndex, mySiteOrigin));
                    }
                }
            }
        }

        @Override
        public boolean tryAdvance(Consumer<? super Byte> action) {
            if (((myTaxaOrigin < myTaxaFence) && (mySiteOrigin < myNumSites))
                    || ((myTaxaOrigin == myTaxaFence) && (mySiteOrigin < mySiteFence))) {
                action.accept(genotype(myTaxaOrigin, mySiteOrigin));
                mySiteOrigin++;
                if (mySiteOrigin >= myNumSites) {
                    mySiteOrigin = 0;
                    myTaxaOrigin++;
                }
                return true;
            } else {
                return false;
            }
        }

        @Override
        public Spliterator<Byte> trySplit() {
            long size = estimateSize();
            if (size > 1) {
                size >>>= 1;
                int loTaxa = myTaxaOrigin;
                int loSite = mySiteOrigin;
                int midTaxa = myTaxaOrigin;
                int midSite = mySiteOrigin;
                midTaxa += size / myNumSites;
                midSite += size % myNumSites;
                if (midSite > myNumSites) {
                    midTaxa++;
                    midSite -= myNumSites;
                }
                myTaxaOrigin = midTaxa;
                mySiteOrigin = midSite;
                return new FilterGenotypeCallTableSpliterator<>(loTaxa, loSite, myNumSites, midTaxa, midSite);
            } else {
                return null;
            }
        }

        @Override
        public long estimateSize() {
            return (long) (myTaxaFence - myTaxaOrigin) * (long) myNumSites - (long) mySiteOrigin + (long) mySiteFence;
        }

        @Override
        public int characteristics() {
            return ORDERED | SIZED | IMMUTABLE | SUBSIZED;
        }
    }
}
