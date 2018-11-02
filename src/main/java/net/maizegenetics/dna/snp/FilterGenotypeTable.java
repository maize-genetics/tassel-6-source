/*
 * FilterGenotypeTable
 */
package net.maizegenetics.dna.snp;

import java.util.*;

import net.maizegenetics.analysis.distance.IBSDistanceMatrixOneByAll;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.AlleleDepthBuilder;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.AlleleProbabilityBuilder;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.DosageBuilder;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbabilityBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

/**
 * Taxa and site filtering of GenotypeTables. The class essentially creates
 * views of the baseGenotypeTable through arrays for indirection.
 *
 * @author Terry Casstevens
 */
public class FilterGenotypeTable {

    private static final Logger myLogger = Logger.getLogger(FilterGenotypeTable.class);

    private FilterGenotypeTable() {
    }

    public static GenotypeTable getInstance(GenotypeTable base, TaxaList taxa, PositionList positions, Translate translate) {

        if (translate == null) {
            throw new IllegalArgumentException("FilterGenotypeTable: getInstance: must specify translation.");
        }

        if (taxa == null) {
            taxa = createTaxa(base, translate);
        }

        if (positions == null) {
            positions = createPositions(base, translate);
        }

        int numTaxa = taxa.numberOfTaxa();
        int numSites = positions.numberOfSites();

        AlleleDepth depth = base.depth();
        if (depth != null) {
            depth = AlleleDepthBuilder.getFilteredInstance(depth, translate);
            if (depth.numTaxa() != numTaxa || depth.numSites() != numSites) {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: depth dimensions doesn't match");
            }
        }

        AlleleProbability alleleProbability = base.alleleProbability();
        if (alleleProbability != null) {
            alleleProbability = AlleleProbabilityBuilder.getFilteredInstance(alleleProbability, translate);
            if (alleleProbability.numTaxa() != numTaxa || alleleProbability.numSites() != numSites) {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: allele probability dimensions doesn't match");
            }
        }

        Dosage dosage = base.dosage();
        if (dosage != null) {
            dosage = DosageBuilder.getFilteredInstance(dosage, translate);
            if (dosage.numTaxa() != numTaxa || dosage.numSites() != numSites) {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: dosage dimensions doesn't match");
            }
        }

        ReferenceProbability referenceProbability = base.referenceProbability();
        if (referenceProbability != null) {
            referenceProbability = ReferenceProbabilityBuilder.getFilteredInstance(referenceProbability, translate);
            if (referenceProbability.numTaxa() != numTaxa || referenceProbability.numSites() != numSites) {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: reference probability dimensions doesn't match");
            }
        }

        GenotypeCallTable genotypeCallTable = base.genotypeMatrix();
        Translate genotypeTranslate = null;
        if (genotypeCallTable != null) {
            Tuple<GenotypeCallTable, Translate> temp = GenotypeCallTableBuilder.getFilteredInstance(genotypeCallTable, translate);
            genotypeCallTable = temp.x;
            genotypeTranslate = temp.y;
            if (genotypeCallTable.numberOfTaxa() != numTaxa || genotypeCallTable.numberOfSites() != numSites) {
                throw new IllegalStateException("FilterGenotypeTable: getInstance: genotype call table dimensions doesn't match");
            }
        }

        return new CoreGenotypeTable(genotypeCallTable, positions, taxa, depth, alleleProbability, referenceProbability, dosage, base.annotations(), genotypeTranslate);

    }

    public static PositionList createPositions(GenotypeTable genotypes, Translate translate) {

        if (translate != null && translate.hasSiteTranslations()) {
            PositionListBuilder builder = new PositionListBuilder();
            PositionList basePositions = genotypes.positions();
            for (int i = 0; i < translate.numSites(); i++) {
                builder.add(basePositions.get(translate.site(i)));
            }
            return builder.build();
        } else {
            return genotypes.positions();
        }

    }

    public static TaxaList createTaxa(GenotypeTable genotypes, Translate translate) {

        if (translate != null && translate.hasTaxaTranslations()) {
            TaxaListBuilder builder = new TaxaListBuilder();
            TaxaList baseTaxa = genotypes.taxa();
            for (int i = 0; i < translate.numTaxa(); i++) {
                builder.add(baseTaxa.get(translate.taxon(i)));
            }
            return builder.build();
        } else {
            return genotypes.taxa();
        }

    }

    /**
     * For use when converting between coordinate systems.
     *
     * @param genotypes base genotypes
     * @param positions new position list
     * @param redirectSites redirected sites from positions to base genotypes
     *
     * @return new reordered genotype table
     */
    public static GenotypeTable getInstance(GenotypeTable genotypes, PositionList positions, int[] redirectSites) {

        if (redirectSites.length != genotypes.numberOfSites() || redirectSites.length != positions.numberOfSites()) {
            throw new IllegalArgumentException("FilterGenotypeTable: getInstance: number of positions should be equal.");
        }

        TranslateIndex translateSite = TranslateIndexBuilder.unorderedTranslation(redirectSites, null);
        Translate translate = TranslateBuilder.getInstance(TranslateIndexBuilder.noTranslation(genotypes.numberOfTaxa()), translateSite);
        return getInstance(genotypes, genotypes.taxa(), positions, translate);

    }

    /**
     * This returns GenotypeTable with only specified subTaxaList. Defaults to
     * retain unknown taxa.
     *
     * @param a original genotype table
     * @param subTaxaList taxa list subset
     *
     * @return genotype table
     */
    public static GenotypeTable getInstance(GenotypeTable a, TaxaList subTaxaList) {
        return getInstance(a, subTaxaList, true);
    }

    public static GenotypeTable getInstanceSortTaxaAlphabetically(GenotypeTable genotypes) {

        TaxaListBuilder builder = new TaxaListBuilder();
        builder.addAll(genotypes.taxa());
        int[] redirect = builder.sortAlphabetically();
        for (int t = 0; t < redirect.length; t++) {
            if (redirect[t] != t) {
                TranslateIndex translateTaxa = TranslateIndexBuilder.unorderedTranslation(redirect, null);
                Translate translate = TranslateBuilder.getInstance(translateTaxa, TranslateIndexBuilder.noTranslation(genotypes.numberOfSites()));
                return getInstance(genotypes, null, null, translate);
            }
        }
        return genotypes;

    }

    public static Tuple<GenotypeTable, double[]> getInstanceTaxaOrderedByGeneticDistance(GenotypeTable genotypes, int taxon) {
        int numTaxa = genotypes.numberOfTaxa();
        double[] distances = IBSDistanceMatrixOneByAll.getInstance(genotypes, taxon);
        GeneticDistance[] distancesForTaxon = new GeneticDistance[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            if (t == taxon) {
                distancesForTaxon[t] = new GeneticDistance(t, -1.0);
            } else if (Double.isNaN(distances[t])) {
                distancesForTaxon[t] = new GeneticDistance(t, 0.0);
            } else {
                distancesForTaxon[t] = new GeneticDistance(t, distances[t]);
            }
        }
        Arrays.sort(distancesForTaxon, (GeneticDistance x, GeneticDistance y) -> (x.myDistance < y.myDistance) ? -1 : ((x.myDistance == y.myDistance) ? 0 : 1));
        int[] taxaRedirect = new int[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            taxaRedirect[t] = distancesForTaxon[t].myIndex;
            distances[t] = distancesForTaxon[t].myDistance;
        }
        TranslateIndex translateTaxa = TranslateIndexBuilder.unorderedTranslation(taxaRedirect, null);
        Translate translate = TranslateBuilder.getInstance(translateTaxa, TranslateIndexBuilder.noTranslation(genotypes.numberOfSites()));
        return new Tuple<>(getInstance(genotypes, null, null, translate), distances);
    }

    private static class GeneticDistance {

        public final int myIndex;
        public final double myDistance;

        public GeneticDistance(int index, double distance) {
            myIndex = index;
            myDistance = distance;
        }

    }

    /**
     * This returns FilterGenotypeTable with only specified subTaxaList. If
     * retainUnknownTaxa is true then Alignment will return unknown values for
     * missing taxa.
     *
     * @param genotype original genotype table
     * @param subTaxaList taxa list subset
     * @param retainUnknownTaxa whether to retain unknown taxa
     *
     * @return genotype table
     */
    public static GenotypeTable getInstance(GenotypeTable genotype, TaxaList subTaxaList, boolean retainUnknownTaxa) {

        int numSites = genotype.numberOfSites();
        int numBaseTaxa = genotype.numberOfTaxa();

        TranslateIndexBuilder builder = TranslateIndexBuilder.getInstance(numBaseTaxa, null).unordered();

        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        boolean noNeedToFilter = true;
        if (subTaxaList.numberOfTaxa() != genotype.numberOfTaxa()) {
            noNeedToFilter = false;
        }

        for (int i = 0, n = subTaxaList.numberOfTaxa(); i < n; i++) {
            int ion = genotype.taxa().indexOf(subTaxaList.get(i));

            if (ion != i) {
                noNeedToFilter = false;
            }

            if (ion == -1) {
                if (retainUnknownTaxa) {
                    builder.keepIndex(-1);
                    taxaBuilder.add(subTaxaList.get(i));
                }
            } else {
                builder.keepIndex(ion);
                taxaBuilder.add(genotype.taxa().get(ion));
            }
        }

        if (noNeedToFilter) {
            return genotype;
        }

        Translate translate = TranslateBuilder.getInstance(builder.build(), TranslateIndexBuilder.noTranslation(numSites));

        return getInstance(genotype, taxaBuilder.build(), null, translate);

    }

    /**
     * Removes specified IDs.
     *
     * @param a alignment to getInstance
     * @param subTaxaList specified IDs
     *
     * @return Filtered Alignment
     */
    public static GenotypeTable getInstanceRemoveIDs(GenotypeTable a, TaxaList subTaxaList) {

        TaxaListBuilder result = new TaxaListBuilder();
        TaxaList current = a.taxa();
        for (int i = 0, n = current.numberOfTaxa(); i < n; i++) {
            if (subTaxaList.indexOf(current.get(i)) < 0) {
                result.add(current.get(i));
            }
        }
        return FilterGenotypeTable.getInstance(a, result.build());

    }

    public static GenotypeTable getInstance(GenotypeTable genotypes, int[] subSites) {

        int numBaseSites = genotypes.numberOfSites();
        int numTaxa = genotypes.numberOfTaxa();

        if (subSites.length > numBaseSites) {
            throw new IllegalArgumentException("FilterGenotypeTable: getInstance: subset of sites: " + subSites.length + " can't be more than original sites: " + genotypes.numberOfSites());
        } else if (subSites.length == numBaseSites) {
            return genotypes;
        }

        TranslateIndexBuilder builder = TranslateIndexBuilder.getInstance(numBaseSites);
        builder.keepIndices(subSites);
        Translate translate = TranslateBuilder.getInstance(TranslateIndexBuilder.noTranslation(numTaxa), builder.build());
        return getInstance(genotypes, genotypes.taxa(), null, translate);

    }

    public static GenotypeTable getInstance(GenotypeTable genotypes, BitSet subSites, boolean includeSites) {

        int numSites = genotypes.numberOfSites();
        int[] newSubSites = null;
        if (includeSites) {
            int numSitesToInclude = (int) subSites.cardinality();
            newSubSites = new int[numSitesToInclude];
            int count = 0;
            for (int s = 0; s < numSites; s++) {
                if (subSites.fastGet(s)) {
                    newSubSites[count++] = s;
                }
            }
        } else {
            int numSitesToInclude = numSites - (int) subSites.cardinality();
            newSubSites = new int[numSitesToInclude];
            int count = 0;
            for (int s = 0; s < numSites; s++) {
                if (!subSites.fastGet(s)) {
                    newSubSites[count++] = s;
                }
            }
        }

        return FilterGenotypeTable.getInstance(genotypes, newSubSites);

    }

    public static GenotypeTable getInstance(GenotypeTable a, List<String> siteNamesToKeep) {
        return getInstance(a, siteNamesToKeep.toArray(new String[siteNamesToKeep.size()]));
    }

    public static GenotypeTable getInstance(GenotypeTable a, String[] siteNamesToKeep) {

        Arrays.sort(siteNamesToKeep);
        int[] temp = new int[siteNamesToKeep.length];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToKeep, a.siteName(i)) >= 0) {
                temp[count++] = i;
                if (count == siteNamesToKeep.length) {
                    break;
                }
            }
        }

        int[] result = null;
        if (count == siteNamesToKeep.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static GenotypeTable getInstanceRemoveSiteNames(GenotypeTable a, List<String> siteNamesToRemove) {
        return getInstanceRemoveSiteNames(a, siteNamesToRemove.toArray(new String[siteNamesToRemove.size()]));
    }

    public static GenotypeTable getInstanceRemoveSiteNames(GenotypeTable a, String[] siteNamesToRemove) {

        Arrays.sort(siteNamesToRemove);
        int[] temp = new int[a.numberOfSites()];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToRemove, a.siteName(i)) < 0) {
                temp[count++] = i;
            }
        }

        int[] result = null;
        if (count == temp.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static GenotypeTable getInstance(GenotypeTable a, PositionList subPositionList) {

        int[] temp = new int[subPositionList.size()];
        int count = 0;
        PositionList positionList = a.positions();
        for (Position position : subPositionList) {
            int index = positionList.indexOf(position);
            if (index >= 0) {
                temp[count++] = index;
            }
        }

        int[] result = null;
        if (count == subPositionList.size()) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static GenotypeTable getInstance(GenotypeTable a, String chromosome, int startPhysicalPos, int endPhysicalPos) {
        return getInstance(a, a.chromosome(chromosome), startPhysicalPos, endPhysicalPos);
    }

    public static GenotypeTable getInstance(GenotypeTable a, Chromosome chromosome, int startPhysicalPos, int endPhysicalPos) {

        int startSite = a.siteOfPhysicalPosition(startPhysicalPos, chromosome);
        if (startSite < 0) {
            startSite = -(startSite + 1);
        }

        int endSite = a.siteOfPhysicalPosition(endPhysicalPos, chromosome);
        if (endSite < 0) {
            endSite = -(endSite + 2);
        }

        if (startSite > endSite) {
            myLogger.warn("getInstance: start site: " + startSite + " from physical pos: " + startPhysicalPos + " is larger than end site: " + endSite + " from physical pos: " + endPhysicalPos);
            return null;
        }

        return getInstance(a, startSite, endSite);

    }

    public static GenotypeTable getInstance(GenotypeTable a, Chromosome chromosome) {
        int[] endStart = a.firstLastSiteOfChromosome(chromosome);
        return getInstance(a, endStart[0], endStart[1]);
    }

    /**
     * Factory method that returns genotypes FilterGenotypeTable viewing sites
     * between start site (inclusive) and end site (inclusive).
     *
     * @param genotypes alignment
     * @param startSite start site
     * @param endSite end site
     *
     * @return Genotype Table
     */
    public static GenotypeTable getInstance(GenotypeTable genotypes, int startSite, int endSite) {

        int numBaseSites = genotypes.numberOfSites();
        int numTaxa = genotypes.numberOfTaxa();

        if ((startSite == 0) && (endSite == genotypes.numberOfSites() - 1)) {
            return genotypes;
        }
        if ((startSite < 0) || (startSite > endSite)) {
            throw new IllegalArgumentException("FilterGenotypeTable: getInstance: startSite: " + startSite + " less than zero or greater than end site.");
        }
        if (endSite >= genotypes.numberOfSites()) {
            throw new IllegalArgumentException("FilterGenotypeTable: getInstance: end site: " + endSite + " greater than or equal to number of sites: " + genotypes.numberOfSites());
        }

        Translate translate = TranslateBuilder.getInstance(TranslateIndexBuilder.noTranslation(numTaxa), TranslateIndexBuilder.range(startSite, endSite, null, numBaseSites));
        return getInstance(genotypes, genotypes.taxa(), null, translate);

    }

    public static GenotypeTable getInstance(GenotypeTable a, int startSite, int endSite, boolean includeSites) {

        if (includeSites) {
            return getInstance(a, startSite, endSite);
        } else {
            int numSites = a.numberOfSites();
            BitSet rangeBits = new OpenBitSet(numSites);
            for (int i = startSite; i <= endSite; i++) {
                rangeBits.fastSet(i);
            }
            return getInstance(a, rangeBits, includeSites);
        }

    }

}
