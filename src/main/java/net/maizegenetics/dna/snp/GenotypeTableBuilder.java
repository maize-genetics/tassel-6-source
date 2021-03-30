/*
 *  GenotypeTableBuilder
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeMergeRule;
import net.maizegenetics.dna.snp.genotypecall.MaskGenotypeCallTable;
import net.maizegenetics.dna.snp.score.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotationStorage;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Builder for GenotypeTables. New genotypeTables are built from a minimum of
 * TaxaList, PositionList, and GenotypeCallTable. Depth and Scores are optional
 * features of GenotypeTables.
 * <p>
 * </p>
 * If you know the taxa,position, and genotypes are known from the beginning
 * use: GenotypeTable a=GenotypeTableBuilder.getInstance(genotype, positionList,
 * taxaList);
 *
 * In many situations only GenotypeTables are built incrementally, either by
 * Taxa or Site.
 * <p>
 * </p>
 * For taxa building:
 * <pre>
 * {@code
 *    GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(gbs.positions(),outFile);
 *    for (int i=0; i<hm2.numberOfTaxa(); i++) {
 *        Taxon taxon=hm2.taxa().get(i);
 *        byte[] geno=hm2.genotypeAllSites(i);
 *        gtb.addTaxon(taxon,geno);
 *        }
 *    GenotypeTable gt=gtb.build();
 * }
 * </pre>
 * <p>
 * </p>
 * In many cases, genotype want to add taxa to an existing genotypeTable. Direct
 * addition is not possible, as GenotypeTables are immutable, but the
 * GenotypeTableBuilder.getTaxaIncremental provides a strategy for creating and
 * merging taxa together. Key to the process is that GenotypeMergeRule defines
 * how the taxa with identical names will be merged.<br></br>
 * <pre>{@code
 * GenotypeTable existingGenotypeTable1, existingGenotypeTable2;
 * GenotypeTableBuilder gtb=GenotypeTableBuilder.getTaxaIncremental(existingGenotypeTable1,
 * new BasicGenotypeMergeRule(0.01));
 * for (int i=0; i<existingGenotypeTable2.numberOfTaxa(); i++) {
 * gtb.addTaxon(existingGenotypeTable2.taxa().get(i), existingGenotypeTable2.genotypeAllSites(i)
 * existingGenotypeTable2.depth().depthAllSitesByte(i));
 * }
 * }</pre>
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class GenotypeTableBuilder {

    //Fields for incremental taxa
    private PositionList positionList = null;
    private TaxaListBuilder taxaListBuilder = null;
    private ArrayList<byte[]> incGeno = null;
    private ArrayList<byte[][]> incDepth = null;
    private AlleleProbabilityBuilder myAlleleProbabilityBuilder = null;
    private ReferenceProbabilityBuilder myReferenceProbabilityBuilder = null;
    private DosageBuilder myDosageBuilder = null;
    private HashMap<Taxon, Integer> incTaxonIndex = null;
    private boolean sortAlphabetically = false;

    //Fields for incremental sites
    private final TaxaList taxaList;
    private PositionListBuilder posListBuilder = null;
    private boolean isTaxaMerge = false; //if in taxa merge mode, this only works with TAXA_INC build type;//, GENO_EDIT}; //GENO_EDIT is not
    private GenotypeMergeRule mergeRule = null;
    private BuildType myBuildType;
    private final GeneralAnnotationStorage.Builder myAnnotationBuilder = GeneralAnnotationStorage.getBuilder();

    /**
     * Builder for in memory taxa incremental
     */
    private GenotypeTableBuilder(PositionList positionList, GenotypeMergeRule mergeRule) {
        this.positionList = positionList;
        this.myBuildType = BuildType.TAXA_INC;
        this.mergeRule = mergeRule;
        if (mergeRule != null) {
            this.isTaxaMerge = true;

        }
        incGeno = new ArrayList<>();
        incDepth = new ArrayList<>();
        incTaxonIndex = new HashMap<>();
        taxaListBuilder = new TaxaListBuilder();
        this.taxaList = null;
    }

    /**
     * Builder for in memory site incremental
     */
    private GenotypeTableBuilder(TaxaList taxaList) {
        this.taxaList = taxaList;
        this.myBuildType = BuildType.SITE_INC;
        incGeno = new ArrayList<>();
        posListBuilder = new PositionListBuilder();
    }

    /**
     * Creates an in memory builder for addition by taxon. Each taxon can only
     * be added once, i.e. merging is not possible
     *
     * @param positionList The positions used for the builder
     *
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncremental(PositionList positionList) {
        return new GenotypeTableBuilder(positionList, (GenotypeMergeRule) null);
    }

    /**
     * Creates an in memory builder for addition by taxon, which permits the
     * merging of taxa.
     *
     * @param positionList The positions used for the builder
     * @param mergeRule rules for merging identically named taxa
     *
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncremental(PositionList positionList, GenotypeMergeRule mergeRule) {
        return new GenotypeTableBuilder(positionList, mergeRule);
    }

    /**
     * Creates a builder initialized with the Genotypes in a existing
     * GenotypeTable. The position list and initial taxa list is derived from
     * the positions, taxa, and genotypes already in the GenotypeTable. The
     * initial GenotypeTable is not changed as it is immutable.
     *
     * @param genotypeTable input genotype table
     * @param mergeRule rules for merging identically named taxa
     *
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncremental(GenotypeTable genotypeTable, GenotypeMergeRule mergeRule) {
        PositionList positionList = genotypeTable.positions();
        GenotypeTableBuilder gtb = new GenotypeTableBuilder(positionList, mergeRule);
        boolean hasDepth = genotypeTable.hasDepth();
        for (int i = 0; i < genotypeTable.numberOfTaxa(); i++) {
            if (hasDepth) {
                gtb.addTaxon(genotypeTable.taxa().get(i), genotypeTable.genotypeAllSites(i), genotypeTable.depth().valuesForTaxonByte(i));
            } else {
                gtb.addTaxon(genotypeTable.taxa().get(i), genotypeTable.genotypeAllSites(i));
            }
        }
        return gtb;
    }

    /**
     * Build an alignment site by site in memory
     *
     * @param taxaList
     *
     * @return builder to add sites to
     */
    public static GenotypeTableBuilder getSiteIncremental(TaxaList taxaList) {
        return new GenotypeTableBuilder(taxaList);
    }

    public static GenotypeTable getInstance(GenotypeTable original, GenotypeCallTable newGenotypes) {
        return getInstance(newGenotypes, original.positions(), original.taxa(), original.depth(), original.alleleProbability(), original.referenceProbability(), original.dosage(), original.annotations());
    }

    /**
     * Standard approach for creating a new Alignment
     *
     * @param genotype
     * @param positionList
     * @param taxaList
     * @param alleleDepth
     * @param alleleProbability
     * @param referenceProbability
     * @param dosage
     * @param annotations
     *
     * @return new genotype table
     */
    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, AlleleDepth alleleDepth, AlleleProbability alleleProbability, ReferenceProbability referenceProbability, Dosage dosage, GeneralAnnotationStorage annotations) {

        if (positionList == null) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: position list is required.");
        }
        int numSites = positionList.numberOfSites();

        if ((genotype != null) && (numSites != genotype.numberOfSites())) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + genotype.numberOfSites() + " doesn't equal number of sites in position list: " + numSites);
        }

        if (taxaList == null) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: taxa list is required.");
        }
        int numTaxa = taxaList.numberOfTaxa();

        if ((genotype != null) && (numTaxa != genotype.numberOfTaxa())) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + genotype.numberOfTaxa() + " doesn't equal number of taxa in taxa list: " + numTaxa);
        }

        if (alleleProbability != null) {
            if (numSites != alleleProbability.numSites()) {
                throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + numSites + " doesn't equal number of sites in allele probability: " + alleleProbability.numSites());
            }
            if (numTaxa != alleleProbability.numTaxa()) {
                throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + numTaxa + " doesn't equal number of taxa in allele probability: " + alleleProbability.numTaxa());
            }
        }

        if (referenceProbability != null) {
            if (numSites != referenceProbability.numSites()) {
                throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + numSites + " doesn't equal number of sites in reference probability: " + referenceProbability.numSites());
            }
            if (numTaxa != referenceProbability.numTaxa()) {
                throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + numTaxa + " doesn't equal number of taxa in reference probability: " + referenceProbability.numTaxa());
            }
        }

        if (dosage != null) {
            if (numSites != dosage.numSites()) {
                throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + numSites + " doesn't equal number of sites in dosage: " + dosage.numSites());
            }
            if (numTaxa != dosage.numTaxa()) {
                throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + numTaxa + " doesn't equal number of taxa in dosage: " + dosage.numTaxa());
            }
        }

        return new CoreGenotypeTable(genotype, positionList, taxaList, alleleDepth, alleleProbability, referenceProbability, dosage, annotations);
    }

    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, AlleleDepth alleleDepth) {
        return getInstance(genotype, positionList, taxaList, alleleDepth, null, null, null, null);
    }

    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList) {
        return getInstance(genotype, positionList, taxaList, null, null, null, null, null);
    }

    public static GenotypeTable getInstance(GenotypeTable base, MaskMatrix mask) {

        AlleleDepth depth = base.depth();
        if (depth != null) {
            depth = AlleleDepthBuilder.getMaskInstance(depth, mask);
        }

        AlleleProbability alleleProbability = base.alleleProbability();
        if (alleleProbability != null) {
            alleleProbability = AlleleProbabilityBuilder.getMaskInstance(alleleProbability, mask);
        }

        Dosage dosage = base.dosage();
        if (dosage != null) {
            dosage = DosageBuilder.getMaskInstance(dosage, mask);
        }

        ReferenceProbability referenceProbability = base.referenceProbability();
        if (referenceProbability != null) {
            referenceProbability = ReferenceProbabilityBuilder.getMaskInstance(referenceProbability, mask);
        }

        return new CoreGenotypeTable(new MaskGenotypeCallTable(base.genotypeMatrix(), mask), base.positions(), base.taxa(), depth, alleleProbability, referenceProbability, dosage, base.annotations());

    }

    public static GenotypeTable getInstanceOnlyMajorMinor(GenotypeTable genotype) {
        return getInstance(genotype, MaskMatrixBuilder.getInstanceRemoveMinorSNPs(genotype.genotypeMatrix()));
    }

    public static GenotypeTable getHomozygousInstance(GenotypeTable genotype) {
        return getInstance(genotype, MaskMatrixBuilder.getInstanceRemoveHeterozygous(genotype.genotypeMatrix()));
    }

    public static GenotypeTable getInstanceMaskIndels(GenotypeTable genotype) {
        return getInstance(genotype, MaskMatrixBuilder.getInstanceRemoveIndels(genotype.genotypeMatrix()));
    }

    /**
     * Creates a GenotypeTable with in-memory instance of GenotypeCallTable.
     * Primarily needed for performance critical situations like imputation.
     *
     * @param genotypeTable genotype table
     *
     * @return alignment backed by a single in-memory SuperByteMatrix
     */
    public static GenotypeTable getGenotypeCopyInstance(GenotypeTable genotypeTable) {
        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getInstanceCopy(genotypeTable.genotypeMatrix());
        return new CoreGenotypeTable(builder.build(), genotypeTable.positions(), genotypeTable.taxa());
    }

    public GenotypeTableBuilder addAlleleProbability(AlleleProbabilityBuilder alleleProbabilityBuilder) {
        myAlleleProbabilityBuilder = alleleProbabilityBuilder;
        return this;
    }

    public GenotypeTableBuilder addReferenceProbability(ReferenceProbabilityBuilder referenceProbabilityBuilder) {
        myReferenceProbabilityBuilder = referenceProbabilityBuilder;
        return this;
    }

    public GenotypeTableBuilder addDosage(DosageBuilder dosageBuilder) {
        myDosageBuilder = dosageBuilder;
        return this;
    }

    public GenotypeTableBuilder addSite(Position pos, byte[] genos) {
        if (myBuildType != BuildType.SITE_INC) {
            throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getSiteIncremental");
        }
        if (genos.length != taxaList.numberOfTaxa()) {
            throw new IndexOutOfBoundsException("Number of taxa and genotypes do not agree");
        }
        synchronized (taxaList) {
            posListBuilder.add(pos);
            incGeno.add(genos);
        }
        return this;
    }

    public GenotypeTableBuilder addTaxon(Taxon taxon, byte[] genos) {
        return addTaxon(taxon, genos, null);
    }

    public GenotypeTableBuilder addTaxon(Taxon taxon, byte[] genos, byte[][] depth) {
        if (myBuildType != BuildType.TAXA_INC) {
            throw new IllegalArgumentException("addTaxon only be used with AlignmentBuilder.getTaxaIncremental");
        }
        if (genos.length != positionList.numberOfSites()) {
            throw new IndexOutOfBoundsException("Number of sites and genotypes do not agree");
        }
        synchronized (taxaListBuilder) {
            if (isTaxaMerge && incTaxonIndex.containsKey(taxon)) {
                mergeTaxonInMemory(taxon, genos, depth);
            } else {
                taxaListBuilder.add(taxon);
                incGeno.add(genos);
                incDepth.add(depth);
                incTaxonIndex.put(taxon, incGeno.size() - 1);
            }
        }

        return this;
    }

    public GenotypeTableBuilder addTaxon(Taxon taxon, int[][] depths, byte[] genos) {  // reversed signature so it does clash with above statement: "return addTaxon(taxon, genos, null);"
        byte[][] byteDepths = AlleleDepthUtil.depthIntToByte(depths);
        return addTaxon(taxon, genos, byteDepths);
    }

    private void mergeTaxonInMemory(Taxon taxon, byte[] genos, byte[][] depth) {
        int taxonIndex = incTaxonIndex.get(taxon);
        byte[] combGenos = new byte[genos.length];
        if (depth != null) {
            byte[][] existingDepth = incDepth.get(taxonIndex);
            byte[][] combDepth = new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][genos.length];
            byte[] currDepths = new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (int site = 0; site < combDepth[0].length; site++) {
                for (int allele = 0; allele < combDepth.length; allele++) {
                    currDepths[allele] = combDepth[allele][site] = AlleleDepthUtil.addByteDepths(depth[allele][site], existingDepth[allele][site]);
                }
                combGenos[site] = mergeRule.callBasedOnDepth(currDepths);
            }
            incGeno.set(taxonIndex, combGenos);
            incDepth.set(taxonIndex, combDepth);
        } else {
            byte[] existingGenos = incGeno.get(taxonIndex);
            for (int site = 0; site < combGenos.length; site++) {
                combGenos[site] = mergeRule.mergeCalls(genos[site], existingGenos[site]);
            }
            incGeno.set(taxonIndex, combGenos);
        }
    }

    /**
     * Set the builder so that when built it will sort the taxa
     */
    public GenotypeTableBuilder sortTaxa() {
        if (myBuildType != BuildType.TAXA_INC) {
            throw new IllegalArgumentException("sortTaxa can only be used with AlignmentBuilder.getTaxaIncremental");
        }
        sortAlphabetically = true;
        return this;
    }

    /**
     * Finishes building the GenotypeTable.
     *
     * @return a genotype table
     */
    public GenotypeTable build() {

        switch (myBuildType) {
            case TAXA_INC: {
                TaxaList tl = (sortAlphabetically) ? taxaListBuilder.sortTaxaAlphabetically().build() : taxaListBuilder.build();
                GenotypeCallTableBuilder gB = GenotypeCallTableBuilder.getInstance(tl.numberOfTaxa(), positionList.numberOfSites());
                boolean hasDepth = (incDepth.size() == tl.numberOfTaxa() && incDepth.get(0) != null);
                AlleleDepthBuilder adb = null;
                if (hasDepth) {
                    adb = AlleleDepthBuilder.getInstance(tl.numberOfTaxa(), positionList.numberOfSites(), taxaList);
                }
                for (int i = 0; i < incGeno.size(); i++) {
                    gB.setBaseRangeForTaxon(i, 0, incGeno.get(incTaxonIndex.get(tl.get(i))));
                    if (hasDepth) {
                        adb.addTaxon(i, incDepth.get(incTaxonIndex.get(tl.get(i))));
                    }
                }
                AlleleDepth ad = (hasDepth) ? adb.build() : null;

                AlleleProbability alleleProbability = null;
                if (myAlleleProbabilityBuilder != null) {
                    alleleProbability = myAlleleProbabilityBuilder.build();
                }

                ReferenceProbability referenceProbability = null;
                if (myReferenceProbabilityBuilder != null) {
                    referenceProbability = myReferenceProbabilityBuilder.build();
                }

                Dosage dosage = null;
                if (myDosageBuilder != null) {
                    dosage = myDosageBuilder.build();
                }

                return getInstance(gB.build(), positionList, tl, ad, alleleProbability, referenceProbability, dosage, myAnnotationBuilder.build());
            }
            case SITE_INC: {
                GenotypeCallTableBuilder gB = GenotypeCallTableBuilder.getInstance(taxaList.numberOfTaxa(), posListBuilder.size());
                for (int s = 0; s < posListBuilder.size(); s++) {
                    byte[] b = incGeno.get(s);
                    for (int t = 0; t < b.length; t++) {
                        gB.setBase(t, s, b[t]);
                    }
                }
                PositionList pl = posListBuilder.build(gB);
                return getInstance(gB.build(), pl, taxaList);
            }
        }
        return null;
    }

    private static enum BuildType {

        TAXA_INC, SITE_INC
    }

    public GenotypeTableBuilder addAnnotation(String key, String value) {
        myAnnotationBuilder.addAnnotation(key, value);
        return this;
    }

    public GenotypeTableBuilder addAnnotation(String key, Number value) {
        myAnnotationBuilder.addAnnotation(key, value);
        return this;
    }

    public GenotypeTableBuilder dataSetName(String dataSetName) {
        myAnnotationBuilder.addAnnotation(GenotypeTable.ANNOTATION_DATA_SET_NAME, dataSetName);
        return this;
    }

    public GenotypeTableBuilder dataSetDescription(String dataSetDescription) {
        myAnnotationBuilder.addAnnotation(GenotypeTable.ANNOTATION_DATA_SET_DESCRIPTION, dataSetDescription);
        return this;
    }
}
