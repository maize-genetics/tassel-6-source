/*
 *  GenotypeTableBuilder
 */
package net.maizegenetics.dna.snp;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
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
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.ArrayList;
import java.util.Arrays;
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
 * Merging is possible with HDF5 files, but only if the closeUnfinished() method
 * was used with the previous building.
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
    private boolean isHDF5 = false;
    private IHDF5Writer writer = null;
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
     * Creates a new HDF5 file if positionList is not null. Opens an existing
     * HDF5 File if positionList is null. Merging is allowed depending on
     * whether a mergeRule is included that can used with TaxaIncremental
     * addition.
     *
     * @param hdf5File
     * @param positionList
     */
    private GenotypeTableBuilder(String hdf5File, PositionList positionList, GenotypeMergeRule mergeRule) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        config.dontUseExtendableDataTypes();
        writer = config.writer();
        if (HDF5Utils.doesGenotypeModuleExist(writer) && HDF5Utils.isHDF5GenotypeLocked(writer)) {
            writer.close();
            throw new UnsupportedOperationException("This file is locked for genotypic additions");
        }
        if (positionList != null) {
            this.positionList = new PositionListBuilder(writer, positionList).build();  //create a new position list
            setupGenotypeTaxaInHDF5(writer);
        } else {
            this.positionList = PositionListBuilder.getInstance(writer);

        }
        this.mergeRule = mergeRule;
        if (mergeRule != null) {
            this.isTaxaMerge = true;
        }
        this.myBuildType = BuildType.TAXA_INC;
        isHDF5 = true;
        this.taxaList = null;
    }

    /**
     * Creates a new HDF5 file if positionList is not null. Opens an existing
     * HDF5 File if positionList is null. Merging is allowed depending on
     * whether a mergeRule is included that can used with TaxaIncremental
     * addition.
     *
     * @param hdf5File
     * @param taxaList
     */
    private GenotypeTableBuilder(String hdf5File, TaxaList taxaList, int numberOfSites) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        config.dontUseExtendableDataTypes();
        writer = config.writer();
        if (HDF5Utils.doesGenotypeModuleExist(writer) && HDF5Utils.isHDF5GenotypeLocked(writer)) {
            writer.close();
            throw new UnsupportedOperationException("This file is locked for genotypic additions");
        }
        this.taxaList = taxaList;
        setupGenotypeTaxaInHDF5(writer);
        posListBuilder = new PositionListBuilder(numberOfSites);
        byte[] missingGenotypes = new byte[numberOfSites];
        Arrays.fill(missingGenotypes, GenotypeTable.UNKNOWN_GENOTYPE);
        for (Taxon taxon : taxaList) {
            HDF5Utils.addTaxon(writer, taxon);
            HDF5Utils.writeHDF5GenotypesCalls(writer, taxon.getName(), missingGenotypes);
        }

        this.myBuildType = BuildType.SITE_INC;
        isHDF5 = true;

    }

    /**
     * Returns a builder to an existing, unfinished HDF5 genotypes file.
     *
     * Can be used if you want to add/modify annotations, etc, and/or call
     * build() to finalize it
     */
    public static GenotypeTableBuilder getBuilder(String existingHDF5File) {
        return new GenotypeTableBuilder(existingHDF5File, null, null);
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
     * Create a new taxa incremental HDF5 GenotypeTableBuilder
     *
     * @param positionList the defined list of positions
     * @param newHDF5File hdf5 file to be created
     *
     * @return the builder to add taxa to
     */
    public static GenotypeTableBuilder getTaxaIncremental(PositionList positionList, String newHDF5File) {
        return new GenotypeTableBuilder(newHDF5File, positionList, null);
    }

    /**
     * Merges taxa to an existing HDF5 file. The position list is derived from
     * the positions already in the existing HDF5 file.
     *
     * @param existingHDF5File
     * @param mergeRule
     *
     * @return builder to merge taxa with
     */
    public static GenotypeTableBuilder mergeTaxaIncremental(String existingHDF5File, GenotypeMergeRule mergeRule) {
        return new GenotypeTableBuilder(existingHDF5File, null, mergeRule);
    }

    /**
     * Creates a new taxa incremental HDF5 GenotypeTableBuilder to which
     * replicate taxa can be added
     *
     * @param newHDF5File
     * @param positionList
     * @param mergeRule
     *
     * @return
     */
    public static GenotypeTableBuilder getTaxaIncrementalWithMerging(String newHDF5File, PositionList positionList, GenotypeMergeRule mergeRule) {
        return new GenotypeTableBuilder(newHDF5File, positionList, mergeRule);
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

    /**
     * Build an GenotypeTable by site block (1<<16 sites). Number of positions
     * (sites) must be known from the beginning. Positions and genotypes must be
     * added by block
     *
     * @param taxaList
     * @param numberOfPositions
     * @param newHDF5File
     *
     * @return builder to add site blocks to
     */
    public static GenotypeTableBuilder getSiteIncremental(TaxaList taxaList, int numberOfPositions, String newHDF5File) {
        return new GenotypeTableBuilder(newHDF5File, taxaList, numberOfPositions);
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

    /**
     * Creates a new HDF5 file alignment based on existing Genotype,
     * PositionList, and TaxaList.
     *
     * @param genotype
     * @param positionList
     * @param taxaList
     * @param hdf5File name of the file
     *
     * @return alignment backed by new HDF5 file
     */
    public static GenotypeTable getInstance(GenotypeCallTable genotype, PositionList positionList, TaxaList taxaList, String hdf5File) {
        if (genotype.numberOfSites() != positionList.numberOfSites()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of sites in genotype: " + genotype.numberOfSites() + " doesn't equal number of sites in position list: " + positionList.numberOfSites());
        }
        if (genotype.numberOfTaxa() != taxaList.numberOfTaxa()) {
            throw new IllegalArgumentException("GenotypeTableBuilder: getInstance: number of taxa in genotype: " + genotype.numberOfTaxa() + " doesn't equal number of taxa in taaxa list: " + taxaList.numberOfTaxa());
        }
        GenotypeTableBuilder aB = GenotypeTableBuilder.getTaxaIncremental(positionList, hdf5File);
        for (int i = 0; i < taxaList.numberOfTaxa(); i++) {
            aB.addTaxon(taxaList.get(i), genotype.genotypeAllSites(i));
        }
        return aB.build();
    }

    /**
     * Creates a new HDF5 file alignment based on an existing alignment.
     *
     * @param a existing alignment
     * @param hdf5File name of the file
     *
     * @return alignment backed by new HDF5 file
     */
    public static GenotypeTable getInstance(GenotypeTable a, String hdf5File) {
        return getInstance(a.genotypeMatrix(), a.positions(), a.taxa(), hdf5File);
    }

    public static GenotypeTable getInstance(String hdf5File) {
        IHDF5Reader reader = HDF5Factory.openForReading(hdf5File);
        TaxaList tL = new TaxaListBuilder().buildFromHDF5Genotypes(reader);
        PositionList pL = PositionListBuilder.getInstance(reader);
        GenotypeCallTable geno = GenotypeCallTableBuilder.buildHDF5(reader);
        AlleleDepth depth = AlleleDepthBuilder.getInstance(reader);
        return GenotypeTableBuilder.getInstance(geno, pL, tL, depth, null, null, null, HDF5Utils.readHDF5Annotation(reader, Tassel5HDF5Constants.ROOT, GenotypeTable.GENOTYPE_TABLE_ANNOTATIONS));
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
        if ((myBuildType != BuildType.SITE_INC) || isHDF5) {
            throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getSiteIncremental and without HDF5");
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

    /**
     * Add TasselHDF5 Block of positions (generally 1<<16 positions). @note
     *
     * This is synchronized, which certainly slows things down but it is needed
     * to prevent the same taxa dataset from being accessed at once. This can
     * probably be rethought with parallelization at this stage across datasets
     *
     * @param startSite start site for positioning blocks correction
     * @param blkPositionList
     * @param blockGenotypes array of genotypes[taxonIndex][siteIndex] true site=startSite+siteIndex
     * @param blockDepths
     */
    public synchronized void addSiteBlock(int startSite, PositionList blkPositionList, byte[][] blockGenotypes, byte[][][] blockDepths) {
        if ((myBuildType != BuildType.SITE_INC) || (isHDF5 == false)) {
            throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getSiteIncremental and with HDF5");
        }
        if (blockGenotypes.length != taxaList.numberOfTaxa()) {
            throw new IndexOutOfBoundsException("Number of taxa and genotypes do not agree");
        }
        int s = startSite;
        System.out.println("startSite = [" + startSite + "], blkPositionList = [" + blkPositionList.size() + "], blockGenotypes = [" + blockGenotypes.length + "], blockDepths = [" + blockDepths + "]");
        for (Position position : blkPositionList) {
            posListBuilder.set(s++, position);
        }
        for (int t = 0; t < taxaList.numberOfTaxa(); t++) {
            HDF5Utils.replaceHDF5GenotypesCalls(writer, taxaList.taxaName(t), startSite, blockGenotypes[t]);
        }
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
        if (isHDF5) {
            if (isTaxaMerge && HDF5Utils.doTaxonCallsExist(writer, taxon)) {
                mergeTaxonInHDF5(writer, taxon, genos, depth);
            } else {
                addTaxon(writer, taxon, genos, depth);
            }
        } else {
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

    public boolean isHDF5() {
        return isHDF5;
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
     * Finishes building the GenotypeTable. For HDF5 files it locks the taxa and
     * genotype modules so that cannot be modified again.
     *
     * @return a genotype table
     */
    public GenotypeTable build() {
        if (isHDF5) {
            switch (myBuildType) {
                case TAXA_INC: {
                    break;
                }
                case SITE_INC: {
                    //copy the in memory position list to the HDF5 file
                    this.positionList = new PositionListBuilder(writer, posListBuilder.build()).build();
                    break;
                }
            }
            HDF5Utils.lockHDF5TaxaModule(writer);
            String name = writer.file().getFile().getAbsolutePath();
            annotateHDF5File(writer);
            HDF5Utils.writeHDF5Annotation(writer, Tassel5HDF5Constants.ROOT, myAnnotationBuilder.build());
            HDF5Utils.lockHDF5GenotypeModule(writer);
            writer.close();
            return getInstance(name);
        }
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

    /**
     * Used to close an HDF5 GenotypeTableBuilder, when it will be reopened
     * later and appended. This file cannot be used for other purposes in this
     * unfinished state.
     */
    public void closeUnfinished() {
        if (isHDF5 == false) {
            throw new UnsupportedOperationException("Only a HDF5 GenotypeTableBuilder can be closed");
        }
        taxaListBuilder = null;
        writer.close();
    }

    /**
     * HDF5 Alignment section.
     */
    private synchronized void setupGenotypeTaxaInHDF5(IHDF5Writer writer) {
        HDF5Utils.createHDF5TaxaModule(writer);
        HDF5Utils.createHDF5GenotypeModule(writer);
        HDF5Utils.writeHDF5GenotypesMaxNumAlleles(writer, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
        HDF5Utils.writeHDF5GenotypesRetainRareAlleles(writer, false);
        HDF5Utils.writeHDF5GenotypesNumTaxa(writer, 0);
        HDF5Utils.writeHDF5GenotypesAlleleStates(writer, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    /**
     * Code needed to add a Taxon to HDF5
     */
    private synchronized void addTaxon(IHDF5Writer myWriter, Taxon id, byte[] genotype, byte[][] depth) {
        boolean goodAdd = HDF5Utils.addTaxon(myWriter, id);
        if (goodAdd == false) {
            throw new IllegalStateException("Taxon [" + id.getName() + "] already exists in the HDF5 file.  Duplicated taxa not allowed.");
        }
        HDF5Utils.writeHDF5GenotypesCalls(myWriter, id.getName(), genotype);
        if (depth != null) {
            if (depth.length != 6) {
                throw new IllegalStateException("Just set A, C, G, T, -, + all at once");
            }
            if (depth[0].length != positionList.numberOfSites()) {
                throw new IllegalStateException("Setting all depth in addTaxon.  Wrong number of sites");
            }
            HDF5Utils.writeHDF5GenotypesDepth(myWriter, id.getName(), depth);
        }
    }

    private synchronized void mergeTaxonInHDF5(IHDF5Writer myWriter, Taxon id, byte[] genotype, byte[][] depth) {
        GeneralAnnotation annotation = HDF5Utils.getTaxon(writer, id.getName()).getAnnotation();
        String[] existingFlowCellLanes;
        if (annotation == null) {
            existingFlowCellLanes = new String[0];
        } else {
            existingFlowCellLanes = annotation.getTextAnnotation("Flowcell_Lane");
        }
        GeneralAnnotation annotation2 = id.getAnnotation();
        String[] newFlowCellLanes;
        if (annotation2 == null) {
            newFlowCellLanes = new String[0];
        } else {
            newFlowCellLanes = annotation2.getTextAnnotation("Flowcell_Lane");
        }
        if (newFlowCellLanes.length > 0) {
            for (String existingFL : existingFlowCellLanes) {
                if (existingFL.equals(newFlowCellLanes[0])) {
                    throw new IllegalStateException("mergeTaxonInHDF5: Reads from flowcell_lane " + id.getAnnotation().getTextAnnotation("Flowcell_Lane")[0]
                            + " previously added to taxon " + id.getName());
                }
            }
            Taxon modifiedTaxon = new Taxon.Builder(HDF5Utils.getTaxon(writer, id.getName())).addAnno("Flowcell_Lane", newFlowCellLanes[0]).build();
            HDF5Utils.replaceTaxonAnnotations(myWriter, modifiedTaxon);
        }
        byte[] combGenos = new byte[genotype.length];
        if (depth != null) {
            byte[][] existingDepth = HDF5Utils.getHDF5GenotypesDepth(myWriter, id.getName());
            if (existingDepth == null) {
                throw new IllegalStateException("mergeTaxonInHDF5: Trying to merge genotypes with and without depth.");
            }
            byte[][] combDepth = new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][genotype.length];
            byte[] currDepths = new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (int site = 0; site < combDepth[0].length; site++) {
                for (int allele = 0; allele < combDepth.length; allele++) {
                    currDepths[allele] = combDepth[allele][site] = AlleleDepthUtil.addByteDepths(depth[allele][site], existingDepth[allele][site]);
                }
                combGenos[site] = mergeRule.callBasedOnDepth(currDepths);
            }
            HDF5Utils.replaceHDF5GenotypesCalls(myWriter, id.getName(), combGenos);
            HDF5Utils.replaceHDF5GenotypesDepth(myWriter, id.getName(), combDepth);
        } else {
            byte[] existingGenos = HDF5Utils.getHDF5GenotypesCalls(myWriter, id.getName());
            for (int site = 0; site < combGenos.length; site++) {
                combGenos[site] = mergeRule.mergeCalls(genotype[site], existingGenos[site]);
            }
            HDF5Utils.replaceHDF5GenotypesCalls(myWriter, id.getName(), combGenos);
        }
    }

    /**
     * Annotates the HDF5 Genotype file with allele frequency information. Can
     * only be called on unlocked HDF5 files. Currently, placed in the
     * GenotypeTableBuilder as it still above genotypes, taxa, and sites.
     *
     * @param writer
     */
    public static void annotateHDF5File(IHDF5Writer writer) {
        if (HDF5Utils.isHDF5GenotypeLocked(writer)) {
            throw new UnsupportedOperationException("This is a locked HDF5 file");
        }
        int hdf5GenoBlock = Tassel5HDF5Constants.BLOCK_SIZE;
        int sites = writer.int32().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH, Tassel5HDF5Constants.POSITION_NUM_SITES);
        TaxaList tL = new TaxaListBuilder().buildFromHDF5Genotypes(writer);
        int taxa = tL.numberOfTaxa();
        writer.int32().setAttr(Tassel5HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, Tassel5HDF5Constants.GENOTYPES_NUM_TAXA, taxa);
        int[][] af = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][sites];
        byte[][] afOrder = new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][sites];
        float[] coverage = new float[taxa];
        float[] hets = new float[taxa];
        for (int taxon = 0; taxon < taxa; taxon++) {
            String basesPath = Tassel5HDF5Constants.getGenotypesCallsPath(tL.taxaName(taxon));
            byte[] genotype = writer.int8().readArray(basesPath);
            int covSum = 0;  //coverage of the taxon
            int hetSum = 0;
            for (int s = 0; s < sites; s++) {
                byte[] b = GenotypeTableUtils.getDiploidValues(genotype[s]);
                if (b[0] < 6) {
                    af[b[0]][s]++;
                }
                if (b[1] < 6) {
                    af[b[1]][s]++;
                }
                if (GenotypeTableUtils.isHeterozygous(genotype[s])) {
                    hetSum++;
                }
                if (genotype[s] != GenotypeTable.UNKNOWN_GENOTYPE) {
                    covSum++;
                }
            }
            coverage[taxon] = (float) covSum / (float) sites;
            hets[taxon] = (float) hetSum / (float) covSum;
        }
        float[] maf = new float[sites];
        float[] paf = new float[sites];
        int baseMask = 0xF;
        for (int s = 0; s < sites; s++) {
            int sum = 0;
            int[] cntAndAllele = new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (byte i = 0; i < 6; i++) {
                cntAndAllele[i] = (af[i][s] << 4) | (5 - i);  //size | allele (the 5-i is to get the sort right, so if case of ties A is first)
                sum += af[i][s];
            }
            Arrays.sort(cntAndAllele);  //ascending quick sort, there are faster ways
            //http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array
            for (byte i = 0; i < 6; i++) {
                afOrder[5 - i][s] = (cntAndAllele[i] > 0xF) ? ((byte) (5 - (baseMask & cntAndAllele[i]))) : GenotypeTable.UNKNOWN_ALLELE;
            }
            if (afOrder[1][s] != GenotypeTable.UNKNOWN_ALLELE) {
                maf[s] = (float) af[afOrder[1][s]][s] / (float) sum;
            }
            paf[s] = (float) sum / (float) (2 * taxa);
        }
        writer.object().createGroup(Tassel5HDF5Constants.GENO_DESC);
        int chunk = (sites < hdf5GenoBlock) ? sites : hdf5GenoBlock;
        writer.int32().createMatrix(Tassel5HDF5Constants.ALLELE_CNT, 6, sites, 1, chunk, Tassel5HDF5Constants.intDeflation);
        writer.int8().createMatrix(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 6, sites, 1, chunk, Tassel5HDF5Constants.intDeflation);
        writer.float32().createArray(Tassel5HDF5Constants.MAF, sites, chunk, Tassel5HDF5Constants.floatDeflation);
        writer.float32().createArray(Tassel5HDF5Constants.SITECOV, sites, chunk, Tassel5HDF5Constants.floatDeflation);
        //  writer.createGroup(Tassel5HDF5Constants.TAXA_DESC);
        chunk = (tL.numberOfTaxa() < hdf5GenoBlock) ? tL.numberOfTaxa() : hdf5GenoBlock;
        writer.float32().createArray(Tassel5HDF5Constants.TAXACOV, tL.numberOfTaxa(), chunk, Tassel5HDF5Constants.floatDeflation);
        writer.float32().createArray(Tassel5HDF5Constants.TAXAHET, tL.numberOfTaxa(), chunk, Tassel5HDF5Constants.floatDeflation);
        if (af[0].length > 0) {
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.ALLELE_CNT, writer, af[0].length, 1 << 16, af);
        }
        if (afOrder[0].length > 0) {
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.ALLELE_FREQ_ORD, writer, afOrder[0].length, 1 << 16, afOrder);
        }
        if (maf.length > 0) {
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.MAF, writer, maf.length, 1 << 16, maf);
        }
        if (paf.length > 0) {
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.SITECOV, writer, paf.length, 1 << 16, paf);
        }

        if (coverage.length > 0) {
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.TAXACOV, writer, coverage.length, 1 << 16, coverage);
        }
        if (hets.length > 0) {
            System.out.println("Number of taxa in HDF5 file:" + hets.length);
            HDF5Utils.writeHDF5EntireArray(Tassel5HDF5Constants.TAXAHET, writer, hets.length, 1 << 16, hets);
        }
    }

    /**
     * Annotates the HDF5 Genotype file with the reference allele and reference
     * genome version. Can only be called on unlocked HDF5 files.
     *
     * @param writer
     * @param refAlleles
     */
    public static void annotateHDF5FileWithRefAllele(IHDF5Writer writer, byte[] refAlleles) {
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
