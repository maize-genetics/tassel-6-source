/*
 * HDF5Constants
 */
package net.maizegenetics.util;

import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;

/**
 * Definition of attributes and paths for Tassel HDF5 file format.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public final class Tassel5HDF5Constants {

    public static final String ROOT = "/";

    // Genotypes Module
    public static final String GENOTYPES_MODULE = "Genotypes";
    public static final String GENOTYPES_ATTRIBUTES_PATH = GENOTYPES_MODULE + "/";
    public static final String GENOTYPES_MAX_NUM_ALLELES = "maxNumAlleles";
    public static final String GENOTYPES_RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final String GENOTYPES_NUM_TAXA = "numTaxa";
    public static final String GENOTYPES_LOCKED = "locked";
    public static final String GENOTYPES_SCORE_TYPE = "scoreType";
    public static final String GENOTYPES_ALLELE_STATES = GENOTYPES_MODULE + "/AlleleStates";
    public static final String GENO_DESC = GENOTYPES_MODULE + "/_Descriptors/";
    public static final String ALLELE_CNT = GENO_DESC + "AlleleCnt";
    public static final String MAF = GENO_DESC + "MAF";
    public static final String SITECOV = GENO_DESC + "SiteCoverage";
    public static final String ALLELE_FREQ_ORD = GENO_DESC + "AlleleFreqOrder";
    public static final String TAXACOV = GENO_DESC + "TaxaCoverage";
    public static final String TAXAHET = GENO_DESC + "TaxaHet";

    public static final int BLOCK_SIZE = 1 << 16;

    public static final String getGenotypesPedigreePath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/pedigree";
    }

    public static final String getGenotypesCallsPath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/calls";
    }

    @Deprecated
    public static final String getGenotypesDepthPath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/depth";
    }

    public static final String getGenotypesSiteScorePath(String taxon, String siteScoreType) {
        return GENOTYPES_MODULE + "/" + taxon + "/" + siteScoreType;
    }

    //Taxa Module
    public static final String TAXA_MODULE = "Taxa";
    public static final String TAXA_ATTRIBUTES_PATH = TAXA_MODULE + "/";
    public static final String TAXA_NUM_TAXA = "numTaxa";
    public static final String TAXA_LOCKED = "locked";
    public static final String TAXA_ORDER = TAXA_ATTRIBUTES_PATH+"TaxaOrder";

    public static final String getTaxonPath(String taxon) {
        return TAXA_MODULE + "/" + taxon;
    }

    //Position Module
    public static final String POSITION_MODULE = "Positions";
    public static final String POSITION_ATTRIBUTES_PATH = POSITION_MODULE + "/";
    public static final String POSITION_NUM_SITES = "numSites";
    public static final String POSITION_HAS_REFEFERENCE = "hasReferenceAlleles";
    public static final String POSITION_GENOME_VERSION = "genomeVersion";
    public static final String POSITIONS = POSITION_ATTRIBUTES_PATH + "Positions";
    public static final String CHROMOSOMES = POSITION_ATTRIBUTES_PATH + "Chromosomes";
    public static final String CHROMOSOME_INDICES = POSITION_ATTRIBUTES_PATH + "ChromosomeIndices";
    public static final String SNP_IDS = POSITION_ATTRIBUTES_PATH + "SnpIds";
    public static final String REF_ALLELES = POSITION_ATTRIBUTES_PATH + "ReferenceAlleles";
    public static final String ANC_ALLELES = POSITION_ATTRIBUTES_PATH + "AncestralAlleles";

    //Standard Compression (deflation) levels
    public static final HDF5IntStorageFeatures intDeflation = HDF5IntStorageFeatures.createDeflation(2);
    public static final HDF5GenericStorageFeatures genDeflation = HDF5GenericStorageFeatures.createDeflation(2);
    public static final HDF5FloatStorageFeatures floatDeflation = HDF5FloatStorageFeatures.createDeflation(2);

    //Tag Module
    public static final String TAG_MODULE = "Tags";
    public static final String TAG_ATTRIBUTES_PATH = TAG_MODULE + "/";
    public static final String TAG_COUNT = "tagCount";
    public static final String TAG_LENGTH_LONG = "tagLengthLong";
    public static final String TAG_LOCKED = "locked";
    public static final String TAGS = TAG_MODULE + "/Tags";
    public static final int TAGS_BIN_NUM = 64;
    public static final int HASH_SHIFT_TO_TAG_BIN=Integer.numberOfLeadingZeros(TAGS_BIN_NUM)+1;

    public static final String TAG_SEQ = "TagSeq";
    public static final String TAG_LENGTHS = "TagLength";
    public static final String TAG_DIST = "TagDist";
    public static final String TAG_DIST_OFFSETS = "TagTaxaDistOffset";
//    public static final String TAG_DIST_CHUNK = "taxaDirection";

    private Tassel5HDF5Constants() {
        // do not instantiate
    }

}
