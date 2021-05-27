package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class SNPSite(feature: GenomicFeature, taxa: TaxaList, private val values: ByteArray, weight: Double? = null, isPhased: Boolean = false) : FeatureSite(feature, taxa, weight, isPhased) {

    init {
        require(taxa.size == values.size) { "Number of taxa: ${taxa.size} should match number of genotypes: ${values.size}." }
    }

    override fun ploidy() = 2

    override fun genotype(taxon: Int): ByteArray = GenotypeTableUtils.getDiploidValues(values[taxon])

    override fun genotypeAsString(taxon: Int): String = NucleotideAlignmentConstants.getNucleotideIUPAC(values[taxon])

}