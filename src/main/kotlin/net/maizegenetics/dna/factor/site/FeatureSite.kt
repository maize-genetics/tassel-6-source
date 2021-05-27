package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.factor.UNKNOWN_ALLELE
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 14, 2018
 */

abstract class FeatureSite(val feature: GenomicFeature, val taxa: TaxaList, val weight: Double? = null, val isPhased: Boolean = false) : Comparable<FeatureSite>, Sequence<ByteArray> {

    val alleleStats by lazy { AlleleStats(this) }

    /**
     * Ploidy of this site. Value of 2 for example is diploid.
     */
    abstract fun ploidy(): Int

    /**
     * This returns an array of allele values for the given taxon index.
     * A diploid site would return and array of size 2
     */
    abstract fun genotype(taxon: Int): ByteArray

    abstract fun genotypeAsString(taxon: Int): String

    fun heterozygousCount(): Int {
        return this
                .map { alleles -> alleles.filter { it != UNKNOWN_ALLELE }.distinct().count() }
                .filter { count -> count > 1 }
                .count()
    }

    override fun compareTo(other: FeatureSite): Int {
        return feature.compareTo(other.feature)
    }

    override fun iterator(): Iterator<ByteArray> {
        return GenotypeIterator(this)
    }

}

private class GenotypeIterator(val site: FeatureSite) : Iterator<ByteArray> {

    private var currentTaxon = 0

    override fun hasNext(): Boolean {
        return currentTaxon < site.taxa.numberOfTaxa()
    }

    override fun next(): ByteArray {
        return site.genotype(currentTaxon++)
    }

}