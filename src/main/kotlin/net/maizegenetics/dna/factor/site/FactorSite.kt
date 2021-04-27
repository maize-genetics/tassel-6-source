package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 14, 2018
 */

abstract class FactorSite(val factor: GenomicFactor, val taxa: TaxaList, val weight: Double? = null, val isPhased: Boolean = false) : Comparable<FactorSite>, Sequence<ByteArray> {

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

    override fun compareTo(other: FactorSite): Int {
        return factor.compareTo(other.factor)
    }

    override fun iterator(): Iterator<ByteArray> {
        return GenotypeIterator(this)
    }

}

private class GenotypeIterator(val site: FactorSite) : Iterator<ByteArray> {

    private var currentTaxon = 0

    override fun hasNext(): Boolean {
        return currentTaxon < site.taxa.numberOfTaxa()
    }

    override fun next(): ByteArray {
        return site.genotype(currentTaxon++)
    }

}