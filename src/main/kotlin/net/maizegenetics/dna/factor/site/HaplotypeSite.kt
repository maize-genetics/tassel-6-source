package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrix
import java.lang.ref.SoftReference

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class HaplotypeSite(factor: GenomicFactor, index: Int, taxa: TaxaList, private val matrices: Array<SuperByteMatrix>, weight: Double? = null, isPhased: Boolean = false) : FactorSite(factor, index, taxa, weight, isPhased) {

    private var values = SoftReference<Array<ByteArray>>(null)

    private fun values(): Array<ByteArray> {
        var result = values.get()
        if (result == null) {
            result = matrices
                    .map { it.getAllRows(index) }
                    .toTypedArray()
            values = SoftReference(result)
        }
        return result!!
    }

    override fun ploidy() = matrices.size

    override fun genotype(taxon: Int): ByteArray {
        return values()
                .map { it[taxon] }
                .toByteArray()
    }

    override fun genotypeAsString(taxon: Int): String {
        TODO("Not yet implemented")
    }

}