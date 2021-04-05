package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrix
import java.lang.ref.SoftReference

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class SNPSite(factor: GenomicFactor, index: Int, taxa: TaxaList, private val matrix: SuperByteMatrix, weight: Double? = null, isPhased: Boolean = false) : FactorSite(factor, index, taxa, weight, isPhased) {

    private var values = SoftReference<ByteArray>(null)

    private fun values(): ByteArray {
        var result = values.get()
        if (result == null) {
            result = matrix.getAllRows(index)
            values = SoftReference(result)
        }
        return result!!
    }

    override fun ploidy() = 2

    override fun genotype(taxon: Int): ByteArray = GenotypeTableUtils.getDiploidValues(values()[taxon])

    override fun genotypeAsString(taxon: Int): String = NucleotideAlignmentConstants.getNucleotideIUPAC(values()[taxon])

}