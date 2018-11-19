package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.HaplotypeSite
import net.maizegenetics.dna.map.GenomicFactorList
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrix
import net.maizegenetics.util.SuperByteMatrixBuilder

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class HaplotypeTable private constructor(taxa: TaxaList, factors: GenomicFactorList, private val matrix: SuperByteMatrix, private val haplotypes: Array<Array<String>>) : FactorTable(taxa, factors) {


    override fun site(index: Int): HaplotypeSite {
        TODO("not implemented")
    }

    class HaplotypeTableBuilder private constructor(val taxa: TaxaList, val factors: GenomicFactorList, private val matrix: SuperByteMatrix) {

        val numTaxa = taxa.numberOfTaxa()
        val haplotypes = ArrayList<Array<String>>(factors.size)

        fun set(site: Int, siteHaplotypes: Array<String>, values: ByteArray): HaplotypeTableBuilder {

            if (values.size != numTaxa) {
                throw IllegalArgumentException("HaplotypeTableBuilder: set: number of values: ${values.size} not equal to number taxa: $numTaxa")
            }
            val numHaplotypes = siteHaplotypes.size
            haplotypes[site] = siteHaplotypes
            values.forEachIndexed { index, value ->
                if (value >= numHaplotypes) {
                    throw IllegalArgumentException("HaplotypeTableBuilder: set: encoded value: $value doesn't have haplotype")
                }
                matrix.set(index, site, value)
            }

            return this
        }

        fun build(): HaplotypeTableBuilder {
            TODO()
        }

        companion object {
            @JvmStatic
            fun instance(taxa: TaxaList, factors: GenomicFactorList): HaplotypeTableBuilder {
                val matrix = SuperByteMatrixBuilder.getInstance(taxa.numberOfTaxa(), factors.size)
                matrix.setAll(GenotypeTable.UNKNOWN_DIPLOID_ALLELE)
                return HaplotypeTableBuilder(taxa, factors, matrix)
            }
        }

    }


}