package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.SNPSite
import net.maizegenetics.dna.map.GenomicFactorList
import net.maizegenetics.dna.map.SNP
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrix
import net.maizegenetics.util.SuperByteMatrixBuilder

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class SNPTable private constructor(taxa: TaxaList, factors: GenomicFactorList, private val matrix: SuperByteMatrix, val isPhased: Boolean = false) : FactorTable(taxa, factors) {

    override fun site(index: Int): SNPSite {
        return SNPSite(factors[index] as SNP, index, taxa, matrix.getAllRows(index))
    }

    class SNPTableBuilder private constructor(val numTaxa: Int, val numFactors: Int, private val matrix: SuperByteMatrix) {

        var taxa: TaxaList? = null
            set(value) {
                if (value?.numberOfTaxa() != numTaxa) {
                    throw IllegalArgumentException("SNPTableBuilder: set: taxa list size: ${value?.numberOfTaxa()} not equal number of taxa: $numTaxa")
                }
                field = value
            }

        var factors: GenomicFactorList? = null
            set(value) {
                if (value?.size != numFactors) {
                    throw IllegalArgumentException("SNPTableBuilder: set: factor list size: ${value?.size} not equal number of factors: $numFactors")
                }
                field = value
            }

        var isPhased = false

        fun set(taxon: Int, site: Int, value: Byte): SNPTableBuilder {
            matrix.set(taxon, site, value)
            return this
        }

        fun setRange(taxon: Int, startSite: Int, value: ByteArray): SNPTableBuilder {
            matrix.arraycopy(taxon, value, startSite)
            return this
        }

        fun set(data: Array<String>): SNPTableBuilder {

            val numTaxa = data.size

            val numSites = data[0].length

            for (site in 0 until numSites) {
                for (taxon in 0 until numTaxa) {
                    set(taxon, site, NucleotideAlignmentConstants.getNucleotideDiploidByte(data[taxon][site]))
                }
            }

            return this

        }

        fun build(): SNPTable {
            if (taxa == null) {
                throw IllegalStateException("SNPTableBuilder: build: taxa list not set")
            }
            if (factors == null) {
                throw IllegalStateException("SNPTableBuilder: build: factor list not set")
            }
            return SNPTable(taxa!!, factors!!, matrix, isPhased)
        }

        companion object {

            /**
             * Get SNP Table Builder given number of taxa and sites. Performance
             * optimized for site loop inside taxon loop. Default is unphased and
             * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
             */
            @JvmStatic
            fun instance(numTaxa: Int, numFactors: Int): SNPTableBuilder {
                val matrix = SuperByteMatrixBuilder.getInstance(numTaxa, numFactors)
                matrix.setAll(GenotypeTable.UNKNOWN_GENOTYPE)
                return SNPTableBuilder(numTaxa, numFactors, matrix)
            }

        }
    }

}