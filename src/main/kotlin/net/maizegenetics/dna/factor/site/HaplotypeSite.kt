package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.factor.UNKNOWN_ALLELE
import net.maizegenetics.dna.factor.UNKNOWN_ALLELE_STR
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrix
import org.apache.logging.log4j.LogManager

/**
 * @author Terry Casstevens
 * Created April 6, 2021
 */

class HaplotypeSite(
    feature: GenomicFeature,
    taxa: TaxaList,
    private val strStates: Array<String>,
    private val genotypes: SuperByteMatrix,
    val ploidy: Int = 2,
    weight: Double? = null,
    isPhased: Boolean = false,
    private val hapAnnotations: Array<HaplotypeAnnotation?>? = null
) : FeatureSite(feature, taxa, weight, isPhased) {

    private val logger = LogManager.getLogger(HaplotypeSite::class.java)

    init {
        require(taxa.size == genotypes.numRows) { "Number of taxa: ${taxa.size} should match number of genotypes: ${genotypes.numRows}." }
    }

    override fun ploidy() = ploidy

    override fun genotype(taxon: Int): ByteArray {
        return genotypes.getAllColumns(taxon)
    }

    override fun genotypeAsString(taxon: Int): String {
        return genotypes.getAllColumns(taxon).joinToString(if (isPhased) "|" else "/") { alleleCode ->
            when (alleleCode) {
                UNKNOWN_ALLELE -> {
                    UNKNOWN_ALLELE_STR
                }
                else -> {
                    strStates[alleleCode.toInt()]
                }
            }
        }
    }

    /**
     * Returns the Haplotype Annotations for given allele code.
     * If no annotation, then returns null
     */
    fun haplotypeAnnotation(allele: Byte): HaplotypeAnnotation? {
        if (allele == UNKNOWN_ALLELE) return null
        try {
            return hapAnnotations?.get(allele.toInt())
        } catch (ie: IndexOutOfBoundsException) {
            logger.warn("haplotypeAnnotation: allele value: $allele doesn't exist for this site")
            return null
        }
    }

}

data class HaplotypeAnnotation(val taxon: String, val asmContig: String, val asmStart: Int, val asmEnd: Int)