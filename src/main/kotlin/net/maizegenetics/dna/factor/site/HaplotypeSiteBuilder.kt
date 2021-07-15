package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.factor.UNKNOWN_ALLELE
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrixBuilder

class HaplotypeSiteBuilder constructor(
    val factor: GenomicFeature,
    val taxa: TaxaList,
    val ploidy: Int,
    val strStates: Array<String>? = null,
    hapAnnotations: Array<HaplotypeAnnotation?>? = null
) {

    private val genotypes = SuperByteMatrixBuilder.getInstance(taxa.size, ploidy)
    private var myHapAnnotations: Array<HaplotypeAnnotation?>? = null

    init {
        require(hapAnnotations == null || strStates != null) { "strStates must be specified when supplying hapAnnotations" }
        require(hapAnnotations == null || (hapAnnotations.size == strStates!!.size)) { "strStates size and hapAnnotations size much equal" }
        myHapAnnotations = if (hapAnnotations?.filterNotNull()?.count() == 0) {
            null
        } else {
            hapAnnotations
        }
        genotypes.setAll(UNKNOWN_ALLELE)
    }

    private val stateMap = strStates?.mapIndexed { index, str -> Pair(str, index.toByte()) }?.toMap() ?: linkedMapOf()

    private var nextCode = 0.toByte()

    private val taxaMap = taxa
            .mapIndexed { index, taxon -> Pair(taxon.name, index) }
            .toMap()

    var isPhased = false

    fun set(taxon: String, values: List<String>): HaplotypeSiteBuilder {
        return set(taxaMap[taxon] ?: throw IllegalArgumentException("$taxon not in taxa list"), values)
    }

    fun set(taxon: Int, values: List<String>): HaplotypeSiteBuilder {

        values
                .map { str ->
                    when (str) {
                        "" -> 0xFF.toByte()
                        "." -> 0xFF.toByte()
                        else -> {
                            if (strStates == null) {
                                if (stateMap[str] == null) (stateMap as MutableMap)[str] = nextCode++
                                stateMap[str] ?: error("This can not be null")
                            } else {
                                stateMap[str]
                                        ?: throw IllegalArgumentException("$str not on allele list: ${strStates.joinToString(",")}")
                            }
                        }
                    }
                }
            .forEachIndexed { index, value ->
                genotypes.set(taxon, index, value)
            }

        return this

    }

    fun build() = HaplotypeSite(
        factor, taxa, strStates ?: stateMap.keys.toTypedArray(),
        genotypes, ploidy, isPhased = isPhased,
        hapAnnotations = myHapAnnotations
    )

}