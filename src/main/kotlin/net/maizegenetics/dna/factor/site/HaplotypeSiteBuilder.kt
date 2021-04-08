package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.util.SuperByteMatrixBuilder

class HaplotypeSiteBuilder constructor(val factor: GenomicFactor, val taxa: TaxaList, val ploidy: Int, val strStates: Array<String>) {

    private val genotypes = SuperByteMatrixBuilder.getInstance(taxa.size, ploidy)

    private val stateMap = strStates
            .mapIndexed { index, str -> Pair(str, index.toByte()) }
            .toMap()

    private val taxaMap = taxa
            .mapIndexed { index, taxon -> Pair(taxon.name, index) }
            .toMap()

    var isPhased = false

    fun set(taxon: Int, values: ByteArray): HaplotypeSiteBuilder {
        values.forEachIndexed { index, value ->
            genotypes.set(taxon, index, value)
        }
        return this
    }

    fun set(taxon: String, values: List<String>): HaplotypeSiteBuilder {
        values.forEachIndexed { index, str ->
            genotypes.set(taxaMap[taxon] ?: throw IllegalArgumentException("$taxon not in taxa list"),
                    index,
                    stateMap[str] ?: throw IllegalArgumentException("$str not on allele list: ${strStates.joinToString(",")}"))
        }
        return this
    }

    fun build() = HaplotypeSite(factor, taxa, strStates, genotypes, ploidy, isPhased = isPhased)

}