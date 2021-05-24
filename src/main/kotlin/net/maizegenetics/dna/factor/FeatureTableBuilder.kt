package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.site.FeatureSite
import net.maizegenetics.dna.factor.site.HaplotypeSite
import net.maizegenetics.dna.factor.site.HaplotypeSiteBuilder
import net.maizegenetics.dna.factor.site.SNPSite
import net.maizegenetics.dna.map.GenomicFeatureList
import net.maizegenetics.taxa.TaxaList
import kotlin.reflect.KClass

class FeatureTableBuilder constructor(val taxa: TaxaList, features: GenomicFeatureList? = null, type: KClass<out FeatureSite> = HaplotypeSite::class, ploidy: Int = 2) {

    init {

        when (type) {
            HaplotypeSite::class -> {
                features
                        ?.forEach { feature ->
                            val builder = HaplotypeSiteBuilder(feature, taxa, ploidy)
                            haplotypeSiteBuilders.add(builder)
                            feature.name?.let { featureMap[it] = builder }
                        }
            }
            SNPSite::class -> {
                TODO()
            }
        }

    }

    private val featureMap = mutableMapOf<String, HaplotypeSiteBuilder>()

    private val sites = mutableListOf<FeatureSite>()

    private val haplotypeSiteBuilders = mutableListOf<HaplotypeSiteBuilder>()

    fun add(site: FeatureSite) {
        sites.add(site)
    }

    fun set(taxon: String, feature: String, values: List<String>) {
        featureMap[feature]?.set(taxon, values) ?: error("FeatureTableBuilder: set: feature: $feature not found")
    }

    fun set(taxon: Int, feature: Int, values: List<String>) {
        haplotypeSiteBuilders[feature].set(taxon, values)
    }

    fun build(): FeatureTable {
        haplotypeSiteBuilders.forEach { sites.add(it.build()) }
        sites.sort()
        return FeatureTable(taxa, sites)
    }

}