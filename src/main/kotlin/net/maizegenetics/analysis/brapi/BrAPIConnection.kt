package net.maizegenetics.analysis.brapi

import khttp.get
import kotlinx.serialization.json.*
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.dna.map.GenomicFeatureList
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.taxa.Taxon
import net.maizegenetics.util.LoggingUtils
import org.apache.log4j.Logger


class BrAPIConnection(val baseURL: String) {

    private val logger = Logger.getLogger(BrAPIConnection::class.java)

    fun getVariantsets(): List<VariantSet>? {

        val url = "$baseURL/variantsets"
        logger.info("getVariantsets: query: $url")
        val json = Json.parseToJsonElement(get(url).text)

        val variantSetsArray = json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray

        return variantSetsArray
                ?.map { variantset(it.jsonObject) }
                ?.toList()

    }

    fun getVariantset(id: String): VariantSet? {

        val url = "$baseURL/variantsets/$id"
        logger.info("getVariantset: query: $url")
        val json = Json.parseToJsonElement(get(url).text)

        return json.jsonObject["result"]?.jsonObject?.let {
            variantset(it)
        }

    }

    // http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid
    private fun variantset(json: JsonObject): VariantSet {
        return VariantSet(json["variantSetDbId"].toString(), json["variantCount"].toString().toInt(), json["callSetCount"].toString().toInt())
    }

    fun getVariants(id: String): GenomicFeatureList {

        val url = "$baseURL/variantsets/$id/variants"
        logger.info("getVariants: query: $url")
        val json = Json.parseToJsonElement(get(url, timeout = 0.0).text)

        val builder = GenomicFeatureList.Builder()
        json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray
                ?.forEach { builder.add(genomicFeature(it.jsonObject)) }

        return builder.build()

    }

    // http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid/variants
    private fun genomicFeature(json: JsonObject): GenomicFeature {
        val chromosome = Chromosome.instance(json["referenceName"].toString())
        val startPos = json["start"].toString().toInt()
        val endPos = json["end"].toString().toInt()
        return GenomicFeature(chromosome, startPos, chromosome, endPos)
    }

    fun getCallsets(id: String): TaxaList {

        val url = "$baseURL/variantsets/$id/callsets"
        logger.info("getCallsets: query: $url")
        val json = Json.parseToJsonElement(get(url, timeout = 0.0).text)

        val builder = TaxaListBuilder()
        json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray
                ?.forEach { builder.add(taxon(it.jsonObject)) }

        return builder.build()

    }

    // http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid/callsets
    private fun taxon(json: JsonObject): Taxon {
        return Taxon(json["callSetName"].toString())
    }

}

data class VariantSet(val id: String, val numVariants: Int, val numCallsets: Int)

fun main() {
    LoggingUtils.setupDebugLogging()
    val connection = BrAPIConnection("http://cbsudc01.biohpc.cornell.edu/brapi/v2")
    println(connection.getVariantsets())

    println(connection.getVariantset("Ames_MergedReadMapping_AllLines_Haploid"))

    //println(connection.getVariants("Ames_MergedReadMapping_AllLines_Haploid").toList())

    println(connection.getCallsets("Ames_MergedReadMapping_AllLines_Haploid").numberOfTaxa())
}