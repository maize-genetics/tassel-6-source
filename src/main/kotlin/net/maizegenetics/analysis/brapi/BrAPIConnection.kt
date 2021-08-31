package net.maizegenetics.analysis.brapi

import khttp.get
import khttp.responses.Response
import kotlinx.serialization.json.*
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.dna.factor.FeatureTableBuilder
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.dna.map.GenomicFeatureList
import net.maizegenetics.gui.basicLoggingInfo
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.taxa.Taxon
import net.maizegenetics.util.setupDebugLogging
import org.apache.logging.log4j.LogManager
import kotlin.system.measureNanoTime


class BrAPIConnection(val baseURL: String) {

    private val logger = LogManager.getLogger(BrAPIConnection::class.java)

    /**
     * Returns list of variant sets available from BrAPI Server.
     *
     * Example: http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets
     */
    fun getVariantsets(): List<VariantSet>? {

        val url = "$baseURL/variantsets"
        logger.info("getVariantsets: query: $url")
        val json = Json.parseToJsonElement(get(url).text)

        val variantSetsArray = json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray

        return variantSetsArray
                ?.map { variantset(it.jsonObject) }
                ?.toList()

    }

    /**
     * Returns the variant set information for the given id.
     *
     * Example: http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid
     */
    fun getVariantset(id: String): VariantSet? {

        val url = "$baseURL/variantsets/$id"
        logger.info("getVariantset: query: $url")
        val json = Json.parseToJsonElement(get(url).text)

        return json.jsonObject["result"]?.jsonObject?.let {
            variantset(it)
        }

    }

    private fun variantset(json: JsonObject): VariantSet {
        return VariantSet(json["variantSetDbId"].toString(), json["variantCount"].toString().toInt(), json["callSetCount"].toString().toInt())
    }

    /**
     * Returns a list of genomic features (i.e. position ranges))
     * for the given variant set id.
     *
     * Example: http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid/variants
     */
    fun getVariants(id: String): GenomicFeatureList {

        val url = "$baseURL/variantsets/$id/variants"
        logger.info("getVariants: query: $url")
        val json = Json.parseToJsonElement(get(url, timeout = 0.0).text)

        val builder = GenomicFeatureList.Builder()
        json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray
                ?.forEach { builder.add(genomicFeature(it.jsonObject)) }

        val result = builder.build()

        logger.info("getVariants: number of features: ${result.size}")

        return result

    }

    private fun genomicFeature(json: JsonElement): GenomicFeature {
        val chromosome = Chromosome.instance(getStringAttribute(json, "referenceName"))
        val startPos = json.jsonObject["start"].toString().toInt()
        val endPos = json.jsonObject["end"].toString().toInt()
        val name = getStringFirstArrayElement(json, "variantNames")
        val id = getStringAttribute(json, "variantDbId")
        return GenomicFeature(chromosome, startPos, chromosome, endPos, name, id)
    }

    /**
     * Returns a list of taxa for the given variant set id.
     *
     * Example: http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid/callsets
     */
    fun getCallsets(id: String): TaxaList {

        val url = "$baseURL/variantsets/$id/callsets"
        logger.info("getCallsets: query: $url")
        val json = Json.parseToJsonElement(get(url, timeout = 0.0).text)

        val builder = TaxaListBuilder()
        json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray
                ?.forEach { builder.add(taxon(it)) }

        val result = builder.build()

        logger.info("getCallsets: number of taxa: ${result.numberOfTaxa()}")

        return result

    }

    private fun taxon(json: JsonElement): Taxon {
        return Taxon.instance(getStringAttribute(json, "callSetName"))
    }

    private val callsPageSize = 1

    /**
     * Returns a feature table for the given variant set id.
     *
     * Example: http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid/calls?pageSize=1
     */
    fun getCalls(id: String): FeatureTable {

        val builder = FeatureTableBuilder(getCallsets(id), getVariants(id))

        var url = "$baseURL/variantsets/$id/calls?pageSize=$callsPageSize"
        logger.info("getCalls: query: $url")
        var json = Json.parseToJsonElement(get(url, timeout = 0.0).text)

        val metadata = metadata(json.jsonObject)

        val numPages = metadata["pagination.totalPages"]?.toInt()
                ?: error("BrAPIConnection: getCalls: can't get number of pages.")

        setGenotypes(json, builder)

        (1 until numPages).forEach {
            url = "$baseURL/variantsets/$id/calls?pageSize=$callsPageSize;page=$it"
            logger.info("getCalls: query: $url")
            lateinit var response: Response
            var notFinished = true
            while (notFinished) {
                try {
                    response = get(url, timeout = 0.0)
                    notFinished = false
                } catch (e: Exception) {
                    logger.warn("getCalls: BrAPI Server query failed: ${e.message}")
                }
            }
            json = Json.parseToJsonElement(response.text)
            setGenotypes(json, builder)
        }

        return builder.build()

    }

    private fun setGenotypes(json: JsonElement, builder: FeatureTableBuilder) {
        json.jsonObject["result"]?.jsonObject?.get("data")?.jsonArray
                ?.forEach { genotype ->
                    val taxon = getStringAttribute(genotype, "callSetName")
                    val id = getStringAttribute(genotype, "variantDbId")
                    val genotypeArray = genotype.jsonObject["genotype"]?.jsonObject?.get("values")?.jsonArray
                    check(genotypeArray != null) { "BrAPIConnection: getCalls: genotype values can't be null" }
                    val genotypeValues = genotypeArray.map { it.toString() }
                    builder.set(taxon, id, genotypeValues)
                }
    }

    private fun getStringAttribute(element: JsonElement, attribute: String): String {
        return element.jsonObject[attribute].toString().removePrefix("\"").removeSuffix("\"").trim()
    }

    private fun getStringFirstArrayElement(element: JsonElement, attribute: String): String {
        return element.jsonObject[attribute]?.jsonArray?.get(0).toString().removePrefix("\"").removeSuffix("\"").trim()
    }

}

data class VariantSet(val id: String, val numVariants: Int, val numCallsets: Int)

fun main() {

    setupDebugLogging()
    basicLoggingInfo()

    val connection = BrAPIConnection("http://cbsudc01.biohpc.cornell.edu/brapi/v2")

    //val taxa = connection.getCallsets("MergedReadMapping_AllNamParents_Haploid")
    //taxa.take(100).forEach { println(it) }

    //taxa.forEachIndexed { index, taxon ->
    //    if (taxon.name == "Z001E0104" || taxon.name == "Z001E0103") println("matched: $index  name: ${taxon.name}")
    //}

    //System.exit(0)

    //println(connection.getVariantsets())

    //println(connection.getVariantset("Ames_MergedReadMapping_AllLines_Haploid"))

    //println(connection.getVariants("Ames_MergedReadMapping_AllLines_Haploid").toList())

    //println(connection.getCallsets("Ames_MergedReadMapping_AllLines_Haploid").numberOfTaxa())

    var featureTable: FeatureTable
    val time = measureNanoTime {
        featureTable = connection.getCalls("MergedReadMapping_AllNamParents_Haploid")

        println("num of features: ${featureTable.numFeatures()}")
        println("num of taxa: ${featureTable.numTaxa()}")
    }

    val taxonIndex = featureTable.taxa.indexOf("Z001E0104")
    featureTable.forEach { site ->
        val genotype = site.genotypeAsString(taxonIndex)
        println("taxon: Z001E0104  id: ${site.feature.id}  name: ${site.feature.name}  chr: ${site.feature.startChr}  pos: ${site.feature.startPos}  genotype: $genotype")
    }

    println("time: ${time / 1e9} secs")

}

fun main1() {

    val connection = BrAPIConnection("http://cbsudc01.biohpc.cornell.edu/brapi/v2")

    //val taxa = connection.getCallsets("MergedReadMapping_AllNamParents_Haploid")
    //taxa.take(100).forEach { println(it) }

    //taxa.forEachIndexed { index, taxon ->
    //    if (taxon.name == "Z001E0104" || taxon.name == "Z001E0103") println("matched: $index  name: ${taxon.name}")
    //}

    //System.exit(0)

    //println(connection.getVariantsets())

    //println(connection.getVariantset("Ames_MergedReadMapping_AllLines_Haploid"))

    //println(connection.getVariants("Ames_MergedReadMapping_AllLines_Haploid").toList())

    //println(connection.getCallsets("Ames_MergedReadMapping_AllLines_Haploid").numberOfTaxa())

    var featureTable: FeatureTable
    val time = measureNanoTime {
        featureTable = connection.getCalls("MergedReadMapping_AllNamParents_Haploid")

        println("num of features: ${featureTable.numFeatures()}")
        println("num of taxa: ${featureTable.numTaxa()}")
    }

    val taxonIndex = featureTable.taxa.indexOf("Z001E0104")
    featureTable.forEach { site ->
        site.forEach { alleles ->
            // TODO
        }
        val genotype = site.genotypeAsString(taxonIndex)
        println("taxon: Z001E0104  id: ${site.feature.id}  name: ${site.feature.name}  chr: ${site.feature.startChr}  pos: ${site.feature.startPos}  genotype: $genotype")
    }

    println("time: ${time / 1e9} secs")

}