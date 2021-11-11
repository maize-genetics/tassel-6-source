package net.maizegenetics.plugindef

import khttp.get
import kotlinx.serialization.json.Json
import kotlinx.serialization.json.JsonElement
import kotlinx.serialization.json.jsonArray
import kotlinx.serialization.json.jsonObject
import net.maizegenetics.analysis.data.ExportPlugin
import net.maizegenetics.analysis.distance.KinshipPlugin
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.dna.factor.FeatureTableBuilder
import net.maizegenetics.dna.factor.io.BuilderFromHaplotypeVCF
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.dna.map.GenomicFeatureList
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.taxa.Taxon
import net.maizegenetics.taxa.distance.DistanceMatrix

class RTasselHelpersTassel6 {
    /**
     * Return 2d array of hap IDs (as converted bytes) from table
     * @author Brandon Monier
     * Created 2021-04-23
     */
    fun getHapArray(ft: FeatureTable): Array<IntArray> {
        val hapArray = Array(ft.numFeatures()) { IntArray(ft.numTaxa()) }

        for (j in 0 until ft.numFeatures()) {
            for (i in 0 until ft.numTaxa()) {
                hapArray[j][i] = ft[j].genotype(i)[0].toInt()
            }
        }
        return hapArray
    }

    /**
     * Return string array of taxa
     * @author Brandon Monier
     * @param ft A FeatureTable object
     * Created 2021-04-30
     */
    fun getTaxaArray(ft: FeatureTable): Array<String?> {
        val taxaArray = Array<String?>(ft.numTaxa()) { null }
        for (i in 0 until ft.numTaxa()) {
            taxaArray[i] = ft.taxa[i].toString()
        }
        return taxaArray
    }

    /**
     * Return 2d array of ref range coordinates
     * @author Brandon Monier
     * @param ft A FeatureTable object
     * Created 2021-04-30
     */
    public fun getRefRanges(ft: FeatureTable): Array<Array<String?>> {
        val obs = 4 // [start chr, end chr, start pos, end pos]
        val rrArray = Array(obs) {Array<String?>(ft.numFeatures()){null} }

        for (i in 0 until ft.numFeatures()) {
            rrArray[0][i] = ft[i].feature.startChr.toString()
            rrArray[1][i] = ft[i].feature.endChr.toString()
            rrArray[2][i] = ft[i].feature.startPos.toString()
            rrArray[3][i] = ft[i].feature.endPos.toString()
        }
        return rrArray
    }

    /**
     * Return 2d array of distance metrics
     * @author Brandon Monier
     * Created 2021-04-30
     */
    public fun getDistArray(dist: DistanceMatrix): Array<FloatArray> {
        val distArray = Array(dist.size) { FloatArray(dist.size) }
        for (i in 0 until dist.size) {
            for (j in 0 until dist.size) {
                distArray[i][j] = dist.getDistance(i, j)
            }
        }
        return distArray
    }

    private fun getUniqueHapIdsFromVariants(urlVariant: String, pageSize: Int): MutableList<IntArray> {
        val varPattern = Regex("/variants$")
        val varConnector = if (varPattern.containsMatchIn(urlVariant)) "?" else "&"
        val varPageSizeSuffix = "${varConnector}pageSize=${pageSize}"

        val jsonVariantArray = Json.parseToJsonElement(get("${urlVariant}${varPageSizeSuffix}", timeout = 0.0).text)
            .jsonObject["result"]
            ?.jsonObject?.get("data")
            ?.jsonArray

        val metadata = Json.parseToJsonElement(get("${urlVariant}${varPageSizeSuffix}", timeout = 0.0).text)
            .jsonObject["metadata"]
            ?.jsonObject
            ?.get("pagination")

        val totalCount = metadata?.jsonObject?.get("totalCount").toString().toInt()
        val totalPages = metadata?.jsonObject?.get("totalPages").toString().toInt()


        val hapIDList = mutableListOf<IntArray>()

        for (page in 0 until totalPages) {
//            if (page % 10 == 0) println("On page: $page")
            val currentURL = "$urlVariant$varPageSizeSuffix&page=$page"
            val currentJsonArray = Json.parseToJsonElement(get(currentURL, timeout = 0.0).text)
                .jsonObject["result"]
                ?.jsonObject?.get("data")
                ?.jsonArray

            currentJsonArray?.forEach {
                val curr = it.jsonObject["alternateBases"]?.jsonArray?.map { j ->
                    j.toString().removeSurrounding("\"").toInt()
                }!!.toIntArray()
                hapIDList.add(curr)
            }
        }

        return hapIDList
    }


    /**
     * Read FeatureTable object from BrAPI `variantTables` endpoint
     * @author Brandon Monier
     * @param urlSample sample URL endpoint
     * @param urlVariant reference range URL endpoint
     * @param urlTable table portion of FeatureTable object as URL endpoint
     * Created 2021-09-30
     */
    public fun readFeatureTableFromBrapi(urlSample: String, urlVariant: String, urlTable: String, pageSize: Int, varPageSize: Int): FeatureTable {
        val taxaBuilder = TaxaListBuilder()
        val refRangeBuilder = GenomicFeatureList.Builder()

        val pattern = Regex("/table$")
        val connector = if (pattern.containsMatchIn(urlTable)) "?" else "&"
        val pageSizeSuffix = "${connector}pageSize=${pageSize}"

        val varPattern = Regex("/variants$")
        val varConnector = if (varPattern.containsMatchIn(urlVariant)) "?" else "&"
        val varPageSizeSuffix = "${varConnector}pageSize=${varPageSize}"

        val jsonSample = getJsonElement(urlSample, "data")
        val jsonRefRange = getJsonElement("$urlVariant$varPageSizeSuffix", "data")
        val jsonTable = Json.parseToJsonElement(get("${urlTable}${pageSizeSuffix}", timeout = 0.0).text)

        jsonSample?.jsonArray?.forEach { taxaBuilder.add(Taxon.instance(it.jsonObject["sampleName"].toString())) }

        val varMetadata = Json.parseToJsonElement(get("${urlVariant}${varPageSizeSuffix}", timeout = 0.0).text)
            .jsonObject["metadata"]
            ?.jsonObject
            ?.get("pagination")
        val varTotalPages = varMetadata?.jsonObject?.get("totalPages").toString().toInt()

        for (page in 0 until varTotalPages) {
//            if (page % 10 == 0) println("On page: $page")
            val currentURL = "$urlVariant$varPageSizeSuffix&page=$page"
            val currentJsonArray = Json.parseToJsonElement(get(currentURL, timeout = 0.0).text)
                .jsonObject["result"]
                ?.jsonObject?.get("data")
                ?.jsonArray
            currentJsonArray?.forEach { refRangeBuilder.add(genomicFeature(it.jsonObject)) }
        }
        val metadata = jsonTable.jsonObject["metadata"]?.jsonObject?.get("pagination")
        val totalTaxa = jsonTable.jsonObject["result"]?.jsonObject?.get("genotypes")?.jsonArray?.get(0)?.jsonArray?.size

        val totalCount = metadata?.jsonObject?.get("totalCount").toString().toInt()
        val totalPages = metadata?.jsonObject?.get("totalPages").toString().toInt()

        val resultTaxa = taxaBuilder.build()
        val resultRefRange = refRangeBuilder.build()
        val ftBuilder = FeatureTableBuilder(resultTaxa, resultRefRange)

        val hapIDList = getUniqueHapIdsFromVariants(urlVariant, 1000)

        var pageSizeInc = 0

        for (page in 0 until totalPages) {
//            if (page % 10 == 0) println("On page: $page")
            val currentURL = "$urlTable$pageSizeSuffix&page=$page"

            val currentJsonArray = Json.parseToJsonElement(get(currentURL, timeout = 0.0).text)
                .jsonObject["result"]
                ?.jsonObject?.get("genotypes")
                ?.jsonArray

            for (i in 0 until currentJsonArray!!.size) {
                for (j in 0 until totalTaxa!!) {
                    val taxon = jsonSample?.jsonArray?.get(j)?.jsonObject?.get("sampleName")
                    val id = jsonRefRange?.jsonArray?.get(i)?.jsonObject?.let { getStringAttribute(it, "variantDbId") }

                    val cell = currentJsonArray[i].jsonArray[j].toString()
                        .removeSurrounding("\"")
                        .split("/")[0]

                    val genotypeValues = if (cell == ".") -1 else hapIDList[i + pageSizeInc][cell.toInt() - 1]

                    ftBuilder.set(taxon.toString(), id.toString(), listOf(genotypeValues.toString(), genotypeValues.toString()))
                }
            }
            pageSizeInc += pageSize
        }

        return ftBuilder.build()
    }

    public fun getVariantTable(url: String): Array<IntArray> {
        val jsonTable = getJsonElement(url, "genotype")

        val nRow = jsonTable?.jsonArray?.size
        val nCol = jsonTable?.jsonArray?.get(0)?.jsonArray?.size
        val hapArray = Array(nCol!!) { IntArray(nRow!!) }

        for (i in 0 until nRow!!) {
            for (j in 0 until nCol!!) {
                hapArray[i][j] = jsonTable.jsonArray[i].jsonArray[j].toString().split("/")[0].toInt()
            }
        }

        return hapArray
    }

    /**
     * Get JSON field from BrAPI endpoints
     */
    private fun getJsonElement(url: String, endPoint: String): JsonElement? {
        return Json.parseToJsonElement(get(url, timeout = 0.0).text).jsonObject["result"]?.jsonObject?.get(endPoint)
    }

    /**
     * Get GenomicFeature object from JSON element
     */
    private fun genomicFeature(json: JsonElement): GenomicFeature {
        val chromosome = Chromosome.instance(getStringAttribute(json, "referenceName"))
        val startPos = json.jsonObject["start"].toString().toInt()
        val endPos = json.jsonObject["end"].toString().toInt()
        val name = getStringFirstArrayElement(json, "variantNames")
        val id = getStringAttribute(json, "variantDbId")
        return GenomicFeature(chromosome, startPos, chromosome, endPos, name, id)
    }

    private fun getStringAttribute(element: JsonElement, attribute: String): String {
        return element.jsonObject[attribute].toString().removePrefix("\"").removeSuffix("\"").trim()
    }

    private fun getStringFirstArrayElement(element: JsonElement, attribute: String): String {
        return element.jsonObject[attribute]?.jsonArray?.get(0).toString().removePrefix("\"").removeSuffix("\"").trim()
    }

    public fun generateKinship(myFT: FeatureTable): DistanceMatrix {
        return KinshipPlugin().run(myFT)
    }
}

fun main() {
    val urlTable = "http://localhost:8080/brapi/v2/variantTables/GATK_PIPELINE/table"
    val urlVariant = "http://localhost:8080/brapi/v2/variantTables/GATK_PIPELINE/variants"
    val urlSample = "http://localhost:8080/brapi/v2/variantTables/GATK_PIPELINE/samples"

    val start = System.currentTimeMillis()

    val myFT = RTasselHelpersTassel6().readFeatureTableFromBrapi(urlSample, urlVariant, urlTable, 1000, 1000)
}