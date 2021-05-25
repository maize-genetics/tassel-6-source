package net.maizegenetics.analysis.brapi

import khttp.get
import kotlinx.serialization.json.*

fun metadata(jsonObject: JsonObject): Map<String, String> {

    val result = mutableMapOf<String, String>()
    metadata(null, jsonObject["metadata"], result)
    return result

}

private fun metadata(path: String?, json: JsonElement?, result: MutableMap<String, String>) {

    when (json) {
        null -> return
        is JsonNull -> return
        is JsonObject -> {
            json.entries.forEach { entry ->
                val newPath = path?.let { "$it.${entry.key}" } ?: entry.key
                metadata(newPath, entry.value, result)
            }
        }
        else -> {
            path?.let { result[path] = json.toString() }
        }
    }

}

fun main() {
    val url = "http://cbsudc01.biohpc.cornell.edu/brapi/v2/variantsets/Ames_MergedReadMapping_AllLines_Haploid/calls"
    val json = Json.parseToJsonElement(get(url, timeout = 0.0).text)
    val metadata = metadata(json.jsonObject)
    println(metadata)
}