package net.maizegenetics.analysis.brapi

import kotlinx.serialization.json.JsonObject
import kotlinx.serialization.json.jsonObject

fun metadata(jsonObject: JsonObject): Map<String, String>? {

    return jsonObject["metadata"]?.jsonObject?.entries
            ?.map { Pair(it.key, it.value.toString()) }
            ?.toMap()

}