@file:JvmName("ChangeLog")

/*
 *  ChangeLog
 *
 *  Created on Oct 8, 2018
 */
package net.maizegenetics.util

import net.maizegenetics.tassel.TASSELMainApp
import java.io.BufferedReader
import java.util.*
import java.util.concurrent.TimeUnit

/**
 * This outputs the change logs messages in HTML for tassel-5-source since this last
 * "New Build Date" message. The results can be added to Tassel5ChangeHistory.html on
 * the website.
 *
 * @author Terry Casstevens
 */
object ChangeLog {

    private val myLogs = LinkedHashSet<String>()

    @JvmStatic
    fun main(args: Array<String>) {

        setupDebugLogging()

        // <h3>(V5.2.14) August 27, 2015</h3>
        // <ul>
        try {
            getLogs().use { reader ->
                var line: String? = reader?.readLine()
                val builder = StringBuilder()
                while (line != null) {
                    line = line.trim { it <= ' ' }
                    if (line.length == 0) {
                        if (builder.length != 0) {
                            builder.append("</li>")
                            if (builder.toString().toLowerCase().equals("<li>new build date</li>") && !myLogs.isEmpty()) {
                                break
                            }
                            addMessage(builder.toString())
                            builder.delete(0, builder.length)
                        }
                    } else {
                        if (builder.length != 0) {
                            builder.append(" ")
                        } else {
                            builder.append("<li>")
                        }
                        builder.append(line)
                    }
                    line = reader?.readLine()
                }
                if (builder.length != 0) {
                    builder.append("</li>")
                    addMessage(builder.toString())
                }
            }
        } catch (e: Exception) {
            e.printStackTrace()
        }

        print("<h3>(V")
        print(TASSELMainApp.version)
        print(") ")
        print(TASSELMainApp.versionDate)
        println("</h3>")
        println("<ul>")
        val itr = myLogs.iterator()
        while (itr.hasNext()) {
            println(itr.next())
        }
        println("</ul>")
    }

    fun getLogs(): BufferedReader? {
        try {
            val proc = ProcessBuilder("/bin/sh", "-c", """git log --since='June 1, 2018' --full-history | grep -v -e '^commit ' -e '^Author:' -e '^Date:'""")
                    .redirectOutput(ProcessBuilder.Redirect.PIPE)
                    .redirectError(ProcessBuilder.Redirect.PIPE)
                    .start()

            proc.waitFor(5, TimeUnit.MINUTES)
            return proc.inputStream.bufferedReader()
        } catch (e: Exception) {
            e.printStackTrace()
            return null
        }
    }

    private fun addMessage(str: String) {
        val lowerCase = str.toLowerCase()
        if (!str.contains("Merge branch")
                && !str.contains("Merge: ")
                && !str.contains("Merged in ")
                && !str.contains("Conflicts: ")
                && !str.contains("Approved-by: ")
                && !str.contains("Merge remote-tracking branch")
                && !lowerCase.contains("removed unused import")
                && !lowerCase.contains("remove unused import")
                && !lowerCase.contains("formatting only")
                && !lowerCase.contains("new build date")) {
            myLogs.add(str)
        }
    }

}
