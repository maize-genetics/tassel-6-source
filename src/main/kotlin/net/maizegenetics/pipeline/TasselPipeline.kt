/*
* TasselPipeline.java
*
* Created on April 29, 2021
*
*/
package net.maizegenetics.pipeline

import net.maizegenetics.analysis.data.*
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType
import net.maizegenetics.gui.basicLoggingInfo
import net.maizegenetics.plugindef.*
import net.maizegenetics.prefs.TasselPrefs
import net.maizegenetics.tassel.TASSELMainApp
import net.maizegenetics.util.ExceptionUtils
import net.maizegenetics.util.LoggingUtils
import net.maizegenetics.util.Utils
import org.apache.log4j.Level
import org.apache.log4j.Logger
import java.io.File
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter
import java.util.*
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import kotlin.collections.HashMap
import kotlin.collections.LinkedHashMap
import kotlin.system.exitProcess

/**
 * @author Terry Casstevens
 */
class TasselPipeline(args: Array<String>, interactive: Boolean = false, name: String? = null) : PluginListener {

    private val myLogger = Logger.getLogger(TasselPipeline::class.java)

    enum class FLAGS {
        t, s, k, q, h, r, plink, fasta, table, vcf, readSerialAlignment, importGuess, projection, convertTOPMtoHDF5, retainRareAlleles, union, intersect, separate, homozygous, synonymizer, mergeGenotypeTables, mergeAlignmentsSameSites, excludeLastTrait, mlm, glm, td_csv, td_tab, td_gui, diversity, ld, ldd, ck, tree, gs, distanceMatrix, distMatrixRanges, genotypeSummary, export, filterAlign, numericalGenoTransform, includeTaxa, includeTaxaInFile, excludeTaxa, excludeTaxaInFile, includeSiteNames, includeSiteNamesInFile, excludeSiteNames, excludeSiteNamesInFile, subsetSites, subsetTaxa, newCoordinates, archaeopteryx, filterTaxaNames, maxThreads, mhd, pca, printGenoSummary, printMemoryUsage;

        override fun toString(): String {
            return "-" + super.toString()
        }
    }

    private val myForks: MutableMap<String, List<Plugin>> = LinkedHashMap()
    private var myCurrentFork: String? = null
    private var myCurrentPipe: MutableList<Plugin>? = null
    private var myFirstPlugin: Plugin? = null
    private val myThreads: MutableList<ThreadedPluginListener> = ArrayList()
    private val myProgressValues: MutableMap<Plugin, Int> = HashMap()
    private val myDeprecatedWarning = StringBuilder()
    private val myIsInteractive: Boolean = interactive
    private val myIsThreaded: Boolean
    private var myDescriptions: Array<String>? = null
    private var myCurrentDescriptionIndex = 0

    fun parseArgs(input: Array<String>) {
        var args = input
        if (args.size >= 1 && args[0].equals("-configFile", ignoreCase = true)) {
            require(args.size >= 2) { "TasselPipeline: parseArgs: a filename must follow -configFile flag." }
            val xmlFilename = args[1].trim { it <= ' ' }
            val tempArgsDesc: Array<Array<String>> = TasselPipelineXMLUtil.readXMLAsArgs(xmlFilename)
            args = tempArgsDesc[0]
            myDescriptions = tempArgsDesc[1]
        } else if (args.isNotEmpty() && args[0].equals("-configResourceFile", ignoreCase = true)) {
            require(args.size >= 2) { "TasselPipeline: parseArgs: a filename must follow -configResourceFile flag." }
            val xmlFilename = args[1].trim { it <= ' ' }
            val tempArgsDesc: Array<Array<String>> = TasselPipelineXMLUtil.readXMLAsArgsFromResource(xmlFilename)
            args = tempArgsDesc[0]
            myDescriptions = tempArgsDesc[1]
        } else {
            args = addForkFlagsIfNeeded(args)
        }
        val argsStr = StringBuilder()
        argsStr.append("[")
        var print = true
        var first = true
        for (current in args!!) {
            if (first) {
                first = false
            } else {
                argsStr.append(", ")
            }
            if (print) {
                argsStr.append(current)
            } else {
                argsStr.append("?????")
                print = true
            }
            if (current!!.toUpperCase().contains("PASSWORD")) {
                print = false
            }
        }
        argsStr.append("]")
        myLogger.info("Tassel Pipeline Arguments: $argsStr")
        var index = 0
        while (index < args.size) {
            myCurrentDescriptionIndex = index
            try {
                var current = args[index++]
                val emDash = "\u2014"
                current = current!!.replaceFirst(emDash.toRegex(), "-")
                require(current.startsWith("-")) { "TasselPipeline: parseArgs: expecting argument beginning with dash: $current" }
                if (current.startsWith("-runfork")) {
                    val key = current.replaceFirst("-runfork".toRegex(), "-fork")
                    val specifiedPipe: List<Plugin>? = myForks[key]
                    requireNotNull(specifiedPipe) { "TasselPipeline: parseArgs: unknown fork: $current" }
                    require(!specifiedPipe.isEmpty()) { "TasselPipeline: parseArgs: empty fork: $current" }
                    require(!(specifiedPipe[0] is AbstractPlugin && !(specifiedPipe[0] as AbstractPlugin).inputs.isEmpty())) { "TasselPipeline: parseArgs: this fork does not need to be explicitly run: it is receiving input from another plugin: $current" }
                    val event = PluginEvent(DataSet(null as Datum?, null))
                    val thread = ThreadedPluginListener(specifiedPipe[0], event)
                    myThreads.add(thread)
                } else if (current.startsWith("-fork")) {
                    if (myCurrentPipe != null && !myCurrentPipe!!.isEmpty()) {
                        myCurrentPipe!![myCurrentPipe!!.size - 1]!!.setThreaded(myIsThreaded)
                    }
                    myCurrentFork = current
                    myCurrentPipe = ArrayList()
                    myForks[myCurrentFork!!] = myCurrentPipe!!
                } else if (current.startsWith("-input")) {
                    val key = current.replaceFirst("-input".toRegex(), "-fork")
                    val specifiedPipe: List<Plugin>? = myForks[key]
                    requireNotNull(specifiedPipe) { "TasselPipeline: parseArgs: unknown input: $current" }
                    var lastCurrentPipe: Plugin? = null
                    lastCurrentPipe = try {
                        myCurrentPipe!![myCurrentPipe!!.size - 1]
                    } catch (e: Exception) {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: -input must come after plugin in current fork.")
                    }
                    val endSpecifiedPipe = specifiedPipe[specifiedPipe.size - 1]
                    lastCurrentPipe!!.receiveInput(endSpecifiedPipe)
                } else if (current.startsWith("-inputOnce")) {
                    val key = current.replaceFirst("-input".toRegex(), "-fork")
                    val specifiedPipe: List<Plugin>? = myForks[key]
                    requireNotNull(specifiedPipe) { "TasselPipeline: parseArgs: unknown input: $current" }
                    val combinePlugin = try {
                        myCurrentPipe!![myCurrentPipe!!.size - 1] as CombineDataSetsPlugin
                    } catch (e: Exception) {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: -inputOnce must follow -combine flag.")
                    }
                    val endSpecifiedPipe = specifiedPipe[specifiedPipe.size - 1]
                    combinePlugin.receiveDataSetOnceFrom(endSpecifiedPipe)
                } else if (current.startsWith("-combine")) {
                    current = current.replaceFirst("-combine".toRegex(), "-fork")
                    if (myCurrentPipe != null && myCurrentPipe!!.isNotEmpty()) {
                        myCurrentPipe!![myCurrentPipe!!.size - 1]!!.setThreaded(myIsThreaded)
                    }
                    myCurrentFork = current
                    myCurrentPipe = ArrayList()
                    myForks[myCurrentFork!!] = myCurrentPipe!!
                    integratePlugin(CombineDataSetsPlugin(), false)
                } else if (current.equals("-maxThreads", ignoreCase = true)) {
                    val str = args[index++]!!.trim { it <= ' ' }
                    val numThreads = try {
                        str.toInt()
                    } catch (e: Exception) {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: Problem with number of max threads: $str")
                    }
                    TasselPrefs.putMaxThreads(numThreads)
                } else if (current.equals("-t", ignoreCase = true)) {
                    val traitFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(traitFile, TasselFileType.Phenotype)
                } else if (current.equals("-s", ignoreCase = true)) {
                    val inputFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(inputFile, TasselFileType.Sequence)
                } else if (current.equals("-k", ignoreCase = true)) {
                    val kinshipFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(kinshipFile, TasselFileType.SqrMatrix)
                } else if (current.equals("-q", ignoreCase = true)) {
                    val populationFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(populationFile, TasselFileType.Phenotype)
                } else if (current.equals("-h", ignoreCase = true)) {
                    val hapFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(hapFile, TasselFileType.Hapmap)
                } else if (current.equals("-r", ignoreCase = true)) {
                    val phenotypeFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(phenotypeFile, TasselFileType.Phenotype)
                } else if (current.equals("-fasta", ignoreCase = true)) {
                    val fastaFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(fastaFile, TasselFileType.Fasta)
                } else if (current.equals("-table", ignoreCase = true)) {
                    val tableFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(tableFile, TasselFileType.Table)
                } else if (current.equals("-vcf", ignoreCase = true)) {
                    val vcfFile = args[index++]!!.trim { it <= ' ' }
                    loadFile(vcfFile, TasselFileType.VCF)
                } else if (current.equals("-readSerialAlignment", ignoreCase = true)) {
                    val file = args[index++]!!.trim { it <= ' ' }
                    loadFile(file, TasselFileType.Serial)
                } else if (current.equals("-importGuess", ignoreCase = true)) {
                    val file = args[index++]!!.trim { it <= ' ' }
                    loadFile(file, TasselFileType.Unknown)
                } else if (current.equals("-sortPositions", ignoreCase = true)) {
                    val plugin = findLastPluginFromCurrentPipe(arrayOf(FileLoadPlugin::class.java)) as FileLoadPlugin?
                    if (plugin != null) {
                        plugin.sortPositions(true)
                    } else {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: No FileLoadPlugin step defined: $current")
                    }
                } else if (current.equals("-noDepth", ignoreCase = true)) {
                    val plugin = findLastPluginFromCurrentPipe(arrayOf(FileLoadPlugin::class.java)) as FileLoadPlugin?
                    if (plugin != null) {
                        plugin.keepDepth(false)
                    } else {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: No FileLoadPlugin step defined: $current")
                    }
                } else if (current.equals("-printGenoSummary", ignoreCase = true)) {
                    val plugin = GenotypeSummaryPlugin(myIsInteractive)
                    integratePlugin(plugin, true)
                    plugin.overview = false
                    plugin.siteSummary = false
                    plugin.taxaSummary = false
                } else if (current.equals("-printMemoryUsage", ignoreCase = true)) {
                    val plugin = MemoryUsagePlugin(myIsInteractive)
                    integratePlugin(plugin, false)
                } else if (current.equals("-retainRareAlleles", ignoreCase = true)) {
                    val temp = args[index++]!!.trim { it <= ' ' }
                    var retain = true
                    retain = if (temp.equals("false", ignoreCase = true)) {
                        false
                    } else if (temp.equals("true", ignoreCase = true)) {
                        true
                    } else {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: -retainRareAlleles parameter must be true or false.")
                    }
                    TasselPrefs.putAlignmentRetainRareAlleles(retain)
                } else if (current.equals("-genotypeSummary", ignoreCase = true)) {
                    val plugin = GenotypeSummaryPlugin(myIsInteractive)
                    integratePlugin(plugin, true)
                    val temp = args[index++]!!.trim { it <= ' ' }
                    val types = temp.split(",".toRegex()).toTypedArray()
                    plugin.overview = false
                    plugin.siteSummary = false
                    plugin.taxaSummary = false
                    for (i in types.indices) {
                        if (types[i].equals("overall", ignoreCase = true)) {
                            plugin.overview = true
                        } else if (types[i].equals("site", ignoreCase = true)) {
                            plugin.siteSummary = true
                        } else if (types[i].equals("taxa", ignoreCase = true)) {
                            plugin.taxaSummary = true
                        } else if (types[i].equals("all", ignoreCase = true)) {
                            plugin.overview = true
                            plugin.siteSummary = true
                            plugin.taxaSummary = true
                        } else {
                            throw IllegalArgumentException("TasselPipeline: parseArgs: -genotypeSummary illegal types: $temp")
                        }
                    }
                } else if (current.equals("-export", ignoreCase = true)) {
                    val plugin = ExportMultiplePlugin()
                    val temp = args[index]!!.trim { it <= ' ' }
                    if (!temp.startsWith("-")) {
                        val filenames = temp.split(",".toRegex()).toTypedArray()
                        plugin.saveFiles = filenames
                        index++
                    }
                    integratePlugin(plugin, false)
                } else if (current.equals("-exportType", ignoreCase = true)) {
                    val plugin: ExportMultiplePlugin? = findLastPluginFromCurrentPipe(arrayOf(ExportMultiplePlugin::class.java)) as ExportMultiplePlugin?
                            ?: throw IllegalArgumentException("TasselPipeline: parseArgs: No Export step defined: $current")
                    val type = args[index++]!!.trim { it <= ' ' }
                    try {
                        plugin!!.fileType(TasselFileType.valueOf(type))
                    } catch (e: Exception) {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: -exportType: Unknown type: " + type + "  Should be: " + TasselFileType.values().contentToString())
                    }
                } else if (current.equals("-exportIncludeAnno", ignoreCase = true)) {
                    val plugin: ExportMultiplePlugin? = findLastPluginFromCurrentPipe(arrayOf(ExportMultiplePlugin::class.java)) as ExportMultiplePlugin?
                            ?: throw IllegalArgumentException("TasselPipeline: parseArgs: No Export step defined: $current")
                    val temp = args[index++]!!.trim { it <= ' ' }
                    val value = if (temp.equals("false", ignoreCase = true)) {
                        false
                    } else if (temp.equals("true", ignoreCase = true)) {
                        true
                    } else {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: -exportIncludeAnno must be true or false: $temp")
                    }
                    plugin!!.includeAnnotations(value)
                } else if (current.equals("-exportIncludeDepth", ignoreCase = true)) {
                    val plugin: ExportMultiplePlugin? = findLastPluginFromCurrentPipe(arrayOf(ExportMultiplePlugin::class.java)) as ExportMultiplePlugin?
                            ?: throw IllegalArgumentException("TasselPipeline: parseArgs: No Export step defined: $current")
                    val temp = args[index++]!!.trim { it <= ' ' }
                    val value = if (temp.equals("false", ignoreCase = true)) {
                        false
                    } else if (temp.equals("true", ignoreCase = true)) {
                        true
                    } else {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: -exportIncludeDepth must be true or false: $temp")
                    }
                    plugin!!.includeDepth(value)
                } else {
                    try {
                        var plugin: Plugin? = null
                        val possibleClassName = current.substring(1)
                        val matches = Utils.getFullyQualifiedClassNames(possibleClassName)
                        for (match in matches) {
                            plugin = Plugin.getPluginInstance(match, myIsInteractive)
                            if (plugin != null) {
                                break
                            }
                        }
                        if (plugin == null) {
                            plugin = Plugin.getPluginInstance(possibleClassName, myIsInteractive)
                        }
                        if (plugin != null) {
                            integratePlugin(plugin, true)
                            val pluginArgs: MutableList<String> = ArrayList()
                            var temp = args[index++]!!.trim { it <= ' ' }
                            while (!temp.equals("-endPlugin", ignoreCase = true)) {
                                if (temp.startsWith("-runfork")) {
                                    index--
                                    break
                                }
                                pluginArgs.add(temp)
                                temp = args[index++]!!.trim { it <= ' ' }
                            }
                            try {
                                plugin.setParameters(pluginArgs.toTypedArray())
                            } catch (e: Exception) {
                                e.printStackTrace()
                                // Self-describing Plugin Should already output Usage and any other error information.
                                ExceptionUtils.logExceptionCauses(e, myLogger, Level.ERROR)
                                exitProcess(1)
                            }
                        } else {
                            throw IllegalArgumentException("TasselPipeline: parseArgs: Unknown parameter: $current")
                        }
                    } catch (usoe: UnsupportedOperationException) {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: this plugin is not self-described: $current")
                    } catch (e: Exception) {
                        throw IllegalArgumentException("TasselPipeline: parseArgs: Unknown parameter: $current")
                    }
                }
            } catch (e: Exception) {
                myLogger.error(e.message)
                myLogger.debug(e.message, e)
                exitProcess(1)
            }
        }
        if (myFirstPlugin != null) {
            tracePipeline()
        } else {
            myLogger.warn("parseArgs: no arguments specified.")
        }
    }

    private fun tracePipeline() {
        for (thread in myThreads) {
            val current = thread.pluginListener as Plugin
            (current as AbstractPlugin).trace(0)
        }
    }

    fun loadFile(filename: String, fileType: TasselFileType?): FileLoadPlugin {
        val plugin = FileLoadPlugin(myIsInteractive)
        integratePlugin(plugin, true)
        if (fileType == null) {
            plugin.theFileType = TasselFileType.Unknown
        } else {
            plugin.theFileType = fileType
        }
        plugin.openFiles = arrayOf(filename)
        return plugin
    }

    private fun integratePlugin(plugin: Plugin, displayDataTree: Boolean) {
        println("plugin: ${plugin.buttonName}  displayDataTree: $displayDataTree")
        if (myFirstPlugin == null) {
            myFirstPlugin = plugin
        }
        if (displayDataTree) {
            plugin.addListener(this)
        }
        if (myCurrentPipe == null) {
            myCurrentPipe = ArrayList()
        }
        if (myCurrentPipe!!.isEmpty()) {
            myCurrentPipe!!.add(plugin)
        } else {
            plugin.receiveInput(myCurrentPipe!![myCurrentPipe!!.size - 1])
            myCurrentPipe!!.add(plugin)
        }
        (plugin as AbstractPlugin?)!!.setConfigParameters()
    }

    private fun findLastPluginFromAll(types: Array<Class<*>>): Plugin? {
        if (myCurrentPipe != null && myCurrentPipe!!.size != 0) {
            for (i in myCurrentPipe!!.indices.reversed()) {
                val current = myCurrentPipe!![i]
                if (matchType(types, current)) {
                    return current
                }
            }
        }
        val keys: List<*> = ArrayList<Any?>(myForks.keys)
        for (i in keys.indices.reversed()) {
            val currentPipe = myForks[keys[i]] as List<*>?
            for (j in currentPipe!!.indices.reversed()) {
                val current = currentPipe[j] as Plugin
                if (matchType(types, current)) {
                    return current
                }
            }
        }
        return null
    }

    private fun findLastPluginFromCurrentPipe(types: Array<Class<*>>): Plugin? {
        if (myCurrentPipe != null && myCurrentPipe!!.size != 0) {
            for (i in myCurrentPipe!!.indices.reversed()) {
                val current = myCurrentPipe!![i]
                if (matchType(types, current)) {
                    return current
                }
            }
        }
        return null
    }

    private fun matchType(types: Array<Class<*>>, test: Any?): Boolean {
        for (i in types.indices) {
            if (types[i].isInstance(test)) {
                return true
            }
        }
        return false
    }

    /**
     * Returns Tassel data set after complete.
     *
     * @param event event
     */
    override fun dataSetReturned(event: PluginEvent) {
        val tds: DataSet? = event.source as DataSet
    }

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    override fun progress(event: PluginEvent) {

        val ds: DataSet? = event.source as DataSet
        if (ds != null) {
            val percentage: List<Datum> = ds.getDataOfType(Int::class.javaObjectType)
            val plugin = ds.creator
            var lastValue = myProgressValues[plugin]
            if (lastValue == null) {
                lastValue = 0
            }
            if (percentage.isNotEmpty()) {
                val datum = percentage[0]
                val percent = datum.data as Int
                if (percent >= lastValue) {
                    val time = LocalDateTime.now()
                    val timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"))
                    myLogger.info(ds.creator.javaClass.name + ": time: " + timeStr + ": progress: " + percent + "%")
                    lastValue += 10
                    myProgressValues[plugin] = lastValue
                }
            }
        }

    }

    companion object {

        @JvmStatic
        fun main(args: Array<String>) {
            val emDash = "\u2014"
            for (i in args.indices) {
                args[i] = args[i].replaceFirst(emDash.toRegex(), "-")
            }
            TasselPrefs.setPersistPreferences(false)
            LoggingUtils.setupLogging()
            if (args.size >= 2 && args[0].equals("-createXML", ignoreCase = true)) {
                val xmlFilename = args[1].trim { it <= ' ' }
                val temp = addForkFlagsIfNeeded(args.copyOfRange(2, args.size))
                TasselPipelineXMLUtil.writeArgsAsXML(xmlFilename, temp)
                return
            }
            if (args.size >= 2 && args[0].equals("-translateXML", ignoreCase = true)) {
                val xmlFilename = args[1].trim { it <= ' ' }
                val result: Array<Array<String>> = TasselPipelineXMLUtil.readXMLAsArgs(xmlFilename)
                for (element in result[0]) {
                    print(element)
                    print(" ")
                }
                println("")
                return
            }
            var currentArgs = args
            var notDone = true
            while (notDone) {
                if (currentArgs.isNotEmpty() && (currentArgs[0].equals("-debug", ignoreCase = true) || currentArgs[0].equals("-log", ignoreCase = true))) {
                    var filename: String? = null
                    if (currentArgs.size >= 2) {
                        filename = currentArgs[1].trim { it <= ' ' }
                    }
                    if (filename != null && !filename.startsWith("-")) {
                        try {
                            if (currentArgs[0].equals("-debug", ignoreCase = true)) {
                                LoggingUtils.setupDebugLogfile(filename)
                            } else {
                                LoggingUtils.setupLogfile(filename)
                            }
                        } catch (e: Exception) {
                            println("Problem with file: " + filename + "\n" + e.message)
                        }
                        currentArgs = currentArgs.copyOfRange(2, currentArgs.size)
                    } else {
                        if (currentArgs[0].equals("-debug", ignoreCase = true)) {
                            LoggingUtils.setupDebugLogging()
                        } else {
                            LoggingUtils.setupLogging()
                        }
                        currentArgs = currentArgs.copyOfRange(1, currentArgs.size)
                    }
                } else if (currentArgs.size >= 2 && currentArgs[0].equals("-configParameters", ignoreCase = true)) {
                    val filename = currentArgs[1].trim { it <= ' ' }
                    if (!File(filename).isFile) {
                        throw IllegalArgumentException("TasselPipeline: main: -configParameters file: $filename doesn't exist or isn't a file.")
                    }
                    ParameterCache.load(filename)
                    currentArgs = currentArgs.copyOfRange(2, currentArgs.size)
                } else {
                    notDone = false
                }
            }
            TasselPipeline(currentArgs)
        }

        fun addForkFlagsIfNeeded(args: Array<String>): Array<String> {

            if (args.isEmpty()) {
                return args
            }
            for (a in args) {
                if (a.toLowerCase().startsWith("-fork") || a.toLowerCase().startsWith("-runfork")) {
                    // If forks included, return arguments unchanged
                    return args
                }
            }

            // If no arguments have "-fork" or "-runfork", add them
            val newArgs = mutableListOf<String>()
            newArgs.add("-fork1")
            newArgs.addAll(args)
            newArgs.add("-runfork1")
            return newArgs.toTypedArray()

        }

    }

    /**
     * Creates a new instance of TasselPipeline
     */
    init {

        myIsThreaded = !myIsInteractive
        if (args.size == 1 && args[0].equals("-versionComment", ignoreCase = true)) {
            System.out.println("Version " + TASSELMainApp.version + " on " + TASSELMainApp.versionDate)
        }
        if (args.size == 1 && args[0].equals("-versionTag", ignoreCase = true)) {
            System.out.println("V" + TASSELMainApp.version)
        }

        basicLoggingInfo()

        val pool: ExecutorService?
        if (myIsInteractive) {
            pool = null
        } else {
            var numThreads = Runtime.getRuntime().availableProcessors() / 2
            numThreads = Math.max(2, numThreads)
            pool = Executors.newFixedThreadPool(numThreads)
        }
        try {
            parseArgs(args)
            for ((_, current) in myForks) {
                if (current != null && !current.isEmpty()) {
                    val first = current[0]
                    if (first is AbstractPlugin && first.inputs.isEmpty()) {
                        var alreadyRun = false
                        for (currentListener in myThreads) {
                            if (currentListener.pluginListener === first) {
                                alreadyRun = true
                                break
                            }
                        }
                        if (!alreadyRun) {
                            val event = PluginEvent(DataSet(null as Datum?, null))
                            val thread = ThreadedPluginListener(first, event)
                            myThreads.add(thread)
                        }
                    }
                }
            }

            val futures: MutableList<Future<*>> = ArrayList()
            myThreads.stream().forEach({ current: ThreadedPluginListener? -> futures.add(pool!!.submit(current)) })
            for (future in futures) {
                future.get()
            }

            if (myDeprecatedWarning.length != 0) {
                myLogger.warn(myDeprecatedWarning.toString())
            }
        } catch (e: Exception) {
            myLogger.error(e.message, e)
            if (myIsInteractive) {
                throw IllegalStateException("TasselPipeline: init: " + e.message)
            } else {
                System.exit(1)
            }
        } finally {
            pool?.shutdown()
        }
    }
}