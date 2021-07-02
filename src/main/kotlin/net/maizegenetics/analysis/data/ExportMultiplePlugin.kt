/*
 * ExportMultiplePlugin.java
 *
 * Created on May 3, 2021
 *
 */
package net.maizegenetics.analysis.data

import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.util.Utils
import org.apache.logging.log4j.LogManager
import java.io.File

/**
 * @author Terry Casstevens
 */
class ExportMultiplePlugin : AbstractPlugin(false) {

    private val logger = LogManager.getLogger(ExportMultiplePlugin::class.java)

    private val myExportPlugin = ExportPlugin()

    private var myFileTypes: Array<TasselFileType>? = null

    var saveFiles: Array<String>? = null

    override fun performFunction(input: DataSet): DataSet? {

        val data = input.dataSet
        var numSaveFiles = 0
        if (saveFiles != null) {
            numSaveFiles = saveFiles!!.size
        }
        check(!(numSaveFiles != 0 && numSaveFiles != 1 && numSaveFiles != data.size)) { "ExportMultiplePlugin: performFunction: number of save files should be either 0, 1 or number of input data sets." }
        if (myFileTypes != null && myFileTypes!!.size != 0) {
            check(!(myFileTypes!!.size != 1 && myFileTypes!!.size != data.size)) { "ExportMultiplePlugin: performFunction: number of files types should be either 0, 1 or number of input data sets." }
        }
        var i = 0
        val n = data.size
        while (i < n) {
            val datum = data[i] as Datum
            val current = DataSet(datum, input.creator)
            if (numSaveFiles == 0) {
                myExportPlugin.saveFile = datum.name
            } else if (numSaveFiles == 1) {
                var temp: String
                if (data.size == 1) {
                    temp = saveFiles!![0]
                } else {
                    val filename = Utils.getFilename(saveFiles!![0])
                    temp = filename.replaceFirst("\\.".toRegex(), (i + 1).toString() + ".")
                    if (temp.length == filename.length) {
                        temp = filename + (i + 1)
                    }
                    val directory = Utils.getDirectory(saveFiles!![0])
                    if (directory != ".") {
                        val dir = File(directory)
                        if (!dir.exists()) {
                            dir.mkdirs()
                        }
                        temp = "$directory/$temp"
                    }
                }
                myExportPlugin.saveFile = temp
            } else {
                myExportPlugin.saveFile = saveFiles!![i]
            }
            if (myFileTypes == null || myFileTypes!!.size == 0) {
                myExportPlugin.fileType = TasselFileType.Unknown
            } else if (myFileTypes!!.size == 1) {
                myExportPlugin.fileType = myFileTypes!![0]
            } else {
                myExportPlugin.fileType = myFileTypes!![i]
            }
            myExportPlugin.performFunction(current)
            i++
        }

        fireProgress(100)

        return null

    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    override fun getButtonName(): String {
        return "Export Multiple"
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    override fun getToolTipText(): String {
        return "Export multiple data sets to files."
    }

    fun fileTypes(types: Array<TasselFileType>?) {
        myFileTypes = types
    }

    fun fileType(type: TasselFileType) {
        myFileTypes = arrayOf(type)
    }

    fun includeAnnotations(include: Boolean) {
        myExportPlugin.includeTaxaAnnotations = include
    }

    fun includeDepth(include: Boolean) {
        myExportPlugin.keepDepth = include
    }

}