/*
 * ExportPlugin.java
 *
 * Created on May 1, 2021
 *
 */
package net.maizegenetics.analysis.data

import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.dna.map.PositionList
import net.maizegenetics.dna.map.PositionListTableReport
import net.maizegenetics.dna.snp.ExportUtils
import net.maizegenetics.dna.snp.FilterList
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.io.FilterJSONUtils
import net.maizegenetics.dna.snp.io.JSONUtils
import net.maizegenetics.dna.snp.io.SiteScoresIO
import net.maizegenetics.phenotype.Phenotype
import net.maizegenetics.phenotype.PhenotypeUtils
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Parameter
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.TaxaListTableReport
import net.maizegenetics.taxa.distance.DistanceMatrix
import net.maizegenetics.taxa.distance.DistanceMatrixUtils
import net.maizegenetics.taxa.distance.WriteDistanceMatrix
import net.maizegenetics.taxa.tree.SimpleTree
import net.maizegenetics.util.TableReport
import net.maizegenetics.util.TableReportUtils
import net.maizegenetics.util.Utils
import org.apache.log4j.Logger
import java.io.File
import java.io.FileWriter
import java.io.PrintWriter
import java.util.*

//import net.maizegenetics.dna.snp.io.FlapjackUtils
//import net.maizegenetics.taxa.tree.NewickUtils


/**
 * @author Terry Casstevens
 */
class ExportPlugin(isInteractive: Boolean = false) : AbstractPlugin(isInteractive) {

    private val myLogger = Logger.getLogger(ExportPlugin::class.java)

    private var mySaveFile = PluginParameter.Builder("saveFile", null, String::class.java)
            .description("Save file as...")
            .outFile()
            .required(true)
            .build()
    var saveFile by Parameter<String>()

    private var myFileType = PluginParameter.Builder("fileType", null, TasselFileType::class.java)
            .description("Export file format (Default format depends on data being exported)")
            .required(true)
            .objectListSingleSelect()
            .range(TasselFileType.values())
            .build()
    var fileType by Parameter<TasselFileType?>()

    private var myKeepDepth = PluginParameter.Builder("keepDepth", true, Boolean::class.java)
            .description("Whether to keep depth if format supports depth.")
            .dependentOnParameter(myFileType, arrayOf(TasselFileType.VCF))
            .build()
    var keepDepth by Parameter<Boolean>()

    private var myIncludeTaxaAnnotations = PluginParameter.Builder("includeTaxaAnnotations", true, Boolean::class.java)
            .description("Whether to include taxa annotations if format supports taxa annotations.")
            .dependentOnParameter(myFileType, arrayOf(TasselFileType.VCF,
                    TasselFileType.Hapmap,
                    TasselFileType.HapmapDiploid,
                    TasselFileType.HapmapLIX))
            .build()
    var includeTaxaAnnotations by Parameter<Boolean>()

    private var myIncludeBranchLengths = PluginParameter.Builder("includeBranchLengths", true, Boolean::class.java)
            .description("Whether to include branch lengths for Newick formatted files.")
            .dependentOnParameter(myFileType, arrayOf(TasselFileType.Newick))
            .build()
    var includeBranchLengths by Parameter<Boolean>()


    override fun preProcessParameters(input: DataSet) {
        require(input.size == 1) { "Please select only one item." }
        val data = input.getData(0).data
        myFileType = when (data) {
            is GenotypeTable -> {
                val genotype = data
                val temp: MutableList<TasselFileType?> = ArrayList()
                if (genotype.hasGenotype()) {
                    temp.add(TasselFileType.HaplotypeVCF)
                    temp.add(TasselFileType.Hapmap)
                    temp.add(TasselFileType.HapmapDiploid)
                    temp.add(TasselFileType.VCF)
                    temp.add(TasselFileType.Plink)
                    temp.add(TasselFileType.Phylip_Seq)
                    temp.add(TasselFileType.Phylip_Inter)
                    temp.add(TasselFileType.Table)
                }
                if (genotype.hasDepth()) {
                    temp.add(TasselFileType.Depth)
                }
                if (genotype.hasReferenceProbablity()) {
                    temp.add(TasselFileType.ReferenceProbability)
                }
                PluginParameter(myFileType, temp)
            }
            is FeatureTable -> {
                PluginParameter(myFileType, listOf(TasselFileType.HaplotypeVCF))
            }
            is Phenotype -> {
                PluginParameter(myFileType,
                        listOf(TasselFileType.Phenotype,
                                TasselFileType.PlinkPhenotype))
            }
            is FilterList -> {
                PluginParameter(myFileType, listOf(TasselFileType.Filter))
            }
            is DistanceMatrix -> {
                PluginParameter(myFileType,
                        listOf(TasselFileType.SqrMatrix,
                                TasselFileType.SqrMatrixBin,
                                TasselFileType.SqrMatrixRaw,
                                TasselFileType.SqrMatrixDARwinDIS))
            }
            is TaxaList -> {
                PluginParameter(myFileType,
                        listOf(TasselFileType.TaxaList,
                                TasselFileType.Table))
            }
            is TaxaListTableReport -> {
                PluginParameter(myFileType,
                        listOf(TasselFileType.TaxaList,
                                TasselFileType.Table))
            }
            is PositionList -> {
                PluginParameter(myFileType,
                        listOf(TasselFileType.PositionList,
                                TasselFileType.Table))
            }
            is PositionListTableReport -> {
                PluginParameter(myFileType,
                        listOf(TasselFileType.PositionList,
                                TasselFileType.Table))
            }
            is TableReport -> {
                PluginParameter(myFileType, listOf(TasselFileType.Table))
            }
            is SimpleTree -> {
                PluginParameter(myFileType,
                        listOf( //FileLoadPlugin.TasselFileType.Newick,
                                TasselFileType.Report))
            }
            else -> {
                throw IllegalStateException("Don't know how to export data type: " + data.javaClass.name)
            }
        }
        if (!isInteractive && (fileType == null || fileType == TasselFileType.Unknown) && myFileType.hasPossibleValues()) {
            fileType = myFileType.possibleValues()[0]
        }
    }

    override fun processData(input: DataSet): DataSet? {
        val filename: String?
        val data = input.getData(0).data
        filename = when (data) {
            is GenotypeTable -> performFunctionForAlignment(data)
            is Phenotype -> performFunctionForPhenotype(data)
            is FilterList -> performFunctionForFilter(data)
            is DistanceMatrix -> performFunctionForDistanceMatrix(data)
            is TaxaList -> performFunctionForTaxaList(data)
            is TaxaListTableReport -> performFunctionForTaxaList(data.taxaList)
            is PositionList -> performFunctionForPositionList(data)
            is PositionListTableReport -> performFunctionForPositionList(data.positionList)
            is TableReport -> performFunctionForTableReport(data)
            is SimpleTree -> performFunctionForSimpleTree(data)
            else -> throw IllegalStateException("Don't know how to export data type: " + data.javaClass.name)
        }
        if (filename != null) {
            myLogger.info("performFunction: wrote dataset: ${input.getData(0).name} to file: $filename")
        }
        return null
    }

    fun write(data: Any, filename: String, type: TasselFileType? = null) {
        saveFile = filename
        fileType = type
        performFunction(DataSet.getDataSet(data))
    }

    private fun performFunctionForDistanceMatrix(input: DistanceMatrix): String {
        return if (fileType == TasselFileType.SqrMatrix) {
            val filename = Utils.addSuffixIfNeeded(saveFile, ".txt", arrayOf(".txt", ".txt.gz"))
            WriteDistanceMatrix.saveDelimitedDistanceMatrix(input, filename)
            filename
        } else if (fileType == TasselFileType.SqrMatrixRaw) {
            val grmFiles = DistanceMatrixUtils.getGRMFilenames(saveFile)
            WriteDistanceMatrix.saveRawMultiBlupMatrix(input, grmFiles[0], grmFiles[3])
            grmFiles[3]
        } else if (fileType == TasselFileType.SqrMatrixBin) {
            val grmFiles = DistanceMatrixUtils.getGRMFilenames(saveFile)
            WriteDistanceMatrix.saveBinMultiBlupMatrix(input, grmFiles[0], grmFiles[1], grmFiles[2])
            grmFiles[1]
        } else if (fileType == TasselFileType.SqrMatrixDARwinDIS) {
            val filename = Utils.addSuffixIfNeeded(saveFile, ".dis")
            WriteDistanceMatrix.saveDARwinMatrix(input, filename)
            filename
        } else {
            throw IllegalArgumentException("ExportPlugin: performFunctionForDistanceMatrix: Unknown file type: $fileType")
        }
    }

    private fun performFunctionForTableReport(input: TableReport): String {
        val theFile = File(Utils.addSuffixIfNeeded(saveFile, ".txt"))
        TableReportUtils.saveDelimitedTableReport(input, "\t", theFile)
        return theFile.absolutePath
    }

    private fun performFunctionForFilter(filter: FilterList): String {
        return FilterJSONUtils.exportFilterToJSON(filter, saveFile)
    }

    private fun performFunctionForPhenotype(input: Phenotype): String {
        val filename = Utils.addSuffixIfNeeded(saveFile, ".txt")
        if (fileType == TasselFileType.Phenotype) {
            PhenotypeUtils.write(input, filename)
        } else if (fileType == TasselFileType.PlinkPhenotype) {
            PhenotypeUtils.writePlink(input, filename)
        }
        return File(filename).absolutePath
    }

    private fun performFunctionForFactorTable(table: FeatureTable): String? {

        var resultFile = saveFile
        when (fileType) {
            TasselFileType.HaplotypeVCF -> TODO()
        }

        return resultFile

    }

    private fun performFunctionForAlignment(inputAlignment: GenotypeTable): String? {
        var resultFile = saveFile
        when (fileType) {
            TasselFileType.ReferenceProbability -> resultFile = SiteScoresIO.writeReferenceProbability(inputAlignment, resultFile)
            TasselFileType.Depth -> resultFile = SiteScoresIO.writeDepth(inputAlignment, resultFile)
            TasselFileType.Hapmap -> resultFile = ExportUtils.writeToHapmap(inputAlignment, false, saveFile, '\t', includeTaxaAnnotations, this)
            TasselFileType.HapmapDiploid -> resultFile = ExportUtils.writeToHapmap(inputAlignment, true, saveFile, '\t', includeTaxaAnnotations, this)
            TasselFileType.Plink -> resultFile = ExportUtils.writeToPlink(inputAlignment, saveFile, '\t')
            //TasselFileType.Flapjack -> resultFile = FlapjackUtils.writeToFlapjack(inputAlignment, saveFile(), '\t')
            TasselFileType.Phylip_Seq -> {
                resultFile = Utils.addSuffixIfNeeded(saveFile, ".phy")
                try {
                    PrintWriter(FileWriter(resultFile)).use { out -> ExportUtils.printSequential(inputAlignment, out) }
                } catch (e: Exception) {
                    myLogger.debug(e.message, e)
                    throw IllegalStateException("ExportPlugin: performFunction: Problem writing file: $resultFile")
                }
            }
            TasselFileType.Phylip_Inter -> {
                resultFile = Utils.addSuffixIfNeeded(saveFile, ".phy")
                try {
                    PrintWriter(FileWriter(resultFile)).use { out -> ExportUtils.printInterleaved(inputAlignment, out) }
                } catch (e: Exception) {
                    myLogger.debug(e.message, e)
                    throw IllegalStateException("ExportPlugin: performFunction: Problem writing file: $resultFile")
                }
            }
            TasselFileType.Table -> resultFile = ExportUtils.saveDelimitedAlignment(inputAlignment, "\t", saveFile)
            TasselFileType.Serial -> resultFile = ExportUtils.writeAlignmentToSerialGZ(inputAlignment, saveFile)
            TasselFileType.VCF -> resultFile = ExportUtils.writeToVCF(inputAlignment, saveFile, keepDepth, this)
            else -> throw IllegalStateException("ExportPlugin: performFunction: Unknown Genotype File Format: $fileType")
        }
        return resultFile
    }

    private fun performFunctionForSimpleTree(input: SimpleTree): String? {
        var resultFile: String? = null
        if (fileType == TasselFileType.Newick) {
            //resultFile = Utils.addSuffixIfNeeded(saveFile(), FileLoadPlugin.FILE_EXT_NEWICK);
            //NewickUtils.write(resultFile, input, includeBranchLengths());
        } else {
            resultFile = Utils.addSuffixIfNeeded(saveFile, ".txt")
            try {
                PrintWriter(resultFile).use { writer -> input.report(writer) }
            } catch (e: Exception) {
                myLogger.debug(e.message, e)
                throw IllegalStateException("ExportPlugin: performFunctionForReport: Problem writing file: $resultFile")
            }
        }
        return resultFile
    }

    private fun performFunctionForTaxaList(input: TaxaList): String {
        return if (fileType == TasselFileType.TaxaList) {
            JSONUtils.exportTaxaListToJSON(input, saveFile)
        } else if (fileType == TasselFileType.Table) {
            val theFile = File(Utils.addSuffixIfNeeded(saveFile, ".txt"))
            TableReportUtils.saveDelimitedTableReport(TaxaListTableReport(input), "\t", theFile)
            theFile.absolutePath
        } else {
            throw IllegalStateException("ExportPlugin: performFunctionForTaxaList: Can't export TaxaList as: " + fileType)
        }
    }

    private fun performFunctionForPositionList(input: PositionList): String {
        return if (fileType == TasselFileType.PositionList) {
            JSONUtils.exportPositionListToJSON(input, saveFile)
        } else if (fileType == TasselFileType.Table) {
            val theFile = File(Utils.addSuffixIfNeeded(saveFile, ".txt"))
            TableReportUtils.saveDelimitedTableReport(PositionListTableReport(input), "\t", theFile)
            theFile.absolutePath
        } else {
            throw IllegalStateException("ExportPlugin: performFunctionForPositionList: Can't export PositionList as: " + fileType)
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    override fun getButtonName(): String {
        return "Save As..."
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    override fun getToolTipText(): String {
        return "Save data to files."
    }

    override fun pluginUserManualURL(): String {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Export/Export"
    }

}