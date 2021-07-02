/*
 * GenotypeSummaryPlugin
 */
package net.maizegenetics.analysis.data

import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.plugindef.*
import net.maizegenetics.util.SimpleTableReport
import net.maizegenetics.util.TableReport
import org.apache.logging.log4j.LogManager
import java.util.*

/**
 * @author Terry Casstevens
 */
class GenotypeSummaryPlugin(isInteractive: Boolean = false) : AbstractPlugin(isInteractive) {
    private var myNumGametesMissing: Long = 0
    private var myNumHeterozygous: Long = 0
    private var myAveMinorAlleleFreq = 0.0

    private var myOverview = PluginParameter.Builder("overview", true, Boolean::class.javaObjectType)
            .description("Get Overview Report").build()
    var overview by Parameter<Boolean>()

    private var mySiteSummary = PluginParameter.Builder("siteSummary", true, Boolean::class.javaObjectType)
            .description("Get Site Summary").build()
    var siteSummary by Parameter<Boolean>()

    private var myTaxaSummary = PluginParameter.Builder("taxaSummary", true, Boolean::class.javaObjectType)
            .description("Get Taxa Summary").build()
    var taxaSummary by Parameter<Boolean>()

    override fun preProcessParameters(input: DataSet) {
        DataSet.data(input, GenotypeTable::class.java)
    }

    override fun processData(input: DataSet): DataSet? {

        if ((!overview)!! && (!siteSummary)!! && (!taxaSummary)!!) {
            printSimpleSummary(input)
            return null
        }

        myNumGametesMissing = 0
        myNumHeterozygous = 0
        myAveMinorAlleleFreq = 0.0

        val data = DataSet.data(input, GenotypeTable::class.java)
        val alignment = data.y
        val name = data.x

        val summaryTables = ArrayList<Datum>()

        var siteSummaryReport: SimpleTableReport? = null
        if (siteSummary) {
            siteSummaryReport = getSiteSummary(alignment)
        }

        var taxaSummaryReport: SimpleTableReport? = null
        if (taxaSummary) {
            taxaSummaryReport = getTaxaSummary(alignment)
        }

        var overallSummaries: Array<SimpleTableReport>? = null
        if (overview) {
            overallSummaries = getOverallSummary(alignment)
            summaryTables.add(Datum(name + "_OverallSummary", overallSummaries[0], "Overall Summary of $name"))
            summaryTables.add(Datum(name + "_AlleleSummary", overallSummaries[1], "Allele Summary of $name"))
        }

        if (siteSummaryReport != null) {
            summaryTables.add(Datum(name + "_SiteSummary", siteSummaryReport, "Site Summary of $name"))
        }
        if (taxaSummaryReport != null) {
            summaryTables.add(Datum(name + "_TaxaSummary", taxaSummaryReport, "Taxa Summary of $name"))
        }

        return when {
            summaryTables.isEmpty() -> null
            else -> DataSet(summaryTables, this)
        }

    }

    private fun getOverallSummary(alignment: GenotypeTable): Array<SimpleTableReport> {

        val firstColumnNames = arrayOf("Stat Type", "Value")

        val numSites = alignment.numberOfSites()
        val numTaxa = alignment.numberOfTaxa().toLong()

        val diploidValueCounts = alignment.genoCounts()
        val numAlleles = diploidValueCounts[0].size

        if (!siteSummary) {
            val totalGametes = numTaxa.toInt() * 2
            for (i in 0 until numSites) {
                val totalGametesNotMissing = alignment.totalGametesNonMissingForSite(i)
                val totalGametesMissing = totalGametes - totalGametesNotMissing
                myNumGametesMissing = myNumGametesMissing + totalGametesMissing.toLong()
                val numHeterozygous = alignment.heterozygousCount(i)
                myNumHeterozygous = myNumHeterozygous + numHeterozygous.toLong()
                myAveMinorAlleleFreq += alignment.minorAlleleFrequency(i)
            }
            myAveMinorAlleleFreq /= numSites.toDouble()
        }

        val totalGametes = numSites * numTaxa * 2L
        val totalGametesNotMissing = totalGametes - myNumGametesMissing

        var numDiploidsMissing: Long = 0
        for (j in 0 until numAlleles) {
            if (diploidValueCounts[0][j] == GenotypeTable.UNKNOWN_ALLELE_STR || diploidValueCounts[0][j] == GenotypeTable.UNKNOWN_GENOTYPE_STR) {
                numDiploidsMissing = diploidValueCounts[1][j] as Long
                break
            }
        }

        val totalDiploids = numSites * numTaxa
        val totalDiploidsNotMissing = totalDiploids - numDiploidsMissing
        var count = 0

        val data = Array(15) { Array<Any?>(firstColumnNames.size) { null } }

        data[count][0] = "Number of Taxa"
        data[count++][1] = numTaxa.toDouble()

        data[count][0] = "Number of Sites"
        data[count++][1] = numSites.toDouble()

        data[count][0] = "Sites x Taxa"
        data[count++][1] = totalDiploids.toDouble()

        data[count][0] = "Number Not Missing"
        data[count++][1] = totalDiploidsNotMissing.toDouble()

        data[count][0] = "Proportion Not Missing"
        data[count++][1] = totalDiploidsNotMissing.toDouble() / totalDiploids.toDouble()

        data[count][0] = "Number Missing"
        data[count++][1] = numDiploidsMissing.toDouble()

        data[count][0] = "Proportion Missing"
        data[count++][1] = numDiploidsMissing.toDouble() / totalDiploids.toDouble()

        data[count][0] = "Number Gametes"
        data[count++][1] = totalGametes.toDouble()

        data[count][0] = "Gametes Not Missing"
        data[count++][1] = totalGametesNotMissing.toDouble()

        data[count][0] = "Proportion Gametes Not Missing"
        data[count++][1] = totalGametesNotMissing.toDouble() / totalGametes.toDouble()

        data[count][0] = "Gametes Missing"
        data[count++][1] = myNumGametesMissing.toDouble()

        data[count][0] = "Proportion Gametes Missing"
        data[count++][1] = myNumGametesMissing.toDouble() / totalGametes.toDouble()

        data[count][0] = "Number Heterozygous"
        data[count++][1] = myNumHeterozygous.toDouble()

        data[count][0] = "Proportion Heterozygous"
        data[count++][1] = myNumHeterozygous.toDouble() / totalDiploids.toDouble()

        data[count][0] = "Average Minor Allele Frequency"
        data[count++][1] = myAveMinorAlleleFreq

        val majorMinorDiploidValueCounts = alignment.majorMinorCounts()
        val numMajorMinorAlleles = majorMinorDiploidValueCounts[0].size

        val alleleColumnNames = arrayOf("Alleles", "Number", "Proportion", "Frequency")
        val data2 = Array(numAlleles + numMajorMinorAlleles) { Array<Any?>(alleleColumnNames.size) { null } }

        count = 0
        for (i in 0 until numAlleles) {
            val value = diploidValueCounts[0][i] as String
            val numValue = diploidValueCounts[1][i] as Long
            data2[count][0] = value
            data2[count][1] = numValue
            data2[count][2] = numValue.toDouble() / totalDiploids.toDouble()
            data2[count++][3] = numValue.toDouble() / totalDiploidsNotMissing.toDouble()
        }

        for (i in 0 until numMajorMinorAlleles) {
            val value = majorMinorDiploidValueCounts[0][i] as String
            val numValue = majorMinorDiploidValueCounts[1][i] as Long
            data2[count][0] = value
            data2[count][1] = numValue
            data2[count++][2] = numValue.toDouble() / numSites.toDouble()
        }

        return arrayOf(SimpleTableReport("Overall Summary", firstColumnNames, data), SimpleTableReport("Allele Summary", alleleColumnNames, data2))
    }

    private fun getSiteSummary(alignment: GenotypeTable): SimpleTableReport {

        val firstColumnNames = arrayOf("Site Number", "Site Name", "Chromosome", "Physical Position", "Number of Taxa", "Ref", "Alt", "Major Allele", "Major Allele Gametes", "Major Allele Proportion", "Major Allele Frequency", "Minor Allele", "Minor Allele Gametes", "Minor Allele Proportion", "Minor Allele Frequency")
        val lastColumnNames = arrayOf("Gametes Missing", "Proportion Missing", "Number Heterozygous", "Proportion Heterozygous", "Inbreeding Coefficient", "Inbreeding Coefficient Scaled by Missing")

        val columnNames = ArrayList(Arrays.asList(*firstColumnNames))

        var maxAlleles = alignment.maxNumAlleles()
        if (alignment.retainsRareAlleles()) {
            maxAlleles++
        }
        for (i in 2 until maxAlleles) {
            val alleleHeading = "Allele " + (i + 1)
            columnNames.add(alleleHeading)
            columnNames.add("$alleleHeading Gametes")
            columnNames.add("$alleleHeading Proportion")
            columnNames.add("$alleleHeading Frequency")
        }

        columnNames.addAll(Arrays.asList(*lastColumnNames))

        val numSites = alignment.numberOfSites()
        val numTaxa = alignment.numberOfTaxa()
        val data = Array(numSites) { Array<Any?>(columnNames.size) { null } }
        val totalGametes = numTaxa * 2

        for (i in 0 until numSites) {

            var count = 0

            data[i][count++] = i
            data[i][count++] = alignment.siteName(i)
            data[i][count++] = alignment.chromosomeName(i)
            data[i][count++] = alignment.chromosomalPosition(i)
            data[i][count++] = numTaxa
            data[i][count++] = alignment.genotypeAsString(i, alignment.referenceAllele(i))
            data[i][count++] = alignment.genotypeAsString(i, alignment.alternateAllele(i))

            val alleles = alignment.allelesSortedByFrequency(i)
            val numAlleles = alleles[0].size
            val totalNotMissing = alignment.totalNonMissingForSite(i)
            val totalGametesNotMissing = AlleleFreqCache.totalGametesNonMissingForSite(alleles)

            for (a in 0 until numAlleles) {
                data[i][count++] = alignment.genotypeAsString(i, alleles[0][a].toByte())
                data[i][count++] = alleles[1][a]
                data[i][count++] = alleles[1][a].toDouble() / totalGametes.toDouble()
                val alleleFreq = alleles[1][a].toDouble() / totalGametesNotMissing.toDouble()
                data[i][count++] = alleleFreq
                if (a == 1) {
                    myAveMinorAlleleFreq += alleleFreq
                }
            }

            for (b in 0 until maxAlleles - numAlleles) {
                data[i][count++] = NA
                data[i][count++] = ZERO_INT
                data[i][count++] = ZERO_DOUBLE
                data[i][count++] = ZERO_DOUBLE
            }

            val totalGametesMissing = totalGametes - totalGametesNotMissing
            myNumGametesMissing += totalGametesMissing.toLong()
            data[i][count++] = totalGametesMissing
            data[i][count++] = totalGametesMissing.toDouble() / totalGametes.toDouble()

            val numHeterozygous = alignment.heterozygousCount(i)
            myNumHeterozygous += numHeterozygous.toLong()
            data[i][count++] = numHeterozygous
            data[i][count++] = numHeterozygous.toDouble() / totalNotMissing.toDouble()

            data[i][count++] = "TBD"
            data[i][count++] = "TBD"

        }

        myAveMinorAlleleFreq /= numSites.toDouble()

        return SimpleTableReport("Site Summary", columnNames.toArray(), data)

    }

    private fun getTaxaSummary(alignment: GenotypeTable): SimpleTableReport {

        val columnNames = arrayOf("Taxa", "Taxa Name", "Number of Sites", "Gametes Missing", "Proportion Missing", "Number Heterozygous", "Proportion Heterozygous", "Inbreeding Coefficient", "Inbreeding Coefficient Scaled by Missing")
        val numSites = alignment.numberOfSites()
        val numTaxa = alignment.numberOfTaxa()
        val totalGametes = numSites * 2

        val gametesNotMissing = IntArray(numTaxa)
        val sitesNotMissing = IntArray(numTaxa)
        val heterozygous = IntArray(numTaxa)
        for (s in 0 until numSites) {
            val genotypes = alignment.genotypeAllTaxa(s)
            for (t in 0 until numTaxa) {
                val alleles = GenotypeTableUtils.getDiploidValues(genotypes[t])
                if (alleles[0] != GenotypeTable.UNKNOWN_ALLELE) {
                    if (alleles[1] != GenotypeTable.UNKNOWN_ALLELE) {
                        gametesNotMissing[t] += 2
                        sitesNotMissing[t]++
                    } else {
                        gametesNotMissing[t]++
                        sitesNotMissing[t]++
                    }
                } else if (alleles[1] != GenotypeTable.UNKNOWN_ALLELE) {
                    gametesNotMissing[t]++
                    sitesNotMissing[t]++
                }

                if (alleles[0] != alleles[1]) {
                    heterozygous[t]++
                }
            }
        }

        val data = Array(numTaxa) { Array<Any?>(columnNames.size) { null } }
        for (i in 0 until numTaxa) {

            val totalGametesMissing = totalGametes - gametesNotMissing[i]

            var count = 0
            data[i][count++] = i
            data[i][count++] = alignment.taxaName(i)
            data[i][count++] = numSites
            data[i][count++] = totalGametesMissing
            data[i][count++] = totalGametesMissing.toDouble() / totalGametes.toDouble()
            data[i][count++] = heterozygous[i]
            data[i][count++] = heterozygous[i].toDouble() / sitesNotMissing[i].toDouble()
            data[i][count++] = "Inbreeding Coefficient"
            data[i][count++] = "ICSBM"
        }

        return SimpleTableReport("Taxa Summary", columnNames, data)

    }

    /**
     * Convenience method to run plugin with one return object.
     */
    fun runPlugin(genotype: GenotypeTable): List<TableReport> {
        val input = DataSet(Datum("Genotype Table", genotype, null), this)
        val dataSet = performFunction(input)
        val result = ArrayList<TableReport>()
        for (i in 0 until dataSet.size) {
            result[i] = dataSet.getData(i)!!.data as TableReport
        }
        return result
    }

    override fun icon(): String? {
        return "/images/summary.gif"
    }

    override fun getButtonName(): String {
        return "Genotype Summary"
    }

    override fun getToolTipText(): String {
        return "Genotype Summary"
    }

    override fun pluginUserManualURL(): String {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/GenoSummary/GenoSummary"
    }

    companion object {

        private val myLogger = LogManager.getLogger(GenotypeSummaryPlugin::class.java)

        private val NA = "NA"
        private val ZERO_DOUBLE = 0.0
        private val ZERO_INT = 0

        @JvmStatic
        fun printSimpleSummary(input: DataSet?) {
            if (input == null) {
                return
            }
            val alignInList = input.getDataOfType(GenotypeTable::class.java)
            if (alignInList!!.isEmpty()) {
                return
            }
            printSimpleSummary(alignInList[0])
        }

        @JvmStatic
        fun printSimpleSummary(current: Datum) {
            val alignment = current.data as GenotypeTable
            val name = current.name
            printSimpleSummary(alignment, name)
        }

        @JvmStatic
        fun printSimpleSummary(alignment: GenotypeTable, name: String) {

            val numSites = alignment.numberOfSites().toLong()
            val numTaxa = alignment.numberOfTaxa().toLong()

            val totalDiploids = numSites * numTaxa

            println("Genotype Table Name: $name")
            println("Number of Taxa: $numTaxa")
            println("Number of Sites: $numSites")
            println("Sites x Taxa: $totalDiploids")

            println("Chromosomes...")
            val chromosomes = alignment.chromosomes()
            val positions = alignment.positions()
            for (i in chromosomes.indices) {
                val startEnd = alignment.firstLastSiteOfChromosome(chromosomes[i])
                println(chromosomes[i].name + ": start site: " + startEnd[0] + " (" + positions[startEnd[0]].position + ") last site: " + startEnd[1] + " (" + positions[startEnd[1]].position + ") total: " + (startEnd[1] - startEnd[0] + 1))
            }
            println()

        }
    }

}

fun main(args: Array<String>) {
    val plugin = GenotypeSummaryPlugin()
    plugin.siteSummary = true
    println(plugin.siteSummary)
}
