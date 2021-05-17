/*
 * KinshipPlugin.java
 *
 * Created on May 3, 2021
 *
 */
package net.maizegenetics.analysis.distance

import com.google.common.collect.Range
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.plugindef.*
import net.maizegenetics.taxa.distance.DistanceMatrix
import org.apache.log4j.Logger
import java.util.*

/**
 * @author Terry Casstevens
 */
class KinshipPlugin(isInteractive: Boolean = false) : AbstractPlugin(isInteractive) {

    private val logger = Logger.getLogger(KinshipPlugin::class.java)

    enum class KINSHIP_METHOD {
        Centered_IBS, Normalized_IBS, Dominance_Centered_IBS, Dominance_Normalized_IBS
    }

    enum class ALGORITHM_VARIATION {
        Observed_Allele_Freq, Proportion_Heterozygous
    }

    private var myMethod = PluginParameter.Builder("method", KINSHIP_METHOD.Centered_IBS, KINSHIP_METHOD::class.java)
            .guiName("Kinship method")
            .range(KINSHIP_METHOD.values())
            .description("The Centered_IBS (Endelman - previously Scaled_IBS) method produces a kinship matrix that is scaled to give a reasonable estimate of additive "
                    + "genetic variance. Uses algorithm http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13. "
                    + "The Normalized_IBS (Previously GCTA) uses the algorithm published here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/pdf/main.pdf.")
            .build()
    var method by Parameter<KINSHIP_METHOD>()

    private var myMaxAlleles = PluginParameter.Builder("maxAlleles", 255, Int::class.java)
            .description("")
            .range(Range.closed<Comparable<Int>>(2, 255))
            .dependentOnParameter(myMethod, arrayOf<Any>(KINSHIP_METHOD.Centered_IBS, KINSHIP_METHOD.Dominance_Centered_IBS))
            .build()
    var maxAlleles by Parameter<Int>()

    private var myAlgorithmVariation = PluginParameter.Builder("algorithmVariation", ALGORITHM_VARIATION.Observed_Allele_Freq, ALGORITHM_VARIATION::class.java)
            .description("")
            .range(ALGORITHM_VARIATION.values())
            .dependentOnParameter(myMethod, arrayOf<Any>(KINSHIP_METHOD.Dominance_Centered_IBS))
            .build()
    var algorithmVariation by Parameter<ALGORITHM_VARIATION>()

    override fun preProcessParameters(input: DataSet) {
        val alignInList = input.getDataOfType(FeatureTable::class.java)
        require(!(alignInList == null || alignInList.isEmpty())) { "KinshipPlugin: Nothing selected. Please select a genotype." }
    }

    override fun processData(input: DataSet): DataSet {
        val alignInList = input.getDataOfType(FeatureTable::class.java)
        val result: MutableList<Datum> = ArrayList()
        val itr: Iterator<Datum> = alignInList.iterator()
        while (itr.hasNext()) {
            val current = itr.next()
            val datasetName = current.name
            var kin: DistanceMatrix? = null
            if (current.data is FeatureTable) {
                val myGenotype = current.data as FeatureTable
                if (method == KINSHIP_METHOD.Centered_IBS) {
                    kin = EndelmanDistanceMatrixBuilder(myGenotype, maxAlleles, this).build()
                } else if (method == KINSHIP_METHOD.Normalized_IBS) {
                    //kin = GCTADistanceMatrix.getInstance(myGenotype, this);
                } else if (method == KINSHIP_METHOD.Dominance_Centered_IBS) {
                    //kin = DominanceRelationshipMatrix.getInstance(myGenotype, maxAlleles(), algorithmVariation(), this);
                } else if (method == KINSHIP_METHOD.Dominance_Normalized_IBS) {
                    //kin = DominanceNormalizedIBSMatrix.getInstance(myGenotype, this);
                } else {
                    throw IllegalArgumentException("Unknown method to calculate kinship: $method")
                }
            } else {
                throw IllegalArgumentException("Invalid selection. Can't create kinship matrix from: $datasetName")
            }
            if (kin != null) {
                val comment = StringBuilder()
                comment.append(method)
                if (method == KINSHIP_METHOD.Dominance_Centered_IBS) {
                    comment.append("(variation: ")
                    comment.append(algorithmVariation)
                    comment.append(")")
                }
                comment.append(" matrix created from ")
                comment.append(datasetName)
                val ds = Datum(method.toString() + "_" + datasetName, kin, comment.toString())
                result.add(ds)
            }
        }
        return DataSet(result, this)
    }

    override fun pluginUserManualURL(): String {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Kinship/Kinship"
    }

    override fun getButtonName(): String {
        return "Kinship"
    }

    override fun getToolTipText(): String {
        return "Calculate kinship from marker data"
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    fun run(input: DataSet): DistanceMatrix {
        return performFunction(input).getData(0).data as DistanceMatrix
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    fun run(input: FeatureTable): DistanceMatrix {
        return performFunction(DataSet.getDataSet(input)).getData(0).data as DistanceMatrix
    }

}