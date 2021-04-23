package net.maizegenetics.analysis.distance;

import com.google.common.collect.Range;
import net.maizegenetics.dna.factor.FactorTable;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import org.apache.log4j.Logger;

import java.awt.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Terry Casstevens
 * @author Zhiwu Zhang
 * @author Peter Bradbury
 */
public class KinshipPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(KinshipPlugin.class);

    public enum KINSHIP_METHOD {

        Centered_IBS,
        Normalized_IBS,
        Dominance_Centered_IBS,
        Dominance_Normalized_IBS
    }

    public enum ALGORITHM_VARIATION {

        Observed_Allele_Freq,
        Proportion_Heterozygous
    }

    private PluginParameter<KINSHIP_METHOD> myMethod = new PluginParameter.Builder<>("method", KINSHIP_METHOD.Centered_IBS, KINSHIP_METHOD.class)
            .guiName("Kinship method")
            .range(KINSHIP_METHOD.values())
            .description("The Centered_IBS (Endelman - previously Scaled_IBS) method produces a kinship matrix that is scaled to give a reasonable estimate of additive "
                    + "genetic variance. Uses algorithm http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13. "
                    + "The Normalized_IBS (Previously GCTA) uses the algorithm published here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/pdf/main.pdf.")
            .build();

    private PluginParameter<Integer> myMaxAlleles = new PluginParameter.Builder<>("maxAlleles", 255, Integer.class)
            .description("")
            .range(Range.closed(2, 6))
            .dependentOnParameter(myMethod, new Object[]{KINSHIP_METHOD.Centered_IBS, KINSHIP_METHOD.Dominance_Centered_IBS})
            .build();

    private PluginParameter<ALGORITHM_VARIATION> myAlgorithmVariation = new PluginParameter.Builder<>("algorithmVariation", ALGORITHM_VARIATION.Observed_Allele_Freq, ALGORITHM_VARIATION.class)
            .description("")
            .range(ALGORITHM_VARIATION.values())
            .dependentOnParameter(myMethod, new Object[]{KINSHIP_METHOD.Dominance_Centered_IBS})
            .build();

    public KinshipPlugin(Frame parentFrame, boolean isInteractive) {
        super(isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(FactorTable.class);
        if ((alignInList == null) || (alignInList.isEmpty())) {
            throw new IllegalArgumentException("KinshipPlugin: Nothing selected. Please select a genotype.");
        }
    }

    @Override
    public void setParameters(String[] args) {
        try {
            for (int i = 0; i < args.length; i++) {
                if (args[i].equals("-method")) {
                    if (args[i + 1].equalsIgnoreCase("GCTA")) {
                        args[i + 1] = KINSHIP_METHOD.Normalized_IBS.name();
                        myLogger.warn("setParameters: Notice GCTA has been changed to Normalized_IBS");
                    } else if (args[i + 1].equalsIgnoreCase("Scaled_IBS")) {
                        args[i + 1] = KINSHIP_METHOD.Centered_IBS.name();
                        myLogger.warn("setParameters: Notice Scaled_IBS has been changed to Centered_IBS");
                    } else if (args[i + 1].equalsIgnoreCase("Dominance")) {
                        args[i + 1] = KINSHIP_METHOD.Dominance_Centered_IBS.name();
                        myLogger.warn("setParameters: Notice Dominance has been changed to Dominance_Centered_IBS");
                    }
                }
                break;
            }
        } catch (Exception e) {
            // do nothing
            myLogger.debug(e.getMessage(), e);
        }
        super.setParameters(args);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> alignInList = input.getDataOfType(FactorTable.class);

        List<Datum> result = new ArrayList<>();
        Iterator<Datum> itr = alignInList.iterator();
        while (itr.hasNext()) {

            Datum current = itr.next();
            String datasetName = current.getName();
            DistanceMatrix kin = null;

            if (current.getData() instanceof FactorTable) {
                FactorTable myGenotype = (FactorTable) current.getData();
                if (kinshipMethod() == KINSHIP_METHOD.Centered_IBS) {
                    System.out.println("KinshipPlugin: maxAlleles: " + maxAlleles());
                    kin = EndelmanDistanceMatrix.getInstance(myGenotype, maxAlleles(), this);
                } else if (kinshipMethod() == KINSHIP_METHOD.Normalized_IBS) {
                    //kin = GCTADistanceMatrix.getInstance(myGenotype, this);
                } else if (kinshipMethod() == KINSHIP_METHOD.Dominance_Centered_IBS) {
                    //kin = DominanceRelationshipMatrix.getInstance(myGenotype, maxAlleles(), algorithmVariation(), this);
                } else if (kinshipMethod() == KINSHIP_METHOD.Dominance_Normalized_IBS) {
                    //kin = DominanceNormalizedIBSMatrix.getInstance(myGenotype, this);
                } else {
                    throw new IllegalArgumentException("Unknown method to calculate kinship: " + kinshipMethod());
                }
            } else {
                throw new IllegalArgumentException("Invalid selection. Can't create kinship matrix from: " + datasetName);
            }

            if (kin != null) {
                StringBuilder comment = new StringBuilder();
                comment.append(kinshipMethod());
                if (kinshipMethod() == KINSHIP_METHOD.Dominance_Centered_IBS) {
                    comment.append("(variation: ");
                    comment.append(algorithmVariation());
                    comment.append(")");
                }
                comment.append(" matrix created from ");
                comment.append(datasetName);
                Datum ds = new Datum(kinshipMethod() + "_" + datasetName, kin, comment.toString());
                result.add(ds);
            }

        }

        return new DataSet(result, this);

    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Kinship/Kinship";
    }

    @Override
    public String getButtonName() {
        return "Kinship";
    }

    @Override
    public String getToolTipText() {
        return "Calculate kinship from marker data";
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public DistanceMatrix runPlugin(DataSet input) {
        return (DistanceMatrix) performFunction(input).getData(0).getData();
    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public DistanceMatrix runPlugin(FactorTable input) {
        return (DistanceMatrix) performFunction(DataSet.getDataSet(input)).getData(0).getData();
    }

    /**
     * The scaled_IBS method produces a kinship matrix that is scaled to give a
     * reasonable estimate of additive genetic variance. The pairwise_IBS
     * method, which is the method used by TASSEL ver.4, may result in an
     * inflated estimate of genetic variance. Either will do a good job of
     * controlling population structure in MLM. The pedigree method is used to
     * calculate a kinship matrix from a pedigree information.
     *
     * @return Kinship method
     */
    public KINSHIP_METHOD kinshipMethod() {
        return myMethod.value();
    }

    /**
     * Set Kinship method. The scaled_IBS method produces a kinship matrix that
     * is scaled to give a reasonable estimate of additive genetic variance. The
     * pairwise_IBS method, which is the method used by TASSEL ver.4, may result
     * in an inflated estimate of genetic variance. Either will do a good job of
     * controlling population structure in MLM. The pedigree method is used to
     * calculate a kinship matrix from a pedigree information.
     *
     * @param value Kinship method
     *
     * @return this plugin
     */
    public KinshipPlugin kinshipMethod(KINSHIP_METHOD value) {
        myMethod = new PluginParameter<>(myMethod, value);
        return this;
    }

    /**
     * Max Alleles
     *
     * @return Max Alleles
     */
    public Integer maxAlleles() {
        return myMaxAlleles.value();
    }

    /**
     * Set Max Alleles. Max Alleles
     *
     * @param value Max Alleles
     *
     * @return this plugin
     */
    public KinshipPlugin maxAlleles(Integer value) {
        myMaxAlleles = new PluginParameter<>(myMaxAlleles, value);
        return this;
    }

    /**
     * Algorithm Variation
     *
     * @return Algorithm Variation
     */
    public ALGORITHM_VARIATION algorithmVariation() {
        return myAlgorithmVariation.value();
    }

    /**
     * Set Algorithm Variation. Algorithm Variation
     *
     * @param value Algorithm Variation
     *
     * @return this plugin
     */
    public KinshipPlugin algorithmVariation(ALGORITHM_VARIATION value) {
        myAlgorithmVariation = new PluginParameter<>(myAlgorithmVariation, value);
        return this;
    }
}
