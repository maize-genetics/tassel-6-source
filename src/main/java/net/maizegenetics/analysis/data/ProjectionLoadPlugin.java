/*
 * ProjectionLoadPlugin
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.util.List;

/**
 * @author Alex Lipka This should enable a used to load a projection alignment using the TASSEL GUI
 */
public class ProjectionLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ProjectionLoadPlugin.class);

    private PluginParameter<String> myRecombinationBreakpoints = new PluginParameter.Builder<>("recombinationBreakpoints", null, String.class).required(true).inFile()
            .description("").build();

    private GenotypeTable myHighDensityMarkersGenotypeTable = null;

    /**
     * Creates a new instance of ProjectionLoadPlugin
     */
    public ProjectionLoadPlugin(boolean isInteractive) {
        super(isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException("ProjectionLoadPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);
        if (genotypeTables.size() == 1) {
            myHighDensityMarkersGenotypeTable = (GenotypeTable) genotypeTables.get(0).getData();
        } else {
            throw new IllegalArgumentException("ProjectionLoadPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        try {
            return loadFile(myRecombinationBreakpoints.value(), myHighDensityMarkersGenotypeTable);
        } catch (Exception e) {
            throw new IllegalStateException("ProjectionLoadPlugin: processData: Problem loading: " + myRecombinationBreakpoints.value() + "\n" + e.getMessage());
        }

    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Load Projection Alignment";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Load Projection Alignments";
    }

    private DataSet loadFile(String theRecombinationBreakpoints, GenotypeTable theHighDensityMarkers) {

        GenotypeTable theAlignmentForGenotype = ProjectionGenotypeIO.getInstance(theRecombinationBreakpoints, theHighDensityMarkers);
        Datum td = new Datum(Utils.getFilename(theRecombinationBreakpoints, FileLoadPlugin.FILE_EXT_HAPMAP), theAlignmentForGenotype, null);
        return new DataSet(td, this);

    }

    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Recombination Breakpoints
     *
     * @return Recombination Breakpoints
     */
    public String recombinationBreakpoints() {
        return myRecombinationBreakpoints.value();
    }

    /**
     * Set Recombination Breakpoints. Recombination Breakpoints
     *
     * @param value Recombination Breakpoints
     *
     * @return this plugin
     */
    public ProjectionLoadPlugin recombinationBreakpoints(String value) {
        myRecombinationBreakpoints = new PluginParameter<>(myRecombinationBreakpoints, value);
        return this;
    }

}
