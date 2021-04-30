/*
 * ExportPlugin.java
 *
 * Created on December 18, 2009
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListTableReport;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.io.FilterJSONUtils;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.dna.snp.io.SiteScoresIO;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListTableReport;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixUtils;
import net.maizegenetics.taxa.distance.WriteDistanceMatrix;
import net.maizegenetics.taxa.tree.SimpleTree;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

//import net.maizegenetics.dna.snp.io.FlapjackUtils;
//import net.maizegenetics.taxa.tree.NewickUtils;

/**
 * @author Terry Casstevens
 */
public class ExportPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ExportPlugin.class);

    private PluginParameter<String> mySaveFile = new PluginParameter.Builder<>("saveAs", null, String.class)
            .description("Save file as...")
            .outFile()
            .required(true)
            .build();

    private PluginParameter<FileLoadPlugin.TasselFileType> myFileType = new PluginParameter.Builder<>("format", null, FileLoadPlugin.TasselFileType.class)
            .description("Export file format (Default format depends on data being exported)")
            .required(true)
            .objectListSingleSelect()
            .range(FileLoadPlugin.TasselFileType.values())
            .build();

    private PluginParameter<Boolean> myKeepDepth = new PluginParameter.Builder<>("keepDepth", true, Boolean.class)
            .description("Whether to keep depth if format supports depth.")
            .dependentOnParameter(myFileType, new FileLoadPlugin.TasselFileType[]{FileLoadPlugin.TasselFileType.VCF})
            .build();

    private PluginParameter<Boolean> myIncludeTaxaAnnotations = new PluginParameter.Builder<>("includeTaxaAnnotations", true, Boolean.class)
            .description("Whether to include taxa annotations if format supports taxa annotations.")
            .dependentOnParameter(myFileType, new FileLoadPlugin.TasselFileType[]{FileLoadPlugin.TasselFileType.VCF,
                    FileLoadPlugin.TasselFileType.Hapmap,
                    FileLoadPlugin.TasselFileType.HapmapDiploid,
                    FileLoadPlugin.TasselFileType.HapmapLIX})
            .build();

    private PluginParameter<Boolean> myIncludeBranchLengths = new PluginParameter.Builder<>("includeBranchLengths", true, Boolean.class)
            .description("Whether to include branch lengths for Newick formatted files.")
            .dependentOnParameter(myFileType, new FileLoadPlugin.TasselFileType[]{FileLoadPlugin.TasselFileType.Newick})
            .build();

    /**
     * Creates a new instance of ExportPlugin
     */
    public ExportPlugin(boolean isInteractive) {
        super(isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        if (input.getSize() != 1) {
            throw new IllegalArgumentException("Please select only one item.");
        }

        Object data = input.getData(0).getData();
        if (data instanceof GenotypeTable) {
            GenotypeTable genotype = (GenotypeTable) data;
            List<FileLoadPlugin.TasselFileType> temp = new ArrayList<>();
            if (genotype.hasGenotype()) {
                temp.add(FileLoadPlugin.TasselFileType.HaplotypeVCF);
                temp.add(FileLoadPlugin.TasselFileType.Hapmap);
                temp.add(FileLoadPlugin.TasselFileType.HapmapDiploid);
                temp.add(FileLoadPlugin.TasselFileType.VCF);
                temp.add(FileLoadPlugin.TasselFileType.Plink);
                temp.add(FileLoadPlugin.TasselFileType.Phylip_Seq);
                temp.add(FileLoadPlugin.TasselFileType.Phylip_Inter);
                temp.add(FileLoadPlugin.TasselFileType.Table);
            }
            if (genotype.hasDepth()) {
                temp.add(FileLoadPlugin.TasselFileType.Depth);
            }
            if (genotype.hasReferenceProbablity()) {
                temp.add(FileLoadPlugin.TasselFileType.ReferenceProbability);
            }
            myFileType = new PluginParameter<>(myFileType, temp);
        } else if (data instanceof Phenotype) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(FileLoadPlugin.TasselFileType.Phenotype,
                            FileLoadPlugin.TasselFileType.PlinkPhenotype));
        } else if (data instanceof FilterList) {
            myFileType = new PluginParameter<>(myFileType,
                    Collections.singletonList(FileLoadPlugin.TasselFileType.Filter));
        } else if (data instanceof DistanceMatrix) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(FileLoadPlugin.TasselFileType.SqrMatrix,
                            FileLoadPlugin.TasselFileType.SqrMatrixBin,
                            FileLoadPlugin.TasselFileType.SqrMatrixRaw,
                            FileLoadPlugin.TasselFileType.SqrMatrixDARwinDIS));
        } else if (data instanceof TaxaList) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(FileLoadPlugin.TasselFileType.TaxaList,
                            FileLoadPlugin.TasselFileType.Table));
        } else if (data instanceof TaxaListTableReport) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(FileLoadPlugin.TasselFileType.TaxaList,
                            FileLoadPlugin.TasselFileType.Table));
        } else if (data instanceof PositionList) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(FileLoadPlugin.TasselFileType.PositionList,
                            FileLoadPlugin.TasselFileType.Table));
        } else if (data instanceof PositionListTableReport) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(FileLoadPlugin.TasselFileType.PositionList,
                            FileLoadPlugin.TasselFileType.Table));
        } else if (data instanceof TableReport) {
            myFileType = new PluginParameter<>(myFileType,
                    Collections.singletonList(FileLoadPlugin.TasselFileType.Table));
        } else if (data instanceof SimpleTree) {
            myFileType = new PluginParameter<>(myFileType,
                    Arrays.asList(//FileLoadPlugin.TasselFileType.Newick,
                            FileLoadPlugin.TasselFileType.Report));
        } else {
            throw new IllegalStateException("Don't know how to export data type: " + data.getClass().getName());
        }

        if (!isInteractive() && myFileType.isEmpty() && myFileType.hasPossibleValues()) {
            fileType(myFileType.possibleValues().get(0));
        }

    }

    @Override
    public DataSet processData(DataSet input) {

        String filename = null;
        Object data = input.getData(0).getData();
        if (data instanceof GenotypeTable) {
            filename = performFunctionForAlignment((GenotypeTable) data);
        } else if (data instanceof Phenotype) {
            filename = performFunctionForPhenotype((Phenotype) data);
        } else if (data instanceof FilterList) {
            filename = performFunctionForFilter((FilterList) data);
        } else if (data instanceof DistanceMatrix) {
            filename = performFunctionForDistanceMatrix((DistanceMatrix) data);
        } else if (data instanceof TaxaList) {
            filename = performFunctionForTaxaList((TaxaList) data);
        } else if (data instanceof TaxaListTableReport) {
            filename = performFunctionForTaxaList(((TaxaListTableReport) data).getTaxaList());
        } else if (data instanceof PositionList) {
            filename = performFunctionForPositionList((PositionList) data);
        } else if (data instanceof PositionListTableReport) {
            filename = performFunctionForPositionList(((PositionListTableReport) data).getPositionList());
        } else if (data instanceof TableReport) {
            filename = performFunctionForTableReport((TableReport) data);
        } else if (data instanceof SimpleTree) {
            filename = performFunctionForSimpleTree((SimpleTree) data);
        } else {
            throw new IllegalStateException("Don't know how to export data type: " + data.getClass().getName());
        }

        if (filename != null) {
            myLogger.info("performFunction: wrote dataset: " + input.getData(0).getName() + " to file: " + filename);
        }

        return null;

    }

    public String performFunctionForDistanceMatrix(DistanceMatrix input) {

        if (fileType() == FileLoadPlugin.TasselFileType.SqrMatrix) {
            String filename = Utils.addSuffixIfNeeded(saveFile(), ".txt", new String[]{".txt", ".txt.gz"});
            WriteDistanceMatrix.saveDelimitedDistanceMatrix(input, filename);
            return filename;
        } else if (fileType() == FileLoadPlugin.TasselFileType.SqrMatrixRaw) {
            String[] grmFiles = DistanceMatrixUtils.getGRMFilenames(saveFile());
            WriteDistanceMatrix.saveRawMultiBlupMatrix(input, grmFiles[0], grmFiles[3]);
            return grmFiles[3];
        } else if (fileType() == FileLoadPlugin.TasselFileType.SqrMatrixBin) {
            String[] grmFiles = DistanceMatrixUtils.getGRMFilenames(saveFile());
            WriteDistanceMatrix.saveBinMultiBlupMatrix(input, grmFiles[0], grmFiles[1], grmFiles[2]);
            return grmFiles[1];
        } else if (fileType() == FileLoadPlugin.TasselFileType.SqrMatrixDARwinDIS) {
            String filename = Utils.addSuffixIfNeeded(saveFile(), ".dis");
            WriteDistanceMatrix.saveDARwinMatrix(input, filename);
            return filename;
        } else {
            throw new IllegalArgumentException("ExportPlugin: performFunctionForDistanceMatrix: Unknown file type: " + fileType());
        }

    }

    public String performFunctionForTableReport(TableReport input) {
        File theFile = new File(Utils.addSuffixIfNeeded(saveFile(), ".txt"));
        TableReportUtils.saveDelimitedTableReport(input, "\t", theFile);
        return theFile.getAbsolutePath();
    }

    public String performFunctionForFilter(FilterList filter) {
        return FilterJSONUtils.exportFilterToJSON(filter, saveFile());
    }

    public String performFunctionForPhenotype(Phenotype input) {
        String filename = Utils.addSuffixIfNeeded(saveFile(), ".txt");
        if (fileType() == FileLoadPlugin.TasselFileType.Phenotype) {
            PhenotypeUtils.write(input, filename);
        } else if (fileType() == FileLoadPlugin.TasselFileType.PlinkPhenotype) {
            PhenotypeUtils.writePlink(input, filename);
        }
        return new File(filename).getAbsolutePath();
    }

    public String performFunctionForAlignment(GenotypeTable inputAlignment) {

        String resultFile = saveFile();

        if (fileType() == FileLoadPlugin.TasselFileType.ReferenceProbability) {
            resultFile = SiteScoresIO.writeReferenceProbability(inputAlignment, resultFile);
        } else if (fileType() == FileLoadPlugin.TasselFileType.Depth) {
            resultFile = SiteScoresIO.writeDepth(inputAlignment, resultFile);
        } else if (fileType() == FileLoadPlugin.TasselFileType.Hapmap) {
            resultFile = ExportUtils.writeToHapmap(inputAlignment, false, saveFile(), '\t', includeTaxaAnnotations(), this);
        } else if (fileType() == FileLoadPlugin.TasselFileType.HapmapDiploid) {
            resultFile = ExportUtils.writeToHapmap(inputAlignment, true, saveFile(), '\t', includeTaxaAnnotations(), this);
        } else if (fileType() == FileLoadPlugin.TasselFileType.Plink) {
            resultFile = ExportUtils.writeToPlink(inputAlignment, saveFile(), '\t');
            //} else if (fileType() == FileLoadPlugin.TasselFileType.Flapjack) {
            //resultFile = FlapjackUtils.writeToFlapjack(inputAlignment, saveFile(), '\t');
        } else if (fileType() == FileLoadPlugin.TasselFileType.Phylip_Seq) {
            resultFile = Utils.addSuffixIfNeeded(saveFile(), ".phy");
            try (PrintWriter out = new PrintWriter(new FileWriter(resultFile))) {
                ExportUtils.printSequential(inputAlignment, out);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("ExportPlugin: performFunction: Problem writing file: " + resultFile);
            }
        } else if (fileType() == FileLoadPlugin.TasselFileType.Phylip_Inter) {
            resultFile = Utils.addSuffixIfNeeded(saveFile(), ".phy");
            try (PrintWriter out = new PrintWriter(new FileWriter(resultFile))) {
                ExportUtils.printInterleaved(inputAlignment, out);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("ExportPlugin: performFunction: Problem writing file: " + resultFile);
            }
        } else if (fileType() == FileLoadPlugin.TasselFileType.Table) {
            resultFile = ExportUtils.saveDelimitedAlignment(inputAlignment, "\t", saveFile());
        } else if (fileType() == FileLoadPlugin.TasselFileType.Serial) {
            resultFile = ExportUtils.writeAlignmentToSerialGZ(inputAlignment, saveFile());
        } else if (fileType() == FileLoadPlugin.TasselFileType.VCF) {
            resultFile = ExportUtils.writeToVCF(inputAlignment, saveFile(), keepDepth(), this);
        } else {
            throw new IllegalStateException("ExportPlugin: performFunction: Unknown Genotype File Format: " + fileType());
        }

        return resultFile;

    }

    public String performFunctionForSimpleTree(SimpleTree input) {

        String resultFile = null;
        if (fileType() == FileLoadPlugin.TasselFileType.Newick) {
            //resultFile = Utils.addSuffixIfNeeded(saveFile(), FileLoadPlugin.FILE_EXT_NEWICK);
            //NewickUtils.write(resultFile, input, includeBranchLengths());
        } else {
            resultFile = Utils.addSuffixIfNeeded(saveFile(), ".txt");
            try (PrintWriter writer = new PrintWriter(resultFile)) {
                input.report(writer);
            } catch (Exception e) {
                myLogger.debug(e.getMessage(), e);
                throw new IllegalStateException("ExportPlugin: performFunctionForReport: Problem writing file: " + resultFile);
            }
        }
        return resultFile;

    }

    public String performFunctionForTaxaList(TaxaList input) {
        if (fileType() == FileLoadPlugin.TasselFileType.TaxaList) {
            return JSONUtils.exportTaxaListToJSON(input, saveFile());
        } else if (fileType() == FileLoadPlugin.TasselFileType.Table) {
            File theFile = new File(Utils.addSuffixIfNeeded(saveFile(), ".txt"));
            TableReportUtils.saveDelimitedTableReport(new TaxaListTableReport(input), "\t", theFile);
            return theFile.getAbsolutePath();
        } else {
            throw new IllegalStateException("ExportPlugin: performFunctionForTaxaList: Can't export TaxaList as: " + fileType());
        }
    }

    public String performFunctionForPositionList(PositionList input) {
        if (fileType() == FileLoadPlugin.TasselFileType.PositionList) {
            return JSONUtils.exportPositionListToJSON(input, saveFile());
        } else if (fileType() == FileLoadPlugin.TasselFileType.Table) {
            File theFile = new File(Utils.addSuffixIfNeeded(saveFile(), ".txt"));
            TableReportUtils.saveDelimitedTableReport(new PositionListTableReport(input), "\t", theFile);
            return theFile.getAbsolutePath();
        } else {
            throw new IllegalStateException("ExportPlugin: performFunctionForPositionList: Can't export PositionList as: " + fileType());
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Save As...";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Save data to files.";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Export/Export";
    }

    /**
     * Save file as...
     *
     * @return Save As
     */
    public String saveFile() {
        return mySaveFile.value();
    }

    /**
     * Save file as...
     *
     * @param value filename
     *
     * @return this plugin
     */
    public ExportPlugin saveFile(String value) {
        mySaveFile = new PluginParameter<>(mySaveFile, value);
        return this;
    }

    @Deprecated
    /**
     * Deprecated use saveFile(value)
     */
    public ExportPlugin setSaveFile(String saveFile) {
        mySaveFile = new PluginParameter<>(mySaveFile, saveFile);
        return this;
    }

    /**
     * Export file format
     *
     * @return Format
     */
    public FileLoadPlugin.TasselFileType fileType() {
        return myFileType.value();
    }

    /**
     * Set Format. Export file format
     *
     * @param value Format
     *
     * @return this plugin
     */
    public ExportPlugin fileType(FileLoadPlugin.TasselFileType value) {
        myFileType = new PluginParameter<>(myFileType, value);
        return this;
    }

    @Deprecated
    /**
     * Deprecated use fileType(value)
     */
    public ExportPlugin setAlignmentFileType(FileLoadPlugin.TasselFileType type) {
        fileType(type);
        return this;
    }

    /**
     * Whether to keep depth if format supports depth.
     *
     * @return Keep Depth
     */
    public Boolean keepDepth() {
        return myKeepDepth.value();
    }

    /**
     * Set Keep Depth. Whether to keep depth if format supports depth.
     *
     * @param value Keep Depth
     *
     * @return this plugin
     */
    public ExportPlugin keepDepth(Boolean value) {
        myKeepDepth = new PluginParameter<>(myKeepDepth, value);
        return this;
    }

    /**
     * Whether to include taxa annotations if format supports taxa annotations.
     *
     * @return Include Taxa Annotations
     */
    public Boolean includeTaxaAnnotations() {
        return myIncludeTaxaAnnotations.value();
    }

    /**
     * Set Include Taxa Annotations. Whether to include taxa annotations if
     * format supports taxa annotations.
     *
     * @param value Include Taxa Annotations
     *
     * @return this plugin
     */
    public ExportPlugin includeTaxaAnnotations(Boolean value) {
        myIncludeTaxaAnnotations = new PluginParameter<>(myIncludeTaxaAnnotations, value);
        return this;
    }

    /**
     * Whether to include branch lengths for Newick formatted
     * files.
     *
     * @return Include Branch Lengths
     */
    public Boolean includeBranchLengths() {
        return myIncludeBranchLengths.value();
    }

    /**
     * Set Include Branch Lengths. Whether to include branch
     * lengths for Newick formatted files.
     *
     * @param value Include Branch Lengths
     *
     * @return this plugin
     */
    public ExportPlugin includeBranchLengths(Boolean value) {
        myIncludeBranchLengths = new PluginParameter<>(myIncludeBranchLengths, value);
        return this;
    }

}
