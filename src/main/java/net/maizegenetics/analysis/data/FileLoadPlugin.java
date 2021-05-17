/*
 * FileLoadPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.data;

import javafx.stage.FileChooser;
import net.maizegenetics.dna.factor.FeatureTable;
import net.maizegenetics.dna.factor.io.BuilderFromHapMap;
import net.maizegenetics.dna.factor.io.BuilderFromHaplotypeVCF;
import net.maizegenetics.dna.snp.io.FilterJSONUtils;
import net.maizegenetics.dna.snp.io.JSONUtils;
import net.maizegenetics.dna.snp.io.LineIndexBuilder;
import net.maizegenetics.dna.snp.io.ReadNumericMarkerUtils;
import net.maizegenetics.gui.FileChooserUtils;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class FileLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FileLoadPlugin.class);

    private PluginParameter<TasselFileType> myFileType = new PluginParameter.Builder<>("format", TasselFileType.Unknown, TasselFileType.class)
            .description("Import file format")
            .objectListSingleSelect()
            .range(TasselFileType.values())
            .build();

    private PluginParameter<Boolean> mySortPositions = new PluginParameter.Builder<>("sortPositions", false, Boolean.class)
            .description("Whether to sort genotype positions if that's possible.")
            .dependentOnParameter(myFileType, new Object[]{TasselFileType.Unknown, TasselFileType.Hapmap, TasselFileType.HapmapDiploid, TasselFileType.VCF, TasselFileType.Plink})
            .build();

    private PluginParameter<Boolean> myKeepDepth = new PluginParameter.Builder<>("keepDepth", true, Boolean.class)
            .description("Whether to keep depth if that's possible.")
            .dependentOnParameter(myFileType, new Object[]{TasselFileType.Unknown, TasselFileType.VCF})
            .build();

    private String[] myOpenFiles = null;
    //private ProjectPcsAndRunModelSelectionPlugin myProjectPcsAndRunModelSelectionPlugin = null;
    //private GOBIIPlugin myGOBIIPlugin = null;
    private final FileChooser myOpenFileChooser;
    private final boolean myHeadless;

    public enum TasselFileType {

        SqrMatrix("Square Matrix"), Sequence("Sequence"), Unknown("Make Best Guess"),
        Fasta("Fasta"), Hapmap("Hapmap"), HapmapLIX("Hapmap LIX"),
        Plink("Plink"), Phenotype("Phenotype"),
        ProjectPCsandRunModelSelection("Project PCs"),
        Phylip_Seq("Phylip (Sequential)"), Phylip_Inter("Phylip (Interleaved)"), Table("Table"),
        Serial("Serial"), HapmapDiploid("Hapmap Diploid"), Text("Text"), VCF("VCF"),
        HaplotypeVCF("Haplotype VCF"), Filter("Filter"), Newick("Newick"),
        NumericGenotype("Numeric Genotype"), TaxaList("Taxa List"), PositionList("Position List"),
        SqrMatrixRaw("Raw MultiBLUP Matrix"), SqrMatrixBin("Binary MultiBLUP Matrix"),
        GOBII("GOBII"), Depth("Depth"), ReferenceProbability("Numeric Genotype"), Report("Report"),
        PlinkPhenotype("Plink Phenotype"), SqrMatrixDARwinDIS("DARwin DIS"), Avro("Avro");

        private final String myText;

        TasselFileType(String text) {
            myText = text;
        }

        @Override
        public String toString() {
            return myText;
        }
    }

    public static final String FILE_EXT_HAPMAP = ".hmp.txt";
    public static final String FILE_EXT_HAPMAP_GZ = ".hmp.txt.gz";
    public static final String FILE_EXT_HAPMAP_GZ_LIX = FILE_EXT_HAPMAP_GZ + LineIndexBuilder.LINE_INDEX_FILE_EXTENSION;
    public static final String FILE_EXT_PLINK_MAP = ".plk.map";
    public static final String FILE_EXT_PLINK_PED = ".plk.ped";
    public static final String FILE_EXT_SERIAL_GZ = ".serial.gz";
    public static final String FILE_EXT_VCF = ".vcf";
    public static final String FILE_EXT_FASTA = ".fasta";
    public static final String FILE_EXT_PHYLIP = ".phy";

    /**
     * Creates a new instance of FileLoadPlugin. This only used by TASSEL GUI to bypass dialog and go straight to file
     * browser. Bypassing the dialog causes it to bypass adding to Data Tree. This constructor tells FileLoadPlugin to
     * add it to the Data Tree Manually.
     */
    public FileLoadPlugin(boolean isInteractive, boolean headless) {
        super(isInteractive);
        if (isInteractive) {
            myOpenFileChooser = new FileChooser();
            myOpenFileChooser.setInitialDirectory(new File(TasselPrefs.getOpenDir()));
        } else {
            myOpenFileChooser = null;
        }
        myHeadless = headless;
    }

    /**
     * Creates a new instance of FileLoadPlugin.
     */
    public FileLoadPlugin(boolean isInteractive) {
        this(isInteractive, false);
    }

    public FileLoadPlugin() {
        this(false, false);
    }

    public Object run(String filename) {
        return runPlugin(filename).getData(0).getData();
    }

    public DataSet runPlugin(String filename) {
        setTheFileType(TasselFileType.Unknown);
        setOpenFiles(filename);
        return performFunction(null);
    }

    public FeatureTable read(String filename) {
        return (FeatureTable) runPlugin(filename).getData(0).getData();
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        List<TasselFileType> temp = new ArrayList<>(Arrays.asList(TasselFileType.Unknown,
                TasselFileType.Hapmap,
                TasselFileType.VCF,
                TasselFileType.Plink,
                TasselFileType.Sequence,
                TasselFileType.Fasta,
                TasselFileType.SqrMatrix,
                TasselFileType.Table));
        myFileType = new PluginParameter<>(myFileType, temp);

        if (!isInteractive() && (fileType() == null || fileType() == TasselFileType.Unknown) && myFileType.hasPossibleValues()) {
            fileType(TasselFileType.Unknown);
        }

    }

    @Override
    public DataSet processData(DataSet input) {

        myWasCancelled = true;

        if (isInteractive()) {

            setOpenFiles(getOpenFilesByChooser());

        }

        if ((myOpenFiles == null) || (myOpenFiles.length == 0)) {
            return null;
        }

        List<DataSet> result = new ArrayList<>();
        ArrayList<String> alreadyLoaded = new ArrayList<>();
        for (String myOpenFile : myOpenFiles) {

            if (alreadyLoaded.contains(myOpenFile)) {
                continue;
            }

            LocalDateTime time = LocalDateTime.now();
            String timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            myLogger.info("Start Loading File: " + myOpenFile + " time: " + timeStr);

            DataSet tds;

            if (fileType() == TasselFileType.Unknown) {
                if (myOpenFile.endsWith(FILE_EXT_HAPMAP_GZ)) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                    alreadyLoaded.add(myOpenFile);
                    tds = processDatum(myOpenFile, TasselFileType.Hapmap);
                } else if (myOpenFile.endsWith(FILE_EXT_HAPMAP)) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                    alreadyLoaded.add(myOpenFile);
                    tds = processDatum(myOpenFile, TasselFileType.Hapmap);
                } else if ((myOpenFile.endsWith(FILE_EXT_VCF) || myOpenFile.endsWith(FILE_EXT_VCF + ".gz")) && myOpenFile.toUpperCase().contains("HAPLOTYPE")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.HaplotypeVCF);
                    alreadyLoaded.add(myOpenFile);
                    tds = processDatum(myOpenFile, TasselFileType.HaplotypeVCF);
                } else if (myOpenFile.endsWith(FILE_EXT_VCF) || myOpenFile.endsWith(FILE_EXT_VCF + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.VCF);
                    alreadyLoaded.add(myOpenFile);
                    tds = processDatum(myOpenFile, TasselFileType.VCF);
                } else if (myOpenFile.endsWith(FILE_EXT_PHYLIP) || myOpenFile.endsWith(FILE_EXT_PHYLIP + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Sequence);
                    alreadyLoaded.add(myOpenFile);
                    tds = processDatum(myOpenFile, TasselFileType.Sequence);
                } else if (myOpenFile.endsWith(FILE_EXT_FASTA) || myOpenFile.endsWith(FILE_EXT_FASTA + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Fasta);
                    alreadyLoaded.add(myOpenFile);
                    tds = processDatum(myOpenFile, TasselFileType.Fasta);
                } else {
                    alreadyLoaded.add(myOpenFile);
                    tds = guessAtUnknowns(myOpenFile);
                }
            } else {
                alreadyLoaded.add(myOpenFile);
                tds = processDatum(myOpenFile, fileType());
            }

            time = LocalDateTime.now();
            timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            if (tds != null) {
                myLogger.info("Finished Loading File: " + myOpenFile + " time: " + timeStr);
                GenotypeSummaryPlugin.printSimpleSummary(tds);
                myWasCancelled = false;
                result.add(tds);
                if (myHeadless) {
                    fireDataSetReturned(new PluginEvent(tds, FileLoadPlugin.class));
                }
            } else {
                myLogger.info("Nothing Loaded for File: " + myOpenFile + " time: " + timeStr);
            }

        }

        return DataSet.getDataSet(result, this);

    }

    private DataSet guessAtUnknowns(String filename) {

        TasselFileType guess = TasselFileType.Table;

        try (BufferedReader br = Utils.getBufferedReader(filename)) {

            String line1 = br.readLine();
            while (line1 != null) {
                line1 = line1.trim();
                if (!line1.isEmpty()) {
                    break;
                }
                line1 = br.readLine();
            }
            if (line1 == null) {
                throw new IllegalArgumentException("FileLoadPlugin: guessAtUnknowns: File is empty: " + filename);
            }
            String[] sval1 = line1.split("\\s");
            String line2 = br.readLine().trim();
            String[] sval2 = line2.split("\\s");
            if (line1.startsWith("{")) {
                String temp;
                if (sval1.length > 1) {
                    temp = sval1[1];
                } else {
                    temp = line2;
                }
                if (temp.startsWith("\"TaxaList\"")) {
                    guess = TasselFileType.TaxaList;
                } else if (temp.startsWith("\"PositionList\"")) {
                    guess = TasselFileType.PositionList;
                } else if (temp.startsWith("\"Filter\"")) {
                    guess = TasselFileType.Filter;
                }
            } else if (line1.startsWith("##")) {
                String matrixStr = "##" + DistanceMatrixBuilder.MATRIX_TYPE;
                if (line1.startsWith(matrixStr) || line2.startsWith(matrixStr)) {
                    guess = TasselFileType.SqrMatrix;
                } else {
                    String line = br.readLine();
                    while ((line != null) && (line.startsWith("##"))) {
                        if (line.startsWith(matrixStr)) {
                            guess = TasselFileType.SqrMatrix;
                            break;
                        }
                        line = br.readLine();
                    }
                }
            } else if (line1.startsWith("<") || line1.startsWith("#")) {
                boolean isTrait = false;
                boolean isMarker = false;
                boolean isNumeric = false;
                Pattern tagPattern = Pattern.compile("[<>\\s]+");
                String[] info1 = tagPattern.split(line1);
                String[] info2 = tagPattern.split(line2);
                if (info1.length > 1) {
                    if (info1[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info1[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info1[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    } else if (info1[1].toUpperCase().startsWith("PHENO")) {
                        isTrait = true;
                    }
                }
                if (info2.length > 1) {
                    if (info2[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info2[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info2[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    }
                } else {
                    guess = null;
                    String inline = br.readLine();
                    while (guess == null && inline != null && (inline.startsWith("#") || inline.startsWith("<"))) {
                        if (inline.startsWith("<")) {
                            String[] info = tagPattern.split(inline);
                            if (info[1].toUpperCase().startsWith("MARKER")) {
                                isMarker = true;
                            } else if (info[1].toUpperCase().startsWith("TRAIT")) {
                                isTrait = true;
                            } else if (info[1].toUpperCase().startsWith("NUMER")) {
                                isNumeric = true;
                            }
                        }
                    }
                }
                if (isTrait) {
                    guess = TasselFileType.Phenotype;
                } else if (isMarker && isNumeric) {
                    guess = TasselFileType.NumericGenotype;
                } else {
                    myLogger.warn("Line1: " + line1);
                    myLogger.warn("Line2: " + line2);
                    throw new IOException("Improperly formatted header. Data will not be imported for file: " + filename);
                }
            } else if ((line1.startsWith(">")) || (line1.startsWith(";"))) {
                guess = TasselFileType.Fasta;
            } else if (sval1.length == 1) {
                guess = TasselFileType.SqrMatrix;
            } else if ((line1.startsWith("#Nexus")) || (line1.startsWith("#NEXUS")) || (line1.startsWith("CLUSTAL"))
                    || ((sval1.length == 2) && (sval2.length == 2))) {
                guess = TasselFileType.Sequence;
            }

            myLogger.info("guessAtUnknowns: type: " + guess);
            return processDatum(filename, guess);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("FileLoadPlugin: Problem loading file: " + filename + ".  Error: " + e.getMessage());
        }

    }

    private DataSet processDatum(String inFile, TasselFileType theFT) {

        Object result;
        String suffix = null;
        try {
            switch (theFT) {
                case Hapmap: {
                    suffix = FILE_EXT_HAPMAP;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_HAPMAP_GZ;
                    }
                    result = BuilderFromHapMap.getBuilder(inFile, this).build();
                    break;
                }
                case HaplotypeVCF: {
                    suffix = FILE_EXT_VCF;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_VCF + ".gz";
                    }
                    result = new BuilderFromHaplotypeVCF().read(inFile);
                    break;
                }
                case VCF: {
                    suffix = FILE_EXT_VCF;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_VCF + ".gz";
                    }
                    //result = ImportUtils.readFromVCF(inFile, this, keepDepth(), sortPositions());
                    result = null;
                    break;
                }
                case Sequence: {
                    //result = ReadSequenceAlignmentUtils.readBasicAlignments(inFile, 40);
                    result = null;
                    break;
                }
                case Fasta: {
                    //result = ImportUtils.readFasta(inFile);
                    result = null;
                    break;
                }
                case SqrMatrix: {
                    result = ReadDistanceMatrix.readDistanceMatrix(inFile);
                    break;
                }
                case SqrMatrixBin: {
                    result = ReadDistanceMatrix.readBinMultiBlupMatrix(inFile);
                    break;
                }
                case Phenotype: {
                    List<Phenotype> phenotypes = new PhenotypeBuilder().fromFile(inFile).build();
                    if (phenotypes.size() != 1) {
                        throw new IllegalStateException("FileLoadPlugin: processDatum: problem loading phenotype file: " + inFile);
                    }
                    result = phenotypes.get(0);
                    break;
                }
                case NumericGenotype: {
                    result = ReadNumericMarkerUtils.readNumericMarkerFile(inFile);
                    break;
                }
                case TaxaList: {
                    result = JSONUtils.importTaxaListFromJSON(inFile);
                    break;
                }
                case PositionList: {
                    result = JSONUtils.importPositionListFromJSON(inFile);
                    break;
                }
                case Table: {
                    result = TableReportUtils.readDelimitedTableReport(inFile, "\t");
                    break;
                }
                case Filter: {
                    result = FilterJSONUtils.importJSONToFilter(inFile);
                    break;
                }
                default: {
                    throw new IllegalStateException("Unknown Format: " + theFT + ".\n  Please check file format or select specific format.");
                }
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("Problem loading file: " + inFile + ".\n  Error: " + e.getMessage());
        }

        if (result != null) {
            Datum td = new Datum(Utils.getFilename(inFile, suffix), result, null);
            return new DataSet(td, this);
        }
        return null;

    }

    /**
     * Provides a open file chooser that remember the last location something was opened from
     */
    private File[] getOpenFilesByChooser() {
        List<File> temp = FileChooserUtils.multipleFiles();
        if (temp == null) return null;
        return temp.toArray(new File[temp.size()]);
    }

    public String[] getOpenFiles() {
        return myOpenFiles;
    }

    public void setOpenFiles(File[] openFiles) {

        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
            return;
        }

        myOpenFiles = new String[openFiles.length];
        for (int i = 0; i < openFiles.length; i++) {
            myOpenFiles[i] = openFiles[i].getPath();
        }

    }

    public void setOpenFiles(String openFile) {
        if ((openFile == null) || openFile.isEmpty()) {
            myOpenFiles = null;
        } else {
            myOpenFiles = new String[]{openFile};
        }
    }

    public void setOpenFiles(String[] openFiles) {
        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
        } else {
            myOpenFiles = openFiles;
        }
    }

    /**
     * Export file format (Default format depends on data being exported)
     *
     * @return Format
     */
    public TasselFileType fileType() {
        return myFileType.value();
    }

    /**
     * Set Format. Export file format (Default format depends on data being exported)
     *
     * @param value Format
     *
     * @return this plugin
     */
    public FileLoadPlugin fileType(TasselFileType value) {
        myFileType = new PluginParameter<>(myFileType, value);
        return this;
    }

    public TasselFileType getTheFileType() {
        return fileType();
    }

    public void setTheFileType(TasselFileType theFileType) {
        fileType(theFileType);
    }

    /**
     * Whether to sort genotype positions if that's possible.
     *
     * @return Sort Positions
     */
    public Boolean sortPositions() {
        return mySortPositions.value();
    }

    /**
     * Set Sort Positions. Whether to sort genotype positions if that's possible.
     *
     * @param value Sort Positions
     *
     * @return this plugin
     */
    public FileLoadPlugin sortPositions(Boolean value) {
        mySortPositions = new PluginParameter<>(mySortPositions, value);
        return this;
    }

    /**
     * Whether to keep depth if that's possible.
     *
     * @return Keep Depth
     */
    public Boolean keepDepth() {
        return myKeepDepth.value();
    }

    /**
     * Set Keep Depth. Whether to keep depth if that's possible.
     *
     * @param value Keep Depth
     *
     * @return this plugin
     */
    public FileLoadPlugin keepDepth(Boolean value) {
        myKeepDepth = new PluginParameter<>(myKeepDepth, value);
        return this;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public String icon() {
        return "/images/LoadFile.gif";
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Open As...";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Open data from filesystem.";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load";
    }
}
