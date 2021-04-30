/*
 * ExportMultiplePlugin.java
 *
 * Created on December 21, 2010
 *
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class ExportMultiplePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(ExportMultiplePlugin.class);
    private FileLoadPlugin.TasselFileType[] myFileTypes = null;
    private String[] mySaveFiles = null;
    private final ExportPlugin myExportPlugin;

    /**
     * Creates a new instance of ExportMultiplePlugin
     */
    public ExportMultiplePlugin() {
        super(false);
        myExportPlugin = new ExportPlugin(false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        List data = input.getDataSet();

        int numSaveFiles = 0;
        if (mySaveFiles != null) {
            numSaveFiles = mySaveFiles.length;
        }

        if ((numSaveFiles != 0) && (numSaveFiles != 1) && (numSaveFiles != data.size())) {
            throw new IllegalStateException("ExportMultiplePlugin: performFunction: number of save files should be either 0, 1 or number of input data sets.");
        }

        if ((myFileTypes != null) && (myFileTypes.length != 0)) {
            if ((myFileTypes.length != 1) && (myFileTypes.length != data.size())) {
                throw new IllegalStateException("ExportMultiplePlugin: performFunction: number of files types should be either 0, 1 or number of input data sets.");
            }
        }

        for (int i = 0, n = data.size(); i < n; i++) {

            Datum datum = (Datum) data.get(i);
            DataSet current = new DataSet(datum, input.getCreator());

            if (numSaveFiles == 0) {
                myExportPlugin.saveFile(datum.getName());
            } else if (numSaveFiles == 1) {
                String temp;
                if (data.size() == 1) {
                    temp = mySaveFiles[0];
                } else {
                    String filename = Utils.getFilename(mySaveFiles[0]);
                    temp = filename.replaceFirst("\\.", (i + 1) + ".");
                    if (temp.length() == filename.length()) {
                        temp = filename + (i + 1);
                    }
                    String directory = Utils.getDirectory(mySaveFiles[0]);
                    if (!directory.equals(".")) {
                        File dir = new File(directory);
                        if (!dir.exists()) {
                            dir.mkdirs();
                        }
                        temp = directory + "/" + temp;
                    }
                }
                myExportPlugin.saveFile(temp);
            } else {
                myExportPlugin.saveFile(mySaveFiles[i]);
            }

            if ((myFileTypes == null) || (myFileTypes.length == 0)) {
                myExportPlugin.fileType(null);
            } else if (myFileTypes.length == 1) {
                myExportPlugin.fileType(myFileTypes[0]);
            } else {
                myExportPlugin.fileType(myFileTypes[i]);
            }

            myExportPlugin.performFunction(current);

        }

        fireProgress(100);
        return null;

    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Export Multiple";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Export multiple data sets to files.";
    }

    public String[] getSaveFiles() {
        return mySaveFiles;
    }

    public void setSaveFiles(String[] saveFiles) {
        mySaveFiles = saveFiles;
    }

    public void setAlignmentFileTypes(FileLoadPlugin.TasselFileType[] types) {
        myFileTypes = types;
    }

    public void setAlignmentFileType(FileLoadPlugin.TasselFileType type) {
        myFileTypes = new FileLoadPlugin.TasselFileType[]{type};
    }

    public void setIncludeAnnotations(boolean include) {
        myExportPlugin.includeTaxaAnnotations(include);
    }

    public void setIncludeDepth(boolean include) {
        myExportPlugin.keepDepth(include);
    }

}
