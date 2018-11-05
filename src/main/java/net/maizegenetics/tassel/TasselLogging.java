/*
 *  TasselLogging
 */
package net.maizegenetics.tassel;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.TextArea;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.FlowPane;
import javafx.stage.FileChooser;
import javafx.stage.Modality;
import javafx.stage.Stage;
import net.maizegenetics.gui.AlertUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class TasselLogging extends AbstractPlugin {

    private static TasselLogging myInstance = null;
    private static final Logger myLogger = Logger.getLogger(TasselLogging.class);

    private final Stage myDialog = new Stage();
    private final TextArea myTextArea = new TextArea();
    private final TextAreaOutputStream myTextAreaOutputStream = new TextAreaOutputStream(myTextArea);
    private final PrintStream myPrintStream = new PrintStream(myTextAreaOutputStream);

    private TasselLogging() {
        super(true);
        myDialog.setResizable(true);
        myDialog.initModality(Modality.NONE);
        myDialog.setTitle("TASSEL Logging");
        createDialog();
        basicLoggingInfo();
        LoggingUtils.setupLogging(myPrintStream);
        myTextAreaOutputStream.clear();
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging();
        }
    }

    public static TasselLogging getInstance() {
        if (myInstance == null) {
            myInstance = new TasselLogging();
        }
        return myInstance;
    }

    public static void closeInstance() {
        if (myInstance != null) {
            myInstance.close();
        }
    }

    private void close() {
        TasselPrefs.putLogXDim((int) myDialog.getWidth());
        TasselPrefs.putLogYDim((int) myDialog.getHeight());
        myDialog.close();
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging();
        }
    }

    public static void updateLoggingLocation() {
        if (myInstance != null) {
            myInstance.updateLogging();
        }
    }

    private void updateLogging() {
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging();
        } else {
            LoggingUtils.setupDebugLogging(myPrintStream);
        }
    }

    private void createDialog() {

        BorderPane main = new BorderPane();
        Scene scene = new Scene(main);

        int x = TasselPrefs.getLogXDim();
        int y = TasselPrefs.getLogYDim();
        if ((x < 50) || (y < 50)) {
            myDialog.setWidth(500.0);
            myDialog.setHeight(400.0);
        } else {
            myDialog.setWidth(x);
            myDialog.setHeight(y);
        }

        myTextArea.setWrapText(true);
        myTextArea.setPadding(new Insets(10.0));
        myTextArea.setEditable(false);

        final CheckBox isDebug = new CheckBox("Debug Level");
        isDebug.setSelected(TasselPrefs.getLogDebug());
        isDebug.setTooltip(new Tooltip("Set to show Debug Logging Messages"));
        isDebug.setOnAction(event -> {
            boolean debugMode = isDebug.isSelected();
            isDebug.setSelected(debugMode);
            TasselPrefs.putLogDebug(debugMode);
            LoggingUtils.setupLogging(myPrintStream);
        });

        Button closeButton = new Button("Close");
        closeButton.setOnAction(event -> close());

        Button clearButton = new Button("Clear");
        clearButton.setOnAction(event -> myTextAreaOutputStream.clear());

        Button saveButton = new Button("Save");
        saveButton.setOnAction(event -> {
            FileChooser chooser = new FileChooser();
            File theFile = chooser.showOpenDialog(myDialog);
            if (theFile != null) {
                try (BufferedWriter writer = Utils.getBufferedWriter(theFile)) {
                    writer.write(myTextArea.getText());
                } catch (Exception ex) {
                    AlertUtils.showError(ex.getMessage());
                }
            }
        });

        FlowPane pnlButtons = new FlowPane();
        pnlButtons.setAlignment(Pos.CENTER);
        pnlButtons.getChildren().add(closeButton);
        pnlButtons.getChildren().add(clearButton);
        pnlButtons.getChildren().add(saveButton);
        pnlButtons.getChildren().add(isDebug);

        main.setCenter(new ScrollPane(myTextArea));
        main.setBottom(pnlButtons);

        myDialog.setResizable(true);

        myDialog.setScene(scene);

    }

    @Override
    public DataSet processData(DataSet input) {
        LoggingUtils.setupLogging(myPrintStream);
        myDialog.show();
        return null;
    }

    @Override
    public String icon() {
        return "/net/maizegenetics/analysis/images/log.gif";
    }

    @Override
    public String getButtonName() {
        return "Logging";
    }

    @Override
    public String getToolTipText() {
        return "Logging";
    }

    public static void basicLoggingInfo() {
        myLogger.info("Tassel Version: " + TASSELMainApp.version + "  Date: " + TASSELMainApp.versionDate);
        myLogger.info("Max Available Memory Reported by JVM: " + Utils.getMaxHeapSizeMB() + " MB");
        myLogger.info("Java Version: " + System.getProperty("java.version"));
        myLogger.info("OS: " + System.getProperty("os.name"));
        myLogger.info("Number of Processors: " + Runtime.getRuntime().availableProcessors());
    }

    class TextAreaOutputStream extends OutputStream {

        private final byte[] myByteArray = new byte[1];
        private TextAppender myTextAppender;

        public TextAreaOutputStream(TextArea textArea) {
            myTextAppender = new TextAppender(textArea);
        }

        public synchronized void clear() {
            if (myTextAppender != null) {
                myTextAppender.clear();
                basicLoggingInfo();
            }
        }

        @Override
        public synchronized void close() {
            myTextAppender = null;
        }

        @Override
        public synchronized void flush() {
        }

        @Override
        public synchronized void write(int val) {
            myByteArray[0] = (byte) val;
            write(myByteArray, 0, 1);
        }

        @Override
        public synchronized void write(byte[] ba) {
            write(ba, 0, ba.length);
        }

        @Override
        public synchronized void write(byte[] ba, int str, int len) {
            if (myTextAppender != null) {
                myTextAppender.append(bytesToString(ba, str, len));
            }
        }

    }

    static private String bytesToString(byte[] ba, int str, int len) {
        try {
            return new String(ba, str, len, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            return new String(ba, str, len);
        }
    }

    static class TextAppender implements Runnable {

        private static final int MAX_NUM_LINES = 1000;

        private final TextArea myTextArea;
        private final LinkedList<Integer> myLineLengths = new LinkedList<>();
        private final List<String> myBufferedText = new ArrayList<>();

        private int myCurrentLineLength = 0;
        private boolean myClearAllText = false;
        private boolean myIfNoTextQueued = true;

        TextAppender(TextArea textArea) {
            myTextArea = textArea;
        }

        synchronized void append(String val) {
            myBufferedText.add(val);
            if (myIfNoTextQueued) {
                myIfNoTextQueued = false;
                //EventQueue.invokeLater(this);
            }
        }

        synchronized void clear() {
            myClearAllText = true;
            myCurrentLineLength = 0;
            myLineLengths.clear();
            myBufferedText.clear();
            if (myIfNoTextQueued) {
                myIfNoTextQueued = false;
                //EventQueue.invokeLater(this);
            }
        }

        @Override
        public synchronized void run() {
            if (myClearAllText) {
                myTextArea.setText(null);
            }
            for (String val : myBufferedText) {
                myCurrentLineLength += val.length();
                if (val.endsWith(END_OF_LINE) || val.endsWith(SYSTEM_END_OF_LINE)) {
                    if (myLineLengths.size() >= MAX_NUM_LINES) {
                        myTextArea.clear();
                    }
                    myLineLengths.addLast(myCurrentLineLength);
                    myCurrentLineLength = 0;
                }
                myTextArea.appendText(val);
            }
            myBufferedText.clear();
            myClearAllText = false;
            myIfNoTextQueued = true;
        }

        static private final String END_OF_LINE = "\n";
        static private final String SYSTEM_END_OF_LINE = System.getProperty("line.separator", END_OF_LINE);
    }

}
