/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
// Title:      TASSELMainApp
// Version:
// Copyright:  Copyright (c) 1998
// Author:     Ed Buckler
package net.maizegenetics.tassel;

import javafx.application.Application;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.LoggingUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

//import net.maizegenetics.pipeline.TasselPipeline;

public class TASSELMainApp {

    private static final Logger myLogger = LogManager.getLogger(TASSELMainApp.class);

    public static final String version = "6.0.0";
    public static final String versionDate = "November 1, 2018";

    private TASSELMainApp() {
    }

    public static void main(String[] args) {
        try {

            TasselPrefs.setPersistPreferences(true);
            LoggingUtils.setupLogging();

            Application.launch(TASSELGUI.class);

            //if (args.length > 0) {
            //    new TasselPipeline(args, frame).run();
            //}

        } catch (OutOfMemoryError e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Out of Memory: \n");
            long heapMaxSize = Runtime.getRuntime().maxMemory() / 1048576l;
            builder.append("Current Max Heap Size: ");
            builder.append(heapMaxSize);
            builder.append(" Mb\n");
            builder.append("Use -Xmx option in start_tassel.pl or start_tassel.bat\n");
            builder.append("to increase heap size.");
            builder.append(" Included with tassel standalone zip.");
            myLogger.error(builder.toString());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
