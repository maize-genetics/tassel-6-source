/*
 * Sizeof.java
 *
 * Created on October 10, 2004, 2:08 AM
 */
package net.maizegenetics.util;

import java.text.NumberFormat;

/**
 *
 * @author Terry Casstevens
 */
public class Sizeof {

    public static void main(String[] args) throws Exception {
        printMemoryUse();
    }

    public static void printMemoryUse() {

        NumberFormat format = NumberFormat.getInstance();

        try {
            System.out.println("-------------------------------");
            System.out.print("Current Heap Size: ");
            long current = Sizeof.getMemoryUse() / 1048576l;
            String currentStr = format.format(current);
            System.out.print(currentStr);
            System.out.println(" MB");
            System.out.print("Max Available Heap: ");
            System.out.print(Utils.getMaxHeapSizeMB());
            System.out.println(" MB");
            System.out.println("-------------------------------");
        } catch (Exception e) {
            System.out.println("Problem getting heap size: " + e.getMessage());
        }

    }

    public static long getMemoryUse() throws Exception {

        // Warm up all classes/methods we will use
        runGC();
        usedMemory();

        runGC();
        return usedMemory();

    }

    private static void runGC() throws Exception {
        // It helps to call Runtime.gc()
        // using several method calls:
        for (int i = 0; i < 4; i++) {
            _runGC();
        }
    }

    private static void _runGC() throws Exception {
        long usedMem1 = usedMemory(), usedMem2 = Long.MAX_VALUE;
        for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++i) {
            s_runtime.runFinalization();
            s_runtime.gc();
            Thread.currentThread().yield();

            usedMem2 = usedMem1;
            usedMem1 = usedMemory();
        }
    }

    private static long usedMemory() {
        return s_runtime.totalMemory() - s_runtime.freeMemory();
    }

    private static final Runtime s_runtime = Runtime.getRuntime();

} // End of class
