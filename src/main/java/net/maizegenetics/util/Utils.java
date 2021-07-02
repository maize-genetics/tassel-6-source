/*
 * Utils.java
 *
 * Created on May 27, 2003, 2:03 AM
 */
package net.maizegenetics.util;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.net.URL;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 * @author terryc
 */
public final class Utils {

    private static final Logger myLogger = LogManager.getLogger(Utils.class);
    private static Collection<String> myJavaPackages = null;

    private Utils() {
        // Utility Class
    }

    /**
     * Returns the base name of a string delimited with periods (i.e. Java
     * Class).
     *
     * @param str string to parse
     *
     * @return base name
     */
    public static String getBasename(String str) {
        int index = str.lastIndexOf('.');
        index++;
        return str.substring(index);
    }

    /**
     * This returns the filename only. Preceding directories are removed and
     * everything after last . is removed.
     *
     * @param str original filename
     *
     * @return trimmed filename
     */
    public static String getFilename(String str) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = 0;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash + 1;
        } else {
            index = indexBackwardSlash + 1;
        }

        String result = str.substring(index);
        if (result.indexOf('.') > 0) {
            result = result.substring(0, result.indexOf('.'));
        }

        return result;

    }

    /**
     * This returns the filename only. Preceding directories are removed and
     * suffix. If suffix not found, then everything after last . is removed.
     *
     * @param str original filename
     * @param suffix suffix
     *
     * @return trimmed filename
     */
    public static String getFilename(String str, String suffix) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = 0;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash + 1;
        } else {
            index = indexBackwardSlash + 1;
        }

        String result = str.substring(index);
        if ((suffix != null) && (result.lastIndexOf(suffix) > 0)) {
            result = result.substring(0, result.lastIndexOf(suffix));
        }
        if (result.endsWith(".gz")) {
            result = result.substring(0, result.lastIndexOf(".gz"));
        }
        if (result.lastIndexOf('.') > 0) {
            result = result.substring(0, result.lastIndexOf('.'));
        }

        return result;

    }

    /**
     * Returns just the directory with the filename removed.
     *
     * @param str original filename
     *
     * @return directory
     */
    public static String getDirectory(String str) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = -1;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash;
        } else {
            index = indexBackwardSlash;
        }

        if (index == -1) {
            return ".";
        } else {
            return str.substring(0, index);
        }

    }

    /**
     * This returns a set of fully qualified resource names that match the
     * specified filename.
     *
     * @param filename filename
     *
     * @return set of resource names
     */
    public static Set<String> getFullyQualifiedResourceNames(String filename) {

        Set<String> result = new LinkedHashSet<>();

        String classpath = System.getProperty("java.class.path");
        String[] paths = classpath.split(File.pathSeparator);
        for (String path : paths) {

            if (path.trim().length() != 0) {
                File file = new File(path);
                if (file.exists()) {

                    try (ZipFile zFile = new ZipFile(file.getAbsolutePath());) {

                        Enumeration<? extends ZipEntry> entries = zFile.entries();
                        while (entries.hasMoreElements()) {
                            ZipEntry entry = entries.nextElement();
                            if (!entry.isDirectory()) {
                                String name = entry.getName();
                                if (name.endsWith(filename)) {
                                    result.add("/" + name);
                                }
                            }
                        }
                    } catch (Exception e) {
                        myLogger.debug(e.getMessage(), e);
                    }

                }
            }
        }

        return result;

    }

    public static List<String> getFullyQualifiedClassNames(String simpleName) {

        if (myJavaPackages == null) {
            myJavaPackages = getJavaPackages();
        }

        List<String> fqns = new ArrayList<String>();
        for (String aPackage : myJavaPackages) {
            try {
                String fqn = aPackage + "." + simpleName;
                Class.forName(fqn);
                fqns.add(fqn);
            } catch (Throwable e) {
                // Do Nothing
            }
        }
        return fqns;

    }

    public static Collection<String> getJavaPackages() {
        String classpath = System.getProperty("java.class.path");
        return getPackagesFromClassPath(classpath);
    }

    public static Set<String> getPackagesFromClassPath(String classpath) {
        Set<String> packages = new HashSet<String>();
        String[] paths = classpath.split(File.pathSeparator);
        for (String path : paths) {
            if (path.trim().length() == 0) {
                continue;
            } else {
                File file = new File(path);
                if (file.exists()) {
                    String childPath = file.getAbsolutePath();
                    if (childPath.endsWith(".jar")) {
                        packages.addAll(readZipFile(childPath));
                    } else {
                        packages.addAll(readDirectory(childPath));
                    }
                }
            }

        }
        return packages;
    }

    public static Set<String> getTasselClasses() {

        String classpath = System.getProperty("java.class.path");
        String[] paths = classpath.split(File.pathSeparator);
        String tasselPath = null;
        for (String path : paths) {
            if (path.trim().length() != 0) {
                File file = new File(path);
                if (file.exists()) {
                    tasselPath = file.getAbsolutePath();
                    if (tasselPath.endsWith("sTASSEL.jar")) {
                        break;
                    }
                }
            }

        }

        Set<String> classes = new LinkedHashSet<>();
        try (ZipFile zFile = new ZipFile(tasselPath);) {

            Enumeration<? extends ZipEntry> entries = zFile.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                if (!entry.isDirectory()) {
                    String name = entry.getName().replace(File.separator, ".");
                    if ((name.endsWith(".class")) && (!name.contains("$"))) {
                        name = name.substring(0, name.lastIndexOf(".class"));
                        classes.add(name);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return classes;
    }

    public static Set<String> readDirectory(String path) {
        Set<String> packages = new HashSet<String>();
        File file = new File(path);
        int startIndex = path.length() + 1;
        for (File child : file.listFiles()) {
            recursiveRead(child, startIndex, packages);
        }
        return packages;
    }

    public static void recursiveRead(File file, int startIndex, Set<String> packages) {
        if (!file.isFile()) {
            packages.add(file.getAbsolutePath().substring(startIndex).replace(File.separator, "."));
            for (File child : file.listFiles()) {
                recursiveRead(child, startIndex, packages);
            }
        }
    }

    public static Set<String> readZipFile(String path) {
        Set<String> packages = new HashSet<String>();
        try {
            ZipFile zFile = new ZipFile(path);
            Enumeration<? extends ZipEntry> entries = zFile.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                if (!entry.isDirectory()) {
                    String dirName = new File(entry.getName()).getParent();
                    if (dirName != null) {
                        String name = dirName.replace(File.separator, ".");
                        if (name.endsWith(".")) {
                            name = name.substring(0, name.length() - 1);
                        }
                        packages.add(name);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return packages;
    }

    public static String shortenStrLineLen(String str, int preferredLen) {
        return shortenStrLineLen(str, preferredLen, -1);
    }

    /**
     *
     */
    public static String shortenStrLineLen(String str, int preferredLen, int preferredLines) {

        StringBuilder buffer = new StringBuilder();

        int startIndex = 0;
        int endIndex = preferredLen;
        int strLen = str.length();
        int numLines = 0;

        while (startIndex < strLen - 1) {

            int place = str.indexOf(' ', endIndex);
            int newLine = str.indexOf('\n', startIndex);

            if ((newLine != -1) && (newLine < place)) {
                place = newLine;
            }

            String part = null;
            if (place == -1) {
                part = str.substring(startIndex);
                buffer.append(part);
                buffer.append("\n");
                break;
            } else {
                place++;
                part = str.substring(startIndex, place);
                buffer.append(part);
                buffer.append("\n");
                startIndex = place;
                endIndex = place + preferredLen;
            }

            numLines++;

            if ((preferredLines > 0) && (numLines >= preferredLines)) {
                return buffer.toString();
            }

        }

        return buffer.toString();

    }

    /**
     * Adds suffix (i.e. .txt) to end of filename if it's not already there.
     *
     * @param filename filename
     * @param suffix suffix
     *
     * @return filename with suffix
     */
    public static String addSuffixIfNeeded(String filename, String suffix) {

        String temp = filename.toLowerCase();

        if (suffix.charAt(0) != '.') {
            suffix = '.' + suffix;
        }

        if (temp.endsWith(suffix)) {
            return filename;
        } else {
            return filename + suffix;
        }

    }

    public static String addGzSuffixIfNeeded(String filename, String suffix) {
        String gzipSuffix = suffix + ".gz";
        String result = addSuffixIfNeeded(filename, gzipSuffix, new String[]{suffix, gzipSuffix});
        return addSuffixIfNeeded(result, ".gz");
    }

    /**
     * Adds default suffix if not already one of the possible suffixes.
     *
     * @param filename filename
     * @param defaultSuffix default suffix
     * @param possible possible suffixes
     *
     * @return filename with suffix
     */
    public static String addSuffixIfNeeded(String filename, String defaultSuffix, String[] possible) {

        String temp = filename.toLowerCase();

        for (String possible1 : possible) {
            String current = possible1.toLowerCase();
            if (current.charAt(0) != '.') {
                current = '.' + current;
            }
            if (temp.endsWith(current)) {
                return filename;
            }
        }

        if (defaultSuffix.charAt(0) != '.') {
            return filename + '.' + defaultSuffix;
        } else {
            return filename + defaultSuffix;
        }

    }

    public static BufferedReader getBufferedReader(String inSourceName) {
        return getBufferedReader(inSourceName, 8192);
    }

    public static BufferedReader getBufferedReader(String inSourceName, int bufSize) {

        try {
            if (bufSize < 1) {
                return getBufferedReader(inSourceName);
            } else if (inSourceName.startsWith("http")) {
                if (inSourceName.endsWith(".gz")) {
                    return new BufferedReader(new InputStreamReader(new GZIPInputStream((new URL(inSourceName)).openStream(), bufSize)), bufSize);
                } else {
                    return new BufferedReader(new InputStreamReader((new URL(inSourceName)).openStream()), bufSize);
                }
            } else if (inSourceName.endsWith(".gz")) {
                return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inSourceName), bufSize)), bufSize);
            } else {
                return new BufferedReader(new InputStreamReader(new FileInputStream(inSourceName)), bufSize);
            }
        } catch (Exception e) {
            myLogger.error("getBufferedReader: Error getting reader for: " + inSourceName);
            e.printStackTrace();
        }
        return null;
    }

    public static BufferedReader getBufferedReader(File file, int bufSize) {
        return getBufferedReader(file.getAbsolutePath(), bufSize);
    }

    /**
     * Read all lines from a file as a {@code Stream}. Unlike Path.readAllLines,
     * this method does not read all lines into a {@code List}, but instead
     * populates lazily as the stream is consumed.
     *
     * <p>
     * Bytes from the file are decoded into characters using the specified
     * charset and the same line terminators as specified by {@code
     * readAllLines} are supported.
     *
     * <p>
     * After this method returns, then any subsequent I/O exception that occurs
     * while reading from the file or when a malformed or unmappable byte
     * sequence is read, is wrapped in an {@link java.io.UncheckedIOException}
     * that will be thrown from the {@link java.util.stream.Stream} method that
     * caused the read to take place. In case an {@code IOException} is thrown
     * when closing the file, it is also wrapped as an
     * {@code UncheckedIOException}.
     *
     * <p>
     * The returned stream encapsulates a {@link java.io.Reader}. If timely
     * disposal of file system resources is required, the try-with-resources
     * construct should be used to ensure that the stream's
     * {@link java.util.stream.Stream#close close} method is invoked after the
     * stream operations are completed.
     *
     * @param path the path to the file
     *
     * @return the lines from the file as a {@code Stream}
     * @throws IOException       if an I/O error occurs opening the file
     * @throws SecurityException In the case of the default provider, and a
     *                           security manager is installed, the
     *                           {@link SecurityManager#checkRead(String) checkRead} method is invoked to
     *                           check read access to the file
     * @see java.io.BufferedReader#lines()
     * @since 1.8
     */
    public static Stream<String> lines(Path path, int bufSize) throws IOException {
        BufferedReader br = getBufferedReader(path.toString(), bufSize);
        try {
            return br.lines().onClose(asUncheckedRunnable(br));
        } catch (Error | RuntimeException e) {
            try {
                br.close();
            } catch (IOException ex) {
                try {
                    e.addSuppressed(ex);
                } catch (Throwable ignore) {
                }
            }
            throw e;
        }
    }

    /**
     * Convert a Closeable to a Runnable by converting checked IOException to
     * UncheckedIOException
     */
    private static Runnable asUncheckedRunnable(Closeable c) {
        return () -> {
            try {
                c.close();
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        };
    }

    public static BufferedWriter getBufferedWriter(String filename) {
        return getBufferedWriter(filename, false);
    }

    public static BufferedWriter getBufferedWriter(String filename, boolean append) {
        return getBufferedWriter(new File(filename), append);
    }

    public static BufferedWriter getBufferedWriter(File file) {
        return getBufferedWriter(file, false);
    }

    public static BufferedWriter getBufferedWriter(File file, boolean append) {

        String filename = null;

        try {
            filename = file.getCanonicalPath();
            if (filename.endsWith(".gz")) {
                return new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file, append))));
            } else {
                return new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file, append)));
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("getBufferedReader: Error getting reader for: " + filename + "\n" + e.getMessage());
        }

    }

    public static DataOutputStream getDataOutputStream(String filename, int bufSize) {

        try {
            if (filename.endsWith(".gz")) {
                return new DataOutputStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(filename), bufSize)));
            } else {
                return new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename), bufSize));
            }
        } catch (Exception e) {
            myLogger.error("getDataOutputStream: Error getting reader for: " + filename);
            e.printStackTrace();
        }
        return null;

    }

    /**
     * Finds index of Nth occurrence of character in string.
     *
     * @param str string
     * @param match character to match
     * @param n Nth occurrence
     *
     * @return index
     */
    public static int findNthOccurrenceInString(String str, char match, int n) {
        int result = str.indexOf(match);
        while (--n > 0 && result != -1) {
            result = str.indexOf(match, result + 1);
        }
        return result;
    }

    /**
     * Returns max heap size in MB.
     *
     * @return max heap size
     */
    public static long getMaxHeapSizeMB() {
        return Runtime.getRuntime().maxMemory() / 1048576l;
    }

    /**
     * Gets input stream for given file.
     *
     * @param filename file name
     *
     * @return input stream
     */
    public static InputStream getInputStream(String filename) {

        try {
            if (filename.startsWith("http")) {
                if (filename.endsWith(".gz")) {
                    return new GZIPInputStream((new URL(filename)).openStream());
                } else {
                    return (new URL(filename)).openStream();
                }
            } else if (filename.endsWith(".gz")) {
                return new GZIPInputStream(new FileInputStream(filename));
            } else {
                return new FileInputStream(filename);
            }
        } catch (Exception e) {
            myLogger.error("getInputStream: Error getting reader for: " + filename);
            e.printStackTrace();
        }
        return null;
    }

    public static BufferedOutputStream getBufferedOutputStream(String filename) {

        try {
            if (filename.endsWith(".gz")) {
                return new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(filename)));
            } else {
                return new BufferedOutputStream(new FileOutputStream(filename));
            }
        } catch (Exception e) {
            myLogger.error("getOutputStream: Error getting output stream for: " + filename);
            myLogger.debug(e.getMessage(), e);
        }
        return null;
    }

    /**
     * Return number of lines in given file.
     *
     * @param filename file name
     *
     * @return number of lines
     */
    public static int getNumberLines(String filename) {

        InputStream input = getInputStream(filename);
        try {

            byte[] buffer = new byte[1024];
            int result = 0;
            int numChrsRead = 0;
            while ((numChrsRead = input.read(buffer)) != -1) {
                for (int i = 0; i < numChrsRead; ++i) {
                    if (buffer[i] == '\n') {
                        ++result;
                    }
                }
            }
            return result;

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("Utils: getNumberLines: Problem getting number lines: " + filename);
        } finally {
            try {
                input.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    /**
     * Return number of lines in given file that doesn't begin with Hash (#) and
     * isn't blank.
     *
     * @param filename file name
     *
     * @return number of lines
     */
    public static int getNumberLinesNotHashOrBlank(String filename) {

        InputStream input = getInputStream(filename);
        try {

            byte[] buffer = new byte[1024];
            int result = 0;
            boolean isHash = false;
            boolean isBlank = true;
            int numChrsRead = input.read(buffer);
            // Check if first char of file is #
            if (numChrsRead > 0) {
                isHash = buffer[0] == '#';
            }
            while (numChrsRead > 0) {

                for (int i = 0; i < numChrsRead - 1; ++i) {
                    if (buffer[i] == '\n') {
                        if (isHash) {
                            isHash = false;
                        } else if (!isBlank) {
                            ++result;
                        }
                        if (buffer[i + 1] == '#') {
                            isHash = true;
                        }
                        isBlank = true;
                    } else {
                        isBlank = false;
                    }
                }

                if (buffer[numChrsRead - 1] == '\n') {
                    if (isHash) {
                        isHash = false;
                    } else if (!isBlank) {
                        ++result;
                    }
                    numChrsRead = input.read(buffer);
                    if (numChrsRead > 0) {
                        isHash = buffer[0] == '#';
                    }
                    isBlank = true;
                } else {
                    numChrsRead = input.read(buffer);
                }

            }

            return result;

        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("Utils: getNumberLinesNotHashOrBlank: Problem getting number lines: " + filename);
        } finally {
            try {
                input.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    public static String readLineSkipComments(BufferedReader br) throws IOException {
        String s = br.readLine();
        while ((s.startsWith("#"))) {
            s = br.readLine();
        }
        return s;
    }
}
