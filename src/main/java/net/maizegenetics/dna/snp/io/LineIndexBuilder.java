/*
 *  LineIndexBuilder
 * 
 *  Created on Aug 29, 2015
 */
package net.maizegenetics.dna.snp.io;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.util.LittleEndianInputStream;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.tribble.util.ParsingUtils;
import java.io.File;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.List;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class LineIndexBuilder {

    private static final Logger myLogger = Logger.getLogger(LineIndexBuilder.class);

    public static final String LINE_INDEX_FILE_EXTENSION = ".lix";
    private static final byte[] MAGIC = {'L', 'I', 'X', 1};
    public static final int MAGIC_NUMBER;

    static {
        final ByteBuffer bb = ByteBuffer.allocate(MAGIC.length);
        bb.put(MAGIC);
        bb.flip();
        MAGIC_NUMBER = bb.order(ByteOrder.LITTLE_ENDIAN).getInt();
    }

    private char myCommentChar = '#';
    private int myNumHeaderLinesToSkip = 1;
    private int myNumColumnsPerRowToKeepInIndex = 0;
    private final String myFileToIndex;
    private final String myIdxFilename;

    public LineIndexBuilder(String fileToIndex) {
        myFileToIndex = fileToIndex;
        myIdxFilename = ParsingUtils.appendToPath(fileToIndex, LINE_INDEX_FILE_EXTENSION);
    }

    public LineIndexBuilder commentChar(char commentChar) {
        myCommentChar = commentChar;
        return this;
    }

    public LineIndexBuilder numHeaderLinesToSkip(int numLines) {
        myNumHeaderLinesToSkip = numLines;
        return this;
    }

    public LineIndexBuilder numColumnsPerRowToKeepInIndex(int numColumns) {
        myNumColumnsPerRowToKeepInIndex = numColumns;
        return this;
    }

    public static Tuple<LineIndex, String[]> readIndex(String idxFilename) {

        try (LittleEndianInputStream input = new LittleEndianInputStream(new BlockCompressedInputStream(new File(idxFilename)))) {

            int magic = input.readInt();
            int commentChar = input.readInt();
            int numHeaderLinesToSkip = input.readInt();
            int numLinesPerInterval = input.readInt();

            int numRowsWithSavedColumns = input.readInt();
            String[] savedColumns = new String[numRowsWithSavedColumns];
            for (int i = 0; i < numRowsWithSavedColumns; i++) {
                savedColumns[i] = input.readString();
            }

            int numVirtualFileOffsets = input.readInt();
            long[] virtualFileOffsets = new long[numVirtualFileOffsets];
            for (int i = 0; i < numVirtualFileOffsets; i++) {
                virtualFileOffsets[i] = input.readLong();
            }

            LineIndex index = new LineIndex(magic, (char) commentChar, numHeaderLinesToSkip, numLinesPerInterval, virtualFileOffsets);

            return new Tuple<>(index, savedColumns);

        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
    }

    public void build() {

        try (LittleEndianOutputStream output = new LittleEndianOutputStream(new BlockCompressedOutputStream(myIdxFilename))) {

            output.writeInt(MAGIC_NUMBER);
            output.writeInt(myCommentChar);
            output.writeInt(myNumHeaderLinesToSkip);
            output.writeInt(LineIndex.NUM_LINES_PER_INTERVAL);

            try (BlockCompressedInputStream input = new BlockCompressedInputStream(new File(myFileToIndex))) {

                boolean notFinished = true;
                int linesSkipped = 0;
                while (notFinished) {
                    String str = input.readLine();
                    if (!(str.charAt(0) == myCommentChar)) {
                        linesSkipped++;
                    }
                    if (myNumHeaderLinesToSkip <= linesSkipped) {
                        notFinished = false;
                    }
                }

                List<Long> virtualFileOffsets = new ArrayList<>();
                List<String> beginningColumnsPerRow = new ArrayList<>();
                notFinished = true;
                while (notFinished) {
                    virtualFileOffsets.add(input.getFilePointer());
                    for (int i = 0; i < LineIndex.NUM_LINES_PER_INTERVAL; i++) {
                        String current = input.readLine();

                        if (current == null) {
                            notFinished = false;
                            break;
                        }

                        if (myNumColumnsPerRowToKeepInIndex > 0) {
                            int n = myNumColumnsPerRowToKeepInIndex - 1;
                            int pos = current.indexOf('\t');
                            while (n-- > 0 && pos != -1) {
                                pos = current.indexOf('\t', pos + 1);
                            }
                            if (pos == -1) {
                                throw new IllegalStateException("LineIndexBuilder: build: " + myNumColumnsPerRowToKeepInIndex + " columns not found.");
                            }
                            beginningColumnsPerRow.add(current.substring(0, pos));
                        }

                    }
                }

                output.writeInt(beginningColumnsPerRow.size());

                for (String current : beginningColumnsPerRow) {
                    output.writeString(current);
                }

                output.writeInt(virtualFileOffsets.size());

                for (Long current : virtualFileOffsets) {
                    output.writeLong(current);
                }

            } catch (SAMFormatException se) {
                myLogger.debug(se.getMessage(), se);
                throw new IllegalStateException("LineIndexBuilder: build: this file is not bgzipped: " + myFileToIndex + ": " + se.getMessage());
            } catch (Exception ex) {
                myLogger.debug(ex.getMessage(), ex);
                throw new IllegalStateException("LineIndexBuilder: build: problem creating index for file: " + myFileToIndex + ": " + ex.getMessage());
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("LineIndexBuilder: build: problem writing index file: " + myIdxFilename + ": " + e.getMessage());
        }

    }

    public static void buildHapmapIndex(String filename) {
        new LineIndexBuilder(filename)
                .commentChar('#')
                .numHeaderLinesToSkip(1)
                .numColumnsPerRowToKeepInIndex(11)
                .build();
    }

}
