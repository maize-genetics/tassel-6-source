/*
 *  ReadBedfile
 * 
 *  Created on Feb 15, 2017
 */
package net.maizegenetics.dna.snp.io;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.collectingAndThen;

/**
 * @author Terry Casstevens
 */
public class ReadBedfile {

    private static final Logger myLogger = Logger.getLogger(ReadBedfile.class);

    private ReadBedfile() {
        // utility
    }

    public static List<BedFileRange> getRanges(String bedFile) {

        List<BedFileRange> result = new ArrayList<>();

        String line = null;
        try (BufferedReader reader = Utils.getBufferedReader(bedFile)) {
            int lineNum = 1;
            line = reader.readLine();
            while (line != null) {
                String[] tokens = line.trim().split("\t");
                if (tokens.length < 3) {
                    throw new IllegalStateException("getRanges: Expecting at least 3 columns on line: " + lineNum);
                }

                // tokens[0] is chromosome
                // tokens[1] is start postion from bed file.
                // plus one because bed files are 0-base
                int startPos = Integer.parseInt(tokens[1]) + 1;

                // tokens[2] is start postion from bed file.
                // plus one because bed files are 0-base
                int endPos = Integer.parseInt(tokens[2]) + 1;

                // tokens[3] is name from bed file
                String name = null;
                if (tokens.length > 3) {
                    if (tokens[3] == null || tokens[3].isEmpty()) {
                        name = null;
                    } else {
                        name = tokens[3];
                    }
                }

                result.add(new BedFileRange(tokens[0], startPos, endPos, name));

                line = reader.readLine();
                lineNum++;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("getRanges: problem reading: " + bedFile + " line: " + line);
        }

        return result;

    }

    public static RangeSet<Position> getRangesAsPositions(String bedfile) {
        return getRanges(bedfile).stream()
                .map(bedFileRange -> Range.closed(Position.of(bedFileRange.myChrInt, bedFileRange.myStartPos),
                        Position.of(bedFileRange.myChrInt, bedFileRange.myEndPos)))
                .collect(collectingAndThen(Collectors.toSet(), TreeRangeSet::create));
    }

    public static RangeMap<Position, String> getRangesAsPositionMap(String bedfile) {
        TreeRangeMap<Position, String> positionNameRangeMap = TreeRangeMap.create();
        for (BedFileRange bedFileRange : getRanges(bedfile)) {
            positionNameRangeMap.put(Range.closed(Position.of(bedFileRange.myChrInt, bedFileRange.myStartPos),
                    Position.of(bedFileRange.myChrInt, bedFileRange.myEndPos)),
                    bedFileRange.myName);
        }
        return positionNameRangeMap;
    }

    public static class BedFileRange implements Comparable<BedFileRange> {

        private final String myChr;
        private final int myChrInt;
        private final int myStartPos;
        private final int myEndPos;
        private final String myName;

        public BedFileRange(String chr, int startPos, int endPos, String name) {
            myChr = chr;
            int temp;
            try {
                temp = Integer.parseInt(chr);
            } catch (Exception e) {
                temp = -1;
            }
            myChrInt = temp;
            myStartPos = startPos;
            myEndPos = endPos;
            myName = name;
        }

        /**
         * Return chromosome
         *
         * @return chromosome
         */
        public String chr() {
            return myChr;
        }

        /**
         * Returns start position (inclusive)
         *
         * @return start position
         */
        public int start() {
            return myStartPos;
        }

        /**
         * Returns end position (exclusive)
         *
         * @return end position
         */
        public int end() {
            return myEndPos;
        }

        public String name() {
            return myName;
        }

        @Override
        public int compareTo(BedFileRange o) {

            if (myChrInt != -1) {
                if (myChrInt < o.myChrInt) {
                    return -1;
                } else if (myChrInt > o.myChrInt) {
                    return 1;
                }
            } else if (!myChr.equals(o.myChr)) {
                return myChr.compareTo(o.myChr);
            }

            if (myStartPos < o.myStartPos) {
                return -1;
            } else if (myStartPos > o.myStartPos) {
                return 1;
            }

            if (myEndPos < o.myEndPos) {
                return -1;
            } else if (myEndPos > o.myEndPos) {
                return 1;
            } else {
                return 0;
            }

        }

    }

}
