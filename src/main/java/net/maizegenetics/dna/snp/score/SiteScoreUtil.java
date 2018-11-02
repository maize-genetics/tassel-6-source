/*
 *  SiteScoreUtil
 */
package net.maizegenetics.dna.snp.score;

import java.math.BigDecimal;

import java.util.Arrays;

import org.apache.log4j.Logger;

/**
 * @author Terry Casstevens
 */
public class SiteScoreUtil {

    private static final Logger myLogger = Logger.getLogger(SiteScoreUtil.class);

    private static final float[] BYTE_TO_FLOAT = new float[256];
    public static final byte BYTE_REPRESENTING_NAN = (byte) 255;

    static {
        Arrays.fill(BYTE_TO_FLOAT, -1);
        for (int i = 0; i < 255; i++) {
            BYTE_TO_FLOAT[i] = decode(i);
        }
        BYTE_TO_FLOAT[255] = Float.NaN;
    }

    private SiteScoreUtil() {
        // utility
    }

    /**
     * Converts float value to byte version.
     */
    public static byte floatToBytePercentage(float value) {
        if ((value < 0.0) || (value > 1.0)) {
            throw new IllegalArgumentException("SiteScoreUtil: floatToBytePercentage: value: " + value + " must be between 0.0 and 1.0");
        }

        if (Float.isNaN(value)) {
            return BYTE_REPRESENTING_NAN;
        } else {
            return (byte) Math.round(254.0f * value);
        }
    }

    public static byte[] floatToBytePercentage(float[] values) {
        byte[] result = new byte[values.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = floatToBytePercentage(values[i]);
        }
        return result;
    }

    public static byte[][] floatToBytePercentage(float[][] values) {
        byte[][] result = new byte[values.length][values[0].length];
        for (int i = 0; i < result.length; i++) {
            result[i] = floatToBytePercentage(values[i]);
        }
        return result;
    }

    /**
     * Converts byte value to float percentage.
     */
    public static float byteToFloatPercentage(byte value) {
        if (value == BYTE_REPRESENTING_NAN) {
            return Float.NaN;
        }
        return BYTE_TO_FLOAT[value & 0xFF];
    }

    public static float[] byteToFloatPercentage(byte[] values) {
        float[] result = new float[values.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = byteToFloatPercentage(values[i]);
        }
        return result;
    }

    private static float decode(int value) {
        float result = value / 254f;
        BigDecimal bd = new BigDecimal(Float.toString(result));
        bd = bd.setScale(3, BigDecimal.ROUND_HALF_UP);
        return bd.floatValue();
    }

}
