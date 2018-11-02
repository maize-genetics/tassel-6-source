// FormattedOutput.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.util;

import java.io.IOException;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.Locale;

/**
 * tools to simplify formatted output to a stream
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class FormattedOutput {

    /**
     * create instance of this class (note that there is no public constructor
     * as this class is a singleton)
     */
    public synchronized static FormattedOutput getInstance() {
        if (singleton == null) {
            singleton = new FormattedOutput();
        }

        return singleton;
    }

    /**
     * print decimal number with a prespecified number of digits after the point
     *
     * @param out output stream
     * @param number to be printed
     * @param width number of fraction digits
     *
     * @return length of the string printed
     */
    public int displayDecimal(Writer out, double number, int width) throws IOException {
        String s = getDecimalString(number, width);

        out.write(s);

        return s.length();
    }

    /**
     * Returns a decimal string representation of a number with constrained
     * width.
     */
    public synchronized String getDecimalString(double number, int width) {
        nf.setMinimumFractionDigits(width);
        nf.setMaximumFractionDigits(width);

        return nf.format(number);
    }

    private static final double round(double number, int sf) {
        double decimals = Math.floor(Math.log(number) / Math.log(10.0));
        double power = Math.pow(10, decimals - sf + 1);
        number /= power;
        number = Math.round(number);
        number *= power;
        return number;
    }

    public String getSFString(double[] numbers, int sf, String delimiter) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < numbers.length; i++) {
            sb.append(getSFString(numbers[i], sf));
            if (i != numbers.length - 1) {
                sb.append(delimiter);
            }
        }
        return sb.toString();
    }

    public String getSFString(double number, int sf) {
        //This seems to give the correct number of sf
        int decimals = (int) Math.ceil(Math.log(number) / Math.log(10));
        double power = Math.pow(10, sf - decimals);
        double result = (Math.round(number * power) / power);
        if (Double.isInfinite(result) || Double.isNaN(result)) {
            return String.valueOf(number);
        }
        return String.valueOf(result);
    }

    /**
     * print label with a prespecified length (label will be shortened or spaces
     * will introduced, if necessary)
     *
     * @param out output stream
     * @param label label to be printed
     * @param width desired length
     */
    public void displayLabel(Writer out, String label, int width) throws IOException {
        int len = label.length();

        if (len == width) {
            // Print as is
            out.write(label);
        } else if (len < width) {
            // fill rest with spaces
            out.write(label);
            multiplePrint(out, ' ', width - len);
        } else {
            // Print first width characters
            for (int i = 0; i < width; i++) {
                out.write(label.charAt(i));
            }
        }
    }

    /**
     * print integer, aligned to a reference number, (introducing space at the
     * left side)
     *
     * @param out output stream
     * @param num number to be printed
     * @param maxNum reference number
     */
    public void displayInteger(Writer out, int num, int maxNum) throws IOException {
        int lenNum = Integer.toString(num).length();
        int lenMaxNum = Integer.toString(maxNum).length();

        if (lenNum < lenMaxNum) {
            multiplePrint(out, ' ', lenMaxNum - lenNum);
        }
        out.write(num);
    }

    /**
     * print whitespace of length of a string displaying a given integer
     *
     * @param out stream
     * @param maxNum number
     */
    public void displayIntegerWhite(Writer out, int maxNum) throws IOException {
        int lenMaxNum = Integer.toString(maxNum).length();

        multiplePrint(out, ' ', lenMaxNum);
    }

    /**
     * repeatedly print a character
     *
     * @param out output stream
     * @param c character
     * @param num number of repeats
     */
    public void multiplePrint(Writer out, char c, int num) throws IOException {
        for (int i = 0; i < num; i++) {
            out.write(c);
        }
    }

    /**
     * returns of string of a given length of a single character.
     *
     * @param size length of the string required
     * @param c character
     */
    public static String space(int size, char c) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < size; i++) {
            sb.append(c);
        }
        return new String(sb);
    }

    //
    // Private stuff
    //
    // private constructor
    private FormattedOutput() {
        nf = NumberFormat.getInstance(Locale.UK);
        nf.setGroupingUsed(false);
    }
    private static FormattedOutput singleton = null;
    private final NumberFormat nf;
}
