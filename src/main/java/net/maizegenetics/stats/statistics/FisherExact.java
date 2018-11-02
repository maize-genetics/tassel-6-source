// ContigencyTable.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.stats.statistics;

import java.util.Arrays;

/**
 * This does a Fisher Exact test.  The Fisher's Exact test procedure calculates an exact probability value
 * for the relationship between two dichotomous variables, as found in a two by two crosstable. The program
 * calculates the difference between the data observed and the data expected, considering the given marginal
 * and the assumptions of the model of independence. It works in exactly the same way as the Chi-square test
 * for independence; however, the Chi-square gives only an estimate of the true probability value, an estimate
 * which might not be very accurate if the marginal is very uneven or if there is a small value (less than five)
 * in one of the cells.
 *
 * It uses an array of factorials initialized at the beginning to provide speed.
 * There could be better ways to do this.
 *
 * @author Ed Buckler
 * @version $Id: FisherExact.java,v 1
 */

public class FisherExact {
	private static final boolean DEBUG = false;
	private static double[] factorialArray;
	private static int  maxSize; // not really size, is highest number

	private static FisherExact myFisherExact;
	
	/**
	 * Static method to get instance
	 */
	public static FisherExact getInstance(int size) {
		if (myFisherExact == null) {
			myFisherExact = new FisherExact(size);
		}
		else if (size > maxSize) {
			factorialArray = resizeArray(size);
			maxSize = size;
		}
		return myFisherExact;
	}

	/**
	 * constructor for FisherExact table
	 *
	 * @param maxSize is the maximum sum that will be encountered by the table (a+b+c+d)
	 */
	private FisherExact(int maxSize) {
		FisherExact.maxSize = maxSize;
		factorialArray = new double[maxSize + 1]; // +1 to account for 0
		factorialArray[0] = 0.0;
		for (int i = 1; i <= FisherExact.maxSize; i++) {
			factorialArray[i] = factorialArray[i - 1] + Math.log(i);
		}
	}

	private static synchronized double[] resizeArray(int size) {
		int flength = factorialArray.length;

		FisherExact.maxSize = size;
		if (flength > size+1) return factorialArray;
		double[] newF = Arrays.copyOf(factorialArray, size+1); //copy old values
		// Calculate new values
		for (int idx = flength; idx <= size; idx++){
			newF[idx] = newF[idx - 1] + Math.log(idx);
		}
		return newF;
	}

	/**
	 * calculates the P-value for this specific state
	 *
	 * @param a     a, b, c, d are the four cells in a 2x2 matrix
	 * @param b
	 * @param c
	 * @param d
	 * @return the P-value
	 */
	public final double getP(int a, int b, int c, int d) {
		int n = a + b + c + d;
		if (n > maxSize) {
			factorialArray = resizeArray(n);
			 //return Double.NaN;
		}
		double p;
		p = (factorialArray[a + b] + factorialArray[c + d] + factorialArray[a + c] + factorialArray[b + d]) - 
				(factorialArray[a] + factorialArray[b] + factorialArray[c] + factorialArray[d] + factorialArray[n]);
		return Math.exp(p);
	}

	/**
	 * Calculates the one-tail P-value for the Fisher Exact test.  Determines whether to calculate the right- or left-
	 * tail, thereby always returning the smallest p-value.
	 *
	 * @param a     a, b, c, d are the four cells in a 2x2 matrix
	 * @param b
	 * @param c
	 * @param d
	 * @return one-tailed P-value (right or left, whichever is smallest)
	 */
	public final double getCumlativeP(int a, int b, int c, int d) {
		int min, i;
		int n = a + b + c + d;
		if (n > maxSize) {
			factorialArray = resizeArray(n);
			//return Double.NaN;
		}
		double p = 0;

		p += getP(a, b, c, d);
		if (DEBUG) {System.out.println("p = " + p);}
		if ((a * d) >= (b * c)) {
			if (DEBUG) {System.out.println("doing R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			min = (c < b) ? c : b;
			for (i = 0; i < min; i++) {
				if (DEBUG) {System.out.print("doing round " + i);}
				p += getP(++a, --b, --c, ++d);
				if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
			}
			//           System.out.println("");
		}
		if ((a * d) < (b * c)) {
			if (DEBUG) {System.out.println("doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
			min = (a < d) ? a : d;
			for (i = 0; i < min; i++) {
				if (DEBUG) {System.out.print("doing round " + i);}
				double pTemp = getP(--a, ++b, ++c, --d);
				if (DEBUG) {System.out.print("\tpTemp = " + pTemp);}
				p += pTemp;
				if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
			}
		}
		return p;
	}

	/**
	 * Calculates the right-tail P-value for the Fisher Exact test.
	 *
	 * @param a     a, b, c, d are the four cells in a 2x2 matrix
	 * @param b
	 * @param c
	 * @param d
	 * @return one-tailed P-value (right-tail)
	 */
	public final double getRightTailedP(int a, int b, int c, int d) {
		int min, i;
		int n = a + b + c + d;
		if (n > maxSize) {
			factorialArray = resizeArray(n);
			//return Double.NaN;
		}
		double p = 0;

		p += getP(a, b, c, d);
		if (DEBUG) {System.out.println("p = " + p);}
		if (DEBUG) {System.out.println("doing R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
		min = (c < b) ? c : b;
		for (i = 0; i < min; i++) {
			p += getP(++a, --b, --c, ++d);
		}
		return p;
	}

	/**
	 * Calculates the right-tail P-value for the Fisher Exact test.
	 *
	 * @param a     a, b, c, d are the four cells in a 2x2 matrix
	 * @param b
	 * @param c
	 * @param d
	 * @return one-tailed P-value (right-tail)
	 */
	public final double getRightTailedPQuick(int a, int b, int c, int d, double maxP) {
		int min, i;
		// int n = a + b + c + d;
		//        if (n > maxSize) {
		//            return Double.NaN;
		//        }
		double p = 0;

		p += getP(a, b, c, d);
		min = (c < b) ? c : b;
		for (i = 0; (i < min) && (p<maxP); i++) {
			p += getP(++a, --b, --c, ++d);

		}
		return p;
	}

	/**
	 * Calculates the left-tail P-value for the Fisher Exact test.
	 *
	 * @param a     a, b, c, d are the four cells in a 2x2 matrix
	 * @param b
	 * @param c
	 * @param d
	 * @return one-tailed P-value (left-tail)
	 */
	public final double getLeftTailedP(int a, int b, int c, int d) {
		int min, i;
		int n = a + b + c + d;
		if (n > maxSize) {
			factorialArray = resizeArray(n);
			//return Double.NaN;
		}
		double p = 0;

		p += getP(a, b, c, d);
		if (DEBUG) {System.out.println("p = " + p);}
		if (DEBUG) {System.out.println("doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
		min = (a < d) ? a : d;
		for (i = 0; i < min; i++) {
			if (DEBUG) {System.out.print("doing round " + i);}
			double pTemp = getP(--a, ++b, ++c, --d);
			if (DEBUG) {System.out.print("\tpTemp = " + pTemp);}
			p += pTemp;
			if (DEBUG) {System.out.println("\ta=" + a + " b=" + b + " c=" + c + " d=" + d);}
		}
		return p;
	}

	/**
	 *   Calculates the two-tailed P-value for the Fisher Exact test.
	 *
	 *   In order for a table under consideration to have its p-value included
	 *   in the final result, it must have a p-value less than the original table's P-value, i.e.
	 *   Fisher's exact test computes the probability, given the observed marginal
	 *   frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
	 *   By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
	 *   occurrence in the same direction (one-tailed) or in both directions (two-tailed).
	 *
	 * @param a     a, b, c, d are the four cells in a 2x2 matrix
	 * @param b
	 * @param c
	 * @param d
	 * @return two-tailed P-value
	 */
	public final double getTwoTailedP(int a, int b, int c, int d) {
		int min, i;
		int n = a + b + c + d;
		if (n > maxSize) {
			System.out.printf("LCJ - FE:getTwoTailedP, resize for a %d, b %d c %d d %d\n", a,b,c,d);
			factorialArray = resizeArray(n);
			//return Double.NaN;
		}
		double p = 0;

		double baseP = getP(a, b, c, d);
		//         in order for a table under consideration to have its p-value included
		//         in the final result, it must have a p-value less than the baseP, i.e.
		//         Fisher's exact test computes the probability, given the observed marginal
		//         frequencies, of obtaining exactly the frequencies observed and any configuration more extreme.
		//         By "more extreme," we mean any configuration (given observed marginals) with a smaller probability of
		//         occurrence in the same direction (one-tailed) or in both directions (two-tailed).

		if (DEBUG) {System.out.println("baseP = " + baseP);}
		int initialA = a, initialB = b, initialC = c, initialD = d;
		p += baseP;
		if (DEBUG) {System.out.println("p = " + p);}
		if (DEBUG) {System.out.println("Starting with R-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
		min = (c < b) ? c : b;
		for (i = 0; i < min; i++) {
			if (DEBUG) {System.out.print("doing round " + i);}
			double tempP = getP(++a, --b, --c, ++d);
			if (tempP <= baseP) {
				if (DEBUG) {System.out.print("\ttempP (" + tempP + ") is less than baseP (" + baseP + ")");}
				p += tempP;
			}
			if (DEBUG) {System.out.println(" a=" + a + " b=" + b + " c=" + c + " d=" + d);}
		}

		// reset the values to their original so we can repeat this process for the other side
		a = initialA;
		b = initialB;
		c = initialC;
		d = initialD;

		if (DEBUG) {System.out.println("Now doing L-tail: a=" + a + " b=" + b + " c=" + c + " d=" + d);}
		min = (a < d) ? a : d;
		if (DEBUG) {System.out.println("min = " + min);}
		for (i = 0; i < min; i++) {
			if (DEBUG) {System.out.print("doing round " + i);}
			double pTemp = getP(--a, ++b, ++c, --d);
			if (DEBUG) {System.out.println("  pTemp = " + pTemp);}
			if (pTemp <= baseP) {
				if (DEBUG) {System.out.print("\ttempP (" + pTemp + ") is less than baseP (" + baseP + ")");}
				p += pTemp;
			}
			if (DEBUG) {System.out.println(" a=" + a + " b=" + b + " c=" + c + " d=" + d);}
		}
		return p;
	}
	 
}