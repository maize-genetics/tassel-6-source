// PenalizedLikelihood.java
//
// (c) 20001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.stats.statistics;




/**
 * Penalized likelihood criteria
 *
 * @version $Id: PenalizedLikelihood.java,v 1.1 2007/01/12 03:26:16 tcasstevens Exp $ 
 *
 * @author Korbinian Strimmer
 */
public class PenalizedLikelihood
{
	//
	// Public stuff
	//
	
	/**
	 * Akaike (AIC) correction (Akaike 1974)
	 *
	 * @param l    log-likelihood
	 * @param k    number of inferred parameters
	 *
	 * @return     l - k   
	 */
	public static double AIC(double l, int k)
	{
		return l - (double) k;
	}
	
	/**
	 * BIC correction (Schwarz 1978)
	 *
	 * @param l    log-likelihood
	 * @param k    number of inferred parameters
	 * @param n    sample size 
	 *
	 * @return     l - k/2 log(n)   
	 */
	public static double BIC(double l, int k, int n)
	{
		return l - (double)k/2.0* Math.log(n);
	}
	
	/**
	 * Second-order Akaike (AICC) correction (Hurvich and Tsai 1989)
	 *
	 * @param l    log-likelihood
	 * @param k    number of inferred parameters
	 * @param n    sample size 
	 *
	 * @return     l - k - (k(k+1))/(n - k - 1)   
	 */
	public static double AICC(double l, int k, int n)
	{
		if (k > n-2) throw new IllegalArgumentException("k must be smaller than n-1");
		
		return  l - k - (double) (k*(k+1.0))/ (double) (n - k - 1.0) ;
	}

}
