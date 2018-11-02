// ParetoDistribution.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.stats.statistics;




/**
 * Pareto distribution
 * (scale-free distribution without characteristic length scale).
 *
 * Parameters: shape parameter k>0, scale parameter m>0 ("minimum income")
 *
 * @version $Id: ParetoDistribution.java,v 1.1 2007/01/12 03:26:16 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 */
public class ParetoDistribution
{
	//
	// Public stuff
	//

	/**
	 * probability density function of the Pareto distribution
	 * 
	 * @param x argument (>=m)
	 * @param k shape parameter (>0)
	 * @param m scale parameter (>0, "minimum income")
	 *
	 * @return pdf value
	 */
	public static double pdf(double x, double k, double m)
	{
		return k*Math.pow(m,k)*Math.pow(x,-(k+1));
	}

	/**
	 * cumulative density function of the Pareto distribution
	 * 
	 * @param x argument (>=m)
	 * @param k shape parameter (>0)
	 * @param m scale parameter (>0, "minimum income")
	 *
	 * @return cdf value
	 */
	public static double cdf(double x, double k, double m)
	{
		return 1.0-Math.pow(m/x, k);
	}


	/**
	 * quantile (inverse cumulative density function) of the Pareto distribution
	 * 
	 * @param p argument (0 < p < 1)
	 * @param k shape parameter (>0)
	 * @param m scale parameter (>0, "minimum income")
	 *
	 * @return icdf value
	 */
	public static double quantile(double p, double k, double m)
	{
		return m/Math.pow(1.0-p,1.0/k);
	}
	
	/**
	 * mean of the Pareto distribution
	 * 
	 * @param k shape parameter (>0)
	 * @param m scale parameter (>0, "minimum income")
	 *
	 * @return mean
	 */
	public static double mean(double k, double m)
	{
		if (k > 1.0)
			return m*k/(k-1.0);
		else
			return Double.POSITIVE_INFINITY;
	}

	/**
	 * variance of the Pareto distribution
	 * 
	 * @param k shape parameter (>0)
	 * @param m scale parameter (>0, "minimum income")
	 *
	 * @return variance
	 */
	public static double variance(double k, double m)
	{
		if (k > 2.0)
			return m*m*k/((k-1.0)*(k-1.0)*(k-2.0));
		else
			return Double.POSITIVE_INFINITY;
	}
	
	/**
	 * moments E(X^n) of the Pareto distribution
	 * 
	 * @param n moment
	 * @param k shape parameter (>0)
	 * @param m scale parameter (>0, "minimum income")
	 *
	 * @return variance
	 */
	public static double moment(int n, double k, double m)
	{
		if (k > n)
			return Math.pow(m,n)*k/(k-n);
		else
			return Double.POSITIVE_INFINITY;
	}
}
