// ExponentialDistribution.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.stats.statistics;




/**
 * exponential distribution.
 *
 * (Parameter: lambda; mean: 1/lambda; variance: 1/lambda^2)
 *
 * The exponential distribution is a special case of the Gamma distribution
 * (shape parameter = 1.0, scale = 1/lambda).
 *
 * @version $Id: ExponentialDistribution.java,v 1.1 2007/01/12 03:26:16 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 */
public class ExponentialDistribution extends GammaDistribution
{
	//
	// Public stuff
	//

	/**
	 * probability density function of the exponential distribution
	 * (mean = 1/lambda)
	 * 
	 * @param x argument
	 * @param lambda parameter of exponential distribution
	 *
	 * @return pdf value
	 */
	public static double pdf(double x, double lambda)
	{
		return lambda*Math.exp(-lambda*x);
		//return pdf(x, 1.0, 1/lambda);
	}

	/**
	 * cumulative density function of the exponential distribution
	 * 
	 * @param x argument
	 * @param lambda parameter of exponential distribution
	 *
	 * @return cdf value
	 */
	public static double cdf(double x, double lambda)
	{
		return 1.0-Math.exp(-lambda*x);
		//return cdf(x, 1.0, 1/lambda);
	}


	/**
	 * quantile (inverse cumulative density function) of the exponential distribution
	 *
	 * @param y argument
	 * @param lambda parameter of exponential distribution
	 *
	 * @return icdf value
	 */
	public static double quantile(double y, double lambda)
	{
		return -(1.0/lambda)*Math.log(1.0-y);
		//return quantile(y, 1.0, 1/lambda);
	}
	
	/**
	 * mean of the exponential distribution
	 *
	 * @param lambda parameter of exponential distribution
	 *
	 * @return mean
	 */
	public static double mean(double lambda)
	{
		return 1.0/(lambda);
	}

	/**
	 * variance of the exponential distribution
	 *
	 * @param lambda parameter of exponential distribution
	 *
	 * @return variance
	 */
	public static double variance(double lambda)
	{
		return 1.0/(lambda*lambda);
	}
}
