// SplitSystem.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.taxa.tree;

import net.maizegenetics.taxa.TaxaList;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * data structure for a set of splits 
 *
 * @version $Id: SplitSystem.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 */
public class SplitSystem
{
	//
	// Public stuff
	//

	/**
	 * @param idGroup  sequence labels
	 * @param size     number of splits
	 */
	public SplitSystem(TaxaList idGroup, int size)
	{
		this.idGroup = idGroup;
		
		labelCount = idGroup.numberOfTaxa();
		splitCount = size;
		
		splits = new boolean[splitCount][labelCount];
	}

	/** get number of splits */
	public int getSplitCount()
	{		
		return splitCount;
	}

	/** get number of labels */
	public int getLabelCount()
	{		
		return labelCount;
	}

	/** get split vector */
	public boolean[][] getSplitVector()
	{		
		return splits;
	}

	/** get split */
	public boolean[] getSplit(int i)
	{		
		return splits[i];
	}


	/** get idGroup */
	public TaxaList getIdGroup()
	{		
		return idGroup;
	}

	/**
	  + test whether a split is contained in this split system
	  * (assuming the same leaf order)
	  *
	  * @param split split
	  */
	public boolean hasSplit(boolean[] split)
	{
		for (int i = 0; i < splitCount; i++)
		{
			if (SplitUtils.isSame(split, splits[i])) return true;
		}
			
		return false;
	}


	/** print split system */
	public String toString()
	{
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		
		for (int i = 0; i < labelCount; i++)
		{
			pw.println(idGroup.get(i));
		}
		pw.println();
		
		
		for (int i = 0; i < splitCount; i++)
		{
			for (int j = 0; j < labelCount; j++)
			{
				if (splits[i][j] == true)
					pw.print('*');
				else
					pw.print('.');
			}
			
			pw.println();
		}

		return sw.toString();
	}

	
	//
	// Private stuff
	//
	
	private int labelCount, splitCount;
	private TaxaList idGroup;
	private boolean[][] splits;
}
