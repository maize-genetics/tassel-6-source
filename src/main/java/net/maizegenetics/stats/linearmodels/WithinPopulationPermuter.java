package net.maizegenetics.stats.linearmodels;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class WithinPopulationPermuter {
	final double[] data;
	final int npops;
	ArrayList<int[]> popIndices = new ArrayList<int[]>();
	static Random randomizer = new Random();
	
	public WithinPopulationPermuter(double[] originalData, FactorModelEffect popEffect) {
		data = originalData;
		npops = popEffect.getNumberOfLevels();
		int[] levels = popEffect.getLevels();
		int[] levelCounts = popEffect.getLevelCounts();

		for(int p = 0; p < npops; p++) {
			popIndices.add(new int[levelCounts[p]]);
		}

		int[] count = new int[npops];

		int n = levels.length;
		for (int i = 0; i < n; i++) {
			int pop = levels[i];
			popIndices.get(pop)[count[pop]++] = i;
		}

	}
	
	public double[] getPermutedData() {
		int ndata = data.length;
		double[] permutedData = Arrays.copyOf(data, ndata);
		
		for (int p = 0; p < npops; p++) {
			int[] ndx = popIndices.get(p);
			
			int n = ndx.length;
			for (int i = n - 1; i >= 1; i--) {
				int j = randomizer.nextInt(i + 1);
				double temp = permutedData[ndx[j]];
				permutedData[ndx[j]] = permutedData[ndx[i]];
				permutedData[ndx[i]] = temp;
			}
		}
		
		return permutedData;
	}
	
	public void permuteData(double[] data) {
		int ndata = data.length;
		
		for (int p = 0; p < npops; p++) {
			int[] ndx = popIndices.get(p);
			
			int n = ndx.length;
			for (int i = n - 1; i >= 1; i--) {
				int j = randomizer.nextInt(i + 1);
				double temp = data[ndx[j]];
				data[ndx[j]] = data[ndx[i]];
				data[ndx[i]] = temp;
			}
		}
		
	}
}