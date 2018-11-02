package net.maizegenetics.stats.linearmodels;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class ModelEffectUtils {
	private static final List<Byte> homGeno = Arrays.asList(
			NucleotideAlignmentConstants.getNucleotideDiploidByte("A"),
			NucleotideAlignmentConstants.getNucleotideDiploidByte("C"),
			NucleotideAlignmentConstants.getNucleotideDiploidByte("G"),
			NucleotideAlignmentConstants.getNucleotideDiploidByte("T"),
			NucleotideAlignmentConstants.getNucleotideDiploidByte("Z"),
			NucleotideAlignmentConstants.getNucleotideDiploidByte("+"),
			NucleotideAlignmentConstants.getNucleotideDiploidByte("-"));
	
	private ModelEffectUtils() {}
	
	public static DoubleMatrix getXtY(ModelEffect X, ModelEffect Y) {
		if (X instanceof FactorModelEffect) {
			FactorModelEffect fme = (FactorModelEffect) X;
			if (Y instanceof FactorModelEffect) {
				return fme.getXtX2((FactorModelEffect) Y);
			} else if (Y instanceof CovariateModelEffect) {
				return fme.getXty(((CovariateModelEffect) Y).getCovariate());
			} else if (Y instanceof NestedCovariateModelEffect) {
				return fme.getX().mult(Y.getX(), true, false);
			}
			
		} else if (X instanceof CovariateModelEffect) {
			double[] cov = ((CovariateModelEffect) X).getCovariate();
			if (Y instanceof FactorModelEffect) {
				return getXtY(Y,X).transpose();
			} else if (Y instanceof CovariateModelEffect) {
				return Y.getXty(cov);
			} else if (Y instanceof NestedCovariateModelEffect) {
				return Y.getXty(cov).transpose();
			}
			
		} else if (X instanceof NestedCovariateModelEffect) {
			if (Y instanceof FactorModelEffect) {
				return X.getX().mult(Y.getX(), true, false);
			} else if (Y instanceof CovariateModelEffect) {
				return X.getXty(((CovariateModelEffect) Y).getCovariate());
			} else if (Y instanceof NestedCovariateModelEffect) {
				return ((NestedCovariateModelEffect) X).getXtX2((NestedCovariateModelEffect) Y);
			}
			
		}
		return null;
	}
	
    public static int[] getIntegerLevels(Object[] originalLevels) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	HashMap<Object, Integer> levelMap = new HashMap<Object, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels[i]);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels[i], ndx);
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	return intLevels;
    }
    
    public static <T> int[] getIntegerLevels(T[] originalLevels, ArrayList<T> ids) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	HashMap<T, Integer> levelMap = new HashMap<T, Integer>();
    	for (int i = 0; i < nLevels; i++) {
    		Integer ndx = levelMap.get(originalLevels[i]);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(originalLevels[i], ndx);
    		}
    		intLevels[i] = ndx.intValue();
    	}
    	
    	if (ids != null) {
        	TreeSet<Entry<T,Integer>> sortedEntries = new TreeSet<Entry<T,Integer>>(new Comparator<Entry<T,Integer>>(){

    			@Override
    			public int compare(Entry<T, Integer> arg0, Entry<T, Integer> arg1) {
    				return arg0.getValue().compareTo(arg1.getValue());
    			}
        		
        	});
        	
        	sortedEntries.addAll(levelMap.entrySet());
        	for (Entry<T, Integer> entry:sortedEntries) {
        		ids.add(entry.getKey());
        	}
    	}
    	
    	return intLevels;
    }
    
    public static int[] getIntegerLevelsSortedByGenotype(String[] originalLevels, ArrayList<String> ids) {
    	int nLevels = originalLevels.length;
    	int[] intLevels = new int[nLevels];
    	
    	String nukes = "ACGT";
    	HashSet<String> idSet = Arrays.stream(originalLevels).collect(Collectors.toCollection(HashSet::new));
    	ids.addAll(idSet);
    	Collections.sort(ids, (a,b) -> {
    		String id1 = a.toString();
    		String id2 = b.toString();
    		if (nukes.contains(id1)) {
    			if (nukes.contains(id2)) return id1.compareTo(id2);
    			else return -1;
    		} else {
    			if (nukes.contains(id2)) return 1;
    			else return id1.compareTo(id2);
    		}
    	});
    	
    	HashMap<String,Integer> levelMap = new HashMap<>();
    	for (int i = 0; i < ids.size(); i++) {
    		levelMap.put(ids.get(i), i);
    	}
    	for (int i = 0; i < nLevels; i++) intLevels[i] = levelMap.get(originalLevels[i]);
    	return intLevels;
    }
    
    /**
     * @param genotypes		a byte array of genotypes
     * @param ids			an empty array list that will hold the Bytes corresponding to each level
     * @return				these genotypes translated into factor levels. On return, ids will hold the Byte labels for each level.
     * 						Levels will be sorted with homozygotes first then heterozygotes.
     */
    public static int[] getIntegerLevels(byte[] genotypes, ArrayList<Byte> ids) {
    	int n = genotypes.length;
    	
    	//create a tree set that sorts genotypes with homozygotes first
    	TreeSet<Byte> genotypeSet = new TreeSet<>((a,b) -> {
    		if (homGeno.contains(a)) {
    			if (homGeno.contains(b)) return a.compareTo(b);
    			return -1;
    		} else {
    			if (homGeno.contains(b)) return 1;
    			else return a.compareTo(b);
    		}
    	});
    	
    	for (byte b : genotypes) genotypeSet.add(b);
    	ids.addAll(genotypeSet);
    	HashMap<Byte,Integer> genoMap = new HashMap<>();
    	for (int i = 0; i < ids.size(); i++) genoMap.put(ids.get(i), i);
    	
    	return IntStream.range(0,n).map(i -> genoMap.get(genotypes[i]).intValue()).toArray();
    }
    
    public static <T> int[] getIntegerLevels(ArrayList<T> originalLevels, ArrayList<T> ids) {
    	int[] intLevels = new int[originalLevels.size()];
    	HashMap<T, Integer> levelMap = new HashMap<T,Integer>();
    	int count = 0;
    	for (T level:originalLevels) {
    		Integer ndx = levelMap.get(level);
    		if (ndx == null) {
    			ndx = new Integer(levelMap.size());
    			levelMap.put(level, ndx);
    		}
    		intLevels[count++] = ndx.intValue();
    	}

    	if (ids != null) {
        	TreeSet<Entry<T,Integer>> sortedEntries = new TreeSet<Entry<T,Integer>>(new Comparator<Entry<T,Integer>>(){

    			@Override
    			public int compare(Entry<T, Integer> arg0, Entry<T, Integer> arg1) {
    				return arg0.getValue().compareTo(arg1.getValue());
    			}
        		
        	});
        	
        	sortedEntries.addAll(levelMap.entrySet());
        	for (Entry<T, Integer> entry:sortedEntries) {
        		ids.add(entry.getKey());
        	}
    	}
    	return intLevels;
    }
    
    public static <T> int[] getIntegerLevels(ArrayList<T> originalLevels) {
    	return getIntegerLevels(originalLevels, null);
    }
    
    public static int[] getIntegerLevels(int[] originalLevels) {
        TreeSet<Integer> originalSet = Arrays.stream(originalLevels).collect(TreeSet<Integer>::new, TreeSet::add, TreeSet::addAll);
        int norig = originalLevels.length;
        if (originalSet.size() == norig) return originalLevels;
        
        int[] index = new int[norig];
        Arrays.fill(index, -1);
        int count = 0;
        for (Integer level : originalSet) index[level] = count++;
        int[] newLevels = new int[norig];
        for (int i = 0; i < norig; i++) newLevels[i] = index[originalLevels[i]];
        return newLevels;
    }
    
    public static double[] getNumericCodingForAdditiveModel(Object[] marker, String allele) {
    	int nmarkers = marker.length;
    	if (allele.equals(GenotypeTable.UNKNOWN_ALLELE_STR)) return new double[nmarkers];
    	
    	String firstMarker = ((String) marker[0]);
    	double[] values = new double[nmarkers];
    	
    	if (firstMarker.contains(":")) {
        	Pattern colon = Pattern.compile(":");
        	for (int m = 0; m < nmarkers; m++) {
        		String markerval = (String) marker[m];
        		String[] markerAlleles = colon.split(markerval);
        		if (markerAlleles[0].equals(allele)) values[m]++;
        		if (markerAlleles[1].equals(allele)) values[m]++;
        	}
    	} else {
    		Pattern nuc = Pattern.compile("[RSYWKM0]");
    		for (int m = 0; m < nmarkers; m++) {
    			String markerval = (String) marker[m];
    			if (markerval.equals(allele)) values[m] = 2;
    			else if (nuc.matcher((String) marker[m]).matches()) values[m] = 1;
    			else values[m] = 0;
    		}
    	}
    	return values;
    }
    
    /**
     * @param marker	the genotypes for a taxon at all sites as bytes
     * @param allele	the allele that should be coded as 1. all others will be 0
     * @return			a vector of values (0,1,2) that encode genotypes equal to the dosage of allele
     */
    public static double[] getNumericCodingForAdditiveModel(byte[] marker, byte allele) {
    	int nmarkers = marker.length;
    	double[] values = new double[nmarkers];
    	
    	for (int m = 0; m < nmarkers; m++) {
    		byte[] markerval = GenotypeTableUtils.getDiploidValues(marker[m]);
    		for (byte alleleval : markerval) if (alleleval == allele) values[m] += 1.0; 
    	}
    	return values;
    }

}
