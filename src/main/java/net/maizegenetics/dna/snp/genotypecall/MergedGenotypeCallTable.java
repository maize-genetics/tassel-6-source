package net.maizegenetics.dna.snp.genotypecall;

import java.util.function.BiFunction;

import net.maizegenetics.dna.snp.GenotypeTableUtils;

public class MergedGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable[] myGenotypeCallTables;
    private int[][] taxonMap;
    private int[][] positionMap;
    private BiFunction<Integer,Integer,Byte> mergeCallFunction;
    /*
    private HashMap<Integer,Integer[]> taxonMap;
    private HashMap<Integer,Integer[]> positionMap;
    */
    //private MergedGenotypeCallTable(GenotypeCallTable[] genotypeCallTables, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings, int maxNumAlleles, HashMap<Integer,Integer[]> taxonMap, HashMap<Integer,Integer[]> positionMap) {
    private MergedGenotypeCallTable(GenotypeCallTable[] genotypeCallTables, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings, int maxNumAlleles, int[][] taxonMap, int[][] positionMap) {
        super(numTaxa, numSites, phased, alleleEncodings, maxNumAlleles);
        
        myGenotypeCallTables = genotypeCallTables;
        this.taxonMap = taxonMap;
        this.positionMap = positionMap;
        this.mergeCallFunction = (taxon,site) -> {
          return myGenotypeCallTables[0].genotype(taxonMap[taxon][0],positionMap[site][0]);  
        };
        /*
        System.out.println("TMAP: ");
        for(int i = 0; i < taxonMap.length; i++) {
            for(int j = 0; j < taxonMap[i].length; j++) {
                System.out.print(taxonMap[i][j]+" ");
            }
            System.out.println();
        }
        */
        
        /*
        System.out.println("PMAP: ");
        for(int i = 0; i < positionMap.length; i++) {
            for(int j = 0; j < positionMap[i].length; j++) {
                System.out.print(positionMap[i][j]+" ");
            }
            System.out.println();
        }
        */
       /*
        System.out.println("Table");
        for(int i = 0; i < taxonMap.length; i++) {
            for(int j = 0; j < positionMap.length; j++) {
                System.out.print(genotype(i,j)+" ");
            }
            System.out.println();
        }
        */
        
    }

    //public static GenotypeCallTable getInstance(GenotypeCallTable[] genotypeCallTables, HashMap<Integer,Integer[]> taxonMap, HashMap<Integer,Integer[]> positionMap) {
    public static GenotypeCallTable getInstance(GenotypeCallTable[] genotypeCallTables, int[][] taxonMap, int[][] positionMap) {

        if ((genotypeCallTables == null) || (genotypeCallTables.length == 0)) {
            throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: must provide genotype call tables.");
        }

        if (genotypeCallTables.length == 1) {
            return genotypeCallTables[0];
        }
        
        boolean allTheSame = true;
        String[][] encodings = genotypeCallTables[0].alleleDefinitions();
        if (encodings.length == 1) {
            for (int i = 1; i < genotypeCallTables.length; i++) {
                String[][] current = genotypeCallTables[i].alleleDefinitions();
                if ((current.length == 1) && (encodings[0].length == current[0].length)) {
                    for (int j = 0; j < encodings[0].length; j++) {
                        if (!current[0][j].equals(encodings[0][j])) {
                            allTheSame = false;
                            break;
                        }
                    }
                } else {
                    allTheSame = false;
                    break;
                }

                if (!allTheSame) {
                    break;
                }
            }
        } else {
            allTheSame = false;
        }

        int totalSites = positionMap.length;
        /*
        int totalSites = 0;
        for (int i = 0; i < genotypeCallTables.length; i++) {
            totalSites += genotypeCallTables[i].numberOfSites();
        }
        */

        String[][] alleleStates;
        if (allTheSame) {
            alleleStates = encodings;
        } else {
            String[][] result = new String[totalSites][];
            int count = 0;
            for (int i = 0; i < genotypeCallTables.length; i++) {
                for (int j = 0, n = genotypeCallTables[i].numberOfSites(); j < n; j++) {
                    result[count++] = genotypeCallTables[i].alleleDefinitions(j);
                }
            }
            alleleStates = result;
        }

        //int numTaxa = genotypeCallTables[0].numberOfTaxa();
        int numTaxa = taxonMap.length;
        boolean phased = genotypeCallTables[0].isPhased();
        int maxNumAlleles = genotypeCallTables[0].maxNumAlleles();
        for (int i = 1; i < genotypeCallTables.length; i++) {
            /*
            if (numTaxa != genotypeCallTables[i].numberOfTaxa()) {
                throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: number of taxa not equal.");
            }
            */
            if (phased != genotypeCallTables[i].isPhased()) {
                throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: phase is different.");
            }
            if (maxNumAlleles != genotypeCallTables[i].maxNumAlleles()) {
                throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: max number of alleles is different.");
            }
        }

        return new MergedGenotypeCallTable(genotypeCallTables, numTaxa, totalSites, phased, alleleStates, maxNumAlleles,taxonMap, positionMap);
    }

    @Override
    public byte genotype(int taxon, int site) {
        //Potential solution
        //return mergeCallFunction.apply(taxon,site);
       
        /* Don't do this if we can help it.
        byte[] genotypeCalls = new byte[myGenotypeCallTables.length];
        for(int i = 0; i < genotypeCalls.length; i++) {
            genotypeCalls[i] = myGenotypeCallTables[i].genotype(taxonMap[taxon][i], positionMap[site][i]);
        }
        */
        //Just Return the value of the first genotypeTable for now. This will change
        return myGenotypeCallTables[0].genotype(taxonMap[taxon][0], positionMap[site][0]);
        
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        //TODO Add in GenotypeCall Logic
        return myGenotypeCallTables[0].genotypeAsString(taxonMap[taxon][0], positionMap[site][0]);
        
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean isSiteOptimized() {
        return myGenotypeCallTables[0].isSiteOptimized();
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return GenotypeTableUtils.getDiploidValues(genotype(taxon, site));
    }

}
