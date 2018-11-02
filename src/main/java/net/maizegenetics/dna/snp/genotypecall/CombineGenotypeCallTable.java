/*
 *  CombineGenotypetypeCallTable
 * 
 *  Created on Dec 9, 2014
 */
package net.maizegenetics.dna.snp.genotypecall;

/**
 *
 * @author Terry Casstevens
 */
public class CombineGenotypeCallTable extends AbstractGenotypeCallTable {

    private final GenotypeCallTable[] myGenotypeCallTables;
    private final int[] mySiteOffsets;

    private CombineGenotypeCallTable(GenotypeCallTable[] genotypeCallTables, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings, int maxNumAlleles) {
        super(numTaxa, numSites, phased, alleleEncodings, maxNumAlleles);
        myGenotypeCallTables = genotypeCallTables;

        mySiteOffsets = new int[genotypeCallTables.length + 1];

        mySiteOffsets[0] = 0;
        int count = 0;
        for (int i = 0; i < genotypeCallTables.length; i++) {
            count = genotypeCallTables[i].numberOfSites() + count;
            mySiteOffsets[i + 1] = count;
        }
    }

    public static GenotypeCallTable getInstance(GenotypeCallTable[] genotypeCallTables) {

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

        int totalSites = 0;
        for (int i = 0; i < genotypeCallTables.length; i++) {
            totalSites += genotypeCallTables[i].numberOfSites();
        }

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

        int numTaxa = genotypeCallTables[0].numberOfTaxa();
        boolean phased = genotypeCallTables[0].isPhased();
        int maxNumAlleles = genotypeCallTables[0].maxNumAlleles();
        for (int i = 1; i < genotypeCallTables.length; i++) {
            if (numTaxa != genotypeCallTables[i].numberOfTaxa()) {
                throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: number of taxa not equal.");
            }
            if (phased != genotypeCallTables[i].isPhased()) {
                throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: phase is different.");
            }
            if (maxNumAlleles != genotypeCallTables[i].maxNumAlleles()) {
                throw new IllegalArgumentException("CombineGenotypeCallTable: getInstance: max number of alleles is different.");
            }
        }

        return new CombineGenotypeCallTable(genotypeCallTables, numTaxa, totalSites, phased, alleleStates, maxNumAlleles);
    }

    @Override
    public byte genotype(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeCallTables[translate].genotype(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeCallTables[translate].genotypeAsString(taxon, site - mySiteOffsets[translate]);
    }
    
    @Override
    public String diploidAsString(int site, byte value) {
        int translate = translateSite(site);
        return myGenotypeCallTables[translate].diploidAsString(site - mySiteOffsets[translate], value);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    private int translateSite(int site) {

        for (int i = 1; i < mySiteOffsets.length; i++) {
            if (mySiteOffsets[i] > site) {
                return i - 1;
            }
        }
        throw new IndexOutOfBoundsException("CombineGenotypeCallTable: translateSite: index out of range: " + site);

    }

    @Override
    public boolean isSiteOptimized() {
        return myGenotypeCallTables[0].isSiteOptimized();
    }

}
