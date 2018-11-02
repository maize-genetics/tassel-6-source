/*
 *  ByteGenotypeCallTable
 */
package net.maizegenetics.dna.snp.genotypecall;

import java.util.stream.Stream;

import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

import org.apache.log4j.Logger;

/**
 * In memory byte implementation of GenotypeCallTable backed the high efficiency
 * SuperByteMatrix class. Although the GenotypeCallTable is accessed as two
 * dimensional array, for efficiency it is actually backed single dimension
 * arrays with either site or taxa as the inner loop.
 *
 * @see SuperByteMatrix
 *
 * @author Terry Casstevens
 */
class ByteGenotypeCallTable extends AbstractGenotypeCallTable {

    private static final Logger myLogger = Logger.getLogger(ByteGenotypeCallTable.class);
    SuperByteMatrix myGenotype;
    private SuperByteMatrix mySiteInnerLoop;
    private SuperByteMatrix myTaxonInnerLoop;

    ByteGenotypeCallTable(byte[][] genotype, boolean phased, String[][] alleleEncodings) {
        super(genotype.length, genotype[0].length, phased, alleleEncodings);
        mySiteInnerLoop = SuperByteMatrixBuilder.getInstance(myTaxaCount, mySiteCount);
        myGenotype = mySiteInnerLoop;
        for (int t = 0; t < myTaxaCount; t++) {
            for (int s = 0; s < mySiteCount; s++) {
                myGenotype.set(t, s, genotype[t][s]);
            }
        }
    }

    ByteGenotypeCallTable(SuperByteMatrix genotype, boolean phased, String[][] alleleEncodings) {
        super(genotype.getNumRows(), genotype.getNumColumns(), phased, alleleEncodings);
        if (genotype.isColumnInnerLoop()) {
            mySiteInnerLoop = genotype;
        } else {
            myTaxonInnerLoop = genotype;
        }
        myGenotype = genotype;
    }

    @Override
    public byte genotype(int taxon, int site) {
        return myGenotype.get(taxon, site);
    }

    @Override
    public byte[] genotypeForAllSites(int taxon) {
        return myGenotype.getAllColumns(taxon);
    }

    @Override
    public byte[] genotypeForSiteRange(int taxon, int start, int end) {
        return myGenotype.getColumnRange(taxon, start, end);
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        return myGenotype.getAllRows(site);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {

        if (siteInnerLoop) {
            if (mySiteInnerLoop == null) {
                mySiteInnerLoop = SuperByteMatrixBuilder.getInstanceTranspose(myTaxonInnerLoop);
            }
            myGenotype = mySiteInnerLoop;
        } else {
            if (myTaxonInnerLoop == null) {
                myTaxonInnerLoop = SuperByteMatrixBuilder.getInstanceTranspose(mySiteInnerLoop);
            }
            myGenotype = myTaxonInnerLoop;
        }

    }

    @Override
    public boolean isSiteOptimized() {
        if (myTaxonInnerLoop != null) {
            return true;
        } else {
            return false;
        }
    }

    @Override
    public Stream<Byte> stream() {
        return myGenotype.stream();
    }

    @Override
    public Stream<Byte> stream(int taxon) {
        return myGenotype.stream(taxon);
    }

}
