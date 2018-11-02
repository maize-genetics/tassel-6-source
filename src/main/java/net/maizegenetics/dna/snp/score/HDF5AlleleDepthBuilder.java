/*
 *  HDF5AlleleDepthBuilder
 */
package net.maizegenetics.dna.snp.score;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;

/**
 * Builder to store information on DNA read depths.
 *
 * @author Terry Casstevens
 */
public class HDF5AlleleDepthBuilder extends AlleleDepthBuilder {

    private IHDF5Writer myHDF5Writer;
    private final int myNumSites;

    private HDF5AlleleDepthBuilder(IHDF5Writer writer, int numSites) {
        super(numSites);
        myHDF5Writer = writer;
        myNumSites = numSites;
    }

    /**
     * AlleleDepthBuilder is created and depths are stored in a HDF5 file.
     * setDepth methods are used to set the depths. Finish the building with
     * build()
     *
     * @param writer
     * @param numSites
     * @return
     */
    public static HDF5AlleleDepthBuilder getHDF5NucleotideInstance(IHDF5Writer writer, int numSites) {
        return new HDF5AlleleDepthBuilder(writer, numSites);
    }

    /**
     * AlleleDepth is returned for an immutable HDF5 file
     *
     * @param reader
     * @return allele depths
     */
    public static AlleleDepth getExistingHDF5Instance(IHDF5Reader reader) {
        //TODO is this the right name for this
        if (HDF5Utils.doesGenotypeDepthExist(reader) == false) {
            return null;
        }
        return new HDF5AlleleDepth(reader);
    }

    /**
     * Add taxon and set values for all sites and alleles for that taxon. First
     * dimension of depths is number of alleles (6 for Nucleotide) and second
     * dimension is sites.
     *
     * @param taxon taxon
     * @param depths depth values
     *
     * @return builder
     */
    public HDF5AlleleDepthBuilder addTaxon(Taxon taxon, byte[][] depths) {
        if ((depths == null) || (depths.length != 6)) {
            throw new IllegalStateException("AlleleDepthBuilder: addTaxon: Set A, C, G, T, -, + at once");
        }
        if (depths[0].length != myNumSites) {
            throw new IllegalStateException("AlleleDepthBuilder: addTaxon: Number of sites: " + depths[0].length + " should be: " + myNumSites);
        }
        synchronized (myHDF5Writer) {
            HDF5Utils.writeHDF5GenotypesDepth(myHDF5Writer, taxon.getName(), depths);
        }
        return this;
    }

    /**
     * Add taxon and set values for all sites and alleles for that taxon. First
     * dimension of depths is number of alleles (6 for Nucleotide) and second
     * dimension is sites.
     *
     * @param taxon taxon
     * @param depths depth values
     *
     * @return builder
     */
    public HDF5AlleleDepthBuilder addTaxon(Taxon taxon, int[][] depths) {
        int numAlleles = depths.length;
        if ((depths == null) || (numAlleles != 6)) {
            throw new IllegalStateException("AlleleDepthBuilder: addTaxon: Set A, C, G, T, -, + at once");
        }
        if (depths[0].length != myNumSites) {
            throw new IllegalStateException("AlleleDepthBuilder: addTaxon: Number of sites: " + depths[0].length + " should be: " + myNumSites);
        }
        byte[][] result = new byte[numAlleles][myNumSites];
        for (int a = 0; a < numAlleles; a++) {
            for (int s = 0; s < myNumSites; s++) {
                result[a][s] = AlleleDepthUtil.depthIntToByte(depths[a][s]);
            }
        }
        return addTaxon(taxon, result);
    }

    @Override
    public AlleleDepth build() {
        IHDF5Reader reader = myHDF5Writer;
        myHDF5Writer = null;
        return new HDF5AlleleDepth(reader);
    }
}
