/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.tag;

import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import java.io.BufferedWriter;
import java.io.FileWriter;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.util.Tassel5HDF5Constants;

/**
 * Basic implementations of HDF5 tags. This is designed to annotate tags with a bunch of attributes, to solve the memory issues
 * @author Fei Lu
 */
public abstract class AbstractTagsHDF5 extends AbstractTags {
    protected IHDF5Writer h5 = null;
    protected int currentBlockIndex = -1;
    protected int currentIndex = -1;
    
    /**
     * Returns the size of HDF5 block
     * @return 
     */
    public int getBlockSize() {
        return Tassel5HDF5Constants.BLOCK_SIZE;
    }

    /**
     * Returns the number of block
     * @return 
     */
    public int getBlockNum() {
        int num = super.getTagCount()/this.getBlockSize();
        if (super.getTagCount()%this.getBlockSize() == 0) return num;
        else return num+1;
    }

    /**
     * Returns the block index based on the tag index in a full list
     * @param tagIndex as the query tag index in the full list
     * @return 
     */
    public int getBlockIndex(int tagIndex) {
        return tagIndex/this.getBlockSize();
    }

    /**
     * Returns current block index
     * @return 
     */
    public int getCurrentBlockIndex() {
        return currentBlockIndex;
    }

    /**
     * Returns current index
     * @return 
     */
    public int getCurrentIndex() {
        return currentIndex;
    } 

    /**
     * Returns current index relative to current block
     * @return 
     */
    public int getCurrentIndexWithinBlock () {
        return this.getCurrentIndex()%this.getBlockSize();
    }
    
    /**
     * Initialize tag matrix
     * @param tagCount
     * @param tagLengthInLong 
     */
    protected void initializeMatrix (int tagCount, int tagLengthInLong) {
        this.tagLengthInLong = tagLengthInLong;
        tags = new long[tagLengthInLong][tagCount];
        tagLength = new byte[tagCount];
    }
    
    /**
     * Return boolean value if an index belongs to current block
     * @param queryIndex
     * @return 
     */
    public boolean isInCurrentBlock(int queryIndex) {
        int queryBlockIndex = queryIndex/this.getBlockSize();
        if (queryBlockIndex == currentBlockIndex) return true;
        else return false;
    }
    
    /**
     * Returns default deflation level of boolean, byte, short, int, long, Enum
     * See {@link net.maizegenetics.util.Tassel5HDF5Constants}
     * @return 
     */
    public HDF5IntStorageFeatures getIntStorageFeatures() {
        return Tassel5HDF5Constants.intDeflation;
    }

    /**
     * Returns default deflation level of String and Class
     * See {@link net.maizegenetics.util.Tassel5HDF5Constants}
     * @return 
     */
    public HDF5GenericStorageFeatures getGenericStorageFeatures() {
        return Tassel5HDF5Constants.genDeflation;
    }

    /**
     * Returns default deflation level of float and double
     * See {@link net.maizegenetics.util.Tassel5HDF5Constants}
     * @return 
     */
    public HDF5FloatStorageFeatures getFloatStorageFeatures() {
        return Tassel5HDF5Constants.floatDeflation;
    }
    
    @Override
    public void sort () {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void swap(int index1, int index2) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public void writeToFasta (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            long[] t;
            for (int i = 0; i < this.getTagCount(); i++) {
                bw.write(">"+String.valueOf(i));
                bw.newLine();
                t = this.getTag(i);
                bw.write(BaseEncoder.getSequenceFromLong(t).substring(0, this.getTagLength(i)));
                bw.newLine();
                if (i%100000 == 0) System.out.println("output " + String.valueOf(i+1) + " tags to Fasta");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Initialize a HDF5 file
     * @param inputFileS 
     */
    public void initializeHDF5 (String inputFileS) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    /**
     * Read in HDF5 file
     * @param hdf5FileS 
     */
    public void readHDF5 (String hdf5FileS) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    /**
     * Populate a block in memory with default values, update current block index at the same time
     * @param blockIndex
     */
    public void populateBlock(int blockIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    /**
     * Read in a block from HDF5 file
     * @param blockIndex 
     */
    public void readBlock(int blockIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    /**
     * Write current block to HDF5 file
     * @param blockIndex 
     */
    public void writeBlock(int blockIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
