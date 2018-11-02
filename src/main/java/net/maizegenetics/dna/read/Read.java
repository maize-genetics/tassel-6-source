/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.read;

import java.util.Arrays;
import static net.maizegenetics.dna.read.ReadUtils.baseCompleMap;

/**
 * Class holding reads of Fastq format.
 * @author Fei Lu
 */
public class Read {
    /**sequence ID*/
    String ID;
    /**sequence*/
    String seq;
    /**description*/
    String des;
    /**quality value string*/
    String qual;
    /**quality value byte*/
    byte[] qualValue = null;
    byte[] sortedQualValue = null;
    
    public Read (String ID, String seq, String des, String qual) {
        this.ID = ID;
        this.seq = seq;
        this.des = des;
        this.qual = qual;
    }
    
    /**
     * Return average quality value of read
     * @return 
     */
    public byte returnAverageQuality (int phredScale) {
        double aver = 0;
        if (qualValue == null) this.buildQualByteArray(phredScale);
        for (int i= 0; i < qualValue.length; i++) {
            aver+=qualValue[i];
        }
        return (byte)(aver/this.getReadLength());
    }
    
    /**
     * Return median quality of read
     * @return 
     */
    public byte getMedianQuality (int phredScale) {
        if (qualValue == null) this.buildQualByteArray(phredScale);
        if (sortedQualValue == null) {
            System.arraycopy(qualValue, 0, sortedQualValue, 0, seq.length());
            Arrays.sort(sortedQualValue);
        }
        return sortedQualValue[sortedQualValue.length/2];
    }
    
    /**
     * Return read length
     * @return 
     */
    public int getReadLength () {
        return this.seq.length();
    }
    
    /**
     * Return quality value of base
     * @param index
     * @return 
     */
    public byte getBaseQuality (int index, int phredScale) {
        if (qualValue == null) this.buildQualByteArray(phredScale);
        return qualValue[index];
    }
    
    private void buildQualByteArray(int phredScale) {
        qualValue = qual.getBytes();
        for (int i = 0; i < qualValue.length; i++) {
            qualValue[i] = (byte)(qualValue[i]-phredScale);
        }
    }
    
    /**
     * Return description of read
     * @return 
     */
    public String getDescription () {
        return this.des;
    }
    
    /**
     * Return ID of read
     * @return 
     */
    public String getID () {
        return this.ID;
    }
    
    /**
     * Return read
     * @return 
     */
    public String getSeq () {
        return this.seq;
    }
    
    /**
     * Return reverse quality string
     * @return 
     */
    public String getReverseQual () {
        return new StringBuilder(qual).reverse().toString();
    }
    
    /**
     * Return quality string
     * @return 
     */
    public String getQual () {
        return this.qual;
    }
    
    /**
     * Return quality string from startIndex to endIndex
     * @param startIndex
     * @param endIndex
     * @return 
     */
    public String getQual (int startIndex, int endIndex) {
        return this.qual.substring(startIndex, endIndex);
    }
    
    /**
     * Return read from startIndex to endIndex
     * @param startIndex
     * @param endIndex
     * @return 
     */
    public String getSeq (int startIndex, int endIndex) {
        return this.seq.substring(startIndex, endIndex);
    }
    
   /**
    * Return reverse complementary sequence
    * @return 
    */
    public String getReverseComplementarySeq () {
        return this.getReverseComplementarySeq(0, this.seq.length());
    }
    
    /**
     * Return reverse complementary sequence
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return 
     */
    public String getReverseComplementarySeq (int startIndex, int endIndex) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < endIndex - startIndex; i++) {
            sb.append(baseCompleMap.get(String.valueOf(this.seq.charAt(i+startIndex))));
        }
        return sb.reverse().toString();
    }
}
