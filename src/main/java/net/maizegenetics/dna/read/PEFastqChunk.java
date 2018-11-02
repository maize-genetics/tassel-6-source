/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.read;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 * Hold PE Fastq file. Providing merging function
 * @author Fei Lu
 */
public class PEFastqChunk {
    PERead[] peReads = null;
    int phredScale = Integer.MIN_VALUE;
    
    /**
     * Constructor to sample PE Fastq files, ignore those bad sequence at the beginning
     * @param fastqR1FileS
     * @param fastqR2FileS
     * @param format
     * @param startIndex
     * @param readNum 
     */
    public PEFastqChunk (String fastqR1FileS, String fastqR2FileS, ReadUtils.ReadFormat format, int startIndex, int readNum) {
        FastqChunk r1c = new FastqChunk (fastqR1FileS, format, startIndex, readNum);
        FastqChunk r2c = new FastqChunk (fastqR2FileS, format, startIndex, readNum);
        this.convert(r1c, r2c);
        this.phredScale = r1c.getPhredScale();
    }
    
    /**
     * Constructor to read in whole PE Fastq files, fastq file should be small for test
     * @param fastqR1FileS
     * @param fastqR2FileS
     * @param format 
     */
    public PEFastqChunk (String fastqR1FileS, String fastqR2FileS, ReadUtils.ReadFormat format) {
        FastqChunk r1c = new FastqChunk (fastqR1FileS, format);
        FastqChunk r2c = new FastqChunk (fastqR2FileS, format);
        this.convert(r1c, r2c);
        this.phredScale = r1c.getPhredScale();
    }
    
    /**
     * Return phred score scale of the fastq, 33 or 64
     * @return 
     */
    public int getPhredScale () {
        return this.phredScale;
    }
    
    private void convert (FastqChunk r1c, FastqChunk r2c) {
        if (r1c.getReadNum() == 0) return;
        peReads = new PERead[r1c.getReadNum()];
        for (int i = 0; i < r1c.getReadNum(); i++) {
            peReads[i] = new PERead(r1c.getRead(i), r2c.getRead(i));
        }
        System.out.println("PEFastqChunk built");
    }
    
    /**
     * Return PE read
     * @param index
     * @return 
     */
    public PERead getPERead (int index) {
        return peReads[index];
    }
    
    /**
     * Return number of PE read
     * @return 
     */
    public int getPEReadNum () {
        if (peReads == null) return 0;
        return peReads.length;
    }
    
    /**
     * Merge PE to Pcontig
     * @param ifPrintAlignment 
     */
    public void merge (boolean ifPrintAlignment) {
        System.out.println("Merging PE reads");
        int cnt = 0;
        for (int i = 0; i < this.getPEReadNum(); i++) {
            if (peReads[i].merge(ifPrintAlignment)) cnt++;
        }
        System.out.println(String.valueOf(cnt) + " out of " + String.valueOf(this.getPEReadNum()) + " (" +String.valueOf((double)cnt/this.getPEReadNum())+") PE reads are merged");
    }
    
    /**
     * Write original PE sequence
     * @param fastqR1FileS
     * @param fastqR2FileS
     * @param format 
     */
    public void writePEFastq (String fastqR1FileS, String fastqR2FileS, ReadUtils.ReadFormat format) {
        Read[] readsR1 = new Read[this.getPEReadNum()];
        Read[] readsR2 = new Read[this.getPEReadNum()];
        for (int i = 0; i < this.getPEReadNum(); i++) {
            readsR1[i] = this.getPERead(i).getForwardRead();
            readsR2[i] = this.getPERead(i).getBackwardRead();
        }
        new FastqChunk(readsR1).writeFastq(fastqR1FileS, format);
        new FastqChunk(readsR2).writeFastq(fastqR2FileS, format);
    }
    
    /**
     * Write PE read and Pcontig read, when there is a Pcontig, PE read will not be output
     * @param fastqR1FileS
     * @param fastqR2FileS
     * @param fastqContigFileS
     * @param format 
     */
    public void writeMergedPEFastq (String fastqR1FileS, String fastqR2FileS, String fastqContigFileS, ReadUtils.ReadFormat format) {
        ArrayList<Read> r1List= new ArrayList();
        ArrayList<Read> r2List= new ArrayList();
        ArrayList<Read> contigList= new ArrayList();
        for (int i = 0; i < this.getPEReadNum(); i++) {
            if (this.getPERead(i).getPContig() == null) {
                r1List.add(this.getPERead(i).getForwardRead());
                r2List.add(this.getPERead(i).getBackwardRead());
            }
            else {
                contigList.add(this.getPERead(i).getPContig());
            }
        }
        new FastqChunk(r1List.toArray(new Read[r1List.size()])).writeFastq(fastqR1FileS, format);
        new FastqChunk(r2List.toArray(new Read[r2List.size()])).writeFastq(fastqR2FileS, format);
        new FastqChunk(contigList.toArray(new Read[contigList.size()])).writeFastq(fastqContigFileS, format);
    }
}
