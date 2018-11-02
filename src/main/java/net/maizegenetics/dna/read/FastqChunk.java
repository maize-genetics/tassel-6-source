/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.read;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;
import net.maizegenetics.dna.read.ReadUtils.ReadFormat;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 * Hold Fastq single end Fastq file
 * @author Fei Lu
 */
public class FastqChunk {
    private int minStartIndex = 100000;
    private int maxReadNum = 1000000;
    Read[] reads = null;
    int phredScale = Integer.MIN_VALUE;
    
    /**
     * Constructor, sample Fastq file, ignore those bad sequence at the beginning
     * @param fastqFileS
     * @param format
     * @param startIndex
     * @param readNum 
     */
    public FastqChunk (String fastqFileS, ReadFormat format, int startIndex, int readNum) {
        if (startIndex < minStartIndex) {
            startIndex = minStartIndex;
            System.out.println("Start index of read was set to " + String.valueOf(startIndex));
        }
        if (readNum > maxReadNum) {
            readNum = maxReadNum;
            System.out.println("Number of read was set to " + String.valueOf(readNum));
        }
        this.readFastq(fastqFileS, format, startIndex, readNum);
        this.setPhredScale();
    }
    
    /**
     * Constructor to read in whole Fastq, fastq file should be small for test
     * @param fastqFileS
     * @param format 
     */
    public FastqChunk (String fastqFileS, ReadFormat format) {
        this.readFastq(fastqFileS, format);
        this.setPhredScale();
    }
    
    public FastqChunk (Read[] reads) {
        this.reads = reads;
        this.setPhredScale();
    }
    
    private void setPhredScale () {
        int size = 10;
        if (this.getReadNum() < 10) size = this.getReadNum();
        for (int i = 0; i < size; i++) {
            byte[] qualB = this.getRead(i).getQual().getBytes();
            for (int j = this.getRead(i).getReadLength()-1; j > -1; j--) {
                if (qualB[j] < 65) {
                    this.phredScale = 33;
                    return;
                }
            }
        }
        this.phredScale = 64;
    }
    
    /**
     * Return phred score scale of the fastq, 33 or 64
     * @return 
     */
    public int getPhredScale () {
        return this.phredScale;
    }
    
    private void readFastq (String fastqFileS, ReadFormat format) {
        System.out.println("Reading fastq file from " + fastqFileS);
        BufferedReader br = null;
        int cnt = 0;
        String temp;
        try {
            if (format == ReadFormat.FastqText) {
                br = new BufferedReader(new FileReader(fastqFileS), 65536);
                while ((temp = br.readLine()) != null) cnt++;
                br = new BufferedReader(new FileReader(fastqFileS), 65536);
            }
            else if (format == ReadFormat.FastqGzip) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fastqFileS))));
                while ((temp = br.readLine()) != null) cnt++;
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fastqFileS))));
            }
            else {}
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.readFastq(br, 0, cnt/4);
    }
    
    private void readFastq (String fastqFileS, ReadFormat format, int startIndex, int readNum) {
        if (readNum <= 0) return;
        BufferedReader br = null;
        try {
            if (format == ReadFormat.FastqText) {
                br = new BufferedReader(new FileReader(fastqFileS), 65536);
            }
            else if (format == ReadFormat.FastqGzip) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(fastqFileS))));
            }
            else {}
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Reading fastq file from " + fastqFileS);
        this.readFastq(br, startIndex, readNum);
    }
    
    private void readFastq (BufferedReader br, int startIndex, int readNum) {
        this.reads = new Read[readNum];
        try {
            Read r = null;
            String temp;
            int index = 0;
            int readCount = 0;
            while ((temp = br.readLine()) != null) {
                if (index >= startIndex) {
                    r = new Read (temp, br.readLine(), br.readLine(), br.readLine());
                    reads[readCount] = r;
                    readCount++;
                    if (readCount == readNum) break;
                }
                index++;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(readNum) + " reads imported");
    }
    
    /**
     * Write Fastq file
     * @param outputFileS
     * @param format 
     */
    public void writeFastq (String outputFileS, ReadFormat format) {
        BufferedWriter bw = null;
        try {
            if (format == ReadFormat.FastqText) {
                bw = new BufferedWriter(new FileWriter(outputFileS), 65536);
            }
            else if (format == ReadFormat.FastqGzip) { 
                bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileS), 65536)), 65536);
            }
            else {}
            for (int i = 0; i < this.getReadNum(); i++) {
                bw.write(reads[i].getID());
                bw.newLine();
                bw.write(reads[i].getSeq());
                bw.newLine();
                bw.write(reads[i].getDescription());
                bw.newLine();
                bw.write(reads[i].getQual());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Fastq file written to " + outputFileS);
    }
    
    public void writeFasta (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(outfileS), 65536);
            for (int i = 0; i < this.getReadNum(); i++) {
                bw.write(">"+this.reads[i].ID);
                bw.newLine();
                bw.write(this.reads[i].seq);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Return Fastq read
     * @param index
     * @return 
     */
    public Read getRead (int index) {
        return reads[index];
    }
    
    /**
     * Return number of Fastq read
     * @return 
     */
    public int getReadNum () {
        if (reads == null) return 0;
        return reads.length;
    }
}
