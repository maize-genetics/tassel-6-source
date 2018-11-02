/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package net.maizegenetics.dna.map;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import net.maizegenetics.dna.tag.AbstractTags;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.BaseEncoder;

/**
 * Class that hold genetic position of tags, genetic positions are the position of markers associated with tags
 * @author Fei Lu
 */
public class TagsOnGeneticMap extends AbstractTags {
    protected int[] gChr;
    protected int[] gPos;
    protected byte[] ifPAV;
    protected float[] prediction;
    
    /**
     * Construct TOGM from a file
     * @param infileS
     *        File name of input TagsOnGeneticMap file
     * @param format
     *        FilePacking format
     */
    public TagsOnGeneticMap (String infileS, FilePacking format) {
        this.readDistFile(infileS, format);
    }
    
    public void mergeTOGM (TagsOnGeneticMap another) {
        int totalCnt = this.getTagCount();
        for (int i = 0; i < another.getTagCount(); i++) {
            long[] t = another.getTag(i);
            if (this.getTagIndex(t) < 0) {
                totalCnt++;
            }
        }
        System.out.println("Identified " + String.valueOf(totalCnt- this.getTagCount()) + " new tags");
        long[][] tags = new long[tagLengthInLong][totalCnt];
        byte[] tagLength = new byte[totalCnt];
        int[] gChr = new int[totalCnt];
        int[] gPos = new int[totalCnt];
        byte[] ifPAV = new byte[totalCnt];
        float[] prediction = new float[totalCnt];
        for (int i = 0; i < this.getTagCount(); i++) {
            long[] t = this.getTag(i);
            for (int j = 0; j < this.tagLengthInLong; j++) {
                tags[j][i] = t[j];
            }
            tagLength[i] = (byte)this.getTagLength(i);
            gChr[i] = this.getGChr(i);
            gPos[i] = this.getGPos(i);
            ifPAV[i] = this.getIfPAV(i);
            prediction[i] = this.getPrediction(i);
        }
        int cnt = this.getTagCount();
        for (int i = 0; i < another.getTagCount(); i++) {
            long[] t = another.getTag(i);
            if (this.getTagIndex(t) < 0) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    tags[j][cnt] = t[j];
                }
                tagLength[cnt] = (byte)another.getTagLength(i);
                gChr[cnt] = another.getGChr(i);
                gPos[cnt] = another.getGPos(i);
                ifPAV[cnt] = another.getIfPAV(i);
                prediction[cnt] = another.getPrediction(i);
                cnt++;
            }
        }
        this.tags = tags;
        this.tagLength = tagLength;
        this.gChr = gChr;
        this.gPos = gPos;
        this.ifPAV = ifPAV;
        this.prediction = prediction;
        this.sort();
    }
    
    /**
     * Return chromosome of genetic position
     * @param index
     * @return Chromosome of genetic position
     */
    public int getGChr (int index) {
        return gChr[index];
    }
    
    /**
     * Return site of genetic position
     * @param index
     * @return Genetic position
     */
    public int getGPos (int index) {
        return gPos[index];
    }
    
    /**
     * Return if this tag is PAV
     * @param index
     * @return 
     */
    public byte getIfPAV (int index) {
        return ifPAV[index];
    }
    
    /**
     * Return if this tag is a PAV
     * @param index
     * @return 
     */
    public boolean isPAV (int index) {
        if (ifPAV[index] == 1) return true;
        return false;
    }
    
    /**
     * Return prediction value from model
     * @param index
     * @return 
     */
    public float getPrediction (int index) {
        return prediction[index];
    }
    
    /**
     * Set PAV value of tag
     * @param index
     * @param value 
     */
    public void setIfPAV (int index, int value) {
        this.ifPAV[index] = (byte)value;
    }
    
    /**
     * Read tagsOnGeneticMap file
     * @param infileS
     * @param format
     */
    public void readDistFile (String infileS, FilePacking format) {
        System.out.println("Reading TOGM file from " + infileS);
        File infile = new File (infileS);
        switch (format) {
            case Text:
                readTextTOGMFile(infile);
                break;
            default:
                readBinaryTOGMFile(infile);
                break;
        }
        System.out.println("TOGM file read. Total: " + this.getTagCount() + " Tags");
    }
    
    /**
     * Read text TOGM file
     * @param infile
     */
    private void readTextTOGMFile (File infile) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(infile), 65536);
            br.readLine();
            tagLengthInLong = br.readLine().split("\t")[0].length()/BaseEncoder.chunkSize;
            int tagNum = 1;
            while ((br.readLine())!=null) {
                tagNum++;
            }
            this.iniMatrix(tagLengthInLong, tagNum);
            br = new BufferedReader (new FileReader(infile), 65536);
            br.readLine();
            for (int i = 0; i < tagNum; i++) {
                String[] temp = br.readLine().split("\t");
                long[] t = BaseEncoder.getLongArrayFromSeq(temp[0]);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = t[j];
                }
                tagLength[i] = Byte.parseByte(temp[1]);
                gChr[i] = Integer.parseInt(temp[2]);
                gPos[i] = Integer.parseInt(temp[3]);
                ifPAV[i] = Byte.parseByte(temp[4]);
                prediction[i] = Float.parseFloat(temp[5]);
                if (i%500000 == 0) System.out.println("Read in " + String.valueOf(i+1)+ " tags");
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Read binary TOGM file
     * @param infile
     */
    private void readBinaryTOGMFile (File infile) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(infile), 65536));
            tagLengthInLong = dis.readInt();
            int tagNum = (int)(infile.length()/(8*tagLengthInLong)+1+4+4+1+4);
            this.iniMatrix(tagLengthInLong, tagNum);
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    this.tags[j][i] = dis.readLong();
                }
                this.tagLength[i] = dis.readByte();
                this.gChr[i] = dis.readInt();
                this.gPos[i] = dis.readInt();
                this.ifPAV[i] = dis.readByte();
                this.prediction[i] = dis.readFloat();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Initialize the matrix of TOGM
     * @param tagLengthInLong
     *        Tag length in Long primitive data type
     * @param tagNum
     *        Total tag number
     */
    protected void iniMatrix (int tagLengthInLong, int tagNum) {
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        gChr = new int[tagNum];
        gPos = new int[tagNum];
        ifPAV = new byte[tagNum];
        prediction = new float[tagNum];
    }
    
    /**
     * Write TagsOnGeneticMap file
     * @param outfileS
     *        File name of output file
     * @param format
     *        FilePacking format
     */
    public void writeDistFile (String outfileS, FilePacking format) {
        System.out.println("Writing TOGM file to " + outfileS);
        switch (format) {
            case Text:
                writeTextTOGMFile(outfileS);
                break;
            default:
                writeBinaryTOGMFile(outfileS);
                break;
        }
        System.out.println("TOGM file written");
    }
    
    /**
     * Write TagsOnGeneticMap file
     * @param outfileS
     * @param ifOut
     * @param format 
     */
    public void writeDistFile (String outfileS, boolean[] ifOut, FilePacking format) {
        System.out.println("Writing TOGM file to " + outfileS);
        switch (format) {
            case Text:
                writeTextTOGMFile(outfileS, ifOut);
                break;
            default:
                writeBinaryTOGMFile(outfileS);
                break;
        }
        System.out.println("TOGM file written");
    }
    
    /**
     * Write text TOGM file
     * @param outfileS
     * @param ifOut 
     */
    private void writeTextTOGMFile (String outfileS, boolean[] ifOut) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Tag\tTagLength\tGChr\tGPos\tIfPAV\tPredictedDistance");
            bw.newLine();
            long[] temp = new long[this.tagLengthInLong];
            for (int i = 0; i < this.getTagCount(); i++) {
                if (!ifOut[i]) continue;
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(temp)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                bw.write(String.valueOf(this.gChr[i])+"\t"+String.valueOf(this.gPos[i])+"\t"+String.valueOf(ifPAV[i])+"\t"+String.valueOf(this.prediction[i]));
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
     * Write text TOGM file
     * @param outfileS
     */
    private void writeTextTOGMFile (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Tag\tTagLength\tGChr\tGPos\tIfPAV\tPredictedDistance");
            bw.newLine();
            long[] temp = new long[this.tagLengthInLong];
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(temp)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                bw.write(String.valueOf(this.gChr[i])+"\t"+String.valueOf(this.gPos[i])+"\t"+String.valueOf(ifPAV[i])+"\t"+String.valueOf(this.prediction[i]));
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
     * Write binary TOGM file
     * @param outfileS
     * @param ifOut 
     */
    private void writeBinaryTOGMFile (String outfileS, boolean[] ifOut) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            for (int i = 0; i < this.getTagCount(); i++) {
                if (!ifOut[i]) continue;
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tags[j][i]);
                }
                dos.writeByte(this.getTagLength(i));
                dos.writeInt(this.getGChr(i));
                dos.writeInt(this.getGPos(i));
                dos.writeByte(this.getIfPAV(i));
                dos.writeFloat(this.getPrediction(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Write binary TOGM file
     * @param outfileS
     */
    private void writeBinaryTOGMFile (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeInt(tagLengthInLong);
            for (int i = 0; i < this.getTagCount(); i++) {
                for (int j = 0; j < this.tagLengthInLong; j++) {
                    dos.writeLong(this.tags[j][i]);
                }
                dos.writeByte(this.getTagLength(i));
                dos.writeInt(this.getGChr(i));
                dos.writeInt(this.getGPos(i));
                dos.writeByte(this.getIfPAV(i));
                dos.writeFloat(this.getPrediction(i));
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Output FASTA file of TOGM
     * @param outfileS 
     */
    public void writeFastA (String outfileS) {
        long[] t = new long[this.getTagSizeInLong()];
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                bw.write(">"+String.valueOf(i));
                bw.newLine();
                for (int j = 0; j < this.getTagSizeInLong(); j++) {
                    t[j] = this.tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(t).substring(0, this.getTagLength(i)));
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
     * Output FASTQ file of TOGM
     * @param outfileS 
     */
    public void writeFastQ (String outfileS) {
        String defaultQualityS = null;
        for (int i = 0; i < this.getTagSizeInLong()*32; i++) {
            defaultQualityS+="f";
        }
        long[] t = new long[this.getTagSizeInLong()];
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            for (int i = 0; i < this.getTagCount(); i++) {
                bw.write("@"+String.valueOf(i));
                bw.newLine();
                for (int j = 0; j < this.getTagSizeInLong(); j++) {
                    t[j] = this.tags[j][i];
                }
                bw.write(BaseEncoder.getSequenceFromLong(t).substring(0, this.getTagLength(i)));
                bw.newLine();
                bw.write("+");
                bw.newLine();
                bw.write(defaultQualityS.substring(0, this.getTagLength(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    @Override
    public void swap(int index1, int index2) {
        super.swap(index1, index2);
        int tInt;
        tInt = gChr[index1];
        gChr[index1] = gChr[index2];
        gChr[index2] = tInt;
        tInt = gPos[index1];
        gPos[index1] = gPos[index2];
        gPos[index2] = tInt;
        byte tByte;
        tByte = ifPAV[index1];
        ifPAV[index1] = ifPAV[index2];
        ifPAV[index2] = tByte;
        float tFloat;
        tFloat = prediction[index1];
        prediction[index1] = prediction[index2];
        prediction[index2] = tFloat;
    }
    
    /**
     * Sort TOGM by genetic position
     */
    public void sortByPosition() {
        System.out.println("Sorting by position");
        GenericSorting.quickSort(0, this.getTagCount(), compPos, this);
        System.out.println("Sorting by position finished");
    }
    
    IntComparator compPos = new IntComparator() {   
        @Override
        public int compare(int index1, int index2) {
            if (gChr[index1] < gChr[index2]) return -1;
            else if (gChr[index1] > gChr[index2]) return 1;
            else {
                if (gPos[index1] < gPos[index2]) return -1;
                else if (gPos[index1] < gPos[index2]) return 1;
                return 0;
            }
        }
    };
}
