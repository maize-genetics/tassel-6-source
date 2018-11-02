/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.map;

import ch.systemsx.cisd.hdf5.HDF5CompoundType;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.map.TagMappingInfoV3.Aligner;
import net.maizegenetics.dna.tag.AbstractTagsHDF5;
import net.maizegenetics.dna.tag.GBSHDF5Constants;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;

/**
 * Holding tag genetic mapping result from GWAS. It includes attributes and methods for machine learning prediction of mapping accuracy 
 * @author Fei Lu
 */
public class TagGWASMap extends AbstractTagsHDF5 {
    HDF5CompoundType<TagGWASMapInfo> tgType = null;
    TagGWASMapInfo[] mapInfo;
    
    public TagGWASMap (String tagGWASMapFileS) {
        this.readHDF5(tagGWASMapFileS);
    }
    
    public TagGWASMap (String gwasMappingResultFileS, String tagCountFileS, String tagGWASMapFileS) {
        this.creatFile(gwasMappingResultFileS, tagCountFileS, tagGWASMapFileS);
    }
    
    public TagGWASMap(String gwasMappingFileS, String topmFileS, String tagCountFileS, String tagGWASMapFileS) {
        this.creatFile(gwasMappingFileS, topmFileS, tagCountFileS, tagGWASMapFileS);
    }
    
    private void creatFile (String gwasMappingResultFileS, String tagCountFileS, String tagGWASMapFileS) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(gwasMappingResultFileS), 65536);
            String temp = null;
            int tagCount = -1;
            tagLengthInLong = 0;
            while ((temp = br.readLine()) != null) {
                if (tagCount == 0) {
                    String[] tem = temp.split("\t");
                    tagLengthInLong = tem[0].length()/BaseEncoder.chunkSize;
                }
                tagCount++;
            }
            this.initializeMatrix(tagCount, tagLengthInLong);
            TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
            br = new BufferedReader (new FileReader(gwasMappingResultFileS), 65536);
            br.readLine();
            int rLen = tagLengthInLong*BaseEncoder.chunkSize;
            long[] t = new long[tagLengthInLong];
            int tagIndex;
            for (int i = 0; i < this.getTagCount(); i++) {
                temp = br.readLine().substring(0, rLen);
                t = BaseEncoder.getLongArrayFromSeq(temp);
                tagIndex = tc.getTagIndex(t);
                for (int j = 0; j < tagLengthInLong; j++) this.tags[j][i] = t[j];
                this.tagLength[i] = (byte)tc.getTagLength(tagIndex);
            }
            this.initializeHDF5(tagGWASMapFileS);
            br = new BufferedReader (new FileReader(gwasMappingResultFileS), 65536);
            br.readLine();
            String[] tem;
            for (int i = 0; i < this.getBlockNum(); i++) {
                this.populateBlock(i);
                int cnt = 0;
                for (int j = 0; j < this.getBlockSize() && (i*this.getBlockSize()+j) < this.getTagCount(); j++) {
                    for (int k = 0; k < this.getTagSizeInLong(); k++)  {
                        t[k] = this.tags[k][cnt];
                    }
                    tem = br.readLine().split("\t");
                    int readCount = Integer.valueOf(tem[1]);
                    int gChr = Integer.valueOf(tem[5]);
                    int gPos = Integer.valueOf(tem[7]);
                    double gwasPValue = Double.valueOf(tem[8]);
                    int numSigSite = Integer.valueOf(tem[9]);
                    int tagTaxaCount = Integer.valueOf(tem[10]);
                    int numSigChr = Integer.valueOf(tem[11]);
                    double lRatioSB = Double.valueOf(tem[12]);
                    if (lRatioSB == Double.POSITIVE_INFINITY) lRatioSB = 305; //305 is the max likelihood observed
                    else if (lRatioSB == Double.NEGATIVE_INFINITY || lRatioSB == 0) lRatioSB = 0.0000001; //Double.NEGATIVE_INFINITY means P-value on the second best chr is also 0, which indicate strong population structure
                    double lRatioMB = Double.valueOf(tem[13]);
                    if (lRatioMB == Double.POSITIVE_INFINITY) lRatioMB = 305;
                    else if (lRatioMB == Double.NEGATIVE_INFINITY || lRatioMB == 0) lRatioMB = 0.0000001;
                    int numSiteOnBestChrThanSecondBest = Integer.valueOf(tem[14]);
                    int sigSiteStart = Integer.valueOf(tem[15]);
                    int sigSiteEnd = Integer.valueOf(tem[16]);
                    mapInfo[cnt] = new TagGWASMapInfo(readCount, gChr, gPos, gwasPValue, numSigSite, tagTaxaCount, numSigChr,
                             lRatioSB, lRatioMB, numSiteOnBestChrThanSecondBest, sigSiteStart, sigSiteEnd);
                    cnt++;
                }
                this.writeBlock(i);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void creatFile (String gwasMappingFileS, String topmFileS, String tagCountFileS, String tagGWASMapFileS) {
        Aligner software = Aligner.Bowtie2;
        try {
            BufferedReader br = new BufferedReader (new FileReader(gwasMappingFileS), 65536);
            String temp = null;
            int tagCount = -1;
            int tagLengthInLong = 0;
            while ((temp = br.readLine()) != null) {
                if (tagCount == 0) {
                    String[] tem = temp.split("\t");
                    tagLengthInLong = tem[0].length()/32;
                }
                tagCount++;
            }
            this.initializeMatrix(tagCount, tagLengthInLong);
            TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
            br = new BufferedReader (new FileReader(gwasMappingFileS), 65536);
            br.readLine();
            int rLen = tagLengthInLong*32;
            long[] t = new long[tagLengthInLong];
            int tagIndex;
            for (int i = 0; i < this.getTagCount(); i++) {
                temp = br.readLine().substring(0, rLen);
                t = BaseEncoder.getLongArrayFromSeq(temp);
                tagIndex = tc.getTagIndex(t);
                for (int j = 0; j < tagLengthInLong; j++) this.tags[j][i] = t[j];
                this.tagLength[i] = (byte)tc.getTagLength(tagIndex);
            }
            this.initializeHDF5(tagGWASMapFileS);
            br = new BufferedReader (new FileReader(gwasMappingFileS), 65536);
            br.readLine();
            String[] tem;
            TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
            int[] mapIndices = topm.getMappingIndicesOfAligner(software);
            TagMappingInfoV3[] pMaps = new TagMappingInfoV3[mapIndices.length];
            boolean ifMap;
            boolean ifRef;
            boolean ifUnique;
            int cnt = 0;
            for (int i = 0; i < this.getBlockNum(); i++) {
                this.populateBlock(i);
                for (int j = 0; j < this.getBlockSize() && (i*this.getBlockSize()+j) < this.getTagCount(); j++) {
                    ifMap = false;
                    ifRef = false;
                    ifUnique = false;
                    for (int k = 0; k < this.getTagSizeInLong(); k++)  {
                        t[k] = this.tags[k][cnt];
                    }
                    tagIndex = topm.getTagIndex(t);
                    for (int k = 0; k < pMaps.length; k++) {
                        pMaps[k] = topm.getMappingInfo(tagIndex, mapIndices[k]);
                    }
                    if (pMaps[0].chromosome != Integer.MIN_VALUE) {
                        ifMap = true;
                        if (pMaps[0].perfectMatch == 1) ifRef = true;
                        if (pMaps[1].chromosome == Integer.MIN_VALUE) ifUnique = true;
                    }
                    tem = br.readLine().split("\t");
                    int readCount = Integer.valueOf(tem[1]);
                    int pChr = pMaps[0].chromosome;
                    int pPos = pMaps[0].startPosition;
                    int gChr = Integer.valueOf(tem[5]);
                    int gPos = Integer.valueOf(tem[7]);
                    double gwasPValue = Double.valueOf(tem[8]);
                    int numSigSite = Integer.valueOf(tem[9]);
                    int tagTaxaCount = Integer.valueOf(tem[10]);
                    int numSigChr = Integer.valueOf(tem[11]);
                    double lRatioSB = Double.valueOf(tem[12]);
                    if (lRatioSB == Double.POSITIVE_INFINITY) lRatioSB = 305; //305 is the max likelihood observed
                    else if (lRatioSB == Double.NEGATIVE_INFINITY) lRatioSB = -305; //Double.NEGATIVE_INFINITY means P-value on the second best chr is also 0, which indicate strong population structure
                    double lRatioMB = Double.valueOf(tem[13]);
                    if (lRatioMB == Double.POSITIVE_INFINITY) lRatioMB = 305;
                    else if (lRatioMB == Double.NEGATIVE_INFINITY) lRatioMB = -305;
                    int numSiteOnBestChrThanSecondBest = Integer.valueOf(tem[14]);
                    int sigSiteStart = Integer.valueOf(tem[15]);
                    int sigSiteEnd = Integer.valueOf(tem[16]);
                    mapInfo[cnt] = new TagGWASMapInfo(readCount, pChr, pPos, ifMap, ifRef, ifUnique, gChr, gPos, gwasPValue, numSigSite, tagTaxaCount, numSigChr,
                             lRatioSB, lRatioMB, numSiteOnBestChrThanSecondBest, sigSiteStart, sigSiteEnd);
                    cnt++;
                }
                this.writeBlock(i);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void addAlignment (TagsOnPhysicalMapV3 topm) {
        Aligner software = Aligner.Bowtie2;
        int[] mapIndices = topm.getMappingIndicesOfAligner(software);
        TagMappingInfoV3[] pMaps = new TagMappingInfoV3[mapIndices.length];
        for (int i = 0; i < this.getTagCount(); i++) {
            boolean ifMap = false;
            boolean ifRef = false;
            boolean ifUnique = false;
            long[] t =new long[this.getTagSizeInLong()];
            for (int j = 0; j < this.getTagSizeInLong(); j++)  {
                t[j] = this.tags[j][i];
            }
            int tagIndex = topm.getTagIndex(t);
            for (int j = 0; j < pMaps.length; j++) {
                pMaps[j] = topm.getMappingInfo(tagIndex, mapIndices[j]);
            }
            int pChr = pMaps[0].chromosome;
            int pPos = pMaps[0].startPosition;
            if (pMaps[0].chromosome != Integer.MIN_VALUE) {
                ifMap = true;
                if (pMaps[0].perfectMatch == 1) ifRef = true;
                if (pMaps[1].chromosome == Integer.MIN_VALUE) ifUnique = true;
            }
            this.getTagGWASMapInfo(i).setAlignment(pChr, pPos, ifMap, ifRef, ifUnique);
            if ((i+1)%this.getBlockSize() == 0) {
                this.writeBlock(this.getCurrentBlockIndex());
                System.out.println("Added alignment to " + String.valueOf(i+1) + " tags");
            }
        }
        this.writeBlock(this.getCurrentBlockIndex());
        System.out.println("Added alignment to " + String.valueOf(this.getTagCount()) + " tags");
        System.out.println("Adding alignment is completed");
    }
    
    @Override
    public void initializeHDF5 (String tagGWASMapFileS) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(tagGWASMapFileS));
        config.overwrite();
        config.useUTF8CharacterEncoding();
        h5 = config.writer();
        h5.int32().setAttr(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG, tagLengthInLong);
        h5.int32().setAttr(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT, this.getTagCount());
        h5.int64().createMatrix(GBSHDF5Constants.TAGS, this.getTagSizeInLong(), this.getTagCount(), this.getTagSizeInLong(), this.getTagCount(), this.getIntStorageFeatures());
        h5.int64().writeMatrix(GBSHDF5Constants.TAGS, tags, this.getIntStorageFeatures());
        System.out.println("...Tags written");
        h5.int8().createArray(GBSHDF5Constants.TAGLENGTH, this.getTagCount(), this.getIntStorageFeatures());
        h5.int8().writeArray(GBSHDF5Constants.TAGLENGTH, tagLength, this.getIntStorageFeatures());
        System.out.println("...Tags lengths written");
        tgType = h5.compound().getInferredType(TagGWASMapInfo.class);
        h5.compound().createArray(GBSHDF5Constants.MAPBASE, tgType, this.getBlockSize()*this.getBlockNum(), this.getBlockSize(), this.getGenericStorageFeatures());
        System.out.println("...MapInfo created");
    }
    
    @Override
    public void readHDF5 (String hdf5FileS) {
        h5 = HDF5Factory.open(hdf5FileS);
        tgType = h5.compound().getInferredType(TagGWASMapInfo.class);
        tagLengthInLong = h5.int32().getAttr(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG);
        int tagCount = h5.int32().getAttr(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT);
        tags = h5.readLongMatrix(GBSHDF5Constants.TAGS);
        tagLength = h5.int8().readArray(GBSHDF5Constants.TAGLENGTH);
        this.getTagGWASMapInfo(0);
    }
    
    @Override
    public void populateBlock(int blockIndex) {
        mapInfo = new TagGWASMapInfo[this.getBlockSize()];
        for (int i = 0; i < this.getBlockSize(); i++) {
            mapInfo[i] = new TagGWASMapInfo();
        }
        this.currentBlockIndex = blockIndex;
    }

    @Override
    public void readBlock(int blockIndex) {
        this.mapInfo = h5.compound().readArrayBlock(GBSHDF5Constants.MAPBASE, tgType, this.getBlockSize(), blockIndex);
        this.currentBlockIndex = blockIndex;
    }

    @Override
    public void writeBlock(int blockIndex) {
        h5.compound().writeArrayBlock(GBSHDF5Constants.MAPBASE, tgType, mapInfo, blockIndex);
    }
    
    public TagGWASMapInfo getTagGWASMapInfo (int tagIndex) {
        if (this.isInCurrentBlock(tagIndex)) {
            this.currentIndex = tagIndex;
        }
        else {
            this.readBlock(this.getBlockIndex(tagIndex));
            this.currentIndex = tagIndex;
        }
        return this.mapInfo[this.getCurrentIndexWithinBlock()];
    }
}
