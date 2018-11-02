/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.map;

/**
 * Stores variables from tag genetic mapping from GWAS. This class is used for I/O of HDF5. 
 * @author Fei Lu
 */
public class TagGWASMapInfo {
    /**Tag count in master tagCount file, unknown = Integer.MIN_VALUE*/
    public int readCount = Integer.MIN_VALUE;
    /**Physical chromosome of tag, unknown = Integer.MIN_VALUE*/
    public int pChr = Integer.MIN_VALUE;
    /**Physical position of tag, unknown = Integer.MIN_VALUE*/
    public int pPos = Integer.MIN_VALUE;
    /**If the tag is mapped by aligner*/
    public boolean ifMap = false;
    /**If the tag is reference tag*/
    public boolean ifRef = false;
    /**If the tag is unique to one position in genome*/
    public boolean ifUnique = false;
    /**Genetic mapping chromosome of a tag, unknown = Integer.MIN_VALUE*/
    public int gChr = Integer.MIN_VALUE;
    /**Genetic mapping position of a tag, unknown  = Integer.MIN_VALUE*/
    public int gPos = Integer.MIN_VALUE;
    /**P-value of genetic mapping of a tag*/
    public double gwasPValue = 1;
    /**Total number of significant site, unknown  = Integer.MIN_VALUE*/
    public int numSigSite = Integer.MIN_VALUE;
    /**Number of taxa where tag exist, unknown  = Integer.MIN_VALUE*/
    public int tagTaxaCount = Integer.MIN_VALUE;
    /**Total number of significant chromosome, unknown  = Integer.MIN_VALUE*/
    public int numSigChr = Integer.MIN_VALUE;
    /**Likelihood ratio of the most significant chromosome Vs the second most significant chromosome, log10 value, unknown  = 0*/
    public double lRatioSB = 0;
    /**Likelihood ratio of the most significant chromosome Vs the median most significant chromosome, log10 value, unknown  = 0*/
    public double lRatioMB = 0;
    /**Number of site on best chromosome having more significant p value than the most significant p value on the second best chromosome, unknown = Integer.MIN_VALUE*/
    public int numSiteOnBestChrThanSecondBest = Integer.MIN_VALUE;
    /**Starting significant site on the best chromosome, unknown  = Integer.MIN_VALUE*/
    public int sigSiteStart = Integer.MIN_VALUE;
    /**Ending significant site on the best chromosome, inclusive, unknown  = Integer.MIN_VALUE*/
    public int sigSiteEnd = Integer.MIN_VALUE;
    /**Predicted distance (log10 value) between mapping position and true position.*/
    public double predictedDistance = Double.NaN;
    
    public TagGWASMapInfo () {}
    
    public TagGWASMapInfo (int readCount, int gChr, int gPos, double gwasPValue, int numSigSite, int tagTaxaCount, int numSigChr,
                             double lRatioSB, double lRatioMB, int numSiteOnBestChrThanSecondBest, int sigSiteStart, int sigSiteEnd) {
        this.readCount = readCount; this.pChr = pChr; this.pPos = pPos; this.ifMap = ifMap; this.ifRef = ifRef; this.ifUnique = ifUnique; this.gChr = gChr; this.gPos = gPos;
        this.gwasPValue = gwasPValue; this.numSigSite = numSigSite; this.tagTaxaCount = tagTaxaCount; this.numSigChr = numSigChr; this.lRatioSB = lRatioSB;
        this.lRatioMB = lRatioMB; this.numSiteOnBestChrThanSecondBest = numSiteOnBestChrThanSecondBest; this.sigSiteStart = sigSiteStart; this.sigSiteEnd = sigSiteEnd;
    }
    
    public TagGWASMapInfo (int readCount, int pChr, int pPos, boolean ifMap, boolean ifRef, boolean ifUnique, int gChr, int gPos, double gwasPValue, int numSigSite, int tagTaxaCount, int numSigChr,
                             double lRatioSB, double lRatioMB, int numSiteOnBestChrThanSecondBest, int sigSiteStart, int sigSiteEnd) {
        this.readCount = readCount; this.pChr = pChr; this.pPos = pPos; this.ifMap = ifMap; this.ifRef = ifRef; this.ifUnique = ifUnique; this.gChr = gChr; this.gPos = gPos;
        this.gwasPValue = gwasPValue; this.numSigSite = numSigSite; this.tagTaxaCount = tagTaxaCount; this.numSigChr = numSigChr; this.lRatioSB = lRatioSB;
        this.lRatioMB = lRatioMB; this.numSiteOnBestChrThanSecondBest = numSiteOnBestChrThanSecondBest; this.sigSiteStart = sigSiteStart; this.sigSiteEnd = sigSiteEnd;
    }
    
    public void setAlignment (int pChr, int pPos, boolean ifMap, boolean ifRef, boolean ifUnique) {
        this.pChr = pChr;
        this.pPos = pPos;
        this.ifMap = ifMap;
        this.ifRef = ifRef;
        this.ifUnique = ifUnique;
    }
    
    public void setPredictedDistance (double predictedDistance) {
        this.predictedDistance = predictedDistance;
    }
    
    public boolean isUniqueRef () {
        if (this.ifRef && this.ifUnique) return true;
        else return false;
    }
    
    public String getBoxcoxAttributesStr (double[] lamdas, String delimiter) {
        StringBuilder sb = new StringBuilder();
        sb.append(this.boxcoxTransform(this.readCount, lamdas[0])).append(delimiter);
        sb.append(this.boxcoxTransform(this.tagTaxaCount, lamdas[1])).append(delimiter);
        sb.append(this.boxcoxTransform(this.getMinusLog10PValue(), lamdas[2])).append(delimiter);
        sb.append(this.boxcoxTransform(this.lRatioSB, lamdas[3])).append(delimiter);
        sb.append(this.boxcoxTransform(this.lRatioMB, lamdas[4])).append(delimiter);
        sb.append(this.boxcoxTransform(this.getAdjustedNumSigChr(), lamdas[5])).append(delimiter);
        sb.append(this.boxcoxTransform(this.getAdjustedNumSigSite(), lamdas[6])).append(delimiter);
        sb.append(this.boxcoxTransform(this.getAdjustedNumSigSiteBC(), lamdas[7])).append(delimiter);
        sb.append(this.boxcoxTransform(this.getSigWidthBC(), lamdas[8])).append(delimiter);
        sb.append(this.getLog10GDist());
        return sb.toString();
    }
    
    public String getAttributesStr (String delimiter) {
        StringBuilder sb = new StringBuilder();
        sb.append(this.readCount).append(delimiter);
        sb.append(this.tagTaxaCount).append(delimiter);
        sb.append(this.getMinusLog10PValue()).append(delimiter);
        sb.append(this.lRatioSB).append(delimiter);
        sb.append(this.lRatioMB).append(delimiter);
        sb.append(this.getAdjustedNumSigChr()).append(delimiter);
        sb.append(this.getAdjustedNumSigSite()).append(delimiter);
        sb.append(this.getAdjustedNumSigSiteBC()).append(delimiter);
        sb.append(this.getSigWidthBC()).append(delimiter);
        sb.append(this.getLog10GDist());
        return sb.toString();
    }
    
    private double getMinusLog10PValue () {
        double value = -Math.log10(this.gwasPValue);
        if (value == Double.POSITIVE_INFINITY) value = 308;
        return value;
    }
    
    private int getAdjustedNumSigSiteBC () {
        if (this.numSiteOnBestChrThanSecondBest == 0) return 1;
        else return this.numSiteOnBestChrThanSecondBest;
    }
    
    private int getAdjustedNumSigSite () {
        if (this.numSigSite == 0) return 1;
        else return this.numSigSite;
    }
    
    private int getAdjustedNumSigChr () {
        if (this.numSigChr == 0) return 1;
        else return this.numSigChr;
    }
    
    private int getSigWidthBC () {
        return Math.abs(this.sigSiteEnd-this.sigSiteStart)+1;
    }
    
    private double getLog10GDist() {
        if (this.gChr == this.pChr) {
            int value = Math.abs(this.gPos-this.pPos);
            if (value < 2) value = 2;
            return Math.log10(value);
        }
        else {
            //return Math.log10((double)Math.abs(this.gChr-this.pChr)*1000000000+Math.abs(this.gPos-this.pPos));
            return Math.log10(Integer.MAX_VALUE);
        }
    }
    
    private double boxcoxTransform (double y, double lambda) {
        if (lambda != 0) {
            return (Math.pow(y, lambda)-1)/lambda;
        }
        else {
            return Math.log(y);
        }
    }
}
