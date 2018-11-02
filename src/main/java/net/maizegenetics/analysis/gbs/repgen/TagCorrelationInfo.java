/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.dna.tag.Tag;

/**
 * @author lcj34
 *
 */
public class TagCorrelationInfo {

    private final Tag tag2;
    private final double t1t2_pearson;
    private final double t1t2_spearman;
    private final double pres_abs_pearson;
    private final double r2;
    
    public TagCorrelationInfo(Tag tag2, double t1t2_pearson, double t1t2_spearman, double pres_abs_pearson, double r2) {
        this.tag2 = tag2;
        this.t1t2_pearson = t1t2_pearson;
        this.t1t2_spearman = t1t2_spearman;
        this.pres_abs_pearson = pres_abs_pearson;
        this.r2 = r2;
    }
    
    public Tag tag2() {
        return tag2;
    }
    
    public double t1t2_pearson() {
        return t1t2_pearson;
    }
    public double t1t2_spearman() {
        return t1t2_spearman;
    }
    public double pres_abs_pearson() {
        return pres_abs_pearson;
    }
    public double r2() {
        return r2;
    }
}
