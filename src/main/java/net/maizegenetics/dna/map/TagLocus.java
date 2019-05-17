/*
 * TagLocus
 */
package net.maizegenetics.dna.map;

import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.dna.snp.io.VCFUtil;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.log4j.Logger;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.Profile;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * Aligns and calls SNPs for all tags at a given locus.
 * 
 * @author jcg233
 */
public class TagLocus {
	private static final Logger myLogger = Logger.getLogger(TagLocus.class);

    ArrayList<SingleTagByTaxa> theTags = new ArrayList<SingleTagByTaxa>();
    private int minStartPosition;
    private int maxStartPosition;
    private int minTagLength;
    private int maxTagLength;
    private int chromosome;
    private byte strand;
    private int indexOfRef = Integer.MIN_VALUE;
    private int[] tagIndices = null;  // redirect from aligned tag indices to index in theTags
    private int[] positionsOfVariableSites;
    private byte[][] allelesAtVariableSitesByTag;
    private byte[] refCallsBySite = null;
    private int nTaxaCovered = Integer.MIN_VALUE;
    private int totalNReads = Integer.MIN_VALUE;
    private String status = "notSet";
    
    // For VCF output with depth (now used for custom SNP report in regular pipeline as well)
    private byte[][] myCommonAlleles = null;
    private byte[][][] myAlleleDepthsInTaxa = null;
        
    private final static int maxSNPsPerLocus = 64;
    private final static int maxAlignmentSize = 10000;
    private final static int maxCountAtGeno = 500;
    private final static int maxNumAlleles = 3;
    private SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
    private SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
    private static int[] likelihoodRatioThreshAlleleCnt = null;  // index = sample size; value = min count of less tagged allele for likelihood ratio > 1
    // if less tagged allele has counts < likelihoodRatioThreshAlleleCnt[totalCount], call it a homozygote
    // where likelihood ratio = (binomial likelihood het) / (binomial likelihood all less tagged alleles are errors)

    static void setLikelihoodThresh(double errorRate) {   // initialize the likelihood ratio cutoffs for quantitative SNP calling
        likelihoodRatioThreshAlleleCnt = new int[maxCountAtGeno];
        System.out.println("\n\nInitializing the cutoffs for quantitative SNP calling likelihood ratio (pHet/pErr) >1\n");
        System.out.println("totalReadsForSNPInIndiv\tminLessTaggedAlleleCountForHet");
        for (int trials = 0; trials < 2; ++trials) {
            likelihoodRatioThreshAlleleCnt[trials] = 1;
        }
        int lastThresh = 1;
        for (int trials = 2; trials < likelihoodRatioThreshAlleleCnt.length; ++trials) {
            BinomialDistribution binomHet = new BinomialDistribution(trials, 0.5);
            BinomialDistribution binomErr = new BinomialDistribution(trials, errorRate);
            double LikeRatio;
            try {
                LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                while (LikeRatio <= 1.0) {
                    ++lastThresh;
                    LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                }
                likelihoodRatioThreshAlleleCnt[trials] = lastThresh;
                System.out.println(trials + "\t" + lastThresh);
            } catch (Exception e) {
                System.err.println("Error in the TagsAtLocus.BinomialDistributionImpl");
            }
        }
        System.out.println("\n");
    }

    public TagLocus(int chromosome, byte strand, int startPosition, int tagLength, 
            boolean includeRefGenome, boolean fuzzyStartPositions, double errorRate) {
        this.chromosome = chromosome;
        this.strand = (includeRefGenome && fuzzyStartPositions) ? 1 : strand;
        this.minStartPosition = startPosition;
        this.maxStartPosition = startPosition;
        this.minTagLength = tagLength;
        this.maxTagLength = tagLength;
        positionsOfVariableSites = null;
        allelesAtVariableSitesByTag = null;
        if (likelihoodRatioThreshAlleleCnt == null) {
            setLikelihoodThresh(errorRate);
        }
    }

    public void addTag(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, 
            boolean includeRefGenome, boolean fuzzyStartPositions) {
        SingleTagByTaxa singleTBT = new SingleTagByTaxa(tagTOPMIndex, theTOPM, theTBT, 
                includeRefGenome, fuzzyStartPositions);
        if (singleTBT.taxaWithTag > 0) {
            theTags.add(singleTBT);
            if (singleTBT.startPosition > minStartPosition) {
                maxStartPosition = singleTBT.startPosition;
            }
            if (singleTBT.tagLength < minTagLength) {
                minTagLength = singleTBT.tagLength;
            }
            if (singleTBT.tagLength > maxTagLength) {
                maxTagLength = singleTBT.tagLength;
            }
        }
    }
    
    public void addRefTag(String refTag, int nLongsPerTag, String nullTag) {
        theTags.add(new SingleTagByTaxa(minStartPosition, strand, refTag, nLongsPerTag, nullTag));
        indexOfRef = theTags.size()-1;
    }

    public int getSize() {
        return theTags.size();
    }

    public int getChromosome() {
        return chromosome;
    }

    public byte getStrand() {
        return strand;
    }

    public int getMinStartPosition() {
        return minStartPosition;
    }

    public int getMaxStartPosition() {
        return maxStartPosition;
    }
    
    public int getMinTagLength() {
        return minTagLength;
    }
    
    public int getMaxTagLength() {
        return maxTagLength;
    }

    public void setMinStartPosition(int newMinStartPosition) {
        minStartPosition = newMinStartPosition;
    }

    public int getTOPMIndexOfTag(int tagIndex) {
        return theTags.get(tagIndex).tagTOPMIndex;
    }

    public int getTBTIndexOfTag(int tagIndex) {
        return theTags.get(tagIndex).tagTBTIndex;
    }

    public int getDivergenceOfTag(int tagIndex) {
        return theTags.get(tagIndex).divergence;
    }

    public byte getCallAtVariableSiteForTag(int site, int tagIndex) {
        return allelesAtVariableSitesByTag[site][tagIndex];
    }
    
    public byte getRefGeno(int site) {
        if (refCallsBySite == null || site > refCallsBySite.length-1) {
            return GenotypeTable.UNKNOWN_GENOTYPE;
        } else {
            return refCallsBySite[site];
        }
    }

    public int getNumberTaxaCovered() {
        if (theTags.size() < 1) {
            return 0;
        }
        if (nTaxaCovered == Integer.MIN_VALUE) {
            nTaxaCovered = 0;
            totalNReads = 0;
            boolean[] covered = new boolean[theTags.get(0).tagDist.length];  // initializes to false
            for (SingleTagByTaxa sTBT : theTags) {
                for (int tx = 0; tx < covered.length; ++tx) {
                    int reads = sTBT.tagDist[tx];
                    totalNReads += reads;
                    if (!covered[tx] && reads > 0) {
                        covered[tx] = true;
                    }
                }
            }
            for (int tx = 0; tx < covered.length; ++tx) {
                if (covered[tx]) {
                    ++nTaxaCovered;
                }
            }
            return nTaxaCovered; 
        } else {
            return nTaxaCovered;
        }
    }

    public int getTotalNReads() {
        if (theTags.size() < 1) {
            return 0;
        }
        if (totalNReads == Integer.MIN_VALUE) {
            getNumberTaxaCovered();
            return totalNReads; 
        } else {
            return totalNReads;
        }
    }

    private void assignRefTag() {
        int lengthOfRef = Integer.MIN_VALUE;
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            if (sTBT.divergence == 0 && sTBT.tagLength > lengthOfRef) {
                indexOfRef = tagIndex;
                lengthOfRef = sTBT.tagLength;
            }
            ++tagIndex;
        }
    }
    
    public byte[][] getCommonAlleles() {
        return myCommonAlleles;
    }
    
    public byte[][][] getAlleleDepthsInTaxa() {
        return myAlleleDepthsInTaxa;
    }

    // Qi's SNP caller for VCF output
    public byte[][] getSNPCallsVCF(boolean callBiallelicSNPsWithGap, boolean includeReferenceTag) {
        if (theTags.size() < 2) {
            status = "invariant";
            return null;
        }
        GenotypeTable tagAlignment = getVariableSites();
        if (tagAlignment == null || tagAlignment.numberOfSites() < 1) {
            status = "invariant";
            return null;
        }
        int nSites = tagAlignment.numberOfSites();
        int nTaxa = theTags.get(0).tagDist.length;
        if (nTaxa < 1) {
            status = "noTaxa";
            return null;
        }
        status = "polymorphic";
        byte[][] callsBySite = new byte[nSites][nTaxa];
        if (includeReferenceTag) refCallsBySite = new byte[nSites];
        populateAllelesAtVariableSitesByTag(tagAlignment, nSites, includeReferenceTag, callBiallelicSNPsWithGap);
        positionsOfVariableSites = new int[nSites];
        myCommonAlleles = new byte[maxNumAlleles][nSites];
        myAlleleDepthsInTaxa = new byte[maxNumAlleles][nSites][nTaxa];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.chromosomalPosition(s);
            byte[] commonAlleles = getCommonAlleles(s, nTaxa, includeReferenceTag);
            int[][] alleleDepthsInTaxa = getAlleleDepthsInTaxa(commonAlleles, s, nTaxa, includeReferenceTag);
            setAlleleDepthsInTaxaForSite(s, alleleDepthsInTaxa, commonAlleles);
            for (int tx = 0; tx < nTaxa; tx++) {
                callsBySite[s][tx] = VCFUtil.resolveVCFGeno(commonAlleles, alleleDepthsInTaxa, tx);
            }
        }
        return callsBySite;
    }

    public byte[][] getSNPCallsQuant(boolean callBiallelicSNPsWithGap, boolean includeReferenceTag) {
        if (theTags.size() < 2) {
            status = "invariant";
            return null;
        }
        GenotypeTable tagAlignment = this.getVariableSites();
        if (tagAlignment == null || tagAlignment.numberOfSites() < 1) {
            status = "invariant";
            return null;
        }
        int nSites = tagAlignment.numberOfSites();
        int nTaxa = theTags.get(0).tagDist.length;
        if (nTaxa < 1) {
            status = "noTaxa";  // this shouldn't happen but is here just as a check
            return null;
        }
        status = "polymorphic";
        byte[][] callsBySite = new byte[nSites][nTaxa];
        if (includeReferenceTag) refCallsBySite = new byte[nSites];
        populateAllelesAtVariableSitesByTag(tagAlignment, nSites, includeReferenceTag, callBiallelicSNPsWithGap);
        positionsOfVariableSites = new int[nSites];
        myCommonAlleles = new byte[maxNumAlleles][nSites];
        myAlleleDepthsInTaxa = new byte[maxNumAlleles][nSites][nTaxa];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.chromosomalPosition(s);
            byte[] commonAlleles = getCommonAlleles(s, nTaxa, includeReferenceTag); // NOTE: gap could be one of the common alleles (even if callBiallelicSNPsWithGap is false)
            int[][] alleleDepthsInTaxa = getAlleleDepthsInTaxa(commonAlleles, s, nTaxa, includeReferenceTag);
            setAlleleDepthsInTaxaForSite(s, alleleDepthsInTaxa, commonAlleles);
            for (int tx = 0; tx < nTaxa; tx++) {
                int count = 0;
                for (int a = 0; a < maxNumAlleles; a++) {
                    count += alleleDepthsInTaxa[a][tx];
                }
                if (count == 0) {
                    callsBySite[s][tx] = GenotypeTable.UNKNOWN_GENOTYPE;
                    continue;
                }
                // check for each possible homozygote
                boolean done = false;
                for (int a = 0; a < maxNumAlleles; a++) {
                    if ((count - alleleDepthsInTaxa[a][tx]) == 0) {
                        callsBySite[s][tx] = (byte) ((commonAlleles[a] << 4) | commonAlleles[a]);
                        done = true;
                        break;
                    }
                }
                if (done) {
                    continue;
                }
                callsBySite[s][tx] = resolveHetGeno(commonAlleles, alleleDepthsInTaxa, tx);
            }
        }
        return callsBySite;
    }

    public byte[][] getSNPCallsQuant(String refSeq, boolean callBiallelicSNPsWithGap) {
        // ToDo: UPDATE THIS FOR Tassel4 allele encoding
        if (theTags.size() < 2) {
            return null;
        }
        GenotypeTable tagAlignment = this.getVariableSites(refSeq);
        if (tagAlignment == null || tagAlignment.numberOfSites() < 1) {
            return null;
        }
        int nSites = tagAlignment.numberOfSites();
        int nTaxa = theTags.get(0).tagDist.length;  // the number of taxa is the same for all tags
        if (nTaxa < 1) {
            return null;
        }
        byte[][] callsBySite = new byte[nSites][nTaxa];
        final int nAlignedTags = tagAlignment.numberOfTaxa();
        tagIndices = new int[nAlignedTags];  // the reference sequence is not included
        allelesAtVariableSitesByTag = new byte[nSites][theTags.size()];
        for (int tg = 0; tg < nAlignedTags; tg++) {
            int indexInTheTags = Integer.parseInt(tagAlignment.taxaName(tg)); // taxaName in tagAlignment is set to indexInTheTags
            tagIndices[tg] = indexInTheTags;
            for (int s = 0; s < nSites; s++) {
                allelesAtVariableSitesByTag[s][tagIndices[tg]] = tagAlignment.genotype(tg, s);
            }
        }
        positionsOfVariableSites = new int[nSites];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.chromosomalPosition(s);
            for (int tx = 0; tx < nTaxa; tx++) {
                int[] alleleCounts = new int[Byte.MAX_VALUE];
                for (int tg = 0; tg < nAlignedTags; tg++) {
                    int tagIndex = tagIndices[tg];
                    byte baseToAdd = allelesAtVariableSitesByTag[s][tagIndex];
                    if (baseToAdd == GenotypeTable.UNKNOWN_GENOTYPE && callBiallelicSNPsWithGap && maxStartPosition == minStartPosition) {
                        baseToAdd = NucleotideAlignmentConstants.GAP_HOMOZYGOUS;
                    }
                    alleleCounts[baseToAdd] += theTags.get(tagIndex).tagDist[tx];
                }
                callsBySite[s][tx] = resolveQuantGeno(alleleCounts);
            }
        }
        return callsBySite;
    }
    
    private void setAlleleDepthsInTaxaForSite(int site, int[][] alleleDepthsInTaxa, byte[] commonAlleles) {
        for (int a = 0; a < commonAlleles.length; a++) {
            myCommonAlleles[a][site] = commonAlleles[a];
        }
        for (int a = 0; a < alleleDepthsInTaxa.length; a++) {
            for (int tx = 0; tx < alleleDepthsInTaxa[a].length; tx++) {
                if (alleleDepthsInTaxa[a][tx] > 127) { // max value is 127
                    alleleDepthsInTaxa[a][tx] = 127;
                }
                myAlleleDepthsInTaxa[a][site][tx] = (byte)alleleDepthsInTaxa[a][tx];
            }
        }
    }

    private void populateAllelesAtVariableSitesByTag(GenotypeTable tagAlignment, int nSites, boolean includeReferenceTag, boolean callBiallelicSNPsWithGap) {
        int nAlignedTags = tagAlignment.numberOfTaxa();
        tagIndices = new int[nAlignedTags];
        allelesAtVariableSitesByTag = new byte[nSites][theTags.size()];
        for (int tg = 0; tg < nAlignedTags; tg++) {
            tagIndices[tg] = Integer.parseInt(tagAlignment.taxaName(tg).split("_")[0]);  // taxaName in tagAlignment is set to indexInTheTags_"refTag"|"no"
            for (int s = 0; s < nSites; s++) {
                if (includeReferenceTag && tagIndices[tg] == theTags.size()-1) {
                    refCallsBySite[s] = tagAlignment.genotype(tg, s); // diploid byte for the reference allele/geno
                } else {
                    byte allele = tagAlignment.genotypeArray(tg, s)[0]; // tags only have one base so the 1st allele (index [0]) sufffices
                    if (callBiallelicSNPsWithGap && allele == GenotypeTable.UNKNOWN_ALLELE) {
                        allele = NucleotideAlignmentConstants.GAP_ALLELE;
                    }
                    allelesAtVariableSitesByTag[s][tagIndices[tg]] = allele;
                }
            }
        }
    }

    public int[] getPositionsOfVariableSites() {
        return positionsOfVariableSites;
    }
    
    public String getLocusReport(int minTaxaWithLocus, boolean[] varSiteKept) {
        int start, end, totalbp, refTag=Integer.MIN_VALUE;
        if (strand == -1) {
            end = minStartPosition;
            start = minStartPosition - maxTagLength + 1;
        } else {
            start = minStartPosition;
            end = minStartPosition + maxTagLength - 1;
        }
        totalbp = end - start + 1;
        int nVarSites=0, nVarSitesKept=0;
        String posVarSites = "", posVarsKept = "";
        if (status.equals("polymorphic")) {
            nVarSites = positionsOfVariableSites.length;
            for (int s = 0; s < nVarSites; s++) {
                posVarSites = s < nVarSites-1 ? posVarSites + positionsOfVariableSites[s] + ":"
                                              : posVarSites + positionsOfVariableSites[s];
                if (varSiteKept[s]) {
                    posVarsKept = posVarsKept + positionsOfVariableSites[s] + ":";
                    nVarSitesKept++;
                }
            }
            if (posVarsKept.length()>0) posVarsKept = posVarsKept.substring(0, posVarsKept.length()-1);
            else posVarsKept = "NA";
        } else {
            posVarSites = "NA";
            posVarsKept = "NA";
            if (status.equals("invariant") || status.contains("tooManyTags")) assignRefTag();
        }
        refTag = (indexOfRef==Integer.MIN_VALUE) ? 0 : 1;
        return
            chromosome +"\t"+
            start +"\t"+
            end +"\t"+
            strand +"\t"+
            totalbp +"\t"+
            theTags.size() +"\t"+
            this.getTotalNReads() +"\t"+
            this.getNumberTaxaCovered() +"\t"+
            minTaxaWithLocus+"\t"+
            status +"\t"+ 
            nVarSites +"\t"+ 
            posVarSites +"\t"+
            nVarSitesKept  +"\t"+
            posVarsKept +"\t"+
            refTag +"\t"+
            maxTagLength +"\t"+
            minTagLength +"\n"
        ;
    }

    private GenotypeTable getVariableSites() {
        if (theTags.size() < 2) {
            status = "invariant";
            return null;
        }
        if (theTags.size() > maxAlignmentSize) {
            status = "tooManyTags(>"+maxAlignmentSize+")";
            return null;
        }
        if (indexOfRef == Integer.MIN_VALUE) this.assignRefTag();
        boolean printOutAlignments = true;
        List<DNASequence> lst = new ArrayList<DNASequence>();
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            try {
                DNASequence ds = new DNASequence(sTBT.tagTrimmed);
                String refMark = (tagIndex == indexOfRef) ? "refTag" : "no";
                ds.setOriginalHeader(tagIndex + "_" + refMark);    // OriginalHeader set to indexInTheTags_'refTag'|'no'
                ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
                lst.add(ds);
               ++tagIndex;
            } catch (CompoundNotFoundException ex) {
                myLogger.error("TagLocus:getVariableSites, compoundNotFound exception from DNASequence call for: " + sTBT.tagTrimmed);
                myLogger.debug(ex.getMessage(), ex);
                return null;
            }
        }
 
        Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        int nSites=profile.getAlignedSequence(1).getSequenceAsString().length();
        String[] alignedSeqs = new String[theTags.size()];
        GenotypeCallTableBuilder gB=GenotypeCallTableBuilder.getInstance(theTags.size(),nSites);
        TaxaListBuilder tlB=new TaxaListBuilder();
        PositionListBuilder pALB=new PositionListBuilder();
        for (int i=0; i<nSites; i++) {pALB.add(new GeneralPosition.Builder(Chromosome.UNKNOWN,i).build());}
        for (int i = 0; i < alignedSeqs.length; i++) {
            alignedSeqs[i] = profile.getAlignedSequence(i + 1).getSequenceAsString();
            String taxonName = profile.getAlignedSequence(i + 1).getOriginalSequence().getOriginalHeader();
            if (taxonName.split("_")[1].equals("refTag")) {  // name was set to indexInTheTags_"refTag"|"no"
                if (alignedSeqs[i].contains("-")) {
                    pALB=new PositionListBuilder();
                    pALB.add(new GeneralPosition.Builder(Chromosome.UNKNOWN,0).build());
                    int prevPosition = 0;
                    for (int site = 1; site < alignedSeqs[i].length(); site++) {
                        int currPosition = (alignedSeqs[i].charAt(site) == '-') ? prevPosition : prevPosition + 1;
                        pALB.add(new GeneralPosition.Builder(Chromosome.UNKNOWN,currPosition).build());
                        prevPosition = currPosition;
                    }
                }
            }
            tlB.add(new Taxon(taxonName));
        }
        profile = null;
        gB.setBases(alignedSeqs);
        GenotypeTable aa = GenotypeTableBuilder.getInstance(gB.build(),pALB.build(),tlB.build());
        GenotypeTable faa = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(aa, 0.000001, 1.0, 2);
        if (printOutAlignments && (minStartPosition % 1000 == 0)) {
            TaxaList tL=tlB.build();
            String tagStr;
            System.out.println("\nHere is an example alignment for a TagLocus (1 out of every 1000 is displayed):");
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  All sites:");
            for (int tg = 0; tg < alignedSeqs.length; tg++) {
                tagStr = aa.genotypeAsStringRow(tg);
                tagStr = tagStr.replaceAll(";", "");
                System.out.println(tagStr + " " + tL.taxaName(tg));
            }
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  Polymorphic sites only:");
            for (int tg = 0; tg < alignedSeqs.length; tg++) {
                tagStr = faa.genotypeAsStringRow(tg);
                tagStr = tagStr.replaceAll(";", "");
                System.out.println(tagStr + " " + tL.taxaName(tg));
            }
            System.out.println();
        }
        if (faa.numberOfSites() > maxSNPsPerLocus) {
            status = "tooManyVariants(>"+maxSNPsPerLocus+")";
            return null;
        }
        if (faa.numberOfSites() < 1) {
            status = "noVarSitesInAlign";
            return null;
        }
        if (faa.numberOfTaxa() < 2) {
            status = "onlyOneTagInAlign";
            return null;
        }
        return faa;
    }

    private GenotypeTable getVariableSites(String refSeqInRegion) {
        if (theTags.size() < 2) {
            return null;
        }
        boolean printOutAlignments = true;
        int startRefGenIndex = 0, endRefGenIndex = 1, startTagIndex = 2, endTagIndex = 3; // relevant indices in alignStats[]  (alignedTagLen=4)
        DNASequence dsRefSeq = null;
        try {
            dsRefSeq = new DNASequence(refSeqInRegion);
        } catch (CompoundNotFoundException ex) {
            // Something's wrong in the code  - this shouldn't happen
            myLogger.error("TagLocus:getVariableSites 2, compoundNotFound exception from DNASequence call for: " + refSeqInRegion);
            myLogger.debug(ex.getMessage(), ex);
            return null;
        }
        dsRefSeq.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
        int minRefGenIndex = Integer.MAX_VALUE, maxRefGenIndex = Integer.MIN_VALUE;
        ArrayList<SequencePair<DNASequence, NucleotideCompound>> pairwiseAligns = new ArrayList<SequencePair<DNASequence, NucleotideCompound>>();
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            DNASequence ds = null;
            try {
                ds = new DNASequence(sTBT.tagTrimmed);
            } catch (CompoundNotFoundException ex) {
                // Something's wrong - this shouldn't happen
                myLogger.error("TagLocus:getVariableSites 3, compoundNotFound exception from DNASequence call for: " + sTBT.tagTrimmed);
                myLogger.debug(ex.getMessage(), ex);
                return null;
            }
            ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
            SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(ds, dsRefSeq, PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
            int[] alignStats = getAlignStats(psa, printOutAlignments, tagIndex, sTBT.tagLength, sTBT.tagStrand);
            minRefGenIndex = adjustMinRefGenIndex(minRefGenIndex, alignStats[startRefGenIndex], alignStats[startTagIndex]);
            maxRefGenIndex = adjustMaxRefGenIndex(maxRefGenIndex, alignStats[endRefGenIndex], alignStats[endTagIndex], sTBT.tagLength, refSeqInRegion);
            pairwiseAligns.add(psa);
            ++tagIndex;
        }
        minStartPosition += minRefGenIndex - 1;
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("minRefGenIndex:" + minRefGenIndex + "  maxRefGenIndex:" + maxRefGenIndex + "  ChrPositionAtMinRefGenIndex:" + minStartPosition + "\n");
        }
        //Todo where are aseqs and names filled out?
        String[] aseqs = new String[theTags.size()];  // omit the reference genome sequence
        String[] names = new String[theTags.size()];
        char[][] myAlign = getAlignment(pairwiseAligns, refSeqInRegion, minRefGenIndex, maxRefGenIndex, aseqs, names);
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            writeAlignment(refSeqInRegion, myAlign, minRefGenIndex, maxRefGenIndex);
        }
        GenotypeTable a = null;
        TaxaList tL=new TaxaListBuilder().addAll(names).build();
        int nSites=aseqs[0].length();
        PositionListBuilder pALB=new PositionListBuilder();
        for (int i=0; i<nSites; i++) {pALB.add(new GeneralPosition.Builder(Chromosome.UNKNOWN,i).build());}
        GenotypeCallTableBuilder gB=GenotypeCallTableBuilder.getInstance(theTags.size(),nSites);
        for (int i=0; i<aseqs.length; i++) {gB.setBaseRangeForTaxon(i,0,aseqs[i].getBytes());}
        a=GenotypeTableBuilder.getInstance(gB.build(),pALB.build(),tL);
//        a = BitAlignment.getNucleotideInstance(tL, aseqs, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Chromosome[]{Chromosome.UNKNOWN}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), true);
        GenotypeTable fa = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(a, 0.000001, 1.0, 2);
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  FA (alignment filtered for polymorphic sites):\n" + fa.toString());
        }
        if (fa.numberOfSites() > maxSNPsPerLocus * 5 || fa.numberOfSites() < 1 || fa.numberOfTaxa() < 2) {
            return null;
        }
        return fa;
    }

    private int[] getAlignStats(SequencePair<DNASequence, NucleotideCompound> psa, boolean printOutAlignments, int tagIndex, int tagLength, byte tagStrand) {
        int[] alignStats = new int[5];
        int startRefGenIndex = 0, endRefGenIndex = 1, startTagIndex = 2, endTagIndex = 3, alignedTagLen = 4; // indices in alignStats[]
        alignStats[startRefGenIndex] = psa.getIndexInTargetAt(1);
        alignStats[endRefGenIndex] = psa.getIndexInTargetAt(psa.getLength());
        alignStats[startTagIndex] = psa.getIndexInQueryAt(1);
        alignStats[endTagIndex] = psa.getIndexInQueryAt(psa.getLength());
        alignStats[alignedTagLen] = alignStats[endTagIndex] - alignStats[startTagIndex] + 1;
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("tagIndex:" + tagIndex
                    + "  startRefGenIndex:" + alignStats[startRefGenIndex]
                    + "  endRefGenIndex:" + alignStats[endRefGenIndex]
                    + "  tagLength:" + tagLength
                    + "  startTagIndex:" + alignStats[startTagIndex]
                    + "  endTagIndex:" + alignStats[endTagIndex]
                    + "  alignedTagLen:" + alignStats[alignedTagLen]
                    + "  originalStrand: " + tagStrand + "\n"
                    + psa);
        }
        return alignStats;
    }

    private int adjustMinRefGenIndex(int minRefGenIndex, int startRefGenIndex, int startTagIndex) {
        if (startRefGenIndex < minRefGenIndex) {
            minRefGenIndex = startRefGenIndex;
        }
        if (startTagIndex > 1 && startTagIndex < 4 && startRefGenIndex - startTagIndex + 1 < minRefGenIndex && startRefGenIndex - startTagIndex + 1 > 0) {
            // extend regional alignment if there was soft clipping of 1 or 2 bases and there is sufficient 5' refSeq available
            minRefGenIndex = startRefGenIndex - startTagIndex + 1;
        }
        return minRefGenIndex;
    }

    private int adjustMaxRefGenIndex(int maxRefGenIndex, int endRefGenIndex, int endTagIndex, int tagLength, String refSeq) {
        if (endRefGenIndex > maxRefGenIndex) {
            maxRefGenIndex = endRefGenIndex;
        }
        if (endTagIndex < tagLength
                && tagLength - endTagIndex < 3
                && endRefGenIndex + tagLength - endTagIndex > maxRefGenIndex
                && endRefGenIndex + tagLength - endTagIndex <= refSeq.length()) {
            // extend regional alignment if there was soft clipping of 1 or 2 bases and there is sufficient 3' refSeq available
            maxRefGenIndex = endRefGenIndex + tagLength - endTagIndex;
        }
        return maxRefGenIndex;
    }

    private char[][] getAlignment(ArrayList<SequencePair<DNASequence, NucleotideCompound>> pairwiseAligns,
            String refSeq, int minRefGenIndex, int maxRefGenIndex, String[] aseqs, String[] names) {
        int totAlignedLen = maxRefGenIndex - minRefGenIndex + 1;
        char[][] myAlign = new char[theTags.size()][totAlignedLen];  // omit the reference genome sequence
        for (int t = 0; t < myAlign.length; t++) {
            Arrays.fill(myAlign[t], 'N');
        }
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            SequencePair<DNASequence, NucleotideCompound> psa = pairwiseAligns.get(tagIndex);
            int tagStart = psa.getIndexInQueryAt(1);
            int tagEnd = psa.getIndexInQueryAt(psa.getLength());
            int refSeqStart = psa.getIndexInTargetAt(1);
            int refSeqEnd = psa.getIndexInTargetAt(psa.getLength());
            if (tagStart > 1 && tagStart < 4 && refSeqStart - minRefGenIndex - tagStart + 1 > -1) {
                // extend tag start if there was soft clipping of 1 or 2 bases and there is sufficient 5' refSeq available
                for (int offset = tagStart - 1; offset > 0; offset--) {
                    myAlign[tagIndex][refSeqStart - minRefGenIndex - offset] = sTBT.tagTrimmed.charAt(tagStart - offset - 1);
                }
            }
            for (int i = 1; i <= psa.getLength(); i++) {
                char refBase = psa.getCompoundInTargetAt(i).getBase().charAt(0);
                if (refBase != '-') {
                    myAlign[tagIndex][psa.getIndexInTargetAt(i) - minRefGenIndex] = psa.getCompoundInQueryAt(i).getBase().charAt(0);
                }
            }
            int extension = sTBT.tagLength - tagEnd;
            if (extension > 0 && extension < 3 && refSeqEnd - minRefGenIndex + extension < refSeq.length()) {
                // extend tag end if there was soft clipping of 1 or 2 bases and there is sufficient 3' refSeq available
                for (int offset = sTBT.tagLength - tagEnd; offset > 0; offset--) {
                    myAlign[tagIndex][refSeqEnd - minRefGenIndex + offset] = sTBT.tagTrimmed.charAt(sTBT.tagLength - offset);
                }
            }
            aseqs[tagIndex] = new String(myAlign[tagIndex]);
            names[tagIndex] = tagIndex + "";
            ++tagIndex;
        }
        return myAlign;
    }

    private void writeAlignment(String refSeq, char[][] myAlign, int minRefGenIndex, int maxRefGenIndex) {
        System.out.println("All tags in the region aligned to the reference sequence (first line) (insertions relative to the reference excluded):");
        System.out.println(refSeq.substring(minRefGenIndex - 1, maxRefGenIndex));
        for (int tagIndex = 0; tagIndex < myAlign.length; tagIndex++) {
            for (int b = 0; b < myAlign[tagIndex].length; b++) {
                System.out.print(myAlign[tagIndex][b]);
            }
            System.out.print("\n");
        }
        System.out.print("\n");
    }

    private byte[] getCommonAlleles(int s, int nTaxa, boolean includeReferenceTag) {
        int[] alleleCounts = new int[16];
        int nTags = includeReferenceTag ? theTags.size()-1 : theTags.size();
        for (int tg = 0; tg < nTags; tg++) {
            byte baseToAdd = allelesAtVariableSitesByTag[s][tg];
            for (int tx = 0; tx < nTaxa; tx++) {
                alleleCounts[baseToAdd] += theTags.get(tg).tagDist[tx];
            }
        }
        byte[] commonAlleles = new byte[maxNumAlleles];
        int[][] sortedAlleleCounts = sortAllelesByCount(alleleCounts);
        // Ties between maxNumAlleles - 1 and maxNumAlleles are not handled
        for (int i = 0; i < maxNumAlleles; i++) {
            commonAlleles[i] = (byte) sortedAlleleCounts[0][i];
        }
        return commonAlleles;
    }

    private int[][] getAlleleDepthsInTaxa(byte[] commonAlleles, int s, int nTaxa, boolean includeReferenceTag) {
        int[][] allelesInTaxa = new int[maxNumAlleles][nTaxa];
        int nTags = includeReferenceTag ? theTags.size()-1 : theTags.size(); // skip the reference tag (=last tag, if present), as it has no tagDist[]
        for (int tg = 0; tg < nTags; tg++) {
            byte baseToAdd = allelesAtVariableSitesByTag[s][tg];
            for (int a = 0; a < maxNumAlleles; a++) {
                if (baseToAdd == commonAlleles[a]) {
                    for (int tx = 0; tx < nTaxa; tx++) {
                        allelesInTaxa[a][tx] += theTags.get(tg).tagDist[tx];
                    }
                }
            }
        }
        return allelesInTaxa;
    }
    
    
    private byte resolveHetGeno(byte[] alleles, int[][] allelesInTaxa, int tx) {
        int max = 0;
        byte maxAllele = GenotypeTable.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = GenotypeTable.UNKNOWN_ALLELE;
        for (int a = 0; a < maxNumAlleles; a++) {
            if (allelesInTaxa[a][tx] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = allelesInTaxa[a][tx];
                maxAllele = alleles[a];
            } else if (allelesInTaxa[a][tx] > nextMax) {
                nextMax = allelesInTaxa[a][tx];
                nextMaxAllele = alleles[a];
            }
        }
        int totCount = max + nextMax;
        if (totCount < maxCountAtGeno) {
            if (nextMax < likelihoodRatioThreshAlleleCnt[totCount]) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        } else {
            if (nextMax / totCount < 0.1) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        }
    }

    private byte resolveQuantGeno(int[] alleleCounts) {
        int[][] sortedAlleleCounts = sortAllelesByCount(alleleCounts);
        int a1Count = sortedAlleleCounts[1][0];
        if (a1Count == 0) {
            return GenotypeTable.UNKNOWN_GENOTYPE;
        }
        int a2Count = sortedAlleleCounts[1][1];  // What if a3Count = a2Count? -- this situation is not dealt with
        byte a1 = (byte) sortedAlleleCounts[0][0];
        if (a2Count == 0) {
            return a1;
        }
        byte a2 = (byte) sortedAlleleCounts[0][1];
        int totCount = a1Count + a2Count;
        if (totCount < maxCountAtGeno) {
            if (a2Count < likelihoodRatioThreshAlleleCnt[totCount]) {
                return a1;  // call it a homozygote
            } else {
                return GenotypeTableUtils.getDiploidValue(a1, a2);  // call it a het
                //return IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(a1, a2); // call it a het
            }
        } else {
            if (a2Count / totCount < 0.1) {
                return a1;  // call it a homozygote
            } else {
                return GenotypeTableUtils.getDiploidValue(a1, a2);  // call it a het
                //return IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(a1, a2); // call it a het
            }
        }
    }

    private int[][] sortAllelesByCount(int[] alleleCounts) {
        // note that 'N' is not included as an allele
        byte[] alleles = {NucleotideAlignmentConstants.A_ALLELE, NucleotideAlignmentConstants.C_ALLELE,
            NucleotideAlignmentConstants.G_ALLELE, NucleotideAlignmentConstants.T_ALLELE, NucleotideAlignmentConstants.GAP_ALLELE};
        int[][] result = new int[2][alleles.length]; // result[0][a]=allele; result[1][a]=count
        for (int i = 0; i < alleles.length; i++) {
            result[0][i] = alleles[i];
            result[1][i] = alleleCounts[alleles[i]];
        }
        boolean change = true;
        while (change) { // sort the alleles by descending frequency
            change = false;
            for (int k = 0; k < alleles.length - 1; k++) {
                if (result[1][k] < result[1][k + 1]) {
                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;
                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;
                    change = true;
                }
            }
        }
        return result;
    }

    private String padTagWithNs(SingleTagByTaxa tag, String refSeq) {
        StringBuilder sb = new StringBuilder();
        char[] nullBases;
        if (tag.tagStrand == -1) {
            nullBases = new char[tag.startPosition - minStartPosition - tag.tagTrimmed.length() + 1];
        } else {
            nullBases = new char[tag.startPosition - minStartPosition];
        }
        Arrays.fill(nullBases, 'N');
        sb.append(nullBases);
        sb.append(tag.tagTrimmed);
        if (tag.tagStrand == -1) {
            nullBases = new char[refSeq.length() - sb.length()];
        } else {
            nullBases = new char[refSeq.length() - sb.length()];
        }

        Arrays.fill(nullBases, 'N');
        sb.append(nullBases);
        return sb.toString();
    }
}
class SingleTagByTaxa {

    int tagTOPMIndex;
    int tagLength;
    int startPosition;
    byte tagStrand;
    int divergence;
    String tagTrimmed;
    int tagTBTIndex; //index in the TBT
    int taxaWithTag;
    byte[] tagDist;  // observed count of the tag for each taxon

    SingleTagByTaxa(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, boolean includeRefGenome, boolean fuzzyStartPositions) {
        tagStrand = Byte.MIN_VALUE;
        this.tagTOPMIndex = tagTOPMIndex;
        long[] tag = theTOPM.getTag(tagTOPMIndex);
        tagTBTIndex = theTBT.getTagIndex(tag);
        taxaWithTag = (tagTBTIndex > -1) ? theTBT.getNumberOfTaxaWithTag(tagTBTIndex) : 0;
        if (taxaWithTag > 0) {  // tags with 0 taxaWithTag will not be added to TagLocus
            startPosition = theTOPM.getStartPosition(tagTOPMIndex);
            tagLength = theTOPM.getTagLength(tagTOPMIndex);
            divergence = theTOPM.getDivergence(tagTOPMIndex);
            tagTrimmed = BaseEncoder.getSequenceFromLong(tag).substring(0, tagLength);
            tagStrand = theTOPM.getStrand(tagTOPMIndex);
            if (includeRefGenome && fuzzyStartPositions) {
                if (tagStrand == -1) {
                    tagTrimmed = BaseEncoder.getReverseComplement(tagTrimmed);
                }
            } else if (tagLength < theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize) {
                tagTrimmed = tagTrimmed
                        + theTOPM.getNullTag().substring(0, theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize - tagLength).replace("A", "N");
            }
            tagDist = theTBT.getTaxaReadCountsForTag(tagTBTIndex);
        }
    }
    
    // this constructor can be used to add a refTag directly from the refererence genome
    SingleTagByTaxa(int startPosition, byte strand, String refTag, int nLongsPerTag, String nullTag) {
        tagTOPMIndex = Integer.MIN_VALUE;
        tagLength = (byte) refTag.length();
        this.startPosition = startPosition;
        tagStrand = strand;
        divergence = 0;
        tagTrimmed = refTag;
        if (tagLength < nLongsPerTag*BaseEncoder.chunkSize) {
            tagTrimmed = tagTrimmed + 
                    nullTag.substring(0,nLongsPerTag*BaseEncoder.chunkSize-tagLength).replace("A","N");
        }
        tagTBTIndex = Integer.MIN_VALUE;
        taxaWithTag = 0;
        tagDist = null;
    }
}
