/*
 *  TOPMGenotypeTable
 */
package net.maizegenetics.dna.snp;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;
import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.map.TOPMInterface;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.GeneralAnnotationStorage;

/**
 *
 * @author Terry Casstevens
 */
public class TOPMGenotypeTable implements GenotypeTable {

    private final TOPMInterface myTOPM;
    private int[] myIndicesOfSortByPosition;
    private int[] mySiteOffsetForEachTag;
    private final int myNumTags;
    private final PositionList myPositionList;
    private final TaxaList myTaxaList;
    private int myBlankCount = 0;

    public TOPMGenotypeTable(TOPMInterface topm) {
        myTOPM = topm;
        myNumTags = myTOPM.getTagCount();
        sortByTag();
        myPositionList = getPositionList();
        myTaxaList = getTaxaList();
    }

    private int getStartPosition(int index) {
        if (myTOPM.getStrand(index) == -1) {
            int tagLength = myTOPM.getTagLength(index);
            if (tagLength < 64) {
                return myTOPM.getStartPosition(index) - (tagLength - 1);
            } else {
                return myTOPM.getEndPosition(index);
            }
        } else {
            return myTOPM.getStartPosition(index);
        }
    }

    private int getEndPosition(int index) {
        if (myTOPM.getStrand(index) == -1) {
            return myTOPM.getStartPosition(index);
        } else {
            return myTOPM.getEndPosition(index);
        }
    }

    private void sortByTag() {

        final int[] indicesOfSortByPosition = new int[myNumTags];
        for (int i = 0; i < indicesOfSortByPosition.length; i++) {
            indicesOfSortByPosition[i] = i;
        }

        Swapper swapTag = new Swapper() {
            @Override
            public void swap(int a, int b) {
                int temp = indicesOfSortByPosition[a];
                indicesOfSortByPosition[a] = indicesOfSortByPosition[b];
                indicesOfSortByPosition[b] = temp;
            }
        };

        IntComparator compTag = new IntComparator() {
            @Override
            public int compare(int a, int b) {
                int chra = myTOPM.getChromosome(indicesOfSortByPosition[a]);
                int chrb = myTOPM.getChromosome(indicesOfSortByPosition[b]);
                if ((chra == Integer.MIN_VALUE) && (chrb == Integer.MIN_VALUE)) {
                    return 0;
                } else if (chra < chrb) {
                    return -1;
                } else if (chra > chrb) {
                    return 1;
                }

                int startPosA = getStartPosition(indicesOfSortByPosition[a]);

                int startPosB = getStartPosition(indicesOfSortByPosition[b]);

                if (startPosA < startPosB) {
                    return -1;
                } else if (startPosB < startPosA) {
                    return 1;
                } else {
                    return 0;
                }
            }
        };

        GenericSorting.quickSort(0, indicesOfSortByPosition.length, compTag, swapTag);

        int i;
        for (i = 0; i < myNumTags; i++) {
            if (myTOPM.getChromosome(indicesOfSortByPosition[i]) != Integer.MIN_VALUE) {
                break;
            }
        }
        int numTagsMapped = myNumTags - i;
        List<Integer> tempIndicesOfSortByPosition = new ArrayList<>(numTagsMapped * 3 / 2);
        List<Integer> siteOffsetForEachTag = new ArrayList<>(numTagsMapped * 3 / 2);
        int previousEndPosition = -1;
        int previousChr = -1;
        int previousIndex = -1;
        for (int j = 0; j < numTagsMapped; j++) {
            int currentIndex = indicesOfSortByPosition[j + i];
            tempIndicesOfSortByPosition.add(currentIndex);

            if (myTOPM.getChromosome(currentIndex) != previousChr) {
                previousEndPosition = -1;
                previousChr = myTOPM.getChromosome(currentIndex);
            }

            int startPosition = getStartPosition(currentIndex);

            if (previousEndPosition >= startPosition) {
                int previousStartPosition = getStartPosition(previousIndex);
                siteOffsetForEachTag.add(startPosition - previousStartPosition + siteOffsetForEachTag.get(siteOffsetForEachTag.size() - 1));
            } else {
                tempIndicesOfSortByPosition.add(tempIndicesOfSortByPosition.size() - 1, -1);
                siteOffsetForEachTag.add(-1);
                siteOffsetForEachTag.add(0);
            }

            previousEndPosition = getEndPosition(currentIndex);

            previousIndex = currentIndex;

        }

        myIndicesOfSortByPosition = new int[tempIndicesOfSortByPosition.size()];
        for (int k = 0; k < tempIndicesOfSortByPosition.size(); k++) {
            myIndicesOfSortByPosition[k] = tempIndicesOfSortByPosition.get(k);
        }

        mySiteOffsetForEachTag = new int[siteOffsetForEachTag.size()];
        for (int k = 0; k < siteOffsetForEachTag.size(); k++) {
            mySiteOffsetForEachTag[k] = siteOffsetForEachTag.get(k);
        }

    }

    private PositionList getPositionList() {

        Map<String, Position> positions = new LinkedHashMap<>();
        Chromosome chrObj = new Chromosome(String.valueOf(0));

        for (int t = 0; t < myIndicesOfSortByPosition.length; t++) {

            if (myIndicesOfSortByPosition[t] != -1) {
                for (int p = mySiteOffsetForEachTag[t]; p < mySiteOffsetForEachTag[t] + 64; p++) {
                    String key = String.valueOf(p);
                    if (!positions.containsKey(key)) {
                        Position currentPos = (new GeneralPosition.Builder(chrObj, p)).build();
                        positions.put(key, currentPos);
                    }
                }

            }

        }

        PositionListBuilder builder = new PositionListBuilder();
        builder.addAll(positions.values());
        return builder.build();

    }

    private TaxaList getTaxaList() {
        TaxaListBuilder builder = new TaxaListBuilder();
        for (int i = 0; i < myIndicesOfSortByPosition.length; i++) {
            if (myIndicesOfSortByPosition[i] != -1) {
                int startPosition = getStartPosition(myIndicesOfSortByPosition[i]);
                String taxaName = String.valueOf(myTOPM.getStrand(myIndicesOfSortByPosition[i])) + String.valueOf(myIndicesOfSortByPosition[i])
                        + "_" + String.valueOf(myTOPM.getChromosome(myIndicesOfSortByPosition[i]))
                        + "_" + String.valueOf(startPosition);
                builder.add(new Taxon.Builder(taxaName).build());
            } else {
                builder.add(new Taxon.Builder("Blank " + String.valueOf(myBlankCount++)).build());
            }
        }
        return builder.build();
    }

    @Override
    public GenotypeCallTable genotypeMatrix() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public boolean isSNP(int taxon, int site) {
        if (myIndicesOfSortByPosition[taxon] == -1) {
            return false;
        }
        byte[] offsets = myTOPM.getVariantPosOffArray(myIndicesOfSortByPosition[taxon]);
        if ((offsets == null) || (offsets.length == 0)) {
            return false;
        }
        int tagLength = myTOPM.getTagLength(myIndicesOfSortByPosition[taxon]);
        if (myTOPM.getStrand(myIndicesOfSortByPosition[taxon]) == -1) {
            // first
            //int index = (site - mySiteOffsetForEachTag[taxon]) - 63;
            // second
            //int startPos = myTOPM.getStartPosition(myIndicesOfSortByPosition[taxon]);
            //int endPos = myTOPM.getEndPosition(myIndicesOfSortByPosition[taxon]);
            //int index = -startPos + (endPos + (site - mySiteOffsetForEachTag[taxon]));
            int index = -myTOPM.getStartPosition(myIndicesOfSortByPosition[taxon]) + getStartPosition(myIndicesOfSortByPosition[taxon]) + site - mySiteOffsetForEachTag[taxon];
            if ((index <= 0) && (index > -tagLength)) {
                for (byte offset : offsets) {
                    if (offset == index) {
                        return true;
                    }
                }
                return false;
            } else {
                return false;
            }
        } else {
            int index = site - mySiteOffsetForEachTag[taxon];
            if ((index >= 0) && (index < tagLength)) {
                for (byte offset : offsets) {
                    if (offset == index) {
                        return true;
                    }
                }
                return false;
            } else {
                return false;
            }
        }
    }

    @Override
    public byte genotype(int taxon, int site) {

        if (myIndicesOfSortByPosition[taxon] == -1) {
            return NucleotideAlignmentConstants.UNDEFINED_DIPLOID_ALLELE;
        }

        int tagLength = myTOPM.getTagLength(myIndicesOfSortByPosition[taxon]);
        if (myTOPM.getStrand(myIndicesOfSortByPosition[taxon]) == -1) {
            int index = tagLength - site + mySiteOffsetForEachTag[taxon] - 1;
            if ((index >= 0) && (index < tagLength)) {
                String tag = BaseEncoder.getSequenceFromLong(myTOPM.getTag(myIndicesOfSortByPosition[taxon]));
                return NucleotideAlignmentConstants.getNucleotideDiploidComplement(NucleotideAlignmentConstants.getNucleotideDiploidByte(tag.charAt(index)));
            } else {
                return NucleotideAlignmentConstants.UNDEFINED_DIPLOID_ALLELE;
            }
        } else {
            int index = site - mySiteOffsetForEachTag[taxon];
            if ((index >= 0) && (index < tagLength)) {
                String tag = BaseEncoder.getSequenceFromLong(myTOPM.getTag(myIndicesOfSortByPosition[taxon]));
                return NucleotideAlignmentConstants.getNucleotideDiploidByte(tag.charAt(index));
            } else {
                return NucleotideAlignmentConstants.UNDEFINED_DIPLOID_ALLELE;
            }
        }

    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return GenotypeTableUtils.getDiploidValues(genotype(taxon, site));
    }

    @Override
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] genotypeAllTaxa(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, WHICH_ALLELE allele, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, WHICH_ALLELE allele, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String genotypeAsString(int taxon, int site) {

        byte genotype = genotype(taxon, site);
        if (genotype == NucleotideAlignmentConstants.UNDEFINED_DIPLOID_ALLELE) {
            return "";
        } else {
            return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype);
        }
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean hasReference() {
        return false;
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int heterozygousCount(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String siteName(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int numberOfSites() {
        return myPositionList.numberOfSites();
    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] firstLastSiteOfChromosome(Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int numberOfTaxa() {
        return myIndicesOfSortByPosition.length;
    }

    @Override
    public PositionList positions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int chromosomalPosition(int site) {
        return myPositionList.chromosomalPosition(site);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
    }

    @Override
    public int[] physicalPositions() {
        return myPositionList.physicalPositions();
    }

    @Override
    public String chromosomeName(int site) {
        return myPositionList.chromosomeName(site);
    }

    @Override
    public Chromosome chromosome(int site) {
        return myPositionList.chromosome(site);
    }

    @Override
    public Chromosome chromosome(String name) {
        return myPositionList.chromosome(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myPositionList.chromosomes();
    }

    @Override
    public int numChromosomes() {
        return myPositionList.numChromosomes();
    }

    @Override
    public int[] chromosomesOffsets() {
        return myPositionList.chromosomesOffsets();
    }

    @Override
    public boolean hasGenotype() {
        return true;
    }

    @Override
    public boolean hasDepth() {
        return false;
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        return null;
    }

    @Override
    public boolean hasAlleleProbabilities() {
        return false;
    }

    @Override
    public boolean hasReferenceProbablity() {
        return false;
    }

    @Override
    public boolean hasDosage() {
        return false;
    }

    @Override
    public int indelSize(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isIndel(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isAllPolymorphic() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPolymorphic(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte majorAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String majorAlleleAsString(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte minorAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String minorAlleleAsString(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] minorAlleles(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] alleles(int site) {
        int[][] alleles = allelesSortedByFrequency(site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    @Override
    public double minorAlleleFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double majorAlleleFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public String taxaName(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String genomeVersion() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPositiveStrand(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public GenotypeTable[] compositeAlignments() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return GenotypeTableUtils.getAllelesSortedByFrequency(this, site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPhased() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean retainsRareAlleles() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[][] alleleDefinitions() {
        return NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0];
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String diploidAsString(int site, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int maxNumAlleles() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalNonMissingForSite(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int minorAlleleCount(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int majorAlleleCount(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] genoCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] majorMinorCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AlleleDepth depth() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitStorage bitStorage(WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte referenceAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte alternateAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] referenceAlleles(int startSite, int endSite) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] referenceAlleleForAllSites() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public AlleleProbability alleleProbability() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public float alleleProbability(int taxon, int site, SITE_SCORE_TYPE type) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Dosage dosage() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte dosage(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public GeneralAnnotationStorage annotations() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public ReferenceProbability referenceProbability() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public float referenceProbability(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Stream<Byte> streamGenotype() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Stream<Byte> streamGenotype(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean hasSiteTranslations() {
        return false;
    }

    @Override
    public int[] siteTranslations() {
        return null;
    }

}
