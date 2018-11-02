/*
 * CombineGenotypeTable
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.score.AlleleProbability;
import net.maizegenetics.dna.snp.score.ReferenceProbability;
import net.maizegenetics.dna.snp.score.AlleleDepth;
import net.maizegenetics.dna.snp.score.Dosage;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.dna.snp.genotypecall.CombineGenotypeCallTable;
import net.maizegenetics.util.GeneralAnnotationStorage;

import java.util.*;
import java.util.stream.Stream;

/**
 * Combines multiple GenotypeTables together.
 *
 * @author Terry Casstevens
 */
public class CombineGenotypeTable implements GenotypeTable {

    private final GenotypeTable[] myGenotypeTables;
    private final GenotypeCallTable myGenotype;
    private final int[] mySiteOffsets;
    private final Map<Chromosome, GenotypeTable> myChromosomes = new HashMap<>();
    private Chromosome[] myChromosomesList;
    private int[] myChromosomesOffsets;
    private final TaxaList myTaxaList;
    private String[][] myAlleleStates;
    private PositionList myPositions = null;

    private CombineGenotypeTable(TaxaList taxaList, GenotypeTable[] genoTables) {

        myTaxaList = taxaList;
        myGenotypeTables = genoTables;
        mySiteOffsets = new int[genoTables.length + 1];

        mySiteOffsets[0] = 0;
        int count = 0;
        for (int i = 0; i < genoTables.length; i++) {
            count = genoTables[i].numberOfSites() + count;
            mySiteOffsets[i + 1] = count;

            Chromosome[] chromosomes = genoTables[i].chromosomes();
            for (int j = 0; j < chromosomes.length; j++) {
                myChromosomes.put(chromosomes[j], genoTables[i]);
            }
        }

        initChromosomes();

        GenotypeCallTable[] genotypeCallTables = new GenotypeCallTable[genoTables.length];
        for (int i = 0; i < genoTables.length; i++) {
            genotypeCallTables[i] = genoTables[i].genotypeMatrix();
        }
        myGenotype = CombineGenotypeCallTable.getInstance(genotypeCallTables);
    }

    /**
     * This factory method combines given genoTables. If only one genotypeTable,
     * then it is returned unchanged. Otherwise, this requires that each
     * genotypeTable has the same Taxa in the same order.
     *
     * @param genoTables
     * @return
     */
    public static GenotypeTable getInstance(GenotypeTable[] genoTables) {

        if ((genoTables == null) || (genoTables.length == 0)) {
            throw new IllegalArgumentException("CombineGenotypeTable: getInstance: must provide genoTables.");
        }

        if (genoTables.length == 1) {
            return genoTables[0];
        }

        TaxaList firstGroup = genoTables[0].taxa();
        for (int i = 1; i < genoTables.length; i++) {
            if (!areTaxaListsEqual(firstGroup, genoTables[i].taxa())) {
                throw new IllegalArgumentException("CombineGenotypeTable: getInstance: TaxaLists do not match.");
            }
        }

        return new CombineGenotypeTable(firstGroup, genoTables);

    }

    /**
     * This factory method combines given genoTables. If only one genotypeTable,
     * then it is returned unchanged. If isUnion equals true, a union join of
     * the Identifiers will be used to construct the combination. Any
     * genotypeTable not containing one of the Identifiers will return unknown
     * value for those locations. If isUnion equals false, a intersect join of
     * the Identifiers will be used.
     *
     * @param genoTables genoTables to combine
     * @param isUnion whether to union or intersect join
     * @return
     */
    public static GenotypeTable getInstance(GenotypeTable[] genoTables, boolean isUnion) {

        if ((genoTables == null) || (genoTables.length == 0)) {
            throw new IllegalArgumentException("CombineGenotypeTable: getInstance: must provide genoTables.");
        }

        if (genoTables.length == 1) {
            return genoTables[0];
        }

        TaxaList[] groups = new TaxaList[genoTables.length];
        for (int i = 0; i < genoTables.length; i++) {
            groups[i] = genoTables[i].taxa();
        }
        TaxaList newTaxa;
        if (isUnion) {
            newTaxa = TaxaListUtils.getAllTaxa(groups, false);
        } else {
            newTaxa = TaxaListUtils.getCommonTaxa(groups, false);
        }

        GenotypeTable[] newAlignmentNews = new GenotypeTable[genoTables.length];
        for (int i = 0; i < genoTables.length; i++) {
            newAlignmentNews[i] = FilterGenotypeTable.getInstance(genoTables[i], newTaxa, true);
        }
        
        for (int i = 0; i < genoTables.length; i++) {
            if (!areTaxaListsEqual(newTaxa, newAlignmentNews[i].taxa())) {
                throw new IllegalArgumentException("CombineGenotypeTable: getInstance: TaxaLists do not match.");
            }
        }

        return new CombineGenotypeTable(newTaxa, newAlignmentNews);

    }

    private static boolean areTaxaListsEqual(TaxaList first, TaxaList second) {

        if (first.numberOfTaxa() != second.numberOfTaxa()) {
            return false;
        }

        for (int i = 0, n = first.numberOfTaxa(); i < n; i++) {
            if (!first.get(i).equals(second.get(i))) {
                return false;
            }
        }

        return true;

    }

    private void initChromosomes() {

        List<Integer> offsets = new ArrayList<>();
        List<Chromosome> chromosomes = new ArrayList<>();
        for (int i = 0; i < myGenotypeTables.length; i++) {
            chromosomes.addAll(Arrays.asList(myGenotypeTables[i].chromosomes()));
            int[] tempOffsets = myGenotypeTables[i].chromosomesOffsets();
            for (int j = 0; j < tempOffsets.length; j++) {
                offsets.add(tempOffsets[j] + mySiteOffsets[i]);
            }
        }

        myChromosomesList = new Chromosome[chromosomes.size()];
        myChromosomesList = chromosomes.toArray(myChromosomesList);

        myChromosomesOffsets = new int[offsets.size()];
        for (int i = 0; i < offsets.size(); i++) {
            myChromosomesOffsets[i] = offsets.get(i);
        }

        if (myChromosomesOffsets.length != myChromosomesList.length) {
            throw new IllegalStateException("CombineGenotypeTable: initChromosomes: number chromosomes offsets should equal number of chromosomes.");
        }

    }

    @Override
    public byte genotype(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotype(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {

        byte[] result = new byte[endSite - startSite];
        int count = 0;
        int firstAlign = translateSite(startSite);
        int secondAlign = translateSite(endSite - 1);
        for (int i = firstAlign; i <= secondAlign; i++) {
            int firstSite = 0;
            if (i == firstAlign) {
                firstSite = startSite - mySiteOffsets[firstAlign];
            }
            int secondSite = 0;
            if (firstAlign == secondAlign) {
                secondSite = endSite - mySiteOffsets[firstAlign];
            } else if (i != secondAlign) {
                secondSite = myGenotypeTables[i].numberOfSites();
            } else {
                secondSite = endSite - mySiteOffsets[secondAlign];
            }
            for (int s = firstSite; s < secondSite; s++) {
                result[count++] = myGenotypeTables[i].genotype(taxon, s);
            }
        }
        return result;

    }

    @Override
    public byte genotype(int taxon, Chromosome locus, int physicalPosition) {
        int site = siteOfPhysicalPosition(physicalPosition, locus);
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotype(taxon, site - mySiteOffsets[translate]);
    }

    /**
     * Returns which genotypeTable to use.
     *
     * @param site
     * @return genotypeTable index.
     */
    public int translateSite(int site) {

        for (int i = 1; i < mySiteOffsets.length; i++) {
            if (mySiteOffsets[i] > site) {
                return i - 1;
            }
        }
        throw new IndexOutOfBoundsException("CombineGenotypeTable: translateSite: index out of range: " + site);

    }

    private int findGenotypeTableIndex(GenotypeTable genotypeTable) {
        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (genotypeTable == myGenotypeTables[i]) {
                return i;
            }
        }
        throw new IllegalArgumentException("CombineGenotypeTable: findGenotypeTableIndex: Genotype Table unknown.");
    }

    @Override
    public boolean hasReference() {

        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (!myGenotypeTables[i].hasReference()) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String siteName(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].siteName(site - mySiteOffsets[translate]);
    }

    @Override
    public int numberOfSites() {
        return mySiteOffsets[mySiteOffsets.length - 1];
    }

    @Override
    public int chromosomeSiteCount(Chromosome locus) {
        return myChromosomes.get(locus).chromosomeSiteCount(locus);
    }

    @Override
    public int chromosomalPosition(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].chromosomalPosition(site - mySiteOffsets[translate]);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome locus) {
        GenotypeTable align = myChromosomes.get(locus);
        int i = -1;
        for (int j = 0; j < myGenotypeTables.length; j++) {
            if (myGenotypeTables[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.siteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome locus, String snpName) {
        GenotypeTable align = myChromosomes.get(locus);
        int i = -1;
        for (int j = 0; j < myGenotypeTables.length; j++) {
            if (myGenotypeTables[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.siteOfPhysicalPosition(physicalPosition, locus, snpName);
    }

    @Override
    public Chromosome chromosome(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].chromosome(site - mySiteOffsets[translate]);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myChromosomesList;
    }

    @Override
    public int numChromosomes() {
        if (myChromosomesList == null) {
            return 0;
        } else {
            return myChromosomesList.length;
        }
    }

    @Override
    public int indelSize(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].indelSize(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isIndel(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].isIndel(site - mySiteOffsets[translate]);
    }

    @Override
    public byte referenceAllele(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].referenceAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte alternateAllele(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].alternateAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public GenotypeTable[] compositeAlignments() {
        return myGenotypeTables;
    }

    @Override
    public byte majorAllele(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].majorAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte minorAllele(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].minorAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] minorAlleles(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].minorAlleles(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] alleles(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].alleles(site - mySiteOffsets[translate]);
    }

    @Override
    public double minorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].minorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].allelesSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotypeArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeAllTaxa(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotypeAllTaxa(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        byte[] result = new byte[numberOfSites()];
        for (int i = 0; i < myGenotypeTables.length; i++) {
            byte[] current = myGenotypeTables[i].genotypeAllSites(taxon);
            System.arraycopy(current, 0, result, myChromosomesOffsets[i], current.length);
        }
        return result;
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("CombineGenotypeTable: getAllelePresenceForAllSites: This operation isn't possible as it spans multiple GenotypeTables.");
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, WHICH_ALLELE allele, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("CombineGenotypeTable: getAllelePresenceForSitesBlock: This operation isn't possible as it spans multiple GenotypeTables.");
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotypeAsString(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotypeAsStringArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] referenceAlleles(int startSite, int endSite) {
        int numSites = endSite - startSite;
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = referenceAllele(startSite + i);
        }
        return result;
    }

    @Override
    public byte[] referenceAlleleForAllSites() {

        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (!myGenotypeTables[i].hasReference()) {
                return null;
            }
        }

        byte[] result = new byte[numberOfSites()];
        int count = 0;
        for (int i = 0; i < myGenotypeTables.length; i++) {
            byte[] current = myGenotypeTables[i].referenceAlleleForAllSites();
            for (int j = 0; j < current.length; j++) {
                result[count++] = current[j];
            }
        }
        return result;

    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].isHeterozygous(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public int[] physicalPositions() {

        boolean allNull = true;
        for (int i = 0; i < myGenotypeTables.length; i++) {
            int[] current = myGenotypeTables[0].physicalPositions();
            if ((current != null) && (current.length != 0)) {
                allNull = false;
                break;
            }
        }

        if (allNull) {
            return null;
        } else {
            int[] result = new int[numberOfSites()];
            int count = 0;
            for (int i = 0; i < myGenotypeTables.length; i++) {
                int[] current = myGenotypeTables[i].physicalPositions();
                for (int j = 0; j < current.length; j++) {
                    result[count++] = current[j];
                }
            }
            return result;
        }
    }

    @Override
    public String chromosomeName(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].chromosomeName(site - mySiteOffsets[translate]);
    }

    @Override
    public int[] chromosomesOffsets() {
        return myChromosomesOffsets;
    }

    @Override
    public Set<SITE_SCORE_TYPE> siteScoreTypes() {
        Set<SITE_SCORE_TYPE> result = new LinkedHashSet<>();
        for (int i = 0; i < myGenotypeTables.length; i++) {
            result.addAll(myGenotypeTables[i].siteScoreTypes());
        }
        return result;
    }

    @Override
    public boolean isAllPolymorphic() {
        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (!myGenotypeTables[i].isAllPolymorphic()) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean isPolymorphic(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].isPolymorphic(site - mySiteOffsets[translate]);
    }

    @Override
    public double majorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].majorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public String genomeVersion() {
        String first = myGenotypeTables[0].genomeVersion();
        if (first == null) {
            return null;
        }
        for (int i = 1; i < myGenotypeTables.length; i++) {
            String current = myGenotypeTables[i].genomeVersion();
            if ((current != null) && (!first.equals(current))) {
                return null;
            }
        }
        return first;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].isPositiveStrand(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isPhased() {
        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (myGenotypeTables[i].isPhased() == false) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean retainsRareAlleles() {
        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (myGenotypeTables[i].retainsRareAlleles() == false) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String[][] alleleDefinitions() {

        if (myAlleleStates != null) {
            return myAlleleStates;
        }

        boolean allTheSame = true;
        String[][] encodings = myGenotypeTables[0].alleleDefinitions();
        if (encodings.length == 1) {
            for (int i = 1; i < myGenotypeTables.length; i++) {
                String[][] current = myGenotypeTables[i].alleleDefinitions();
                if ((current.length == 1) && (encodings[0].length == current[0].length)) {
                    for (int j = 0; j < encodings[0].length; j++) {
                        if (!current[0][j].equals(encodings[0][j])) {
                            allTheSame = false;
                            break;
                        }
                    }
                } else {
                    allTheSame = false;
                    break;
                }

                if (!allTheSame) {
                    break;
                }
            }
        } else {
            allTheSame = false;
        }

        if (allTheSame) {
            myAlleleStates = encodings;
        } else {
            String[][] result = new String[numberOfSites()][];
            int count = 0;
            for (int i = 0; i < myGenotypeTables.length; i++) {
                for (int j = 0, n = myGenotypeTables[i].numberOfSites(); j < n; j++) {
                    result[count++] = myGenotypeTables[i].alleleDefinitions(j);
                }
            }
            myAlleleStates = result;
        }

        return myAlleleStates;

    }

    @Override
    public String[] alleleDefinitions(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].alleleDefinitions(site - mySiteOffsets[translate]);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genotypeAsString(site - mySiteOffsets[translate], value);
    }

    @Override
    public int maxNumAlleles() {
        int result = 999999;
        for (int i = 0; i < myGenotypeTables.length; i++) {
            if (myGenotypeTables[i].maxNumAlleles() < result) {
                result = myGenotypeTables[i].maxNumAlleles();
            }
        }
        return result;
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].totalGametesNonMissingForSite(site - mySiteOffsets[translate]);
    }

    @Override
    public int heterozygousCount(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].heterozygousCount(site - mySiteOffsets[translate]);
    }

    @Override
    public int minorAlleleCount(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].minorAlleleCount(site - mySiteOffsets[translate]);
    }

    @Override
    public int majorAlleleCount(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].majorAlleleCount(site - mySiteOffsets[translate]);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].genosSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].allelesBySortType(scope, site - mySiteOffsets[translate]);
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, WHICH_ALLELE allele) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].allelePresenceForAllTaxa(site - mySiteOffsets[translate], allele);
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
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        int firstGenotype = translateSite(startSite);
        int secondGenotype = translateSite(endSite);
        if (firstGenotype == secondGenotype) {
            return myGenotypeTables[firstGenotype].genotypeAsStringRange(taxon, startSite - mySiteOffsets[firstGenotype], endSite - mySiteOffsets[firstGenotype]);
        } else if (secondGenotype - firstGenotype == 1) {
            StringBuilder builder = new StringBuilder();
            builder.append(myGenotypeTables[firstGenotype].genotypeAsStringRange(taxon, startSite - mySiteOffsets[firstGenotype], myGenotypeTables[firstGenotype].numberOfSites()));
            builder.append(";");
            builder.append(myGenotypeTables[secondGenotype].genotypeAsStringRange(taxon, 0, endSite - mySiteOffsets[secondGenotype]));
            return builder.toString();
        } else {
            StringBuilder builder = new StringBuilder();
            builder.append(myGenotypeTables[firstGenotype].genotypeAsStringRange(taxon, startSite - mySiteOffsets[firstGenotype], myGenotypeTables[firstGenotype].numberOfSites()));
            for (int i = firstGenotype + 1; i < secondGenotype; i++) {
                builder.append(";");
                builder.append(myGenotypeTables[i].genotypeAsStringRow(taxon));
            }
            builder.append(";");
            builder.append(myGenotypeTables[secondGenotype].genotypeAsStringRange(taxon, 0, endSite - mySiteOffsets[secondGenotype]));
            return builder.toString();
        }
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        StringBuilder builder = new StringBuilder();
        boolean first = true;
        for (GenotypeTable current : myGenotypeTables) {
            if (first) {
                first = false;
            } else {
                builder.append(";");
            }
            builder.append(current.genotypeAsStringRow(taxon));
        }
        return builder.toString();
    }

    @Override
    public int[] firstLastSiteOfChromosome(Chromosome chromosome) {
        GenotypeTable genotypeTable = myChromosomes.get(chromosome);
        int index = findGenotypeTableIndex(genotypeTable);
        int[] result = genotypeTable.firstLastSiteOfChromosome(chromosome);
        result[0] += myChromosomesOffsets[index];
        result[1] += myChromosomesOffsets[index];
        return result;
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaList.size();
    }

    @Override
    public Chromosome chromosome(String name) {
        for (Chromosome current : myChromosomesList) {
            if (current.getName().equals(name)) {
                return current;
            }
        }
        return null;
    }

    @Override
    public String majorAlleleAsString(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].majorAlleleAsString(site - mySiteOffsets[translate]);
    }

    @Override
    public String minorAlleleAsString(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].minorAlleleAsString(site - mySiteOffsets[translate]);
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.taxaName(index);
    }

    @Override
    public String diploidAsString(int site, byte value) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].diploidAsString(site - mySiteOffsets[translate], value);
    }

    @Override
    public int totalNonMissingForSite(int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].totalNonMissingForSite(site - mySiteOffsets[translate]);
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
        int result = 0;
        for (GenotypeTable current : myGenotypeTables) {
            result += current.totalGametesNonMissingForTaxon(taxon);
        }
        return result;
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        int result = 0;
        for (GenotypeTable current : myGenotypeTables) {
            result += current.heterozygousCountForTaxon(taxon);
        }
        return result;
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        int result = 0;
        for (GenotypeTable current : myGenotypeTables) {
            result += current.totalNonMissingForTaxon(taxon);
        }
        return result;
    }

    @Override
    public boolean hasGenotype() {
        boolean result = true;
        for (GenotypeTable current : myGenotypeTables) {
            if (!current.hasGenotype()) {
                result = false;
            }
        }
        return result;
    }

    @Override
    public boolean hasDepth() {
        boolean result = true;
        for (GenotypeTable current : myGenotypeTables) {
            if (!current.hasDepth()) {
                result = false;
            }
        }
        return result;
    }

    @Override
    public boolean hasAlleleProbabilities() {
        boolean result = true;
        for (GenotypeTable current : myGenotypeTables) {
            if (!current.hasAlleleProbabilities()) {
                result = false;
            }
        }
        return result;
    }

    @Override
    public boolean hasReferenceProbablity() {
        boolean result = true;
        for (GenotypeTable current : myGenotypeTables) {
            if (!current.hasReferenceProbablity()) {
                result = false;
            }
        }
        return result;
    }

    @Override
    public boolean hasDosage() {
        boolean result = true;
        for (GenotypeTable current : myGenotypeTables) {
            if (!current.hasDosage()) {
                result = false;
            }
        }
        return result;
    }

    @Override
    public AlleleDepth depth() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].depthForAlleles(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public BitStorage bitStorage(WHICH_ALLELE allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public PositionList positions() {
        if (myPositions == null) {
            PositionListBuilder builder = new PositionListBuilder();
            for (GenotypeTable current : myGenotypeTables) {
                builder.addAll(current.positions());
            }
            myPositions = builder.build();
        }
        return myPositions;
    }

    @Override
    public GenotypeCallTable genotypeMatrix() {
        return myGenotype;
    }

    @Override
    public AlleleProbability alleleProbability() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public float alleleProbability(int taxon, int site, SITE_SCORE_TYPE type) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].alleleProbability(taxon, site - mySiteOffsets[translate], type);
    }

    @Override
    public ReferenceProbability referenceProbability() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public float referenceProbability(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].referenceProbability(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public Dosage dosage() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte dosage(int taxon, int site) {
        int translate = translateSite(site);
        return myGenotypeTables[translate].dosage(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public GeneralAnnotationStorage annotations() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Stream<Byte> streamGenotype() {
        return myGenotype.stream();
    }

    @Override
    public Stream<Byte> streamGenotype(int taxon) {
        return myGenotype.stream(taxon);
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
