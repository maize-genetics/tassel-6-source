/*
 * GenotypeSummaryPlugin
 */
package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.SimpleTableReport;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Tuple;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Terry Casstevens
 */
public class GenotypeSummaryPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(GenotypeSummaryPlugin.class);
    private static final String NA = "NA";
    private static final Double ZERO_DOUBLE = 0.0;
    private static final int ZERO_INT = 0;
    private long myNumGametesMissing = 0;
    private long myNumHeterozygous = 0;
    private double myAveMinorAlleleFreq = 0.0;

    private PluginParameter<Boolean> myOverview = new PluginParameter.Builder<>("overview", true, Boolean.class)
            .description("Get Overview Report").build();
    private PluginParameter<Boolean> mySiteSummary = new PluginParameter.Builder<>("siteSummary", true, Boolean.class)
            .description("Get Site Summary").build();
    private PluginParameter<Boolean> myTaxaSummary = new PluginParameter.Builder<>("taxaSummary", true, Boolean.class)
            .description("Get Taxa Summary").build();

    public GenotypeSummaryPlugin() {
        super(false);
        overview(false);
        siteSummary(false);
        taxaSummary(false);
    }

    public GenotypeSummaryPlugin(boolean isInteractive) {
        super(isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        DataSet.data(input, GenotypeTable.class);
    }

    @Override
    public DataSet processData(DataSet input) {

        if (!overview() && !siteSummary() && !taxaSummary()) {
            printSimpleSummary(input);
            return null;
        }

        myNumGametesMissing = 0;
        myNumHeterozygous = 0;
        myAveMinorAlleleFreq = 0.0;

        Tuple<String, GenotypeTable> data = DataSet.data(input, GenotypeTable.class);
        GenotypeTable alignment = data.y;
        String name = data.x;

        List<Datum> summaryTables = new ArrayList<>();

        SimpleTableReport siteSummary = null;
        if (siteSummary()) {
            siteSummary = getSiteSummary(alignment);
        }

        SimpleTableReport taxaSummary = null;
        if (taxaSummary()) {
            taxaSummary = getTaxaSummary(alignment);
        }

        SimpleTableReport[] overallSummaries = null;
        if (overview()) {
            overallSummaries = getOverallSummary(alignment);
            summaryTables.add(new Datum(name + "_OverallSummary", overallSummaries[0], "Overall Summary of " + name));
            summaryTables.add(new Datum(name + "_AlleleSummary", overallSummaries[1], "Allele Summary of " + name));
        }

        if (siteSummary != null) {
            summaryTables.add(new Datum(name + "_SiteSummary", siteSummary, "Site Summary of " + name));
        }
        if (taxaSummary != null) {
            summaryTables.add(new Datum(name + "_TaxaSummary", taxaSummary, "Taxa Summary of " + name));
        }

        if (summaryTables.isEmpty()) {
            return null;
        }

        return new DataSet(summaryTables, this);

    }

    private SimpleTableReport[] getOverallSummary(GenotypeTable alignment) {

        Object[] firstColumnNames = new String[]{"Stat Type", "Value"};

        long numSites = alignment.numberOfSites();
        long numTaxa = alignment.numberOfTaxa();

        Object[][] diploidValueCounts = alignment.genoCounts();
        int numAlleles = diploidValueCounts[0].length;

        if (!siteSummary()) {
            int totalGametes = (int) numTaxa * 2;
            for (int i = 0; i < numSites; i++) {
                int totalGametesNotMissing = alignment.totalGametesNonMissingForSite(i);
                int totalGametesMissing = totalGametes - totalGametesNotMissing;
                myNumGametesMissing = myNumGametesMissing + (long) totalGametesMissing;
                int numHeterozygous = alignment.heterozygousCount(i);
                myNumHeterozygous = myNumHeterozygous + (long) numHeterozygous;
                myAveMinorAlleleFreq += alignment.minorAlleleFrequency(i);
            }
            myAveMinorAlleleFreq /= numSites;
        }

        long totalGametes = numSites * numTaxa * 2L;
        long totalGametesNotMissing = totalGametes - myNumGametesMissing;

        long numDiploidsMissing = 0;
        for (int j = 0; j < numAlleles; j++) {
            if ((diploidValueCounts[0][j].equals(GenotypeTable.UNKNOWN_ALLELE_STR)) || (diploidValueCounts[0][j].equals(GenotypeTable.UNKNOWN_DIPLOID_ALLELE_STR))) {
                numDiploidsMissing = (Long) diploidValueCounts[1][j];
                break;
            }
        }

        long totalDiploids = numSites * numTaxa;
        long totalDiploidsNotMissing = totalDiploids - numDiploidsMissing;
        int count = 0;

        Object[][] data = new Object[15][firstColumnNames.length];

        data[count][0] = "Number of Taxa";
        data[count++][1] = (double) numTaxa;

        data[count][0] = "Number of Sites";
        data[count++][1] = (double) numSites;

        data[count][0] = "Sites x Taxa";
        data[count++][1] = (double) totalDiploids;

        data[count][0] = "Number Not Missing";
        data[count++][1] = (double) totalDiploidsNotMissing;

        data[count][0] = "Proportion Not Missing";
        data[count++][1] = (double) totalDiploidsNotMissing / (double) totalDiploids;

        data[count][0] = "Number Missing";
        data[count++][1] = (double) numDiploidsMissing;

        data[count][0] = "Proportion Missing";
        data[count++][1] = (double) numDiploidsMissing / (double) totalDiploids;

        data[count][0] = "Number Gametes";
        data[count++][1] = (double) totalGametes;

        data[count][0] = "Gametes Not Missing";
        data[count++][1] = (double) totalGametesNotMissing;

        data[count][0] = "Proportion Gametes Not Missing";
        data[count++][1] = (double) totalGametesNotMissing / (double) totalGametes;

        data[count][0] = "Gametes Missing";
        data[count++][1] = (double) myNumGametesMissing;

        data[count][0] = "Proportion Gametes Missing";
        data[count++][1] = (double) myNumGametesMissing / (double) totalGametes;

        data[count][0] = "Number Heterozygous";
        data[count++][1] = (double) myNumHeterozygous;

        data[count][0] = "Proportion Heterozygous";
        data[count++][1] = (double) myNumHeterozygous / (double) totalDiploids;

        data[count][0] = "Average Minor Allele Frequency";
        data[count++][1] = myAveMinorAlleleFreq;

        Object[][] majorMinorDiploidValueCounts = alignment.majorMinorCounts();
        int numMajorMinorAlleles = majorMinorDiploidValueCounts[0].length;

        Object[] alleleColumnNames = new String[]{"Alleles", "Number", "Proportion", "Frequency"};
        Object[][] data2 = new Object[numAlleles + numMajorMinorAlleles][alleleColumnNames.length];

        count = 0;
        for (int i = 0; i < numAlleles; i++) {
            String value = (String) diploidValueCounts[0][i];
            Long numValue = (Long) diploidValueCounts[1][i];
            data2[count][0] = value;
            data2[count][1] = numValue;
            data2[count][2] = numValue.doubleValue() / (double) totalDiploids;
            data2[count++][3] = numValue.doubleValue() / (double) totalDiploidsNotMissing;
        }

        for (int i = 0; i < numMajorMinorAlleles; i++) {
            String value = (String) majorMinorDiploidValueCounts[0][i];
            Long numValue = (Long) majorMinorDiploidValueCounts[1][i];
            data2[count][0] = value;
            data2[count][1] = numValue;
            data2[count++][2] = numValue.doubleValue() / (double) numSites;
        }

        return new SimpleTableReport[]{new SimpleTableReport("Overall Summary", firstColumnNames, data),
                new SimpleTableReport("Allele Summary", alleleColumnNames, data2)};
    }

    public static void printSimpleSummary(DataSet input) {
        if (input == null) {
            return;
        }
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if (alignInList.isEmpty()) {
            return;
        }
        printSimpleSummary(alignInList.get(0));
    }

    public static void printSimpleSummary(Datum current) {
        GenotypeTable alignment = (GenotypeTable) current.getData();
        String name = current.getName();
        printSimpleSummary(alignment, name);
    }

    public static void printSimpleSummary(GenotypeTable alignment, String name) {

        long numSites = alignment.numberOfSites();
        long numTaxa = alignment.numberOfTaxa();

        long totalDiploids = numSites * numTaxa;

        System.out.println("Genotype Table Name: " + name);
        System.out.println("Number of Taxa: " + numTaxa);
        System.out.println("Number of Sites: " + numSites);
        System.out.println("Sites x Taxa: " + totalDiploids);

        System.out.println("Chromosomes...");
        Chromosome[] chromosomes = alignment.chromosomes();
        PositionList positions = alignment.positions();
        for (int i = 0; i < chromosomes.length; i++) {
            int[] startEnd = alignment.firstLastSiteOfChromosome(chromosomes[i]);
            System.out.println(chromosomes[i].getName() + ": start site: " + startEnd[0] + " (" + positions.get(startEnd[0]).getPosition() + ") last site: " + startEnd[1] + " (" + positions.get(startEnd[1]).getPosition() + ") total: " + (startEnd[1] - startEnd[0] + 1));
        }
        System.out.println();

    }

    private SimpleTableReport getSiteSummary(GenotypeTable alignment) {

        String[] firstColumnNames = new String[]{"Site Number", "Site Name", "Chromosome", "Physical Position", "Number of Taxa", "Ref", "Alt", "Major Allele", "Major Allele Gametes", "Major Allele Proportion", "Major Allele Frequency",
                "Minor Allele", "Minor Allele Gametes", "Minor Allele Proportion", "Minor Allele Frequency"};
        String[] lastColumnNames = new String[]{"Gametes Missing", "Proportion Missing", "Number Heterozygous", "Proportion Heterozygous",
                "Inbreeding Coefficient", "Inbreeding Coefficient Scaled by Missing"};

        List<String> columnNames = new ArrayList<>(Arrays.asList(firstColumnNames));

        int maxAlleles = alignment.maxNumAlleles();
        if (alignment.retainsRareAlleles()) {
            maxAlleles++;
        }
        for (int i = 2; i < maxAlleles; i++) {
            String alleleHeading = "Allele " + (i + 1);
            columnNames.add(alleleHeading);
            columnNames.add(alleleHeading + " Gametes");
            columnNames.add(alleleHeading + " Proportion");
            columnNames.add(alleleHeading + " Frequency");
        }

        columnNames.addAll(Arrays.asList(lastColumnNames));

        int numSites = alignment.numberOfSites();
        int numTaxa = alignment.numberOfTaxa();
        Object[][] data = new Object[numSites][columnNames.size()];
        int totalGametes = numTaxa * 2;

        for (int i = 0; i < numSites; i++) {

            int count = 0;

            data[i][count++] = i;
            data[i][count++] = alignment.siteName(i);
            data[i][count++] = alignment.chromosomeName(i);
            data[i][count++] = alignment.chromosomalPosition(i);
            data[i][count++] = numTaxa;
            data[i][count++] = alignment.genotypeAsString(i, alignment.referenceAllele(i));
            data[i][count++] = alignment.genotypeAsString(i, alignment.alternateAllele(i));

            int[][] alleles = alignment.allelesSortedByFrequency(i);
            int numAlleles = alleles[0].length;
            int totalNotMissing = alignment.totalNonMissingForSite(i);
            int totalGametesNotMissing = AlleleFreqCache.totalGametesNonMissingForSite(alleles);

            for (int a = 0; a < numAlleles; a++) {
                data[i][count++] = alignment.genotypeAsString(i, (byte) alleles[0][a]);
                data[i][count++] = alleles[1][a];
                data[i][count++] = (double) alleles[1][a] / (double) totalGametes;
                double alleleFreq = (double) alleles[1][a] / (double) totalGametesNotMissing;
                data[i][count++] = alleleFreq;
                if (a == 1) {
                    myAveMinorAlleleFreq += alleleFreq;
                }
            }

            for (int b = 0; b < (maxAlleles - numAlleles); b++) {
                data[i][count++] = NA;
                data[i][count++] = ZERO_INT;
                data[i][count++] = ZERO_DOUBLE;
                data[i][count++] = ZERO_DOUBLE;
            }

            int totalGametesMissing = totalGametes - totalGametesNotMissing;
            myNumGametesMissing = myNumGametesMissing + (long) totalGametesMissing;
            data[i][count++] = totalGametesMissing;
            data[i][count++] = (double) totalGametesMissing / (double) totalGametes;

            int numHeterozygous = alignment.heterozygousCount(i);
            myNumHeterozygous = myNumHeterozygous + (long) numHeterozygous;
            data[i][count++] = numHeterozygous;
            data[i][count++] = (double) numHeterozygous / (double) totalNotMissing;

            data[i][count++] = "TBD";
            data[i][count++] = "TBD";

        }

        myAveMinorAlleleFreq /= numSites;

        String[] columnNameStrings = new String[columnNames.size()];
        columnNames.toArray(columnNameStrings);
        return new SimpleTableReport("Site Summary", columnNameStrings, data);

    }

    private SimpleTableReport getTaxaSummary(GenotypeTable alignment) {

        Object[] columnNames = new String[]{"Taxa", "Taxa Name", "Number of Sites", "Gametes Missing", "Proportion Missing",
                "Number Heterozygous", "Proportion Heterozygous", "Inbreeding Coefficient",
                "Inbreeding Coefficient Scaled by Missing"};
        int numSites = alignment.numberOfSites();
        int numTaxa = alignment.numberOfTaxa();
        int totalGametes = numSites * 2;

        int[] gametesNotMissing = new int[numTaxa];
        int[] sitesNotMissing = new int[numTaxa];
        int[] heterozygous = new int[numTaxa];
        for (int s = 0; s < numSites; s++) {
            byte[] genotypes = alignment.genotypeAllTaxa(s);
            for (int t = 0; t < numTaxa; t++) {
                byte[] alleles = GenotypeTableUtils.getDiploidValues(genotypes[t]);
                if (alleles[0] != GenotypeTable.UNKNOWN_ALLELE) {
                    if (alleles[1] != GenotypeTable.UNKNOWN_ALLELE) {
                        gametesNotMissing[t] += 2;
                        sitesNotMissing[t]++;
                    } else {
                        gametesNotMissing[t]++;
                        sitesNotMissing[t]++;
                    }
                } else if (alleles[1] != GenotypeTable.UNKNOWN_ALLELE) {
                    gametesNotMissing[t]++;
                    sitesNotMissing[t]++;
                }

                if (alleles[0] != alleles[1]) {
                    heterozygous[t]++;
                }
            }
        }

        Object[][] data = new Object[numTaxa][columnNames.length];
        for (int i = 0; i < numTaxa; i++) {

            int totalGametesMissing = totalGametes - gametesNotMissing[i];

            int count = 0;
            data[i][count++] = i;
            data[i][count++] = alignment.taxaName(i);
            data[i][count++] = numSites;
            data[i][count++] = totalGametesMissing;
            data[i][count++] = (double) totalGametesMissing / (double) totalGametes;
            data[i][count++] = heterozygous[i];
            data[i][count++] = (double) heterozygous[i] / (double) sitesNotMissing[i];
            data[i][count++] = "Inbreeding Coefficient";
            data[i][count++] = "ICSBM";
        }

        return new SimpleTableReport("Taxa Summary", columnNames, data);

    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(GenotypeSummaryPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TableReport[] runPlugin(GenotypeTable genotype) {
        DataSet input = new DataSet(new Datum("Genotype Table", genotype, null), this);
        DataSet dataSet = performFunction(input);
        TableReport[] result = new TableReport[dataSet.getSize()];
        for (int i = 0; i < dataSet.getSize(); i++) {
            result[i] = (TableReport) dataSet.getData(i).getData();
        }
        return result;
    }

    /**
     * Get Overview Report
     *
     * @return Overview
     */
    public Boolean overview() {
        return myOverview.value();
    }

    /**
     * Set Overview. Get Overview Report
     *
     * @param value Overview
     *
     * @return this plugin
     */
    public GenotypeSummaryPlugin overview(Boolean value) {
        myOverview = new PluginParameter<>(myOverview, value);
        return this;
    }

    /**
     * Get Site Summary
     *
     * @return Site Summary
     */
    public Boolean siteSummary() {
        return mySiteSummary.value();
    }

    /**
     * Set Site Summary. Get Site Summary
     *
     * @param value Site Summary
     *
     * @return this plugin
     */
    public GenotypeSummaryPlugin siteSummary(Boolean value) {
        mySiteSummary = new PluginParameter<>(mySiteSummary, value);
        return this;
    }

    /**
     * Get Taxa Summary
     *
     * @return Taxa Summary
     */
    public Boolean taxaSummary() {
        return myTaxaSummary.value();
    }

    /**
     * Set Taxa Summary. Get Taxa Summary
     *
     * @param value Taxa Summary
     *
     * @return this plugin
     */
    public GenotypeSummaryPlugin taxaSummary(Boolean value) {
        myTaxaSummary = new PluginParameter<>(myTaxaSummary, value);
        return this;
    }

    @Override
    public String icon() {
        return "/images/summary.gif";
    }

    @Override
    public String getButtonName() {
        return "Genotype Summary";
    }

    @Override
    public String getToolTipText() {
        return "Genotype Summary";
    }

    @Override
    public String pluginUserManualURL() {
        return "https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/GenoSummary/GenoSummary";
    }

}
