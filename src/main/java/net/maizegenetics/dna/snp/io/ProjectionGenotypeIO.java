package net.maizegenetics.dna.snp.io;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.ProjectionBuilder;
import net.maizegenetics.dna.snp.genotypecall.ProjectionGenotypeCallTable;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.dna.map.DonorHaplotypes;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.*;

/**
 * Methods for reading and writing ProjectionGenotypes to files. ProjectionGenotypes have two parts - one is the
 * Projection file that has the names of high density genotyped taxa and the recombination breakpoints for
 * each of the low density taxa that point to the high density taxa.
 *
 * @author Ed Buckler
 */
public class ProjectionGenotypeIO {

    /**
     * Returns a genotypeTable based on a projection genotype file and high density genotype table
     * @param paFile file name for projection file
     * @param baseHighDensityAlignmentFile file name for high density file
     * @return GenotypeTable based on both
     */
    public static GenotypeTable getInstance(String paFile, String baseHighDensityAlignmentFile) {
        if(baseHighDensityAlignmentFile.endsWith(".h5")) return getInstance(paFile, GenotypeTableBuilder.getInstance(baseHighDensityAlignmentFile));
        return getInstance(paFile, ImportUtils.readFromHapmap(baseHighDensityAlignmentFile, null));
    }

    /**
     * Returns a genotypeTable based on a projection genotype file and high density genotype table
     * @param paFile file name for projection file
     * @param baseHighDensityAlignment GenotypeTable of high density taxa
     * @return Projection GenotypeTable based on both
     */
    public static GenotypeTable getInstance(String paFile, GenotypeTable baseHighDensityAlignment) {
        BufferedReader br = null;
        try {
            br = Utils.getBufferedReader(paFile);
            String[] sl=Utils.readLineSkipComments(br).split("\t");
            int baseTaxaCnt=Integer.parseInt(sl[0]);
            if(baseTaxaCnt!=baseHighDensityAlignment.numberOfTaxa()) {
                System.err.println("Error in number of base taxa"); return null;
            }
            int taxaCnt=Integer.parseInt(sl[1]);
            //Reading the map of high density taxa index to high density taxa name, e.g. 1 B73\n
            Map<Integer,Integer> paIndexToBaseIndex=new HashMap<>();
            for (int i = 0; i < baseTaxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                int index=Integer.parseInt(sl[0]);
                Taxon taxon=new Taxon(sl[1]);
                int matches=baseHighDensityAlignment.taxa().indexOf(taxon);
                if(matches<0) {
                    throw new NoSuchElementException("Taxon "+sl[1]+" not found within base taxa");
                }
                paIndexToBaseIndex.put(index, matches);
            }
            //Reading the recombination breakpoints for each of the projected taxa
            ProjectionBuilder pb=new ProjectionBuilder(baseHighDensityAlignment);
            for (int i = 0; i < taxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                int breakTotal=sl.length-1;
                if(breakTotal==0) continue;  //no data
                NavigableSet<DonorHaplotypes> breakForTaxon=new TreeSet<>();
                for (int bp = 0; bp < breakTotal; bp++) {
                    String[] bptext=sl[bp+1].split(":");
                    Chromosome chr=new Chromosome(bptext[0]);
                    int baseParent1=paIndexToBaseIndex.get(Integer.parseInt(bptext[3]));
                    int baseParent2=paIndexToBaseIndex.get(Integer.parseInt(bptext[4]));
                    DonorHaplotypes dh=new DonorHaplotypes(chr, Integer.parseInt(bptext[1]), Integer.parseInt(bptext[2]),
                            baseParent1, baseParent2);
                    breakForTaxon.add(dh);
                }
                pb.addTaxon(new Taxon(sl[0]), breakForTaxon);
            }
            return pb.build();
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error reading Projection file: " + paFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Exports the ProjectionGenotypes to a file
     * @param outfile the path and name of the projection file
     * @param pa GenotypeTable with a ProjectionGenotypeCallTable
     */
    public static void writeToFile(String outfile, GenotypeTable pa) {
        if(!(pa.genotypeMatrix() instanceof ProjectionGenotypeCallTable)) {
            throw new UnsupportedOperationException("Save only works for Alignments with projection genotypes");
        }
        ProjectionGenotypeCallTable pg=(ProjectionGenotypeCallTable)pa.genotypeMatrix();
        GenotypeTable baseAlignment=pg.getBaseGenotypeTable();
        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(outfile, ".pa.txt.gz", new String[]{".pa.txt", ".pa.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
            bw.write(baseAlignment.numberOfTaxa()+"\t"+pa.numberOfTaxa()+"\n");
            bw.write("#Donor Haplotypes\n");
            for (int i = 0; i < baseAlignment.numberOfTaxa(); i++) {
                bw.write(i+"\t"+baseAlignment.taxaName(i)+"\n");
            }
            bw.write("#Taxa Breakpoints\n");
            bw.write("#Block are defined chr:startPos:endPos:donor1:donor2 (-1 means no hypothesis)\n");
            for (int i = 0; i < pa.numberOfTaxa(); i++) {
                bw.write(pa.taxaName(i)+"\t");
                NavigableSet<DonorHaplotypes> theDH=pg.getDonorHaplotypes(i);
                for (DonorHaplotypes dh : theDH) {
                    bw.write(dh.getChromosome().getName()+":"+dh.getStartPosition()+":"+dh.getEndPosition()+":"+
                            dh.getParent1index()+":"+dh.getParent2index()+"\t");
                }
                bw.write("\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing Projection file: " + outfile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

}
