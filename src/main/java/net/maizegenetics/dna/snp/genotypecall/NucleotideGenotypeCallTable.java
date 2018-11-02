/*
 *  NucleotideGenotype
 */
package net.maizegenetics.dna.snp.genotypecall;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;

/**
 * In memory GenotypeCallTable for nucleotide data.
 *
 * @author Terry Casstevens
 */
class NucleotideGenotypeCallTable extends ByteGenotypeCallTable {

    NucleotideGenotypeCallTable(SuperByteMatrix genotype, boolean phased) {
        super(genotype, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype(taxon, site));
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(genotypeAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public int maxNumAlleles() {
        return NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES;
    }

    @Override
    public boolean retainsRareAlleles() {
        return false;
    }
}
