/**
 * Genotype calls for all the taxa and genomic positions scored.
 * <p></p>
 * Genotype calls are stored as one byte per genotype, and they are accessed through their taxa and site index.  Generally,
 * the only class the API user need is only needs the {@link GenotypeTableBuilder},
 * which provides mechanisms for creating basic GenotypeTables.  For more complex manipulation of the GenotypeCalls,
 * {@link GenotypeCallTableBuilder} can be used.
 * <p></p>
 * A DNA nucleotide example are covered is covered below, but other genotype encoding systems
 * with less than 15 states can be supported with this code base.
 *
 Alleles [A, C, G, T, -(gap), +(insertion)], are stored in 4 bits. The class NucleotideAlignmentConstants stores this information.

 <p></p>
 Alleles
 Allele	Binary	Hex
 A	0000	0x0
 C	0001	0x1
 G	0010	0x2
 T	0011	0x3
 INS(+)	0100	0x4
 GAP(-)	0101	0x5
 N	1111	0xF

 <p></p>
 Genotypes
 Genotype	Binary	Hex
 A/A	00000000	0x0
 A/C	00000001	0x1
 C/A	00010000	0x10
 T/G	00110010	0x32
 -/-	01010101	0x55
 -/+=0	01010100	0x54
 N/N	11111111	0xFF
 <p></p>
 Genotypes are passed around using a single byte (8-bits), and are combinations of two alleles,
 with allele 1 in the first 4-bits and allele 2 in the last 4-bits.
 * @see net.maizegenetics.dna.snp.GenotypeTable
 * @see net.maizegenetics.dna.snp.GenotypeTableBuilder
 * @see net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable
 * @see net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder
 *
 *
 */

package net.maizegenetics.dna.snp.genotypecall;