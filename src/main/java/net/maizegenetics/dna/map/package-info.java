/**
 * Supports genetic and physical map representations of genomes including their annotation.
 * <p>
 *     Definitions:
 *     <p></p>
 *     Sites are the list of positions recorded in a genotype table.  They are organized into a list.
 *     <p></p>
 *     Positions are any defined point in the genome.  Generally they are polymorphic, but they do not need to
 *     be.
 *     <p></p>
 *     Chromosome can be viewed as physical chromosome, but it can also be used as an order set of positions.
 *     For example, chromosome
 *     could be used to reference a large series of contigs.
 *     <p></p>
 *     Annotations are descriptors about a position, e.g. reference allele, minor allele frequency, strand, etc.
 */
package net.maizegenetics.dna.map;
