/**
 * Data structures for holding raw sequence data (only Fastq currently), including Illumina, PacBio, Nanapore, etc.
 * This package provides methods for quality control, reformatting and contiging pair end sequence.
 *<p>
 * Definitions:
 * <li>Read:  a single sequence from a sequencer
 * <Li>PERead a paired end sequence, including ReadF and ReadB
 * <li>ReadF:  a single sequence from paired end sequence of R1
 * <li>ReadB:  a single sequence from paired end sequence of R2
 * <li>Pcontig:  a Pcontig is formed by a PE Tag whose forward tag and backward tag have overlap 
 *<p>
 * This package provides methods for quality control, reformatting and contiging pair end sequence.
 * All sequences are stored in string. Large chunk of reads do not go into memory
 */
package net.maizegenetics.dna.read;
