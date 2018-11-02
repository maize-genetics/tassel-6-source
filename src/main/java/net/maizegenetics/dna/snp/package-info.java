/**
 * Data structures for holding SNP and indel variation across taxa (samples) and genomic positions.
 * For information on other DNA related data structures see {@link net.maizegenetics.dna.package-info.java}
 * <p></p>
 * {@link net.maizegenetics.dna.snp.GenotypeTable} is the key interface to accessing SNP and indel variation.  GenotypeTables
 * are immutable, but they can be built through the {@link net.maizegenetics.dna.snp.GenotypeTableBuilder}.
 * GenotypeBuilder at a basic level are the product of a TaxaList (describing the germplasm or samples),
 * a PositionList (describing genomic locations), and GenotypeCalls (the SNP or indels calls for each taxon and position).
 * <p></p>
 *
 * @see net.maizegenetics.dna.map.PositionList
 * @see net.maizegenetics.taxa.TaxaList
 */
package net.maizegenetics.dna.snp;