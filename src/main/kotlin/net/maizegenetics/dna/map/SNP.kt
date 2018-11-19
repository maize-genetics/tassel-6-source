package net.maizegenetics.dna.map

/**
 * @author Terry Casstevens
 * Created November 12, 2018
 */

class SNP(val chr: Chromosome, val position: Int, val snpName: String) : GenomicFactor(startChr = chr, startPos = position) {

    // TODO("Variants")
    // TODO("Builder?")

    override fun compareTo(other: GenomicFactor): Int {

        val result = super.compareTo(other)
        if (other !is SNP || result != 0) return result
        return snpName.compareTo(other.snpName)

    }

}