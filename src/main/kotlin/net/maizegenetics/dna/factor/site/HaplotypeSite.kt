package net.maizegenetics.dna.factor.site

import net.maizegenetics.dna.map.SNP
import net.maizegenetics.taxa.TaxaList

/**
 * @author Terry Casstevens
 * Created November 16, 2018
 */

class HaplotypeSite(factor: SNP, index: Int, taxa: TaxaList, private val values: ByteArray) : FactorSite(factor, index, taxa) {

}