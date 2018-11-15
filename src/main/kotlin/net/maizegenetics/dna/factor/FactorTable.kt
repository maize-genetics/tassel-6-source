package net.maizegenetics.dna.factor

import net.maizegenetics.dna.factor.factorvalues.FactorSite
import net.maizegenetics.dna.factor.factorvalues.FactorValues
import net.maizegenetics.dna.map.GenomicFactorList
import net.maizegenetics.taxa.TaxaList
import java.util.stream.Stream

/**
 * @author Terry Casstevens
 * Created November 13, 2018
 */

class FactorTable(val taxa: TaxaList, val factors: GenomicFactorList, val values: FactorValues) {

    fun sites(): Stream<FactorSite> {
        TODO()
    }

}