package net.maizegenetics.dna.factor.io

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.dna.factor.FactorTable
import net.maizegenetics.dna.factor.FactorTableBuilder
import net.maizegenetics.dna.factor.site.HaplotypeSite
import net.maizegenetics.dna.factor.site.HaplotypeSiteBuilder
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.taxa.Taxon
import java.io.File

class BuilderFromHaplotypeVCF {

    private val processingChannel = Channel<HaplotypeSite>(5)
    private lateinit var taxa: TaxaList

    fun read(filename: String): FactorTable {

        VCFFileReader(File(filename), false).use { reader ->

            val header = reader.fileHeader
            val samples = header.sampleNamesInOrder
            val taxaListBuilder = TaxaListBuilder()
            samples
                    .map { Taxon(it) }
                    .forEach { taxaListBuilder.add(it) }
            taxa = taxaListBuilder.build()

            CoroutineScope(Dispatchers.IO).launch {
                processPositions(reader)
            }

            return runBlocking { addSitesToTable() }

        }

    }

    private suspend fun processPositions(reader: VCFFileReader) = withContext(Dispatchers.IO) {

        var index = 0
        reader.forEach { context ->
            processingChannel.send(contextToSite(context, index++))
        }

        processingChannel.close()

    }

    private fun contextToSite(context: VariantContext, index: Int): HaplotypeSite {
        println("${index}_1 $context")
        val factor = GenomicFactor(Chromosome.instance(context.contig), context.start, endPos = context.end)
        val strStates = context.alleles
                .map { it.displayString }
                .toTypedArray()
        val builder = HaplotypeSiteBuilder(factor, taxa, context.getMaxPloidy(2), strStates)
        println("${index}_2 $context")
        try {
            context.genotypes.forEach { genotype ->
                builder.set(genotype.sampleName, genotype.alleles.map { it.displayString })
                println(genotype.alleles)
            }
        } catch (e: Exception) {
            println("${index}_3 $context")
            throw e
        }
        println("${index}_4 $context")
        return builder.build()
    }

    private suspend fun addSitesToTable(): FactorTable {
        val builder = FactorTableBuilder(taxa)
        for (site in processingChannel) {
            builder.add(site)
        }
        return builder.build()
    }

}

fun main() {
    val filename = "/Users/tmc46/git/tassel-5-standalone/small_seq_haplotypes.vcf"
    val factorTable = BuilderFromHaplotypeVCF().read(filename)
    println(factorTable.numFactors())
}