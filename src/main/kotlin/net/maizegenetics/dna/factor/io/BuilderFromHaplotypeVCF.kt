package net.maizegenetics.dna.factor.io

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.dna.factor.FeatureTableBuilder
import net.maizegenetics.dna.factor.site.HaplotypeSite
import net.maizegenetics.dna.factor.site.HaplotypeSiteBuilder
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.dna.map.GenomicFeature
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.taxa.Taxon
import java.io.File

class BuilderFromHaplotypeVCF {

    private val processingChannel = Channel<HaplotypeSite>(5)
    private lateinit var taxa: TaxaList

    fun read(filename: String): FeatureTable {

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

        reader.forEach { context ->
            processingChannel.send(contextToSite(context))
        }

        processingChannel.close()

    }

    private fun contextToSite(context: VariantContext): HaplotypeSite {
        val factor = GenomicFeature(Chromosome.instance(context.contig), context.start, endPos = context.end)
        val strStates = context.alleles
                .map { it.displayString }
                .toTypedArray()
        val builder = HaplotypeSiteBuilder(factor, taxa, context.getMaxPloidy(2), strStates)
        try {
            context.genotypes.forEach { genotype ->
                builder.set(genotype.sampleName, genotype.alleles.map { it.displayString })
            }
        } catch (e: Exception) {
            throw e
        }
        return builder.build()
    }

    private suspend fun addSitesToTable(): FeatureTable {
        val builder = FeatureTableBuilder(taxa)
        for (site in processingChannel) {
            builder.add(site)
        }
        return builder.build()
    }

}

fun main() {
    val filename = "/Users/tmc46/git/tassel-5-standalone/small_seq_haplotypes.vcf"
    val factorTable = BuilderFromHaplotypeVCF().read(filename)
    println(factorTable.numFeatures())
}