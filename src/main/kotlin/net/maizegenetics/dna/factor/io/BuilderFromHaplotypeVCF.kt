package net.maizegenetics.dna.factor.io

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFAltHeaderLine
import htsjdk.variant.vcf.VCFFileReader
import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.dna.factor.FeatureTableBuilder
import net.maizegenetics.dna.factor.site.HaplotypeAnnotation
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

            val altHeaderLines = header.idHeaderLines
                .filter { it is VCFAltHeaderLine }
                .map { Pair(it.id, it.toString().substringAfter("Description=\"").substringBefore("\"")) }
                .toMap()

            val samples = header.sampleNamesInOrder

            val taxaListBuilder = TaxaListBuilder()
            samples
                .map { Taxon.instance(it) }
                .forEach { taxaListBuilder.add(it) }
            taxa = taxaListBuilder.build()

            CoroutineScope(Dispatchers.IO).launch {
                processPositions(reader, altHeaderLines)
            }

            return runBlocking { addSitesToTable() }

        }

    }

    private suspend fun processPositions(reader: VCFFileReader, altHeaderLines: Map<String, String>) =
        withContext(Dispatchers.IO) {

            reader.forEach { context ->
                processingChannel.send(contextToSite(context, altHeaderLines))
            }

            processingChannel.close()

        }

    private fun contextToSite(context: VariantContext, altHeaderLines: Map<String, String>): HaplotypeSite {

        val factor = GenomicFeature(Chromosome.instance(context.contig), context.start, endPos = context.end)
        val strStates = context.alleles
            .map { it.displayString.substringAfter("<").substringBefore(">") }
            .toTypedArray()

        val haplotypeAnnotations = strStates
            .map { altHeaderLines[it] }
            .map {
                it?.let {
                    val taxon = it.substringBefore(":").substringBefore(",")
                    val range = it.substringAfter(":")
                    val asmContig = range.substringBefore(":")
                    val positions = range.substringAfterLast(":")
                    val asmStart = positions.substringBefore("-").toInt()
                    val asmEnd = positions.substringAfter("-").toInt()
                    HaplotypeAnnotation(taxon, asmContig, asmStart, asmEnd)
                }
            }
            .toTypedArray()

        val builder = HaplotypeSiteBuilder(factor, taxa, context.getMaxPloidy(2), strStates, haplotypeAnnotations)
        try {
            context.genotypes.forEach { genotype ->
                builder.set(
                    genotype.sampleName,
                    genotype.alleles.map { it.displayString.substringAfter("<").substringBefore(">") })
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
    println((factorTable.get(0) as HaplotypeSite).haplotypeAnnotation(0x2))
    println((factorTable.get(0) as HaplotypeSite).haplotypeAnnotation(0x2)?.asmContig)
}