package net.maizegenetics.analysis.distance

import kotlinx.coroutines.*
import kotlinx.coroutines.channels.Channel
import net.maizegenetics.dna.factor.FeatureTable
import net.maizegenetics.dna.factor.UNKNOWN_ALLELE
import net.maizegenetics.dna.factor.site.FeatureSite
import net.maizegenetics.taxa.distance.DistanceMatrix
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder
import net.maizegenetics.util.GeneralAnnotationStorage
import net.maizegenetics.util.ProgressListener
import org.apache.log4j.Logger
import java.util.*
import kotlin.math.roundToLong
import kotlin.system.measureNanoTime

class EndelmanDistanceMatrixBuilder(val table: FeatureTable, val maxAlleles: Int = 255, private val listener: ProgressListener? = null) {

    private val logger = Logger.getLogger(EndelmanDistanceMatrixBuilder::class.java)

    private val psuedoSiteChannel = Channel<List<List<PsuedoSite>>>(1000)

    private val resultsChannel = Channel<CountersDistances>(30)

    private var numProcessingThreads: Int = 1

    data class PsuedoSite(val site: FeatureSite, val allele: Byte, val alleleFreq: Float)

    fun build(): DistanceMatrix {

        val time = measureNanoTime {

            logger.debug("EndelmanDistanceMatrixBuilder: factor table num taxa: ${table.numTaxa()}  num factors: ${table.numFeatures()}")
            numProcessingThreads = (Runtime.getRuntime().availableProcessors() - 2).coerceAtLeast(1)
            logger.debug("EndelmanDistanceMatrixBuilder: numProcessingThreads: $numProcessingThreads")

            CoroutineScope(Dispatchers.IO).launch { createPsuedoSites() }

            runBlocking {
                val jobs = List(numProcessingThreads) {
                    launch(Dispatchers.Default) {
                        processPsuedoSites()
                    }
                }
                jobs.forEach { it.join() }

                resultsChannel.close()
            }

        }

        val estimatedNumMinutesToRun = (time.toDouble() / 1e9 / 60.0).roundToLong()
        if (estimatedNumMinutesToRun < 60L) {
            logger.info("EndelmanDistanceMatrixBuilder: actual time to complete: $estimatedNumMinutesToRun minutes")
        } else {
            logger.info("EndelmanDistanceMatrixBuilder: actual time to complete: ${estimatedNumMinutesToRun / 60} hours ${estimatedNumMinutesToRun % 60} minutes")
        }

        return runBlocking { accumulateResults() }

    }

    private val numPsuedoSitesPerBlock = 15
    private val numBlocksPerChunk = 200
    private var aveAllelesPerSite = 0.0

    private suspend fun createPsuedoSites() = withContext(Dispatchers.IO) {

        var numSitesProcessed = 0
        var totalNumAllelesToEvaluate = 0

        table
                .map { site ->
                    numSitesProcessed++
                    val siteStats = site.alleleStats
                    totalNumAllelesToEvaluate += siteStats.numAlleles - 1
                    aveAllelesPerSite = totalNumAllelesToEvaluate.toDouble() / numSitesProcessed.toDouble()
                    siteStats.alleleCounts
                            .dropLast(1)
                            .map { alleleCount ->
                                PsuedoSite(site, alleleCount.allele, alleleCount.count.toFloat() / siteStats.totalNonMissingAlleles.toFloat())
                            }
                }
                .flatten()
                .chunked(numPsuedoSitesPerBlock)
                .chunked(numBlocksPerChunk)
                .forEach { psuedoSites ->
                    psuedoSiteChannel.send(psuedoSites)
                }

        logger.debug("Number Factors Processed: $numSitesProcessed")
        logger.debug("Total Number of Psuedo Sites: $totalNumAllelesToEvaluate")
        logger.debug("Average Alleles Evaluation Per Site: $aveAllelesPerSite")

        psuedoSiteChannel.close()

        val estimatedNumMinutesToRun: Long = (table.numTaxa() * (table.numTaxa() + 1.0) / 2.0 * totalNumAllelesToEvaluate.toDouble() / numProcessingThreads.toDouble() / 1.02e11).roundToLong()
        if (estimatedNumMinutesToRun < 60L) {
            logger.info("EndelmanDistanceMatrixBuilder: estimated time to complete: $estimatedNumMinutesToRun minutes")
        } else {
            logger.info("EndelmanDistanceMatrixBuilder: estimated time to complete: ${estimatedNumMinutesToRun / 60} hours ${estimatedNumMinutesToRun % 60} minutes")
        }

    }

    private var numPsuedoSitesProcessed = 0

    private suspend fun processPsuedoSites() {

        val result = CountersDistances(table.numTaxa())
        val distances: FloatArray = result.distances
        val sumpi = DoubleArray(1)

        val answer1 = FloatArray(32768)
        val answer2 = FloatArray(32768)
        val answer3 = FloatArray(32768)

        for (psuedoSitesBlock in psuedoSiteChannel) {

            for (psuedoSites in psuedoSitesBlock) {

                //
                // Pre-calculates possible terms and gets counts for
                // three blocks for five (pseudo-)sites.
                //
                val blocksOfSites = getBlocksOfSites(psuedoSites, sumpi, table.numTaxa())

                val possibleTerms = blocksOfSites.second[0]
                val alleleCount1 = blocksOfSites.first[0]

                val possibleTerms2 = blocksOfSites.second[1]
                val alleleCount2 = blocksOfSites.first[1]

                val possibleTerms3 = blocksOfSites.second[2]
                val alleleCount3 = blocksOfSites.first[2]

                //
                // Using possible terms, calculates all possible answers
                // for each site block.
                //
                for (i in 0..32767) {
                    answer1[i] = possibleTerms[i and 0x7000 ushr 12] + possibleTerms[i and 0xE00 ushr 9 or 0x8] + possibleTerms[i and 0x1C0 ushr 6 or 0x10] + possibleTerms[i and 0x38 ushr 3 or 0x18] + possibleTerms[i and 0x7 or 0x20]
                    answer2[i] = possibleTerms2[i and 0x7000 ushr 12] + possibleTerms2[i and 0xE00 ushr 9 or 0x8] + possibleTerms2[i and 0x1C0 ushr 6 or 0x10] + possibleTerms2[i and 0x38 ushr 3 or 0x18] + possibleTerms2[i and 0x7 or 0x20]
                    answer3[i] = possibleTerms3[i and 0x7000 ushr 12] + possibleTerms3[i and 0xE00 ushr 9 or 0x8] + possibleTerms3[i and 0x1C0 ushr 6 or 0x10] + possibleTerms3[i and 0x38 ushr 3 or 0x18] + possibleTerms3[i and 0x7 or 0x20]
                }

                //
                // Iterates through all pair-wise combinations of taxa adding
                // distance comparisons and site counts.
                //
                var index = 0
                for (firstTaxa in 0 until table.numTaxa()) {
                    //
                    // Can skip inter-loop if all fifteen sites for first
                    // taxon is Unknown diploid allele values
                    //
                    if (alleleCount1[firstTaxa] != 0x7FFF.toShort() || alleleCount2[firstTaxa] != 0x7FFF.toShort() || alleleCount3[firstTaxa] != 0x7FFF.toShort()) {
                        for (secondTaxa in firstTaxa until table.numTaxa()) {
                            //
                            // Combine first taxon's allele counts with
                            // second taxon's major allele counts to
                            // create index into pre-calculated answers
                            //
                            distances[index] += answer1[(alleleCount1[firstTaxa].toInt() or alleleCount1[secondTaxa].toInt()) and 0xFFFF] + answer2[(alleleCount2[firstTaxa].toInt() or alleleCount2[secondTaxa].toInt()) and 0xFFFF] + answer3[(alleleCount3[firstTaxa].toInt() or alleleCount3[secondTaxa].toInt()) and 0xFFFF]
                            index++
                        }
                    } else {
                        index += table.numTaxa() - firstTaxa
                    }
                }
            }

            numPsuedoSitesProcessed += numPsuedoSitesPerBlock * numBlocksPerChunk
            val percent = (numPsuedoSitesProcessed.toDouble() / aveAllelesPerSite / table.numFeatures().toDouble() * 100.0).toInt()
            fireProgress(percent, listener)

        }

        result.sumPi = sumpi[0]

        resultsChannel.send(result)

    }

    private suspend fun accumulateResults(): DistanceMatrix {

        val result = resultsChannel.receive()

        for (intermediateResult in resultsChannel) {
            result.addAll(intermediateResult)
        }

        var sumpk = result.sumPi
        val distances = result.distances

        //
        // This does the final division of the frequency sum into
        // the distance sums.
        //
        sumpk *= 2.0

        val annotations = GeneralAnnotationStorage.getBuilder()
        annotations.addAnnotation(DistanceMatrixBuilder.MATRIX_TYPE, KinshipPlugin.KINSHIP_METHOD.Centered_IBS.toString())
        annotations.addAnnotation(DistanceMatrixBuilder.CENTERED_IBS_SUMPK, sumpk)

        val builder: DistanceMatrixBuilder = DistanceMatrixBuilder.getInstance(table.taxa)
        builder.annotation(annotations.build())
        var index = 0
        for (t in 0 until table.numTaxa()) {
            var i = 0
            val n: Int = table.numTaxa() - t
            while (i < n) {
                builder[t, t + i] = distances[index] / sumpk
                index++
                i++
            }
        }

        return builder.build()

    }

    private fun getBlocksOfSites(psuedoSites: List<PsuedoSite>, sumpi: DoubleArray, numTaxa: Int): Pair<Array<ShortArray>, Array<FloatArray>> {

        val numBlocks = 3
        val numSitesPerBlock = 5
        var currentBlock = 0
        var currentSiteNum = 0

        //
        // This hold possible terms for the Endelman summation given
        // site's allele frequency.  First three bits
        // identifies relative site (0, 1, 2, 3, 4).  Remaining three bits
        // the allele counts encoding.
        //
        val possibleTerms = Array(numBlocks) { FloatArray(40) }

        //
        // This holds count of allele for each taxa.
        // Each short holds count (0, 1, 2, 3) for all four sites
        // at given taxon.  The count encodings are stored in three
        // bits each.
        //
        val alleleCount = Array(numBlocks) { ShortArray(numTaxa) }

        //
        // This initializes the counts to 0x7FFF.  That means
        // diploid allele values for the five pseudo-sites are Unknown.
        //
        for (i in 0 until numBlocks) {
            Arrays.fill(alleleCount[i], 0x7FFF.toShort())
        }

        psuedoSites
                .forEach { psuedoSite ->

                    val allele = psuedoSite.allele
                    val alleleFreq = psuedoSite.alleleFreq
                    val alleleFreqTimes2 = alleleFreq * 2.0f
                    sumpi[0] += alleleFreq * (1.0 - alleleFreq)

                    val term0 = 0.0f - alleleFreqTimes2
                    val term1 = 1.0f - alleleFreqTimes2
                    val term2 = 2.0f - alleleFreqTimes2

                    //
                    // Pre-calculates all possible terms of the summation
                    // for this current (pseudo-) site.
                    // Counts (0,0; 0,1; 0,2; 1,1; 1,2; 2,2)
                    //
                    val siteNumIncrement = currentSiteNum * 8
                    possibleTerms[currentBlock][siteNumIncrement + 1] = term0 * term0
                    possibleTerms[currentBlock][siteNumIncrement + 3] = term0 * term1
                    possibleTerms[currentBlock][siteNumIncrement + 5] = term0 * term2
                    possibleTerms[currentBlock][siteNumIncrement + 2] = term1 * term1
                    possibleTerms[currentBlock][siteNumIncrement + 6] = term1 * term2
                    possibleTerms[currentBlock][siteNumIncrement + 4] = term2 * term2

                    //
                    // Records allele counts for current site in
                    // three bits.
                    //
                    val shift = (numSitesPerBlock - currentSiteNum - 1) * 3
                    val mask = (0x7 shl shift).inv() and 0x7FFF
                    for (i in 0 until numTaxa) {
                        val taxonAlleles = psuedoSite.site.genotype(i)
                        alleleCount[currentBlock][i] = (alleleCount[currentBlock][i].toInt() and (mask or (calculateCount(allele, taxonAlleles[0], taxonAlleles[1]) shl shift))).toShort()
                    }

                    currentSiteNum++
                    if (currentSiteNum == numSitesPerBlock) {
                        currentSiteNum = 0
                        currentBlock++
                    }

                }

        return Pair(alleleCount, possibleTerms)

    }

    private fun calculateCount(allele: Byte, value1: Byte, value2: Byte): Int {

        if (allele == UNKNOWN_ALLELE || (value1 == UNKNOWN_ALLELE && value2 == UNKNOWN_ALLELE)) return 7
        var result = 0
        if (value1 == allele) {
            result = 2
        }
        if (value2 == allele) {
            result += 2
        }
        if (result == 0) {
            result = 1
        }
        return result

    }

    private fun fireProgress(percent: Int, listener: ProgressListener?) {
        listener?.progress(if (percent > 100) 100 else percent, null)
    }

    //
    // Each CPU thread (process) creates an instance of this class
    // to acculate terms of the Endelman equation. These are
    // combined with addAll() to result in one instance at the end.
    //
    private class CountersDistances(numTaxa: Int) {

        var sumPi = 0.0
        val distances: FloatArray = FloatArray(numTaxa * (numTaxa + 1) / 2)

        fun addAll(counters: CountersDistances) {
            val otherDistances = counters.distances
            var t = 0
            val n = distances.size
            while (t < n) {
                distances[t] += otherDistances[t]
                t++
            }
            sumPi += counters.sumPi
        }

    }

}
