package net.maizegenetics.dna.factor.io

import com.google.common.collect.SetMultimap
import net.maizegenetics.dna.factor.FactorTable
import net.maizegenetics.dna.factor.site.SNPSiteBuilder
import net.maizegenetics.dna.map.Chromosome
import net.maizegenetics.dna.map.GenomicFactor
import net.maizegenetics.dna.map.GenomicFactorList
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants
import net.maizegenetics.taxa.TaxaListBuilder
import net.maizegenetics.taxa.TaxaListIOUtils
import net.maizegenetics.taxa.Taxon
import net.maizegenetics.util.ProgressListener
import net.maizegenetics.util.SuperByteMatrix
import net.maizegenetics.util.SuperByteMatrixBuilder
import net.maizegenetics.util.Utils
import org.apache.log4j.Logger
import java.util.*
import java.util.concurrent.*
import java.util.regex.Pattern

/*
 *  BuilderFromHapMap
 */

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
class BuilderFromHapMap private constructor(private val myHapmapFile: String, private val myProgressListener: ProgressListener?) {

    private var mySortPositions = false

    fun build(): FactorTable {

        var pool: ExecutorService? = null
        try {
            Utils.getBufferedReader(myHapmapFile, 1 shl 20)!!.use { reader ->

                val sampAnnoBuild = TreeMap<String, SetMultimap<String, String>>()

                var currLine: String? = reader.readLine()
                while (currLine != null && currLine.startsWith("##")) {

                    val cat = currLine.split("=".toRegex(), 2).toTypedArray()
                    if (cat.size < 2) {
                        continue
                    }
                    if (cat[0].startsWith("##SAMPLE")) {

                        val mapOfAnno = TaxaListIOUtils.parseVCFHeadersIntoMap(cat[1])
                        val taxaID = mapOfAnno.get("ID").iterator().next()
                        if (taxaID != null) {
                            sampAnnoBuild[taxaID] = mapOfAnno
                        }

                    }

                    currLine = reader.readLine()
                }

                val taxaList = processTaxa(currLine, sampAnnoBuild)
                val taxa = taxaList.build()
                val numTaxa = taxa.numberOfTaxa()

                val chromosomeLookup = ConcurrentHashMap<String, Chromosome>()

                currLine = reader.readLine()

                var isOneLetter = false
                val tokens = WHITESPACE_PATTERN.split(currLine, NUM_HAPMAP_NON_TAXA_HEADERS + 1)
                if (tokens.size <= NUM_HAPMAP_NON_TAXA_HEADERS) {
                    throw IllegalStateException("BuilderFromHapMap: Header Incorrectly Formatted: See:\nhttps://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load#markdown-header-hapmap")
                }
                val avg = (tokens[NUM_HAPMAP_NON_TAXA_HEADERS].length + 1).toDouble() / numTaxa.toDouble()
                if (avg > 1.99 && avg < 2.01) {
                    isOneLetter = true
                } else if (avg > 2.99 && avg < 3.01) {
                    isOneLetter = false
                } else {
                    throw IllegalStateException("BuilderFromHapMap: Genotype coded wrong use 1 or 2 letters per genotype. Average chars including tab: $avg  Or first site has incorrect number of values. Number of taxa: $numTaxa")
                }

                val numThreads = Runtime.getRuntime().availableProcessors()
                pool = Executors.newFixedThreadPool(numThreads)
                val futures = ArrayList<Future<ProcessHapmapBlock>>()

                var numSitesToProcessTogether = NUM_VALUES_PROCESSED_TOGETHER / numTaxa
                numSitesToProcessTogether = Math.min(1 shl 16, numSitesToProcessTogether)
                numSitesToProcessTogether = Math.max(512, numSitesToProcessTogether)

                var textLines = ArrayList<String>(numSitesToProcessTogether)
                var numLines = 0
                while (currLine != null) {
                    textLines.add(currLine)
                    numLines++
                    if (numLines % numSitesToProcessTogether == 0) {
                        val processBlock = ProcessHapmapBlock(textLines, numTaxa, chromosomeLookup, isOneLetter)
                        futures.add(pool!!.submit(processBlock))
                        textLines = ArrayList(numSitesToProcessTogether)
                    }
                    currLine = reader.readLine()
                }

                if (textLines.size > 0) {
                    val processBlock = ProcessHapmapBlock(textLines, numTaxa, chromosomeLookup, isOneLetter)
                    futures.add(pool!!.submit(processBlock))
                }

                var currentSite = 0
                val positions = GenomicFactorList.Builder()
                val genotypes = SNPSiteBuilder.instance(numTaxa, numLines)

                val numFutures = futures.size
                var count = 0
                for (future in futures) {
                    val pb = future.get()
                    positions.addAll(pb.positions)
                    val bgTS = pb.genotypes
                    for (t in 0 until bgTS!!.numRows) {
                        for (s in 0 until bgTS.numColumns) {
                            genotypes.set(t, currentSite + s, bgTS.get(t, s))
                        }
                    }
                    currentSite += pb.numberSitesProcessed
                    if (myProgressListener != null) {
                        count++
                        myProgressListener.progress(count * 100 / numFutures, null)
                    }
                }
                pool!!.shutdown()

                //if (mySortPositions) {
                //    positions.sortPositions(genotypes)
                //}

                //if (positions.validateOrdering() == false) {
                //    throw IllegalStateException("BuilderFromHapMap: Ordering incorrect. HapMap must be ordered by position. Please first use SortGenotypeFilePlugin to correctly order the file.")
                //}

                genotypes.taxa = taxa
                genotypes.factors = positions.build()
                return genotypes.build()

            }
        } catch (e: Exception) {
            if (pool != null) {
                pool!!.shutdown()
            }
            myLogger.debug(e.message, e)
            throw IllegalStateException(e.message)
        }

    }

    private inner class ProcessHapmapBlock(private var myInputLines: List<String>?, private val myNumTaxa: Int, private val myChromosomeLookup: MutableMap<String, Chromosome>, private val myIsOneLetter: Boolean) : Callable<ProcessHapmapBlock> {
        private val myPositionList: MutableList<GenomicFactor>
        val numberSitesProcessed: Int
        var genotypes: SuperByteMatrix? = null
            private set

        val positions: List<GenomicFactor>
            get() = myPositionList

        init {
            numberSitesProcessed = myInputLines!!.size
            myPositionList = ArrayList(numberSitesProcessed)
        }

        @Throws(Exception::class)
        override fun call(): ProcessHapmapBlock {

            genotypes = SuperByteMatrixBuilder.getInstance(myNumTaxa, numberSitesProcessed)
            for (site in 0 until numberSitesProcessed) {
                val input = myInputLines!![site]
                try {
                    val tabPos = IntArray(NUM_HAPMAP_NON_TAXA_HEADERS)
                    var tabIndex = 0
                    val len = input.length
                    run {
                        var i = 0
                        while (tabIndex < NUM_HAPMAP_NON_TAXA_HEADERS && i < len) {
                            if (input[i] == '\t') {
                                tabPos[tabIndex++] = i
                            }
                            i++
                        }
                    }
                    val chrName = input.substring(tabPos[CHROMOSOME_INDEX - 1] + 1, tabPos[CHROMOSOME_INDEX])

                    val tempChr = myChromosomeLookup[chrName]
                    val currChr = when (tempChr) {
                        null -> {
                            Chromosome.instance(chrName).also { myChromosomeLookup[chrName] = it }
                        }
                        else -> tempChr
                    }

                    val variants = input.substring(tabPos[VARIANT_INDEX - 1] + 1, tabPos[VARIANT_INDEX])
                    val physicalPos: Int
                    try {
                        physicalPos = Integer.parseInt(input.substring(tabPos[POSITION_INDEX - 1] + 1, tabPos[POSITION_INDEX]))
                    } catch (ex: Exception) {
                        throw IllegalArgumentException("BuilderFromHapMap: Position must be an integer: " + input.substring(tabPos[POSITION_INDEX - 1] + 1, tabPos[POSITION_INDEX]).trim { it <= ' ' })
                    }

                    val apb = GenomicFactor(input.substring(0, tabPos[SNPID_INDEX]), currChr, physicalPos)
                    // TODO("knownVariants(variants) //TODO   strand, variants,")

                    //val glbMajor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants[0])
                    //apb.allele(WHICH_ALLELE.GlobalMajor, glbMajor)
                    //if (variants.length == 3) {
                    //    val glbMinor = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants[2])
                    //    apb.allele(WHICH_ALLELE.GlobalMinor, glbMinor)
                    //}

                    myPositionList.add(apb)
                    val offset = tabPos[NUM_HAPMAP_NON_TAXA_HEADERS - 1] + 1

                    var taxon = 0
                    if (myIsOneLetter) {
                        var i = offset
                        while (i < len) {
                            if (taxon >= myNumTaxa) {
                                throw IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList[myPositionList.size - 1].factorName + " has too many values.")
                            }
                            val value = NucleotideAlignmentConstants.getNucleotideDiploidByte(input[i])
                            if (value == NucleotideAlignmentConstants.UNDEFINED_HOMOZYGOUS) {
                                throw IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList[myPositionList.size - 1].factorName + " has illegal value: " + input[i])
                            }
                            genotypes!!.set(taxon++, site, value)
                            i += 2
                        }
                    } else {
                        var i = offset
                        while (i < len) {
                            if (taxon >= myNumTaxa) {
                                throw IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList[myPositionList.size - 1].factorName + " has too many values.")
                            }
                            // there is a phasing conflict with the existing import approach
                            val value = GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.getNucleotideDiploidByte(input[i + 1]),
                                    NucleotideAlignmentConstants.getNucleotideDiploidByte(input[i]))
                            if (value == NucleotideAlignmentConstants.UNDEFINED_HOMOZYGOUS) {
                                throw IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList[myPositionList.size - 1].factorName + " has illegal value: " + input[i] + input[i + 1])
                            }
                            genotypes!!.set(taxon++, site, value)
                            i += 3
                        }
                    }
                    if (taxon != myNumTaxa) {
                        throw IllegalStateException("BuilderFromHapMap: SNP Named: " + myPositionList[myPositionList.size - 1].factorName + " has too few values.")
                    }

                    swapSitesIfOutOfOrder(site)

                } catch (e: Exception) {
                    myLogger.error("Error parsing this row $input")
                    myLogger.debug(e.message, e)
                    throw IllegalStateException("BuilderFromHapMap: Error Parsing Line: " + input.substring(0, Math.min(25, input.length)) + "...\n" + e.message)
                }

            }
            myInputLines = null

            return this

        }

        // Swap adjacent misordered sites, often caused by two sites at the same positions with a different name order
        private fun swapSitesIfOutOfOrder(site: Int) {
            if (site < 1) {
                return
            }
            if (myPositionList[site - 1].compareTo(myPositionList[site]) > 0) {
                //swap
                val tempP = myPositionList[site - 1]
                myLogger.warn("Swapping:" + tempP.toString() + " <-> " + myPositionList[site].toString())
                myPositionList[site - 1] = myPositionList[site]
                myPositionList[site] = tempP
                for (t in 0 until genotypes!!.numRows) {
                    val tempG = genotypes!!.get(t, site - 1)
                    genotypes!!.set(t, site - 1, genotypes!!.get(t, site))
                    genotypes!!.set(t, site, tempG)
                }
            }

        }

    }

    /**
     * Set the builder so that when built it will sort positions.
     */
    fun sortPositions(): BuilderFromHapMap {
        mySortPositions = true
        return this
    }

    companion object {

        private val myLogger = Logger.getLogger(BuilderFromHapMap::class.java)
        private val WHITESPACE_PATTERN = Pattern.compile("\\t")
        private val NUM_HAPMAP_NON_TAXA_HEADERS = 11
        private val SNPID_INDEX = 0
        private val VARIANT_INDEX = 1
        private val CHROMOSOME_INDEX = 2
        private val POSITION_INDEX = 3
        private val NUM_VALUES_PROCESSED_TOGETHER = 7 shl 20

        @JvmStatic
        fun getBuilder(hapmapFile: String): BuilderFromHapMap {
            return BuilderFromHapMap(hapmapFile, null)
        }

        @JvmStatic
        fun getBuilder(hapmapFile: String, listener: ProgressListener): BuilderFromHapMap {
            return BuilderFromHapMap(hapmapFile, listener)
        }

        internal fun processTaxa(readLn: String?, taxaAnnotation: Map<String, SetMultimap<String, String>>): TaxaListBuilder {

            val header = WHITESPACE_PATTERN.split(readLn)
            val numTaxa = header.size - NUM_HAPMAP_NON_TAXA_HEADERS
            val tlb = TaxaListBuilder()
            for (i in 0 until numTaxa) {
                val taxonID = header[i + NUM_HAPMAP_NON_TAXA_HEADERS]
                if (taxonID == null || taxonID.isEmpty()) {
                    throw IllegalStateException("BuilderFromHapMap: processTaxa: Taxa names should be separated by a single tab and contain no spaces.")
                }
                val at = Taxon.Builder(taxonID)
                val taMap = taxaAnnotation[taxonID]
                if (taMap != null) {
                    for ((key, value) in taMap.entries()) {
                        if (key == "ID") {
                            continue //skip the IDs as these became the name
                        }
                        val s = value.replace("\"", "")
                        at.addAnno(key, s)
                    }
                }
                tlb.add(at.build())
            }
            return tlb

        }
    }
}
