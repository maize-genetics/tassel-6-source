/**
 * 
 */
package net.maizegenetics.dna.map;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.LongAdder;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;

import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;
import net.maizegenetics.util.Tuple;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;


/**
 * Builder for a chromosome genome sequence that hides the ReferenceGenomeSequence implementation.
 * The sequence is read from a fasta file and stored in a byte array, 2 alleles packed per byte.
 * 
 * @author Lynn Johnson
 *
 */
public class GenomeSequenceBuilder {
    private static final Logger myLogger = Logger.getLogger(GenomeSequenceBuilder.class);

    /**
     * Builds GenomeSequence from a fasta file.
     * @param fastaFileName full path to fasta file
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(String fastaFileName) {
        Function<Character, Character> charConversion= (c) -> c;
        return instance(fastaFileName,charConversion);
    }

    /**
     * Builds GenomeSequence from a fasta file.  The char conversion provide a mechanism to convert upper and lower case
     * or convert one case to N.  This is useful if a case if used to define a certain class of bases
     * @param fastaFileName full path to fasta file
     * @param charConversion lambda Function to convert characters
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(String fastaFileName, Function<Character, Character> charConversion) {
        Map<Chromosome, byte[]> chromPositionMap = readReferenceGenomeChr(fastaFileName, charConversion);
        return new HalfByteGenomeSequence(chromPositionMap);
    }

    /**
     * Builds GenomeSequence from a String with one Chromosome.
     * @param chromosome chromosome object
     * @param sequence in upper case character
     * @return GenomeSequence object
     */
    public static GenomeSequence instance(Chromosome chromosome, String sequence) {
        Function<Character, Character> charConversion= (c) -> c;
        Map<Chromosome, byte[]> chromPositionMap = new HashMap<>();
        Chromosome currChr=new Chromosome(chromosome.getName(),sequence.length(),chromosome.getAnnotation());
        chromPositionMap.put(currChr,halfByteCompression(sequence.getBytes(),charConversion));
        return new HalfByteGenomeSequence(chromPositionMap);
    }

    public static void main(String[] args) {
        Chromosome chr = new Chromosome("3");
        String sequence = "ACGT";
        GenomeSequence temp = instance(chr, sequence);
        System.out.println(temp.genomeSequenceAsString(0, sequence.length()));
    }

    protected static  Map<Chromosome, byte[]> readReferenceGenomeChr(String fastaFileName, Function<Character, Character> charConversion) {
        // Read specified file, return entire sequence for requested chromosome
        String base="ACGTNacgtn";
        String conv=base.chars().mapToObj(ci -> charConversion.apply((char)ci).toString()).collect(Collectors.joining());
        System.out.println("Genome FASTA character conversion: "+base+" to "+conv);
        Map<Chromosome, byte[]> chromPositionMap = new HashMap<Chromosome, byte[]>();
        Chromosome currChr = null;	
        ByteArrayOutputStream currSeq = new ByteArrayOutputStream();
        String line = null;
        try {
            boolean found = false;
            BufferedReader br = Utils.getBufferedReader(fastaFileName); // this takes care of .gz

            while ((line = br.readLine()) != null && !found) {
                line = line.trim();

                if (line.startsWith(">")) {
                    if (currChr != null) {
                        // end processing current chromosome sequence
                        currChr=new Chromosome(currChr.getName(),currSeq.size(),currChr.getAnnotation());
                        chromPositionMap.put(currChr, halfByteCompression(currSeq.toByteArray(),charConversion));
                    }
                    currChr = parseChromosome(line); 
                    currSeq = new ByteArrayOutputStream();
                } else {
                    currSeq.write(line.getBytes());
                }
            }
            // reached end of file - write last bytes
            if (currSeq.size() > 0) {
                currChr=new Chromosome(currChr.getName(),currSeq.size(),currChr.getAnnotation());
                chromPositionMap.put(currChr, halfByteCompression(currSeq.toByteArray(),charConversion));
            }
            br.close();
        } catch (IOException ioe) {
            System.out.println("ReferenceGenomeSequence: caught buffered read exception");
        }
        return chromPositionMap;
    }

    private static byte[] halfByteCompression(byte[] unpkSequence, Function<Character, Character> charConversion){
        // Take byte array, turn bytes into NucleotideAlignmentConstant
        // allele values, store as half bytes
        int nBytes = (unpkSequence.length+1)/2;
        byte[] packedSequence = new byte[nBytes];
        for (int i = 0; i < unpkSequence.length; i++) {
            byte halfByte = NucleotideAlignmentConstants.getNucleotideAlleleByte(charConversion.apply((char)unpkSequence[i]));
            if(i%2==0) halfByte<<=4;
            packedSequence[i/2]|=halfByte;
        }
        return packedSequence;
    }

    private static  Chromosome parseChromosome (String chromString) {
        String chrS = chromString.replace(">","");
        chrS = chrS.toUpperCase();
        chrS = chrS.replace("CHROMOSOME", ""); 
        chrS = chrS.replace("CHR", ""); // keep chromosome string, minus any leading "chr" or "chromosome"
        GeneralAnnotation myAnnotations = null;

        String currChrDesc = null;
        int spaceIndex = chrS.indexOf(" ");
        if (spaceIndex > 0) {			
            currChrDesc = chrS.substring(chrS.indexOf(" ") + 1);
            myAnnotations = GeneralAnnotationStorage.getBuilder().addAnnotation("Description", currChrDesc).build();
            chrS = chrS.substring(0,chrS.indexOf(" "));
        } 
        return new Chromosome(chrS, -1, myAnnotations);
    }

}

/**
 * ReferenceGenomeSequence class.  This class is used to read chromosome sequences
 * from fasta files.  Data is stored as half-bytes packed into a byte array.
 * This byte array comprises the "value" for a hash map whose key is a
 * Chromosome object.
 *
 * The class also contains methods to obtain a full or partial genome sequence for a
 * specified stored chromosome.
 *
 * @author Lynn Johnson
 *
 */
class HalfByteGenomeSequence implements GenomeSequence{
    private Map<Chromosome, byte[]> chromPositionMap;
    private Map<Chromosome, Integer> chromLengthLookup=new HashMap<>();
    private RangeMap<Long,Chromosome> wholeGenomeIndexMap= TreeRangeMap.create();
    private final long genomeSize;


    protected HalfByteGenomeSequence(Map<Chromosome, byte[]>chromPositionMap) {
        this.chromPositionMap = chromPositionMap;
        chromPositionMap.entrySet().stream()
        .forEach(e -> chromLengthLookup.put(e.getKey(),e.getKey().getLength()));
        LongAdder genomeIndex=new LongAdder();
        chromosomes().stream().sorted()
        .forEach(chrom -> {
            int length=chromLengthLookup.get(chrom);
            wholeGenomeIndexMap.put(Range.closed(genomeIndex.longValue(),
                    genomeIndex.longValue()+length-1),chrom);
            genomeIndex.add(length);}
                );
        genomeSize=genomeIndex.longValue();
    }
    @Override
    public Set<Chromosome> chromosomes() {
        return chromPositionMap.keySet();
    }

    @Override
    public byte[] chromosomeSequence(Chromosome chrom) {
        return chromosomeSequence(chrom,1,chromLengthLookup.get(chrom));
    }

    @Override
    // Code assumes 1-based coordinates have been passed.  It will catch and return
    // null if the startSite is 0.  Otherwise, the user is on their own to ensure
    // input is 1-based.
    public byte[] chromosomeSequence(Chromosome chrom, int startSite, int lastSite) {
        startSite--;  //shift over to zero base
        lastSite--;   //shift over to zero base
        if (startSite < 0) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: starting parameter is less than 1 for 1-based method");; // method needs 1-based coordinates.
        byte[] packedBytes = chromPositionMap.get(chrom);
        if (packedBytes == null) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); // chromosome not found
        if (startSite > packedBytes.length*2 || lastSite > packedBytes.length*2 ) {
            System.out.println("GenomeSequenceBuilder, chrom: " + chrom.getName() + ", startSite: " 
              + startSite + ", lastSite: " + lastSite + ", chromSize: " + chromosomeSize(chrom)); 
            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested sequence is out of range"); // requested sequence is out of range
        }
        byte[] fullBytes = new byte[lastSite - startSite + 1];
        for (int i = startSite; i <= lastSite; i++) {
            fullBytes[i - startSite] = (byte) ((i % 2 == 0) ? ((packedBytes[i / 2] & 0xF0) >> 4) : (packedBytes[i / 2] & 0x0F));
        }
        return fullBytes;
    }

    @Override
    public byte genotype(Chromosome chrom, int position) {
        position--; // shift to 0-based
        byte[] packedBytes = chromPositionMap.get(chrom);
        if (packedBytes == null) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); 
        if (position > packedBytes.length*2 ) {
            System.out.println("GenomeSequenceBuilder, chrom: " + chrom.getName() + ", position: " 
              + position + ", chromSize: " + chromosomeSize(chrom)); 
            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested position is out of range");            
        }
        byte alleleAtPosition = (byte) ((position % 2 == 0) ? ((packedBytes[position / 2] & 0xF0) >> 4) : (packedBytes[position / 2] & 0x0F));
        return alleleAtPosition;
    }
    
    @Override
    public byte genotype(Chromosome chrom, Position positionObject) {
        int position = positionObject.getPosition();
        position--; // shift to 0-based
        byte[] packedBytes = chromPositionMap.get(chrom);
        if (packedBytes == null) throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: chromosome not found"); 
        if (position > packedBytes.length*2 ) {
            System.out.println("GenomeSequenceBuilder, chrom: " + chrom.getName() + ", position: " 
              + position + ", chromSize: " + chromosomeSize(chrom)); 
            throw new IllegalArgumentException("GenomeSequenceBuilder.chromosomeSequence: requested position is out of range");            
        }
        byte valueAtPosition = (byte) ((position % 2 == 0) ? ((packedBytes[position / 2] & 0xF0) >> 4) : (packedBytes[position / 2] & 0x0F));
        return valueAtPosition;
    }
    
    @Override
    public String genotypeAsString(Chromosome chrom, int position) {
        // get the bytes, convert to string, return the string
        byte valueAtPosition = genotype(chrom,position);
        return NucleotideAlignmentConstants.getHaplotypeNucleotide(valueAtPosition);
    }
    
    @Override
    public String genotypeAsString(Chromosome chrom, Position positionObject) {
        // get the bytes, convert to string, return the string
        // position is flipped to 0 based in call to genotype(chrom,position)
        int position = positionObject.getPosition();
        byte valueAtPosition = genotype(chrom,position);
        return NucleotideAlignmentConstants.getHaplotypeNucleotide(valueAtPosition);
    }
    
    @Override
    public String genotypeAsString(Chromosome chrom, int startSite, int endSite) {
        // This method takes physical positions.  chromosomeSequence flips the start/end to 0 based. 
        byte[] chromBytes = chromosomeSequence( chrom,  startSite,  endSite);
        StringBuilder builder = new StringBuilder();
        for (int idx = 0; idx < chromBytes.length; idx++) {
            builder.append(NucleotideAlignmentConstants.getHaplotypeNucleotide(chromBytes[idx]));
        }
        return builder.toString();
    }
    
    @Override
    public byte[] genomeSequence(long startSite, long lastSite) {
        if(lastSite-startSite>Integer.MAX_VALUE) throw
            new IllegalArgumentException("Less than "+Integer.MAX_VALUE+" sites must be requested at a time");
        byte[] fullBytes=new byte[(int)(lastSite-startSite+1)];
        long currentSiteToGet=startSite;
        while(currentSiteToGet<lastSite) {
            Map.Entry<Range<Long>,Chromosome> rangeChromEntry=wholeGenomeIndexMap.getEntry(currentSiteToGet);
            int chrStart=(int)(currentSiteToGet-rangeChromEntry.getKey().lowerEndpoint());
            int chrLast=(int)Math.min(rangeChromEntry.getKey().upperEndpoint()-rangeChromEntry.getKey().lowerEndpoint(),lastSite-rangeChromEntry.getKey().lowerEndpoint());
            byte[] chromoSeq=chromosomeSequence(rangeChromEntry.getValue(), chrStart+1,chrLast+1);  //+1 for 0 based genome, 1 based chromosomes
            System.arraycopy(chromoSeq,0,fullBytes,(int)(currentSiteToGet-startSite),chromoSeq.length);
            currentSiteToGet+=chromoSeq.length;
        }
        return fullBytes;
    }

    @Override
    public int chromosomeSize(Chromosome chromosome) {
        return chromLengthLookup.get(chromosome);
    }

    @Override
    public long genomeSize() {
        return genomeSize;
    }

    @Override
    public int numberOfChromosomes() {
        return chromPositionMap.size();
    }
    
    @Override
    public Map<Long, Tuple<Chromosome, Integer>> fullRefCoordinateToChromCoordinate(ArrayList<Long> coordinates) {
        // Returns 0-based value from Chromosome array (values are stored as 0-based) 
        Map<Long, Tuple<Chromosome, Integer>> mappedCoordinates = new ConcurrentHashMap<Long, Tuple<Chromosome, Integer>>();
        coordinates.stream().parallel().forEach(coordinate -> {
            Map.Entry<Range<Long>,Chromosome> rangeChromEntry=wholeGenomeIndexMap.getEntry(coordinate);
            Chromosome curChrom = rangeChromEntry.getValue();
            long chromCoordinate = coordinate - rangeChromEntry.getKey().lowerEndpoint();
            Tuple<Chromosome, Integer> chromWithCoordinate = new Tuple<>(curChrom, (int)chromCoordinate);
            mappedCoordinates.put(coordinate, chromWithCoordinate);
        });
        return mappedCoordinates;
    }
}