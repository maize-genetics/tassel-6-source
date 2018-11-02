package net.maizegenetics.dna.tag;

import cern.colt.list.IntArrayList;
import cern.colt.list.ShortArrayList;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.Multiset;
import com.google.common.primitives.Ints;
import com.google.common.primitives.UnsignedBytes;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TShortIntIterator;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TShortIntHashMap;
import org.xerial.snappy.Snappy;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;


/**
 * Builder for TaxaDistribution.  Deals with the
 */
public class TaxaDistBuilder {
    private TaxaDistBuilder() {}
    
    /**
     * Create a TaxaDistribution with only a single taxa with a tag.  Very memory efficient, but needs conversion add
     * an additional taxon with a tag
     * @param maxTaxa
     * @param taxaWithTag
     * @return
     */
    public static TaxaDistribution create(int maxTaxa, int taxaWithTag) {
        TaxaDistribution td=create(maxTaxa);
        return td.increment(taxaWithTag);
    }
    
    /**
     * Create an expandable TaxaDistribution with no initial values
     * @param maxTaxa
     * @return
     */
    public static TaxaDistribution create(int maxTaxa) {
        if(maxTaxa<Short.MAX_VALUE) return new TaxaDistShortExpandable(maxTaxa);
        return new TaxaDistIntExpandable(maxTaxa);
    }
    
    /**
     * Create an fixed TaxaDistribution with set values
     * @param maxTaxa
     * @return
     */
    public static TaxaDistribution create(int maxTaxa, int[] taxaWithTags, int[] depthOfTags) {
        return new TaxaDistFixed(maxTaxa,taxaWithTags,depthOfTags);
        
    }
    
    /**
     * Create an fixed TaxaDistribution with set values
     * @param encodedTaxaDistribution byte array of encoded TaxaDistribution
     * @return
     */
    public static TaxaDistribution create(byte[] encodedTaxaDistribution) {
        int[][] decodedTD=getDepthMatrixForEncodedDepths(encodedTaxaDistribution);
        return new TaxaDistFixed(decodedTD[2][0],decodedTD[0], decodedTD[1]);
        
    }
    
    /**
     * Copies a TaxaDistribution to an expandable TaxaDistribution.  Can be used to convert a single TaxaDistribution
     * to an expandable version
     * @param srcTaxaDist
     * @return
     */
    public static TaxaDistribution create(TaxaDistribution srcTaxaDist) {
        TaxaDistribution dstTD=create(srcTaxaDist.maxTaxa());
        int[] depths=srcTaxaDist.depths();
        for (int taxaIndex = 0; taxaIndex < depths.length; taxaIndex++) {
            for (int j = 0; j < depths[taxaIndex]; j++) {
                dstTD.increment(taxaIndex);
            }
        }
        return dstTD;
    }
    
    /**
     * Combines a TaxaDistribution to an expandable TaxaDistribution.  Can be used to convert a single TaxaDistribution
     * to an expandable version
     * @param taxaDist1
     * @return
     */
    public static TaxaDistribution combine(TaxaDistribution taxaDist1, TaxaDistribution taxaDist2) {
        if(taxaDist1.maxTaxa()!=taxaDist2.maxTaxa()) throw new IllegalStateException("TaxaDistributions not of same size");
        int[] depths1=taxaDist1.depths();
        int[] depths2=taxaDist2.depths();
        int taxaWithDepth=0;
        for (int i = 0; i < depths1.length; i++) {
            depths1[i]+=depths2[i];
            if(depths1[i]>0) taxaWithDepth++;
        }
        int[] taxaWithTags=new int[taxaWithDepth];
        int[] depthOfTags=new int[taxaWithDepth];
        taxaWithDepth=0;
        for (int i = 0; i < depths1.length; i++) {
            if(depths1[i]>0) {
                taxaWithTags[taxaWithDepth]=i;
                depthOfTags[taxaWithDepth]=depths1[i];
                taxaWithDepth++;
            }
        }
        return create(taxaDist1.maxTaxa(),taxaWithTags,depthOfTags);
    }
    
    /**
     * Create expandable TaxaDistribution from encoded taxa distribution.  The int[] encoding use first three bytes
     * for taxa index, and last byte for depth as unsigned byte.  Depths greater than 256 just increase.
     * @param maxTaxa
     * @param encodeTaxaDepths
     * @return
     */
    public static TaxaDistribution create(int maxTaxa, int[] encodeTaxaDepths) {
        TaxaDistribution dstTD=create(maxTaxa);
        for (int taxaDepth : encodeTaxaDepths) {
            int taxa=taxaDepth>>>8;
            int depth=UnsignedBytes.toInt((byte)taxaDepth);
            for (int i = 0; i < depth; i++) {
                dstTD.increment(taxa);
            }
        }
        return dstTD;
    }
    
    public static int[][] getDepthMatrixForEncodedDepths(byte[] input) {
        try{
            final int maxValueInInt=UnsignedBytes.toInt(UnsignedBytes.MAX_VALUE);
            ByteBuffer bb=ByteBuffer.wrap(Snappy.uncompress(input));
            int maxTaxa=bb.getInt();
            int taxaWithDepth=bb.getInt();
            int[][] result=new int[3][];
            result[0]=new int[taxaWithDepth];
            result[1]=new int[taxaWithDepth];
            result[2]=new int[]{maxTaxa};
            for (int i = 0; i < taxaWithDepth; i++) {
                int space=(i>0)?result[0][i-1]:0;
                int inc=UnsignedBytes.toInt(bb.get());
                while(inc==maxValueInInt) {
                    space+=inc;
                    inc=UnsignedBytes.toInt(bb.get());
                }
                space+=inc;
                result[0][i]=space;
            }
            for (int i = 0; i < taxaWithDepth; i++) {
                int depth=0;
                int inc=UnsignedBytes.toInt(bb.get());
                while(inc==maxValueInInt) {
                    depth+=inc;
                    inc=UnsignedBytes.toInt(bb.get());
                }
                depth+=inc;
                result[1][i]=depth;
            }
            return result;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        
        
        
    }
    
    
}

abstract class AbstractTaxaDistribution implements TaxaDistribution {
    
    @Override
    public int[][] taxaWithDepths() {
        int[] depths=depths();
        int countNoZero=0;
        for (int depth : depths) if(depth>0) countNoZero++;
        int[][] taxaDepth=new int[2][countNoZero];
        countNoZero=0;
        for (int i = 0; i <depths.length ; i++) {
            if(depths[i]>0) {
                taxaDepth[0][countNoZero]=i;
                taxaDepth[1][countNoZero]=depths[i];
                countNoZero++;
            }
        }
        return taxaDepth;
    }
    
    @Override
    public byte[] encodeTaxaDepth() {
        int[][] tds=taxaWithDepths();
        ByteBuffer bb=ByteBuffer.allocate(8+(maxTaxa()/64)+(2*tds[0].length)+(totalDepth()/64));
        bb.putInt(maxTaxa());  //maximum number of taxa with depth
        bb.putInt(tds[0].length);  //number of taxa with depth
        int priorTaxa=0;
        for (int i = 0; i < tds[0].length; i++) {
            int space=tds[0][i]-priorTaxa;
            while(space>=0) {
                bb.put(UnsignedBytes.saturatedCast(space));
                space-=255;
            }
            priorTaxa=tds[0][i];
        }
        for (int i = 0; i < tds[1].length; i++) {
            int depth=tds[1][i];
            while(depth>=0) {
                bb.put(UnsignedBytes.saturatedCast(depth));
                depth-=255;
            }
        }
        try{
            return Snappy.compress(Arrays.copyOf(bb.array(), bb.position()));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
    
    @Override
    public Multiset<Integer> taxaDepthMap() {
        int[][] tds=taxaWithDepths();
        ImmutableMultiset.Builder<Integer> taxaCnts = new ImmutableMultiset.Builder<>();
        for (int i = 0; i < tds[0].length; i++) {
            taxaCnts.setCount(tds[0][i],tds[1][i]);
        }
        return taxaCnts.build();
    }
    
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null) return false;
        if (!(o instanceof TaxaDistribution)) return false;
        
        TaxaDistribution that = (TaxaDistribution) o;
        if (maxTaxa() != that.maxTaxa()) return false;
        int[][] thisTDS=taxaWithDepths();
        int[][] thatTDS=that.taxaWithDepths();
        if (!Arrays.equals(thisTDS[0], thatTDS[0])) return false;
        if (!Arrays.equals(thisTDS[1], thatTDS[1])) return false;
        return true;
    }
    
    @Override
    public String toString() {
        return "TaxaDist{" +
                "taxaWithRead=" + numberOfTaxaWithTag() +
                ", totalDepth=" + totalDepth() +
                ", maxTaxa=" + maxTaxa() +
                ", "+taxaDepthMap().toString()+
                '}';
    }
}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistShortExpandable extends AbstractTaxaDistribution  {
    
    private ShortArrayList taxaWithTag;
    private TShortIntHashMap taxaTagMap= null;
    private int totalDepth;
    private final int maxTaxa;
    
    public TaxaDistShortExpandable(int maxTaxa) {
        this.maxTaxa=maxTaxa;
        taxaWithTag = new ShortArrayList(1); // defaults to 10 items, we only want size of 1 initially
    }
    
    @Override
    public synchronized TaxaDistribution increment(int taxaNum) {
        if (true) {
            taxaWithTag.add((short)taxaNum);
        } else { // add previous values to map
            if(taxaTagMap==null) convertListToMap();
            taxaTagMap.increment((short)taxaNum);
        }
        totalDepth++;
        return this;
    }
    
    private void convertListToMap() {
        taxaTagMap= new TShortIntHashMap(maxTaxa);
        for (short taxaIndex : taxaWithTag.elements()) {
            taxaTagMap.adjustOrPutValue(taxaIndex, (short) 1, (short) 1);
        }
        taxaWithTag=null;
    }
    
    @Override
    public int[] depths() {
        int[] depths=new int[maxTaxa];
        if (taxaWithTag!=null) {
            for (int i=0; i<taxaWithTag.size(); i++) {
                depths[taxaWithTag.getQuick(i)]++;
            }
        } else {
            for (TShortIntIterator sst = taxaTagMap.iterator(); sst.hasNext();) {
                sst.advance();
                depths[sst.key()]=sst.value();
            }
        }
        return depths;
    }
    
    @Override
    public int totalDepth(){
        return totalDepth;
    }
    @Override
    public int maxTaxa() {
        return maxTaxa;
    }
    
    @Override
    public int numberOfTaxaWithTag() {
        return taxaWithDepths()[0].length;
    }
    
    @Override
    public int memorySize() {
        //minimal size 8 (object) + 12 (outer short array) + 12 (sizeArray) + 4+ 4 = 40
        int size=40;
        if (taxaWithTag != null) {
            size+= taxaWithTag.elements().length*2;
            //size+= taxaWithTag.size() * 2; // 8 + 4
        }
        return size;
    }
    
    private int unSignShort(short v) {
        if(v<0) return -v+1;
        return v;
    }
}

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistIntExpandable extends AbstractTaxaDistribution  {
    
    private IntArrayList taxaWithTag;
    private TIntIntHashMap taxaTagMap= null;
    private int totalDepth;
    private final int maxTaxa;
    
    public TaxaDistIntExpandable(int maxTaxa) {
        this.maxTaxa=maxTaxa;
        taxaWithTag = new IntArrayList(1); // defaults to 10 items, we only want size of 1 initially
    }
    
    @Override
    public synchronized TaxaDistribution increment(int taxaNum) {
        if (true) {
            taxaWithTag.add(taxaNum);
        } else { // add previous values to map
            if(taxaTagMap==null) convertListToMap();
            taxaTagMap.increment(taxaNum);
        }
        totalDepth++;
        return this;
    }
    
    private void convertListToMap() {
        taxaTagMap= new TIntIntHashMap(maxTaxa);
        for (int taxaIndex : taxaWithTag.elements()) {
            taxaTagMap.adjustOrPutValue(taxaIndex, 1, 1);
        }
        taxaWithTag=null;
    }
    
    @Override
    public int[] depths() {
        int[] depths=new int[maxTaxa];
        if (taxaWithTag!=null) {
            for (int i=0; i<taxaWithTag.size(); i++) {
                depths[taxaWithTag.getQuick(i)]++;
            }
        } else {
            for (TIntIntIterator sst = taxaTagMap.iterator(); sst.hasNext();) {
                sst.advance();
                depths[sst.key()]=sst.value();
            }
        }
        return depths;
    }
    
    @Override
    public int totalDepth(){
        return totalDepth;
    }
    @Override
    public int maxTaxa() {
        return maxTaxa;
    }
    
    @Override
    public int numberOfTaxaWithTag() {
        return taxaWithDepths()[0].length;
    }
    
    @Override
    public int memorySize() {
        //minimal size 8 (object) + 12 (outer short array) + 12 (sizeArray) + 4+ 4 = 40
        int size=40;
        if (taxaWithTag != null) {
            size+= taxaWithTag.size() * 8; // 8 + 4
        }
        if (taxaTagMap!= null) {
            size+= taxaTagMap.size() * 12; // 4 + 4
        }
        return size;
    }

}
    

/**
 * This is a specialized multiset for recording the distribution of a single across taxa.
 *
 * HashMap or Multiset or large arrays could be reasonable approaches by they do not scale well with hundreds of
 * taxa scored out of the thousands.
 * @author Ed Buckler
 */
class TaxaDistFixed extends AbstractTaxaDistribution  {
    //minimal size 8 + 12 + 12 + 10 + 12 + 4+ 4 = 66
    //This could be changed for the singletons by just making a new class
    private byte[] compTaxaSize;
    private int numTaxaWithTags;
    private int totalDepth;
    private final int maxTaxa;
    
    public TaxaDistFixed(int maxTaxa, int[] taxaWithTags, int[] depthOfTags) {
        this.maxTaxa=maxTaxa;
        numTaxaWithTags =taxaWithTags.length;
        totalDepth=0;
        for (int depthOfTag : depthOfTags) {totalDepth+=depthOfTag;}
        try{
            compTaxaSize=Snappy.compress(Ints.concat(taxaWithTags,depthOfTags));
            //System.out.printf("tn:%d\tdepth:%d\ts:%d%n", numTaxaWithTags,totalDepth,compTaxaSize.length);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    @Override
    public synchronized TaxaDistribution increment(int taxaNum) {
        throw new UnsupportedOperationException("TaxaDistFixed cannot be increment.  Change to expandable first first.");
    }
    
    @Override
    public int[] depths() {
        try {
            int[] depths = new int[maxTaxa];
            int[] taxaSize = Snappy.uncompressIntArray(compTaxaSize);
            for (int i = 0; i < numTaxaWithTags; i++) {
                depths[taxaSize[i]]=taxaSize[i+ numTaxaWithTags];
            }
            return depths;
        }  catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }
    
    @Override
    public int[][] taxaWithDepths() {
        try {
            int[][] taxaDepth=new int[2][numTaxaWithTags];
            int[] taxaSize = Snappy.uncompressIntArray(compTaxaSize);
            taxaDepth[0]=Arrays.copyOf(taxaSize, numTaxaWithTags);
            taxaDepth[1]=Arrays.copyOfRange(taxaSize, numTaxaWithTags,taxaSize.length);
            return taxaDepth;
        }  catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }
    
    @Override
    public int totalDepth(){
        return totalDepth;
    }
    
    @Override
    public int maxTaxa() {
        return maxTaxa;
    }
    
    @Override
    public int numberOfTaxaWithTag() {
        return numTaxaWithTags;
    }
    
    @Override
    public int memorySize() {
        //minimal size 8 (object) + 12 (sizeArray) + 4+ 4 = 40
        return 28+compTaxaSize.length;
    }
}
