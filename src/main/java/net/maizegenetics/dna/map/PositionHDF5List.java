package net.maizegenetics.dna.map;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.lang.reflect.Array;
import java.nio.IntBuffer;
import java.util.*;
import java.util.concurrent.ExecutionException;

/**
 * HDF5 immutable instance of {@link PositionList}.  Use the {@link PositionListBuilder}
 * to create the list.
 *
 * @author Ed Buckler
 */
final class PositionHDF5List implements PositionList {
    private final IHDF5Reader reader;
    private final int numPositions;
    private final Map<Chromosome,ChrOffPos> myChrOffPosTree;
    private final Map<String,Chromosome> myChrNameHash;
    private final int[] chrOffsets;  //starting site for each chromosome
    private final Chromosome[] chrIndex;
    private final byte[][] alleles;  //store reference, major, etc. alleles only fully initialized if requested once.
    private final String genomeVersion;

    /*Byte representations of DNA sequences are stored in blocks of 65536 sites*/
    public static final int BLOCKSIZE=1<<16;
    public static final int blockMask=BLOCKSIZE-1;
    public static final int siteMask=~(BLOCKSIZE-1);

    private LoadingCache<Integer,Position> mySiteList; //key site > AnnoPos
    private CacheLoader<Integer,Position> annoPosLoader = new CacheLoader<Integer,Position>()  {
        @Override
        public Position load(Integer key) {
            List<Integer> toFill=new ArrayList<>();
            toFill.add(key);
            try {
                mySiteList.putAll(loadAll(toFill));
                return get(key);
            } catch (Exception e) {
                e.printStackTrace();
                return null;
            }
        }

        @Override
        public Map<Integer, Position> loadAll(Iterable<? extends Integer> keys) throws Exception {
            int key=keys.iterator().next();
            HashMap<Integer, Position> result=new HashMap<Integer, Position>(BLOCKSIZE);
            byte[][] afOrder=new byte[2][];
            byte[] ref, anc;
            float[] maf;
            float[] paf;
            String[] snpIDs;
            int startSite=key&siteMask;
            int length=((numPositions-startSite)<BLOCKSIZE)?numPositions-startSite:BLOCKSIZE;
            synchronized(reader) {
                afOrder = reader.int8().readMatrixBlockWithOffset(Tassel5HDF5Constants.ALLELE_FREQ_ORD, 2, length, 0l, startSite);
                ref=HDF5Utils.getHDF5ReferenceAlleles(reader,startSite,length);
                anc=HDF5Utils.getHDF5AncestralAlleles(reader,startSite,length);
                maf= reader.float32().readArrayBlockWithOffset(Tassel5HDF5Constants.MAF,length, startSite);
                paf= reader.float32().readArrayBlockWithOffset(Tassel5HDF5Constants.SITECOV,length, startSite);
                snpIDs=reader.string().readArrayBlockWithOffset(Tassel5HDF5Constants.SNP_IDS,length, startSite);
            }
            for (int i=0; i<length; i++) {
                int site=i+startSite;
                Chromosome chr=chromosome(site);
                ChrOffPos cop=myChrOffPosTree.get(chr);
                int pos=cop.position[site-cop.startSiteOff];
                Position p=new GeneralPosition.Builder(chr,pos)
                        .snpName(snpIDs[i])
                        .allele(WHICH_ALLELE.GlobalMajor,afOrder[0][i])
                        .allele(WHICH_ALLELE.GlobalMinor,afOrder[1][i])
                        .allele(WHICH_ALLELE.Reference, ref[i])
                        .allele(WHICH_ALLELE.Ancestral, anc[i])
                        .maf(maf[i])
                        .siteCoverage(paf[i])
                        .build();
                result.put(site,p);
            }
            return result;
        }
    };

    private class ChrOffPos {
        final int startSiteOff;
        final int endSiteOff;
        final int[] position;
        private ChrOffPos(int startSiteOff, int endSiteOff, int[] position) {
            this.startSiteOff=startSiteOff;
            this.endSiteOff=endSiteOff;
            this.position=position;
        }
    }

    PositionHDF5List(IHDF5Reader reader) {
        this.reader=reader;
        if (reader.hasAttribute(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH,Tassel5HDF5Constants.POSITION_GENOME_VERSION)) {
            genomeVersion = reader.string().getAttr(Tassel5HDF5Constants.POSITION_ATTRIBUTES_PATH,Tassel5HDF5Constants.POSITION_GENOME_VERSION);
        } else {
            genomeVersion = null;
        }
        int[] variableSites = reader.readIntArray(Tassel5HDF5Constants.POSITIONS);
        this.numPositions=variableSites.length;
        alleles=new byte[WHICH_ALLELE.COUNT][];//only fully initialized if requested
        String[] lociStrings = reader.readStringArray(Tassel5HDF5Constants.CHROMOSOMES);
        ArrayList<Chromosome> chrs=new ArrayList<>();
        for (String ls : lociStrings) {
            chrs.add(new Chromosome(ls));
        }

        int[] locusIndices = reader.readIntArray(Tassel5HDF5Constants.CHROMOSOME_INDICES);
        myChrOffPosTree=new TreeMap<>();
        myChrNameHash=new HashMap<>();
        int currStart=0;
        int currLocusIndex=locusIndices[0];
        chrOffsets=new int[chrs.size()];
        chrIndex=new Chromosome[chrs.size()];
        int cI=0;
        for (int i=0; i<locusIndices.length; i++) {
            // Check for next loci.  if loci has changed, store
            // the beginning/end position indices for the previous loci
            if((i==(locusIndices.length-1))||currLocusIndex!=locusIndices[i]) {  
                //int end=(i==locusIndices.length-1)?i:i-1;
                int end=((i==(locusIndices.length-1 ))&& (currLocusIndex==locusIndices[i]))?i:i-1;
                
                // Using endRange variable instead of "end/end+i" as the value 
                // varies based on whether we are at the last element in the array.  
                // Found that without this, when both the last and the second-to-last
                // positions each belong to a single chromosome, the last position ended
                // up being counted with both the second-to-last chromosome AND the last
                // chromosome.  See TAS-915
                int endRange = (i==locusIndices.length-1)?i+1:i; 
                int[] cPos=Arrays.copyOfRange(variableSites,currStart,endRange); 
                Chromosome currChr=chrs.get(currLocusIndex);
                myChrOffPosTree.put(currChr, new ChrOffPos(currStart, end, cPos));
                myChrNameHash.put(currChr.getName(),currChr);
                chrOffsets[cI]=currStart;
                chrIndex[cI]=currChr;
                cI++;
                if (i==(locusIndices.length-1) && currLocusIndex !=locusIndices[i]){
                    // This is the case where the last chromosome has only 1 position
                    // Could be case of some chromosomes and lots of small contigs
                    // not mapped onto the chromosomes.
                    // We finished off the previous one above, now handle the last one
                    // Without this check, the last position, if it belongs to a new loci,
                    // is skipped.
                    currLocusIndex=locusIndices[i];
                    currStart = i;
                    end = i; // start and end are the same - only 1 element
                    int[] cPosLast=Arrays.copyOfRange(variableSites,currStart, end);
                    Chromosome lastChr=chrs.get(currLocusIndex); 
                    myChrOffPosTree.put(lastChr, new ChrOffPos(currStart, end, cPosLast));
                    myChrNameHash.put(lastChr.getName(),lastChr);
                    chrOffsets[cI]=currStart;
                    chrIndex[cI]=currChr;
                }                   
                currStart=i;
                currLocusIndex=locusIndices[i];
            }
        }
       // rangeMap=rangeMapBuild.build();
        mySiteList= CacheBuilder.newBuilder()
                .maximumSize(1000000)
                .build(annoPosLoader);
    }

    private void loadAllele(WHICH_ALLELE alleleType) {
        if(alleles[alleleType.index()]==null) {
            switch (alleleType) {
                case Reference:
                    alleles[alleleType.index()]=HDF5Utils.getHDF5ReferenceAlleles(reader);
                    break;
                case GlobalMajor:
                    alleles[alleleType.index()]=HDF5Utils.getHDF5Alleles(reader, WHICH_ALLELE.Major);
                    break;
                case GlobalMinor:
                    alleles[alleleType.index()]=HDF5Utils.getHDF5Alleles(reader, WHICH_ALLELE.Minor);
                    break;
                case Ancestral:
                    alleles[alleleType.index()]=HDF5Utils.getHDF5AncestralAlleles(reader);
                    break;
                case HighCoverage:
                    break;
            }
        }

    }

    @Override
    public byte allele(WHICH_ALLELE alleleType, int site) {
        try {
            return mySiteList.get(site).getAllele(alleleType);
        } catch (ExecutionException e) {
            e.printStackTrace();
            return GenotypeTable.UNKNOWN_ALLELE;
        }
    }

    @Override
    public byte[] alleles(WHICH_ALLELE alleleType, int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        System.arraycopy(alleleForAllSites(alleleType),startSite,result,0, result.length);
        return result;
    }

    @Override
    public byte[] alleleForAllSites(WHICH_ALLELE alleleType) {
        if(alleles[alleleType.index()]==null) {loadAllele(alleleType);}
        return alleles[alleleType.index()];
    }

    @Override
    public boolean hasReference() {
        if (genomeVersion == null) {
            return false;
        }
        return true;
    }

    @Override
    public String siteName(int site) {
        try {
            return mySiteList.get(site).getSNPID();
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public int numberOfSites() {
        return numPositions;
    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        return myChrOffPosTree.get(chromosome).position.length;
    }

    @Override
    public int[] startAndEndOfChromosome(Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return null;
        return new int[]{cop.startSiteOff,cop.endSiteOff};
    }

    @Override
    public int chromosomalPosition(int site) {
        int i=Arrays.binarySearch(chrOffsets,site);
        if(i<0) i=-(i+1)-1;
        Chromosome chr=chrIndex[i];
//        Chromosome chr=rangeMap.get(site);
        ChrOffPos cop=myChrOffPosTree.get(chr);
        return cop.position[site-cop.startSiteOff];
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return Integer.MIN_VALUE;
        int i=Arrays.binarySearch(cop.position, physicalPosition); //AvgPerObj:227.5715ns  for 2million positions
        while((i>0)&&(physicalPosition==cop.position[i-1])) {i--;} //backup to the first position if there are duplicates
        i+=(i<0)?-cop.startSiteOff:cop.startSiteOff;
        return i;
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        int result=siteOfPhysicalPosition(physicalPosition, chromosome);
        if (result < 0) {return result;}
        else {
            if (snpName.equals(siteName(result))) {return result;
            } else {
                int index=result;
                while ((index < numPositions) && (chromosomalPosition(index) == physicalPosition)) {
                    if (snpName.equals(siteName(index))) {return index;}
                    result++;
                }
                return -result - 1;
            }
        }
    }

    @Override
    public int[] physicalPositions() {
        int[] result=new int[numPositions];
        IntBuffer ib=IntBuffer.wrap(result);
        for (ChrOffPos cop: myChrOffPosTree.values()) {
            ib.put(cop.position);
        }
        return result;
    }

    @Override
    public String chromosomeName(int site) {
        return chromosome(site).getName();
       // return rangeMap.get(site).getName();
    }

    @Override
    public Chromosome chromosome(int site) {
        int i=Arrays.binarySearch(chrOffsets,site);
        if(i<0) i=-(i+1)-1;
        Chromosome chr=chrIndex[i];
        return chr;
       // return rangeMap.get(site);
    }

    @Override
    public Chromosome chromosome(String name) {
        return myChrNameHash.get(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myChrOffPosTree.keySet().toArray(new Chromosome[0]);
    }

    @Override
    public int numChromosomes() {
        return myChrOffPosTree.size();
    }

    @Override
    public int[] chromosomesOffsets() {
        int[] result=new int[myChrOffPosTree.size()];
        int index=0;
        for (ChrOffPos cop: myChrOffPosTree.values()) {
            result[index++]=cop.startSiteOff;
        }
        return result;
    }

    @Override
    public int indelSize(int site) {
        try{return mySiteList.get(site).getKnownVariants()[1].length();}
    catch (ExecutionException e) {
        e.printStackTrace();
        return -1;
    }
    }

    @Override
    public boolean isIndel(int site) {
        try{return mySiteList.get(site).isIndel();}
        catch (ExecutionException e) {
            e.printStackTrace();
            return false;
        }
    }

    @Override
    public String genomeVersion() {
        return genomeVersion;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return true;
    }

    // List methods

    @Override
    public int size() {
        return numPositions;
    }

    @Override
    public boolean isEmpty() {
        return (numPositions==0);
    }

    @Override
    public boolean contains(Object o) {
        if(o instanceof Position) {
            Position p=(Position)o;
            int site=siteOfPhysicalPosition(p.getPosition(),p.getChromosome());
            //test for SNP ID also?
            if(site>=0) return true;
        }
        return false;
    }

    @Override
    public Iterator<Position> iterator() {
        Iterator<Position> it = new Iterator<Position>() {
            private int currentIndex = 0;
            @Override
            public boolean hasNext() {
                return currentIndex < numPositions;
            }
            @Override
            public Position next() {
                return get(currentIndex++);
            }
            @Override
            public void remove() {
                throw new UnsupportedOperationException("This Class is Immutable.");
            }
        };
        return it;
       // return mySiteList.iterator();
    }

    @Override
    public Object[] toArray() {
        Position[] aps=new Position[numPositions];
        for (int i=0; i<numPositions; i++) {
            aps[i]=get(i);
        }
        return aps;
    }

    @Override
    public <Position> Position[] toArray(Position[] a) {
        if (a.length < numPositions) {
            // If array is too small, allocate the new one with the same component type
            a = (Position[])Array.newInstance(a.getClass().getComponentType(), numPositions);
        } else if (a.length > numPositions) {
            // If array is to large, set the first unassigned element to null
            a[numPositions] = null;
        }
        for (int i=0; i<numPositions; i++) {
            a[i]=(Position)get(i);
        }
        return a;
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean add(Position e) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not implemented yet.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(Collection<? extends Position> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(int index, Collection<? extends Position> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void clear() {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public Position get(int index) {
        try {
            return mySiteList.get(index);
        } catch (ExecutionException e) {
            return null;
        }
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public Position set(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void add(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public Position remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        if(o instanceof Position) {
            Position p=(Position)o;
            int site=siteOfPhysicalPosition(p.getPosition(),p.getChromosome());
            if(site>=0) return site;
        }
        return -1;
    }

    @Override
    public int lastIndexOf(Object o) {
        throw new UnsupportedOperationException("Not implemented yet.");
       // return mySiteList.lastIndexOf(o);
    }

    @Override
    public ListIterator<Position> listIterator() {
        throw new UnsupportedOperationException("Not implemented yet.");
     //   return mySiteList.listIterator();
    }

    @Override
    public ListIterator<Position> listIterator(int index) {
        throw new UnsupportedOperationException("Not implemented yet.");
       // return mySiteList.listIterator(index);
    }

    @Override
    public List<Position> subList(int fromIndex, int toIndex) {
        throw new UnsupportedOperationException("Not implemented yet.");
        //return mySiteList.subList(fromIndex, toIndex);
    }
}


