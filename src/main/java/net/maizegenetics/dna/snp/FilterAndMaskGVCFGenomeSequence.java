package net.maizegenetics.dna.snp;

import com.google.common.collect.SetMultimap;
import net.maizegenetics.dna.map.GVCFGenomeSequence;
import net.maizegenetics.dna.map.GVCFGenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.util.BitSet;
import org.apache.log4j.Logger;

/**
 * Created by zrm22 on 4/3/17.
 */
public class FilterAndMaskGVCFGenomeSequence {
    private static final Logger myLogger = Logger.getLogger(FilterAndMaskGVCFGenomeSequence.class);

    private FilterAndMaskGVCFGenomeSequence() {
    }

    public static GVCFGenomeSequence getInstance(GVCFGenomeSequence base, String annoName, int symbol, int threshold,boolean mode) {
        //mode = false for filter
        //mode = true for masking with Ns
        PositionList positions = base.getGVCFPositions();
        BitSet currentBitSet = null;

        if(mode) {
            currentBitSet = base.getMaskBitSet();
        }
        else {
            currentBitSet = base.getFilterBitSet();

        }
        for(int i = 0; i < positions.size(); i++) {
            Position position = positions.get(i);
            //check the position
            boolean alterPosition = filterPosition(position,annoName, symbol,threshold);
            boolean currentBitValue = (mode)?base.getMaskBitSet().fastGet(i):base.getFilterBitSet().fastGet(i);
            //Check if bit has already been flipped
            if(!currentBitValue) {
                //if not check if we need to filter or mask the Position
                if(alterPosition) {
                    //if so flip bit
                    currentBitSet.fastFlip(i);
                }
            }
        }
        try {
            if(mode) {
                return (GVCFGenomeSequence) GVCFGenomeSequenceBuilder.instance(base, currentBitSet, base.getFilterBitSet());
            }
            else {
                return (GVCFGenomeSequence) GVCFGenomeSequenceBuilder.instance(base, base.getMaskBitSet(), currentBitSet);
            }
        }
        catch(Exception e) {
            System.out.println("ERROR");
            e.printStackTrace();
            return null;
        }

    }

    //Return false to keep position
    //Return true to filter out
    private static boolean filterPosition(Position position, String annoName,int symbol, int threshold) {
        //0,  1,  2,  3, 4
        //<, <=, ==, >=, >
        SetMultimap<String,String> annos = position.getAnnotation().getAnnotationAsMap();

        switch(symbol) {
            case 0:
                //less than
                if(!annos.containsKey(annoName)) {
                    return false;
                }
                else {
                    if(Integer.parseInt((String)annos.get(annoName).toArray()[0])<threshold) {
                        return false;
                    }
                    else {
                        return true;
                    }
                }

            case 1:
                //less than or equal
                if(!annos.containsKey(annoName)) {
                    return false;
                }
                else {
                    if(Integer.parseInt((String)annos.get(annoName).toArray()[0])<=threshold) {
                        return false;
                    }
                    else {
                        return true;
                    }
                }

            case 2:
                //equal
                if(!annos.containsKey(annoName)) {
                    return true;
                }
                else {
                    if(Integer.parseInt((String)annos.get(annoName).toArray()[0])==threshold) {
                        return false;
                    }
                    else {
                        return true;
                    }
                }

            case 3:
                //greater than or equal
                if(!annos.containsKey(annoName)) {
                    return true;
                }
                else {
                    if(Integer.parseInt((String)annos.get(annoName).toArray()[0])>=threshold) {
                        return false;
                    }
                    else {
                        return true;
                    }
                }

            case 4:
                //greater than
                if(!annos.containsKey(annoName)) {
                    return true;
                }
                else {
                    if(Integer.parseInt((String)annos.get(annoName).toArray()[0])>threshold) {
                        return false;
                    }
                    else {
                        return true;
                    }
                }
        }
        return true;
    }
}
