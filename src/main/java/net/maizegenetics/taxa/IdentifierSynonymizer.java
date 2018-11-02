package net.maizegenetics.taxa;


import net.maizegenetics.taxa.Taxon.Builder;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.AbstractTableReport;

import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;

import org.apache.commons.codec.language.DoubleMetaphone;
import org.apache.commons.codec.language.Metaphone;
import org.apache.commons.codec.language.RefinedSoundex;
import org.apache.commons.codec.language.Soundex;

import com.google.common.collect.ImmutableMultimap;


/**
 * User: Ed
 * Date: Mar 30, 2005
 * Time: 1:39:47 PM
 */
public class IdentifierSynonymizer extends AbstractTableReport implements Serializable, TableReport {

    ArrayList<TaxaList> taxaListSynonymized = new ArrayList<>();
    TaxaList tempTaxaList = null;
    private TaxaList referenceIDGroup;
    private int matchCount = 0;
    private int unmatchCount = 0;
    private int technique = 0;
    private String delimiter = "";
    private double globalMin = Double.POSITIVE_INFINITY;
    private double globalMax = Double.NEGATIVE_INFINITY;

    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets) {
        init2(preferredTaxa, alternateTaxaSets);
    }
    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets,int technique) {
        this.technique = technique;
        init2(preferredTaxa, alternateTaxaSets);
    }
    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList[] alternateTaxaSets,int technique,String delimiter) {
        this.technique = technique;
        this.delimiter = delimiter;
        init2(preferredTaxa, alternateTaxaSets);
    }

    public IdentifierSynonymizer(TaxaList preferredTaxa, TaxaList alternateTaxa) {
        TaxaList[] alternateTaxaSets = new TaxaList[1];
        alternateTaxaSets[0] = alternateTaxa;
        init2(preferredTaxa, alternateTaxaSets);
    }


    //Taxa List implementation This will become init when finished
    //TODO Stream implementation
    private void init2(TaxaList referenceTaxa, TaxaList[] alternateTaxaSets) {
        referenceIDGroup = referenceTaxa;
        ImmutableMultimap.Builder<String,String> changeSynBuild = new ImmutableMultimap.Builder();
        
        //Go Through each alternateTaxaSet and compute the similarity to the referenceTaxa entries
        for(TaxaList altTaxaList:alternateTaxaSets) { 
            for(Taxon altTaxon:altTaxaList){
                ArrayList<String> theBest = findBestMatch(altTaxon.getName(),referenceTaxa);
                //If the score is better than the current score, put it on the change Syn
                if (theBest.size() == 1) {
                    String bs = theBest.get(0);
                    if(!bs.equals(altTaxon.getName())) {
                        changeSynBuild.put(altTaxon.getName(), bs);
                    }
                    matchCount++;
                } else {
                    for(String score:theBest) {
                        changeSynBuild.put(altTaxon.getName(),score);
                    }
                    unmatchCount++;
                }
            } 
        }
        ImmutableMultimap<String,String> changeSyn = changeSynBuild.build();
        //Then go through and build a new taxaList adding the Synonyms to this
        for(TaxaList altTaxaList:alternateTaxaSets) {
            //Create a taxaList builder to make the new list
            TaxaListBuilder tlb = new TaxaListBuilder();
            
            //Go through previous list(Taxon level)
            for(Taxon taxon:altTaxaList) {
                //Copy all previous datafrom altTaxaList
                GeneralAnnotation ga = taxon.getAnnotation();
                Taxon.Builder tb = new Taxon.Builder(taxon.getName());
                Set<String> keys = ga.getAnnotationKeys();
                //Copy Keys
                for(String key:keys) {
                    if(!key.equals(Taxon.SynonymKey)) {
                        String[] values = ga.getTextAnnotation(key);
                        for(String value:values) {
                            tb.addAnno(key, value);
                        }
                    }
                }
                //If changeSyn has key for a given Taxon
                if(changeSyn.keySet().contains(taxon.getName())) {
                    //Add synonym to Taxon object
                    for(String entry:changeSyn.get(taxon.getName())) {
                        tb.addAnno(Taxon.SynonymKey,entry);
                    }
                }
                //Build Taxon and Add Taxon to new list builder tlb
                tlb.add(tb.build());
            }
            //Build tlb and Add to taxaListToBeSyn
            taxaListSynonymized.add(tlb.build());
        }
        resetTempTaxaList();
        System.out.println(toString());
    }
   
    public TaxaList getTaxaList() {
        return taxaListSynonymized.get(0);
    }

    public int getTechnique() {
        return technique;
    }
    public void setGlobalMax(double max) {
        this.globalMax = max;
    }
    private ArrayList<String> findBestMatch(String unmatchedString) {
        ArrayList<String> bestMatches = new ArrayList<>();
        double maxScore = -1;
        double minScore = Double.POSITIVE_INFINITY;
        double sm;
        int levelOfRestriction = 0;
        boolean ignoreCase = true, ignoreWhite = false, ignorePunc = false;
        while ((bestMatches.size() != 1) && (levelOfRestriction < 4)) {
            switch (levelOfRestriction) {
                case 1:
                    ignoreCase = true;
                    break;
                case 2:
                    ignoreWhite = true;
                    break;
                case 3:
                    ignorePunc = true;
                    break;
            }
            /*
            for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
                sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique);
                //sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
                if (sm > maxScore) {
                    bestMatches.clear();
                    bestMatches.add(referenceIDGroup.taxaName(i));
                    maxScore = sm;
                } else if (sm == maxScore) {
                    bestMatches.add(referenceIDGroup.taxaName(i));
                }
            }*/
            
            for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
                if(technique==7) {
                    sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique,delimiter);
                    
                }
                else {
                    sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique);
                }
                
                //sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
                if (sm < minScore) {
                    bestMatches.clear();
                    bestMatches.add(referenceIDGroup.taxaName(i));
                    minScore = sm;
                    if(minScore<globalMin) {
                        System.out.println(minScore);
                        globalMin = minScore;
                    } 
                } else if (sm == minScore) {
                    bestMatches.add(referenceIDGroup.taxaName(i));
                }
            }
            if(minScore>globalMax) {
                globalMax = minScore;
            }
            
            levelOfRestriction++;
        }
        return bestMatches;
    }
    public ArrayList<String> findBestMatch(String taxaName, TaxaList referenceTaxa) {
        ArrayList<String> bestMatches = new ArrayList<>();
        double maxScore = -1;
        double minScore = Double.POSITIVE_INFINITY;
        double sm;
        int levelOfRestriction = 0;
        boolean ignoreCase = true, ignoreWhite = false, ignorePunc = false;
        while ((bestMatches.size() != 1) && (levelOfRestriction < 4)) {
            switch (levelOfRestriction) {
                case 1:
                    ignoreCase = true;
                    break;
                case 2:
                    ignoreWhite = true;
                    break;
                case 3:
                    ignorePunc = true;
                    break;
            }
            for(Taxon refTaxa:referenceTaxa) {
                if(technique==7) {
                    sm = getScore(refTaxa.getName(),taxaName,ignoreCase,ignoreWhite,ignorePunc,technique,delimiter);
                }
                else {
                    sm = getScore(refTaxa.getName(), taxaName, ignoreCase, ignoreWhite, ignorePunc,technique);
                }
                if (sm < minScore) {
                    bestMatches.clear();
                    bestMatches.add(refTaxa.getName());
                    minScore = sm;
                    if(minScore<globalMin) {
                        globalMin = minScore;
                    } 
                } else if (sm == minScore) {
                    if(!bestMatches.contains(refTaxa.getName())) {
                        bestMatches.add(refTaxa.getName());
                    }
                }
            }
            if(minScore>globalMax) {
                globalMax = minScore;
            }
            
            levelOfRestriction++;
        }   
        System.out.println("GlobalMin"+globalMin);
        System.out.println("GlobalMax"+globalMax);
        return bestMatches;
    }
    public ArrayList<String> findOrderedMatches(String unmatchedString, int levelOfRestriction) {
        SortedMap<Double,String> theSortMap = new TreeMap<>();
        double sm;
        boolean ignoreCase = false, ignoreWhite = false, ignorePunc = false;
        if (levelOfRestriction > 0) {
            ignoreCase = true;
        }
        if (levelOfRestriction > 1) {
            ignoreWhite = true;
        }
        if (levelOfRestriction > 2) {
            ignorePunc = true;
        }
        for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
            //sm = scoreMatch(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc);
            if(technique==7) {
                sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique,delimiter);
            }
            else {
                sm = getScore(referenceIDGroup.taxaName(i), unmatchedString, ignoreCase, ignoreWhite, ignorePunc,technique);
            }
            sm = 1.0-((sm - globalMin)/(globalMax-globalMin));
            theSortMap.put(1 - sm - ((double) i / 100000.0), referenceIDGroup.taxaName(i));
            //theSortMap.put(sm - ((double) i / 100000.0), referenceIDGroup.taxaName(i));
        }
        return new ArrayList<>(theSortMap.values());
    }

    public static double getScore(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc, int technique) {
        double score = 0.0;
        if(s1.equals(s2)) {
            return score;
        }
      
        //dice need to do a 1- as high similarity = low distance
        if(technique == 0) {
            score = 1.0 - scoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //String edit
        else if(technique == 1) {
            score = editDistanceScoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //DTW with hamming
        else if(technique == 2) {
            score = dtwDist(s1,s2,"hamming",ignoreCase,true,ignorePunc);
        }
        //DTW with keyboard dist
        else if(technique == 3) {
            score = dtwDist(s1,s2,"key",ignoreCase,true,ignorePunc);
        }
        //Hamming with soundex
        else if(technique == 4) {
            score = hammingDistSoundex(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //Dice with metaphone  need to do a 1- as high similarity = low distance
        else if(technique == 5) {
            score = 1 - diceWithMetaphone(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //Edit Distance with metaphone
        else if(technique == 6) {
            score=editWithMetaphone(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        
        return score;
    }
    public static double getScore(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc, int technique,String delimiter) {
        double score = 0.0;
        if(s1.equals(s2)) {
            return score;
        }
      
        //dice need to do a 1- as high similarity = low distance
        if(technique == 0) {
            score = 1.0 - scoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //String edit
        else if(technique == 1) {
            score = editDistanceScoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //DTW with hamming
        else if(technique == 2) {
            score = dtwDist(s1,s2,"hamming",ignoreCase,true,ignorePunc);
        }
        //DTW with keyboard dist
        else if(technique == 3) {
            score = dtwDist(s1,s2,"key",ignoreCase,true,ignorePunc);
        }
        //Hamming with soundex
        else if(technique == 4) {
            score = hammingDistSoundex(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //Dice with metaphone  need to do a 1- as high similarity = low distance
        else if(technique == 5) {
            score = 1 - diceWithMetaphone(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        //Edit Distance with metaphone
        else if(technique == 6) {
            score=editWithMetaphone(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
        }
        else if(technique == 7) {
            score= 1- delimiterDistance(s1,s2,ignoreCase,ignoreWhite,ignorePunc,delimiter);
        }
        return score;
    }
    
    public static double hammingDistSoundex(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = soundex2(s1, true, true, true);
        s2 = soundex2(s2, true, true, true);
        int sum = 0;
        for(int i = 0; i<s1.length();i++) {
            sum += hammingDist(s1.charAt(i), s2.charAt(i));
        }
        return (double)sum;
    }
    public static double diceWithMetaphone(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = metaphone2(s1,true,true,true);
        s2 = metaphone2(s2,true,true,true);
        return scoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
    }
    
    public static double editWithMetaphone(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = metaphone2(s1,true,true,true);
        s2 = metaphone2(s2,true,true,true);
        return editDistanceScoreMatch(s1,s2,ignoreCase,ignoreWhite,ignorePunc);
    }
    private double scoreMatch2(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        //idea from http://www.catalysoft.com/articles/StrikeAMatch.html?article=How_to_Strike_a_Match_15
        //this is faster but it can be tricked if there are long runs of characters in s1
        int score = 0;
        double sm;
        s1 = cleanName(s1, ignoreCase, ignoreWhite, ignorePunc);
        s2 = cleanName(s2, ignoreCase, ignoreWhite, ignorePunc);
//        System.out.println("s1="+s1+"  s2="+s2);
        for (int c1 = 0; c1 < (s1.length() - 1); c1++) {
            for (int c2 = 0; c2 < (s2.length() - 1); c2++) {
                if ((s1.charAt(c1) == s2.charAt(c2)) && (s1.charAt(c1 + 1) == s2.charAt(c2 + 1))) {
                    score++;
                    break;
                }
            }
        }
        sm = (2.0 * (double) score) / (s1.length() + s2.length() - 2);
        return sm;
    }

    /** @return lexical similarity value in the range [0,1] */
    public static double scoreMatch(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        //idea from http://www.catalysoft.com/articles/StrikeAMatch.html?article=How_to_Strike_a_Match_15
        //this is slower but it will not be tricked if there are long runs of characters in s1
        s1 = cleanName(s1, ignoreCase, ignoreWhite, ignorePunc);
        s2 = cleanName(s2, ignoreCase, ignoreWhite, ignorePunc);
        ArrayList<String> pairs1 = letterPairs(s1);
        ArrayList<String> pairs2 = letterPairs(s2);

        int intersection = 0;
        int union = pairs1.size() + pairs2.size();
        for (int i = 0; i < pairs1.size(); i++) {
            Object pair1 = pairs1.get(i);
            for (int j = 0; j < pairs2.size(); j++) {
                Object pair2 = pairs2.get(j);
                if (pair1.equals(pair2)) {
                    intersection++;
                    pairs2.remove(j);
                    break;
                }
            }
        }
        return (2.0 * intersection) / union;
    }
    
    public static double editDistanceScoreMatch(String s1, String s2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        
        s1 = cleanName(s1, ignoreCase, ignoreWhite, ignorePunc);
        s2 = cleanName(s2, ignoreCase, ignoreWhite, ignorePunc);
        
        if(s1.equals("")) {
            return s2.length();
        }
        if(s2.equals("")) {
            return s1.length();
        }
        
        double[][] editMatrix = new double[s1.length()][s2.length()];
        
        //Init first row and column of editMatrix
        editMatrix[0][0] = (s1.charAt(0) == s2.charAt(0))?0.0:1.0;
        for(int i = 1; i<editMatrix.length; i++) {
            editMatrix[i][0] = (s1.charAt(i)==s2.charAt(0))?editMatrix[i-1][0]: editMatrix[i-1][0]+1;
        }
        for(int i = 1; i<editMatrix[0].length;i++) {
            editMatrix[0][i] = (s1.charAt(0) == s2.charAt(i))? editMatrix[0][i-1]:editMatrix[0][i-1]+1;
        }
        
        //Fill in the rest of the matrix
        for(int i = 1; i < editMatrix.length; i++) {
            for(int j = 1; j< editMatrix[i].length;j++) {
                if(s1.charAt(i) == s2.charAt(j)) {
                    editMatrix[i][j] = editMatrix[i-1][j-1];
                }
                else {
                    double diagCost = editMatrix[i-1][j-1];
                    double upCost = editMatrix[i-1][j];
                    double leftCost = editMatrix[i][j-1];
                    
                    if(diagCost <= upCost && diagCost <= leftCost) {
                        //Substitution
                        editMatrix[i][j] = diagCost + 1;
                    }
                    else if(upCost <= leftCost && upCost <= diagCost) {
                        //Deletion
                        editMatrix[i][j] = upCost + 1;
                    }
                    else {
                        //Insertion
                        editMatrix[i][j] = leftCost + 1;
                    }
                }
            }
        }
        return editMatrix[editMatrix.length-1][editMatrix[editMatrix.length-1].length-1];
    }

    public static String metaphone2(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        if(s1.equals("")) {
            return "";
        }
        //Parse out numbers i.e. ab12cd3e becomes ["ab","12","cd","3","e"] then metaphone all letter only strings
        ArrayList<String> parsed = new ArrayList<String>();
        String current = "";
        boolean digitMode = false;
        if(Character.isDigit(s1.charAt(0))) {
            current+=s1.charAt(0);
            digitMode = true;
        }
        for(int i = 0; i<s1.length();i++) {
            if(Character.isDigit(s1.charAt(i)) && !digitMode) {
                parsed.add(current);
                current = ""+s1.charAt(i);
                digitMode = true;
            }
            else if(!Character.isDigit(s1.charAt(i)) && digitMode) {
                parsed.add(current);
                current = ""+s1.charAt(i);
                digitMode = false;
            }
            else {
                current += s1.charAt(i);
            }
        }
        parsed.add(current);
        Metaphone metaphone = new Metaphone();
        String encodedString = "";
        for(int i = 0; i<parsed.size();i++) {
            if(!Character.isDigit(parsed.get(i).charAt(0))) {
               encodedString += metaphone.encode(parsed.get(i)); 
            }
            else {
                encodedString+=parsed.get(i);
            }
        }
        return encodedString;
        //return metaphone.metaphone(s1);
    }
    public static String soundex2(String s1, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        s1 = cleanName(s1, true, true, true);
        if(s1.equals("")) {
            return "1abcd";
        }
        Soundex soundex = new Soundex();
        return soundex.soundex(s1);
    }
    
    /*
     * Method to determine keyboard distance.
     * Idea borrowed from http://search.cpan.org/~krburton/String-KeyboardDistance-1.01/KeyboardDistance.pm
     * This needs refractoring it is messy
     */
    private static int keyboardDist(char firstChar, char secondChar) {
        /*    | 0   1   2   3   4   5   6   7   8   9   10  11  12  13
            --+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
            0 | ` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 0 | - | = |   |
            1 |   | q | w | e | r | t | y | u | i | o | p | [ | ] | \ |
            2 |   | a | s | d | f | g | h | j | k | l | ; | ' |   |   |
            3 |   | z | x | c | v | b | n | m | , | . | / |   |   |   |
            --+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
        */
        HashMap<Character,Integer[]> map = new HashMap<>();
        map.put('`', new Integer[]{0,0});
        map.put('~', new Integer[]{0,0});
        map.put('1', new Integer[]{0,1});
        map.put('!', new Integer[]{0,1});
        map.put('2', new Integer[]{0,2});
        map.put('@', new Integer[]{0,2});
        map.put('3', new Integer[]{0,3});
        map.put('#', new Integer[]{0,3});
        map.put('4', new Integer[]{0,4});
        map.put('$', new Integer[]{0,4});
        map.put('5', new Integer[]{0,5});
        map.put('%', new Integer[]{0,5});
        map.put('6', new Integer[]{0,6});
        map.put('^', new Integer[]{0,6});
        map.put('7', new Integer[]{0,7});
        map.put('&', new Integer[]{0,7});
        map.put('8', new Integer[]{0,8});
        map.put('*', new Integer[]{0,8});
        map.put('9', new Integer[]{0,9});
        map.put('(', new Integer[]{0,9});
        map.put('0', new Integer[]{0,10});
        map.put(')', new Integer[]{0,10});
        map.put('-', new Integer[]{0,11});
        map.put('_', new Integer[]{0,11});
        map.put('=', new Integer[]{0,12});
        map.put('+', new Integer[]{0,12});
        
        map.put('q', new Integer[]{1,1});
        map.put('Q', new Integer[]{1,1});
        map.put('w', new Integer[]{1,2});
        map.put('W', new Integer[]{1,2});
        map.put('e', new Integer[]{1,3});
        map.put('E', new Integer[]{1,3});
        map.put('r', new Integer[]{1,4});
        map.put('R', new Integer[]{1,4});
        map.put('t', new Integer[]{1,5});
        map.put('T', new Integer[]{1,5});
        map.put('y', new Integer[]{1,6});
        map.put('Y', new Integer[]{1,6});
        map.put('u', new Integer[]{1,7});
        map.put('U', new Integer[]{1,7});
        map.put('i', new Integer[]{1,8});
        map.put('I', new Integer[]{1,8});
        map.put('o', new Integer[]{1,9});
        map.put('O', new Integer[]{1,9});
        map.put('p', new Integer[]{1,10});
        map.put('P', new Integer[]{1,10});
        map.put('[', new Integer[]{1,11});
        map.put('{', new Integer[]{1,11});
        map.put(']', new Integer[]{1,12});
        map.put('}', new Integer[]{1,12});
        map.put('\\', new Integer[]{1,13});
        map.put('|', new Integer[]{1,13});
        
        map.put('a', new Integer[]{2,1});
        map.put('A', new Integer[]{2,1});
        map.put('s', new Integer[]{2,2});
        map.put('S', new Integer[]{2,2});
        map.put('d', new Integer[]{2,3});
        map.put('D', new Integer[]{2,3});
        map.put('f', new Integer[]{2,4});
        map.put('F', new Integer[]{2,4});
        map.put('g', new Integer[]{2,5});
        map.put('G', new Integer[]{2,5});
        map.put('h', new Integer[]{2,6});
        map.put('H', new Integer[]{2,6});
        map.put('j', new Integer[]{2,7});
        map.put('J', new Integer[]{2,7});
        map.put('k', new Integer[]{2,8});
        map.put('K', new Integer[]{2,8});
        map.put('l', new Integer[]{2,9});
        map.put('L', new Integer[]{2,9});
        map.put(';', new Integer[]{2,10});
        map.put(':', new Integer[]{2,10});
        map.put('\'', new Integer[]{2,11});
        map.put('"', new Integer[]{2,11});
        
        map.put('z', new Integer[]{3,1});
        map.put('Z', new Integer[]{3,1});
        map.put('x', new Integer[]{3,2});
        map.put('X', new Integer[]{3,2});
        map.put('c', new Integer[]{3,3});
        map.put('C', new Integer[]{3,3});
        map.put('v', new Integer[]{3,4});
        map.put('V', new Integer[]{3,4});
        map.put('b', new Integer[]{3,5});
        map.put('B', new Integer[]{3,5});
        map.put('n', new Integer[]{3,6});
        map.put('N', new Integer[]{3,6});
        map.put('m', new Integer[]{3,7});
        map.put('M', new Integer[]{3,7});
        map.put(',', new Integer[]{3,8});
        map.put('<', new Integer[]{3,8});
        map.put('.', new Integer[]{3,9});
        map.put('>', new Integer[]{3,9});
        map.put('/', new Integer[]{3,10});
        map.put('?', new Integer[]{3,10});
        
       Integer[] coords1 = map.get(firstChar);
       Integer[] coords2 = map.get(secondChar);
       
       //calculate manhattan distance between the characters
       return Math.abs(coords1[0] - coords2[0]) + Math.abs(coords1[1] - coords2[1]);
        
    }
    
    private static int hammingDist(char firstChar, char secondChar) {
        if(firstChar == secondChar) {
            return 0;
        }
        else{
            return 1;
        }   
    }
   
    /*
     * Simple implementation of Dynamic Time Warping distance
     * 
     * Currently uses KeyboardDistance as the distance measurement
     */
    private static double dtwDist(String str1, String str2, String distMeas,boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        str1 = cleanName(str1,ignoreCase,ignoreWhite,ignorePunc);
        str2 = cleanName(str2,ignoreCase,ignoreWhite,ignorePunc);
        double[][] costMat = new double[str1.length()+1][str2.length()+1];
        //Initialize arrays
        for(int i = 0; i<costMat.length;i++) {
            costMat[i][0] = Double.POSITIVE_INFINITY;
        }
        for(int i = 0; i<costMat[0].length;i++) {
            costMat[0][i] = Double.POSITIVE_INFINITY;
        }
        
        double currDist = 0.0;
        if(distMeas.equals("key")) {
            costMat[1][1] = (double)keyboardDist(str1.charAt(0),str2.charAt(0));
        }
        else if(distMeas.equals("hamming")) {
            costMat[1][1] = (double)hammingDist(str1.charAt(0), str2.charAt(0));
        }
        for(int i = 2;i<costMat[1].length;i++) {
            if(distMeas.equals("key")) {
                costMat[1][i] = costMat[1][i-1]+(double)keyboardDist(str1.charAt(0),str2.charAt(i-1));
            }
            else if(distMeas.equals("hamming")) {
                costMat[1][i] = costMat[1][i-1]+(double)hammingDist(str1.charAt(0), str2.charAt(i-1));
            }        
        }
        for(int i = 2; i<costMat.length;i++) {
            if(distMeas.equals("key")) {
                costMat[i][1] = costMat[i-1][1]+(double)keyboardDist(str1.charAt(i-1),str2.charAt(0));
            }
            else if(distMeas.equals("hamming")) {
                costMat[i][1] = costMat[i-1][1]+(double)hammingDist(str1.charAt(i-1), str2.charAt(0));
            }
        }
        for(int i = 2; i<costMat.length;i++) {
            for(int j = 2; j<costMat[i].length;j++) {
                currDist = 0.0;
                //get distance
                if(distMeas.equals("key")) {
                    currDist = (double)keyboardDist(str1.charAt(i-1),str2.charAt(j-1));
                    //System.out.println(currDist);
                }
                else if(distMeas.equals("hamming")) {
                    currDist = (double)hammingDist(str1.charAt(i-1), str2.charAt(j-1));
                }
                
                if(costMat[i-1][j-1] < costMat[i-1][j] && costMat[i-1][j-1] < costMat[i][j-1]) {
                    costMat[i][j] = costMat[i-1][j-1] + currDist;
                }
                else if(costMat[i-1][j]<costMat[i-1][j-1] && costMat[i-1][j] < costMat[i][j-1]) {
                    costMat[i][j] = costMat[i-1][j] + currDist;
                }
                else {
                    costMat[i][j] = costMat[i][j-1] + currDist;
                }
            }
        }
        //return end point
        return costMat[costMat.length-1][costMat[costMat.length-1].length-1];
    }
    
    private static double subSetDist(String str1, String str2,boolean ignoreCase,boolean ignoreWhite, boolean ignorePunc) {
        str1 = cleanName(str1,ignoreCase,ignoreWhite,ignorePunc);
        str2 = cleanName(str2,ignoreCase,ignoreWhite,ignorePunc);
        
        double distance = 0.0;
        int length = (str1.length()>str2.length())?str2.length():str1.length();
        int errorCount = 0;
        for(int i = 0; i<length;i++) {
            if(str1.charAt(i)!=str2.charAt(i)) {
                errorCount++;
            }
        }
        distance = (double)errorCount/length;
        return distance;
    }
    
    private static double delimiterDistance(String str1, String str2, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc,String delimiter) {
        double distance = 0.0;
       
        String[] str1Split = str1.split(delimiter);
        String[] str2Split = str2.split(delimiter);
        
        if(str1Split.length == str2Split.length) {
            for(int i = 0; i<str1Split.length;i++) {
                distance += scoreMatch(str1Split[i],str2Split[i],ignoreCase,ignoreWhite,ignorePunc);
            }
            distance/=str1Split.length;
            return distance;
        }
        //Make sure Str1Split is the longer of the two
        if(str1Split.length<str2Split.length) {
            String[] temp = str1Split;
            str1Split = str2Split;
            str2Split = temp;
        }
        //Use a sliding window
        for(int i = 0; i<str1Split.length-str2Split.length;i++) {
            double currDist = 0.0;
            for(int j = 0; j<str2Split.length;j++) {
                currDist += scoreMatch(str1Split[i+j],str2Split[j],ignoreCase,ignoreWhite,ignorePunc);
            }
            currDist/=str2Split.length;
            if(currDist>distance) {
                distance = currDist;
            }
        }
        
        return distance;
    }
    
    /** @return an array of adjacent letter pairs contained in the input string */
    private static ArrayList<String> letterPairs(String str) {
        ArrayList<String> allPairs = new ArrayList<>();
        //int numPairs = str.length()-1;
        //String[] pairs = new String[numPairs];
        for (int i = 0; i < (str.length() - 1); i++) {
            allPairs.add(str.substring(i, i + 2));
        }
        return allPairs;
    }

    private static String cleanName(String s, boolean ignoreCase, boolean ignoreWhite, boolean ignorePunc) {
        if (ignoreCase) {
            s = s.toUpperCase();
        }
        //StringBuffer sb=new StringBuffer(s);
        //int x;
        if (ignoreWhite) {
            s.replaceAll("\\s", "");
        // while((x=sb.indexOf(" "))>=0) {sb.deleteCharAt(x);}
        }
        if (ignorePunc) {
            //           s=s.replaceAll("\\W","");
            s = s.replaceAll("[^a-zA-Z0-9]", "");
        }
        // sb=new StringBuffer(s);
        return s;
    }

    public void changeAlignmentIdentifiers(TaxaList alternateIdGroups) {
        TaxaList[] aidg = new TaxaList[1];
        aidg[0] = alternateIdGroups;
        changeAlignmentIdentifiers(aidg[0]);
    }

    
    public String toString() {
        String s = "Synonym Table\n";
        int counter = 0;
        for(TaxaList tl:taxaListSynonymized) {
            s+="SynonymFile "+counter+":\n";
            for(Taxon tx:tl) {
                s+="Name: " + tx.getName()+", Synonyms: ";
                for(String syn:tx.getAnnotation().getTextAnnotation(Taxon.SynonymKey)) {
                    s+=syn+", ";
                }
                s+="\n";
            }
            counter++;
        }
        return s;    //To change body of overridden methods use File | Settings | File Templates.
    }

    public void deleteByThreshold(double threshold) {
        ImmutableMultimap.Builder<String,String> removeBuilder = new ImmutableMultimap.Builder();
        //Go through taxa list
        for(Taxon tx:tempTaxaList) {
            String taxonName = tx.getName();
            //Go through the list of synonyms
            for(String synName:tx.getAnnotation().getTextAnnotation(Taxon.SynonymKey)) {
                double score = 0.0;
                if(technique==7) {
                    score = getScore(taxonName,synName,true,false,false,technique,delimiter);
                }
                else {
                    score = getScore(taxonName,synName,true,false,false,technique);
                }
                //double score = getScore(taxonName,synName,true,false,false,technique);
                
                if(technique==4) {
                    globalMax = 4.0;
                }
                score = 1.0-((score - globalMin)/(globalMax-globalMin));
                //If Score is less than thresh add to remove list
                if(score<threshold) {
                    removeBuilder.put(taxonName,synName);
                }
            }
        }
        ImmutableMultimap<String,String> removeMap = removeBuilder.build();
        TaxaListBuilder tlb = new TaxaListBuilder();
        for(Taxon tx:tempTaxaList) {
            GeneralAnnotation ga = tx.getAnnotation();
            Taxon.Builder tb = new Taxon.Builder(tx.getName());
            Set<String> keys = ga.getAnnotationKeys();
            //Copy Keys
            for(String key:keys) {
                if(!key.equals(Taxon.SynonymKey)) {
                    String[] values = ga.getTextAnnotation(key);
                    for(String value:values) {
                        tb.addAnno(key, value);
                    }
                }
            }
            //If removeMap has key for a given Taxon
            if(removeMap.keySet().contains(tx.getName())) {
                //Loop through and ignore any that are in removeMap
                for(String value:ga.getTextAnnotation(Taxon.SynonymKey)) {
                    if(!removeMap.get(tx.getName()).contains(value)) {
                        tb.addAnno(Taxon.SynonymKey, value);
                    }
                }
            }
            else {
                //Loop through add add all the synonyms
                for(String value:ga.getTextAnnotation(Taxon.SynonymKey)) {
                    tb.addAnno(Taxon.SynonymKey, value);
                }
            }
            //Build Taxon and Add Taxon to new list builder tlb
            tlb.add(tb.build());
        }
        tempTaxaList = tlb.build();
    }

    public Object[] getRealNames() {
        Object[] idArray = new Object[referenceIDGroup.numberOfTaxa()];
        for (int i = 0; i < referenceIDGroup.numberOfTaxa(); i++) {
            idArray[i] = referenceIDGroup.get(i).toString();
        }
        return idArray;
    }

    public void report(PrintWriter out) {
        //String s="Synonym Table\n"+idSynonyms.toString()+"\n\n"+"Unmatched\n"+unmatchedIDs.toString();
        out.println("Synonym Table");
        //out.println(idSynonyms.size() + " unique matches");
        out.println(matchCount + " unique matches");
        out.println(unmatchCount + " multiple matches:");
    }

    public Object[] getTableColumnNames() {
        String[] cn = new String[3];
        cn[0] = "TaxaName";
        cn[1] = "Synonyms";
        cn[2] = "MatchScore";
        return cn;
    }

    
    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(long rowLong) {
        int row = (int) rowLong;
        Object[] data = new Object[3];
        
        //TaxaList tl = taxaListSynonymized.get(0);
        TaxaList tl = tempTaxaList;
        Taxon tx = tl.get(row);
        
        //Object[] keyArray = idSynonyms.keySet().toArray();
        data[0] = tx.getName();
        String synString = "";
        String firstMatch = "";
        boolean first = true;
        for(String syn:tx.getAnnotation().getTextAnnotation(Taxon.SynonymKey)){
            synString+=syn+",";
            if(first) {
                firstMatch = syn;
                first = false;
            }
        }
        
        data[2] = "";
        if(firstMatch.equals("")) {
            data[1] = "";
            data[2] = "1.0";
        }
        else {
            data[1] = synString.substring(0,synString.length()-1);
            if(technique == 4) {
                globalMax = 4.0;
            }
            if(technique==0) {
                data[2] = "" + scoreMatch("" + data[0], "" + firstMatch, true, false, false);
            }
            else if(technique==5) {
                data[2] = "" + (1.0-getScore("" + data[0], "" + firstMatch, true, false, false,technique));  
            }
            else if(technique==7) {
                data[2] = "" + (1.0-getScore("" + data[0], "" + firstMatch, true, false, false,technique,delimiter));  
            }
            else {
                //To fix the - number bug
                if((1.0-((getScore("" + data[0], "" + firstMatch, true, false, false,technique) - globalMin)/(globalMax-globalMin)))<0.0) {
                    if((1.0-((getScore("" + data[0], "" + firstMatch, true, true, false,technique) - globalMin)/(globalMax-globalMin)))<0.0) {
                        if ((1.0-((getScore("" + data[0], "" + firstMatch, true, true, true,technique) - globalMin)/(globalMax-globalMin)))<0.0) {
                            data[2] = ""+0.0;
                        }
                        else {
                            data[2] = "" + (1.0-((getScore("" + data[0], "" + firstMatch, true, true, true,technique) - globalMin)/(globalMax-globalMin)));
                        }
                    }
                    else {
                        data[2] = "" + (1.0-((getScore("" + data[0], "" + firstMatch, true, true, false,technique) - globalMin)/(globalMax-globalMin)));
                    }
                }
                else {
                    data[2] = "" + (1.0-((getScore("" + data[0], "" + firstMatch, true, false, false,technique) - globalMin)/(globalMax-globalMin)));
                }
            }
        }
        return data;
    }
    
    

    public void resetTempTaxaList() {
        TaxaListBuilder tlb = new TaxaListBuilder();
        for(Taxon tx:taxaListSynonymized.get(0)) {
            tlb.add(tx);
        }
        tempTaxaList = tlb.build();
    }
    //Method to save changes from TempTaxaList and reset the TempTaxaList
    public void saveTempTaxaList() {
       taxaListSynonymized.set(0, tempTaxaList);
       resetTempTaxaList();
    }

    public void removeSynonyms(int rowNumber) {
        TaxaListBuilder tlb = new TaxaListBuilder();
        int rowCounter = 0;
        for(Taxon tx:tempTaxaList) {
            if(rowNumber == rowCounter) {
                //Copy all annotations except for Synonyms
                GeneralAnnotation ga = tx.getAnnotation();
                Taxon.Builder tb = new Taxon.Builder(tx.getName());
                Set<String> keys = ga.getAnnotationKeys();
                //Copy Keys
                for(String key:keys) {
                    if(!key.equals(Taxon.SynonymKey)) {
                        String[] values = ga.getTextAnnotation(key);
                        for(String value:values) {
                            tb.addAnno(key, value);
                        }
                    }
                }
                tlb.add(tb.build());
            }
            else {
                //Copy taxon
                tlb.add(tx);
            }
            rowCounter++;
        }
        tempTaxaList = tlb.build();
    }
    public void updateSynonym(int rowNumber, String newName) {
        TaxaListBuilder tlb = new TaxaListBuilder();
        int rowCounter = 0;
        for(Taxon tx:tempTaxaList) {
            if(rowNumber == rowCounter) {
                //Copy all annotations except for Synonyms
                GeneralAnnotation ga = tx.getAnnotation();
                Taxon.Builder tb = new Taxon.Builder(tx.getName());
                Set<String> keys = ga.getAnnotationKeys();
                //Copy Keys
                for(String key:keys) {
                    if(!key.equals(Taxon.SynonymKey)) {
                        String[] values = ga.getTextAnnotation(key);
                        for(String value:values) {
                            tb.addAnno(key, value);
                        }
                    }  
                }
                tb.addAnno(Taxon.SynonymKey,newName);
                tlb.add(tb.build());
            }
            else {
                //Copy taxon
                tlb.add(tx);
            }
            rowCounter++;
        }
        tempTaxaList = tlb.build();
    }
    
    
    public void deleteElements(String name) {
        TaxaListBuilder tlb = new TaxaListBuilder();
        for(Taxon tx:tempTaxaList) {
            if(!tx.getName().equals(name)) {
                tlb.add(tx);
            }
        }
        tempTaxaList = tlb.build();
    }

    public boolean checkSynForDups() {
        for(TaxaList tl:taxaListSynonymized) {
            ArrayList<String> viewedTaxa = new ArrayList<String>();
            for(Taxon tx:tl) {
                GeneralAnnotation ga = tx.getAnnotation();
                String[] values = ga.getTextAnnotation(Taxon.SynonymKey);
                
                if(values.length==0||values==null){
                    viewedTaxa.add(tx.getName());
                }
                else if(viewedTaxa.contains(values[0])) {
                    return true;
                }
                else {
                    viewedTaxa.add(values[0]);
                }
            }
                
        }
        return false;
    }
    
    public ArrayList<TaxaList> swapSynonyms() {
        ArrayList<TaxaList> newTaxaListArray = new ArrayList<>();
        for(TaxaList tl:taxaListSynonymized) {
            TaxaListBuilder tlb = new TaxaListBuilder();
            for(Taxon tx:tl) {
                GeneralAnnotation ga = tx.getAnnotation();
                Taxon.Builder tb = null;
                Set<String> keys = ga.getAnnotationKeys();
                if(ga.getTextAnnotation(Taxon.SynonymKey).length==0) {
                    tb = new Taxon.Builder(tx);
                }
                else {
                    tb = new Taxon.Builder(ga.getTextAnnotation(Taxon.SynonymKey)[0]);
                    //Copy Keys
                    for(String key:keys) {
                        if(!key.equals(Taxon.SynonymKey)) {
                            String[] values = ga.getTextAnnotation(key);
                            for(String value:values) {
                                tb.addAnno(key, value);
                            }
                        }
                    }
                    String[] synVals = ga.getTextAnnotation(Taxon.SynonymKey);
                    tb.addAnno(Taxon.SynonymKey, tx.getName());
                    for(int i = 1; i<synVals.length;i++) {
                        tb.addAnno(Taxon.SynonymKey, synVals[i]);
                    }
                }
                tlb.add(tb.build());
            }
            newTaxaListArray.add(tlb.build());
        }
        return newTaxaListArray;
    }
    
    //Table methods
    public String getTableTitle() {
        return "Taxa Synonym Table";
    }

    // interface ExtendedTableReport
    public int getColumnCount() {
        return 3;
    }

    public long getRowCount() {
        return tempTaxaList.size();
        //return idSynonyms.size();
    }

    public long getElementCount() {
        return getColumnCount() * getRowCount();
    }
}
