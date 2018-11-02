package net.maizegenetics.taxa;

public class TaxaTissueDist {
    private final int maxTaxa;
    private final int maxTissue;
    int[][] tissueTaxa;       
    public TaxaTissueDist(int mtissue, int mtaxa) {
        maxTaxa = mtaxa;
        maxTissue = mtissue;
        tissueTaxa = new int[maxTissue][maxTaxa];
    }
    
   public synchronized void increment(int tissueNum, int taxaNum) {
       tissueTaxa[tissueNum][taxaNum]++;
   }

   public synchronized int getTissueDepths(int tissueNum, int taxaNum) {
       return tissueTaxa[tissueNum][taxaNum];
   }
}
