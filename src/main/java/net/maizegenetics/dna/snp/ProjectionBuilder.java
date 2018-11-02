package net.maizegenetics.dna.snp;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.genotypecall.ProjectionGenotypeCallTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.dna.map.DonorHaplotypes;

import java.util.Collection;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * Builder for a projection alignment.  Projection alignments use defined haplotypes and breakpoints
 * that point to high density genotypes (base GenotypeTable).  These are used to efficiently store
 * and connect low density maps with imputed high density genotypes.
 * <p></p>
 * The alignment built by this builder is a CoreGenotypeTable with a ProjectionGenotypeCallTable.  The
 * taxa come from the projection alignment file, while the sites and positions are the
 * same as the base alignment.
 *
 * @author Ed Buckler
 */
public class ProjectionBuilder {
    private final GenotypeTable myBaseAlignment;  //high density marker alignment that is being projected.
    private ImmutableMap.Builder<Taxon,NavigableSet<DonorHaplotypes>> allBreakPoints;

    public static GenotypeTable getInstance(GenotypeTable baseAlignment, ImmutableMap<Taxon, NavigableSet<DonorHaplotypes>> allBreakPoints) {
        TaxaList tl=new TaxaListBuilder().addAll(allBreakPoints.keySet()).build();
        ImmutableList breakList=ImmutableList.builder().addAll(allBreakPoints.values()).build();
        return GenotypeTableBuilder.getInstance(new ProjectionGenotypeCallTable(baseAlignment, breakList),
                baseAlignment.positions(),tl);
    }

    public ProjectionBuilder(GenotypeTable myBaseAlignment) {
        this.myBaseAlignment=myBaseAlignment;
        allBreakPoints=new ImmutableMap.Builder<>();
    }

    public synchronized ProjectionBuilder addTaxon(Taxon taxon, Map<Position,Taxon[]> breakPoints) {
        NavigableSet<DonorHaplotypes> intBreak=convertToIndexBreakPoints(breakPoints);
        allBreakPoints.put(taxon,intBreak);
        return this;
    }

    public synchronized ProjectionBuilder addTaxon(Taxon taxon, NavigableSet<DonorHaplotypes> breakPoints) {
        allBreakPoints.put(taxon,breakPoints);
        return this;
    }

    public synchronized ProjectionBuilder addTaxon(Taxon taxon, Collection<DonorHaplotypes> breakPoints) {
        NavigableSet<DonorHaplotypes> tBreak=new TreeSet<>();
        for (DonorHaplotypes breakPoint : breakPoints) {
            tBreak.add(breakPoint);
        }
        return addTaxon(taxon,tBreak);
    }

    public GenotypeTable build() {
        ImmutableMap<Taxon,NavigableSet<DonorHaplotypes>> immBreak=allBreakPoints.build();
        System.out.println(immBreak.size());
        return getInstance(myBaseAlignment,immBreak);
    }

    /**
     * Conversion and validation routine.  Converts Taxon to indices in the BaseAlignment, and it removes redundant breaks.
     * This also ensures the resulting map is unconnected with the outside (essentially a defensive copy)
     * @param breakPoints
     * @return
     */
    private NavigableSet<DonorHaplotypes> convertToIndexBreakPoints(Map<Position, Taxon[]> breakPoints) {
        DonorHaplotypes lastP=new DonorHaplotypes(null, -1, -1, -1, -1);
        NavigableSet<DonorHaplotypes> intBreak=new TreeSet<>();
        TaxaList tl=myBaseAlignment.taxa();
        for (Map.Entry<Position,Taxon[]> bp : breakPoints.entrySet()) {
            Taxon[] ts=bp.getValue();
            if(ts.length!=2) throw new IllegalArgumentException("Two parents required for DonorHaplotypes");
            int[] iT=new int[ts.length];
            for (int i=0; i<iT.length; i++) {
                int r=tl.indexOf(ts[i]);
                if(r>=0) {iT[i]=r;}
                else {throw new IllegalArgumentException("Taxa not found or duplicated:"+ts[i].getName());}
            }
            DonorHaplotypes dh=new DonorHaplotypes(bp.getKey().getChromosome(), bp.getKey().getPosition(),
                    Integer.MAX_VALUE, iT[0], iT[1]);
            //TODO need to set the end Position based on next start position
            if(lastP!=dh) {
                intBreak.add(dh);
                lastP=dh;
            }
        }
        return intBreak;
    }

}
