package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;

/**
 * Basic implementation of allele.
 * TODO: Consider whether this should be kept.
 * TODO: Should there be a builder?
 * TODO: Should there simple nucleotide one versus other states.  Memory could be an important consideration.
 *
 *
 * @author Ed Buckler
 */
public class SimpleAllele implements Allele {
    private final byte myAllele;
    private final Position myPosition;
    private final GeneralAnnotationStorage myAnnotations;

    public SimpleAllele(byte myAllele, Position myPosition) {
        this.myAllele=myAllele;
        this.myPosition=myPosition;
        myAnnotations= GeneralAnnotationStorage.getBuilder().build();
    }

    @Override
    public Position position() {
        return myPosition;
    }

    @Override
    public byte allele() {
        return myAllele;
    }

    @Override
    public String alleleAsString() {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(myAllele);
    }

    @Override
    public GeneralAnnotation annotations() {
        return myAnnotations;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SimpleAllele that = (SimpleAllele) o;

        if (myAllele != that.myAllele) return false;
        if (!myPosition.equals(that.myPosition)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = (int) myAllele;
        result = 31 * result + myPosition.hashCode();
        return result;
    }
}
