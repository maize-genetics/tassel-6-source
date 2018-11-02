package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.util.GeneralAnnotation;

/**
 * Defines Allele present at a genomic Position
 *
 * @author Ed Buckler
 */
public interface Allele {

    /*
    Returns Allele state
     */
    byte allele();

    /*
    Returns Allele state as String
     */
    String alleleAsString();


    /*Returns of the Position associated with the allele*/
    Position position();

    /*Returns the annotations for allele*/
    GeneralAnnotation annotations();

}
