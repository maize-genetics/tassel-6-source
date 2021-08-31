// Taxon.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.taxa;

import net.maizegenetics.util.GeneralAnnotation;
import net.maizegenetics.util.GeneralAnnotationStorage;

import java.util.HashMap;
import java.util.Map;

/**
 * An identifier for some sampled data. This will most often be for example, the
 * accession number of a DNA sequence, or the taxonomic name that the sequence
 * represents, et cetera.
 *
 * The generally used class for defining a taxon. Contains its name, plus a
 * series of optional annotations. Use the builder to create immutable
 * instances.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class Taxon implements Comparable<Taxon> {

    private static final Map<String, Taxon> TAXON_CACHE = new HashMap();

    public static final String DELIMITER = ":";
    public static Taxon ANONYMOUS = Taxon.instance("");
    /**
     * Standard key for the father of the taxon
     */
    public static final String FatherKey = "FATHER";
    /**
     * Standard key for the mother of the taxon
     */
    public static final String MotherKey = "MOTHER";
    /**
     * Standard key for the pedigree of the taxon
     */
    public static final String PedigreeKey = "PEDIGREE";
    /**
     * Standard key for the sex of the taxon
     */
    public static final String SexKey = "SEX";
    /**
     * Standard key for inbreeding coefficient of the taxon
     */
    public static final String InbreedFKey = "INBREEDF";
    /**
     * Standard key for a synonym of the taxon
     */
    public static final String SynonymKey = "SYNONYM";
    /**
     * Standard key for the latitude the taxon was sampled
     */
    public static final String LatitudeKey = "LATITUDE";
    /**
     * Standard key for the longitude the taxon was sampled
     */
    public static final String LongitudeKey = "LONGITUDE";
    /**
     * Standard key for altitude the taxon was sampled
     */
    public static final String AltitudeKey = "ALTITUDE";
    /**
     * Standard key for a genus of the taxon
     */
    public static final String GenusKey = "GENUS";
    /**
     * Standard key for a species of the taxon
     */
    public static final String SpeciesKey = "SPECIES";
    private final GeneralAnnotation myAnno;
    private final String myName;
    private final int hashCode;

    // Use instance(String name)
    private Taxon(String name) {
        this(name, GeneralAnnotationStorage.EMPTY_ANNOTATION_STORAGE);
    }

    public Taxon(String name, GeneralAnnotation anno) {
        myName = name.trim();
        hashCode = myName.hashCode();
        myAnno = anno;
    }

    public static Taxon instance(String name) {
        Taxon result = TAXON_CACHE.get(name);
        if (result == null) {
            result = new Taxon(name);
            TAXON_CACHE.put(name, result);
        }
        return result;
    }

    @Override
    public String toString() {
        return getName();
    }

    public String toStringWithVCFAnnotation() {
        StringBuilder sb = new StringBuilder("<");
        sb.append("ID=" + getName() + ",");
        for (Map.Entry<String, String> en : myAnno.getAllAnnotationEntries()) {
            sb.append(en.getKey() + "=" + en.getValue() + ",");
        }
        if (myAnno.numAnnotations() > 0) {
            sb.deleteCharAt(sb.length() - 1);
        }
        sb.append(">");
        return sb.toString();
    }

    @Override
    public int compareTo(Taxon c) {
        if (this == c) {
            return 0;
        } else {
            return myName.compareTo(c.getName());
        }
    }

    @Override
    public boolean equals(Object c) {

        if (this == c) {
            return true;
        }

        if (c instanceof Taxon) {
            return myName.equals(((Taxon) c).getName());
        } else if (c instanceof String) {
            return myName.equals((String) c);
        } else {
            return false;
        }

    }

    public String getName() {
        return myName;
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    public GeneralAnnotation getAnnotation() {
        return myAnno;
    }

    /**
     * A builder for creating immutable Taxon instances.
     * <p>
     * Example:
     * <pre>   {@code
     * Taxon cp= new Taxon.Builder("Z001E0001:Line:mays:Zea")
     *   .inbreedF(0.99)
     *   .parents("B73","B97")
     *   .pedigree("(B73xB97)S6I1")
     *   .addAnno("Group","Dent")
     *   .build();}</pre>
     * <p>
     * This would create an Taxon.
     */
    public static class Builder {

        // Required parameters
        private String myTaxonName;
        private GeneralAnnotationStorage.Builder myAnnotations = GeneralAnnotationStorage.getBuilder();

        /**
         * Constructor for Builder, requires a Taxon object
         *
         * @param aTaxon taxon object
         */
        public Builder(Taxon aTaxon) {
            myTaxonName = aTaxon.getName();
            myAnnotations.addAnnotations(aTaxon.getAnnotation());
        }

        /**
         * Constructor for Builder, requires a Taxon name
         *
         * @param aTaxonName name of the taxon
         */
        public Builder(String aTaxonName) {
            myTaxonName = aTaxonName;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnno(String key, String value) {
            myAnnotations.addAnnotation(key, value);
            return this;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnno(String key, Number value) {
            myAnnotations.addAnnotation(key, value);
            return this;
        }

        /**
         * Change the name
         */
        /**
         * Set sex: 0=both, 1=female, 2=male (default=0 Both)
         */
        public Builder name(String newName) {
            myTaxonName = newName;
            return this;
        }

        /**
         * Set sex: 0=both, 1=female, 2=male (default=0 Both)
         */
        public Builder sex(byte val) {
            return addAnno(SexKey, val);
        }

        /**
         * Set inbreeding coefficient (default=Float.NaN)
         */
        public Builder inbreedF(float val) {
            return addAnno(InbreedFKey, val);
        }

        /**
         * Set text definition of parents (default=null)
         */
        public Builder parents(String mom, String dad) {
            addAnno(MotherKey, mom);
            return addAnno(FatherKey, dad);
        }

        /**
         * Set text definition of pedigree (default=null)
         */
        public Builder pedigree(String val) {
            return addAnno(PedigreeKey, val);
        }

        /**
         * Set text definition of pedigree (default=null)
         */
        public Builder synonym(String val) {
            return addAnno(SynonymKey, val);
        }

        public Taxon build() {
            return new Taxon(myTaxonName, myAnnotations.build());
        }
    }
}
