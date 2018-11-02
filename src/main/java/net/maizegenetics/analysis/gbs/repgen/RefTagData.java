/**
 * 
 */
package net.maizegenetics.analysis.gbs.repgen;

import net.maizegenetics.dna.tag.Tag;

/**
 *  Class needed for storing reference tags into RepGen SQLite
 *  tables
 * @author lcj34
 *
 */
public class RefTagData {
    private final Tag myTag;
    private final String chromosome;
    private final int position;
    private final String refGenome;

    public RefTagData(Tag myTag, String chromosome,int position, String refGenome) {
        this.myTag = myTag;
        this.chromosome = chromosome;
        this.position = position;
        this.refGenome = refGenome;
    }

    public Tag tag() {
        return myTag;
    }

    public String chromosome() {
        return chromosome;
    }
    public int position() {
        return position;
    }
    public String refGenome() {
        return refGenome;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        RefTagData that = (RefTagData) obj;

        if (!(tag().equals(that.tag()))) return false;
        if (!(chromosome().equals(that.chromosome()))) return false;
        if (position() != that.position()) return false;
        if (!(refGenome().equals(that.refGenome()))) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 37 * hash + this.myTag.hashCode();
        hash = 37 * hash + this.chromosome().hashCode();
        hash = 37 * hash + this.position;
        hash = 37 * hash + this.refGenome.hashCode();
        return hash;
    }
}
