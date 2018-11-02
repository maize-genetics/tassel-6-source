// LableMapping.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.taxa.tree;

/**
 * Title:        LabelMapping
 * Description:  Allows for the substitution of one label for another
 * @author			 Matthew Goode
 * @version 1.0
 */

import net.maizegenetics.taxa.Taxon;

import java.util.Hashtable;
public class LabelMapping implements java.io.Serializable {
	Hashtable<String,String> mappings_ = new Hashtable<>();

	//
	// Serialization code
	//
	private static final long serialVersionUID=-9217142171228146380L;

	//serialver -classpath ./classes net.maizegenetics.pal.tree.LabelMapping
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeObject(mappings_);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				mappings_ = (Hashtable<String,String>)in.readObject();
				break;
			}
		}
	}

	public String getLabel(String id, String defaultLabel) {
		if(id==null||!mappings_.containsKey(id)) {
			return defaultLabel;
		}
		return mappings_.get(id);
	}

	public Taxon getLabelIdentifier(Taxon id) {
		if(id==null) {
			return null;
		}
		return new Taxon(getLabel(id.getName(),id.getName()));
	}

}