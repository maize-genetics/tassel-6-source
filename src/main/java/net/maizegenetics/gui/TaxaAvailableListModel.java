/*
 *  TaxaAvailableListModel
 * 
 *  Created on Jul 2, 2014
 */
package net.maizegenetics.gui;

import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class TaxaAvailableListModel extends AbstractAvailableListModel {

    private final TaxaList myTaxa;

    public TaxaAvailableListModel(TaxaList taxa) {
        super(taxa.numberOfTaxa());
        myTaxa = taxa;
    }

    @Override
    public String getRealElementAt(int index) {
        return myTaxa.taxaName(index);
    }

}
