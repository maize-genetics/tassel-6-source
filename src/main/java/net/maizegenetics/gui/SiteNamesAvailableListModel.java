/*
 *  SiteNamesAvailableListModel
 * 
 *  Created on Jul 3, 2014
 */
package net.maizegenetics.gui;

import net.maizegenetics.dna.map.PositionList;

/**
 *
 * @author Terry Casstevens
 */
public class SiteNamesAvailableListModel extends AbstractAvailableListModel {
    
    private final PositionList myPositions;
    
    public SiteNamesAvailableListModel(PositionList positions) {
        super(positions.numberOfSites());
        myPositions = positions;
    }
    
    @Override
    public String getRealElementAt(int index) {
        return myPositions.siteName(index);
    }
    
}
