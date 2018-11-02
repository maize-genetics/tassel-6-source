/*
 * TOPMTableReport
 */
package net.maizegenetics.dna.map;

import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.TableReport;

/**
 *
 * @author Terry Casstevens
 */
public class TOPMTableReport implements TableReport {

    private static final String[] DEFAULT_COLUMN_HEADINGS = new String[]{"Tag", "Tag Length", "Multi Maps", "Chromosome", "Strand", "Start Pos", "End Pos", "Divergence"};

    private final TOPMInterface myTOPM;
    private final String[] myColumnHeadings;

    public TOPMTableReport(TOPMInterface topm) {
        myTOPM = topm;
        int numVariants = myTOPM.getMaxNumVariants();
        int totalHeadings = DEFAULT_COLUMN_HEADINGS.length + numVariants * 2;
        myColumnHeadings = new String[totalHeadings];
        for (int i = 0; i < DEFAULT_COLUMN_HEADINGS.length; i++) {
            myColumnHeadings[i] = DEFAULT_COLUMN_HEADINGS[i];
        }
        for (int j = DEFAULT_COLUMN_HEADINGS.length; j < totalHeadings; j += 2) {
            myColumnHeadings[j] = "Offset";
            myColumnHeadings[j + 1] = "Variant";
        }
    }

    @Override
    public Object[] getTableColumnNames() {
        return myColumnHeadings;
    }

    @Override
    public String getTableTitle() {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public int getColumnCount() {
        return myColumnHeadings.length;
    }

    @Override
    public long getRowCount() {
        return myTOPM.getSize();
    }

    @Override
    public long getElementCount() {
        return getColumnCount() * getRowCount();
    }

    @Override
    public Object[] getRow(long row) {
        throw new UnsupportedOperationException("Not supported.");
    }

    @Override
    public Object getValueAt(long rowLong, int col) {
        int row = (int) rowLong;
        switch (col) {
            case 0:
                return BaseEncoder.getSequenceFromLong(myTOPM.getTag(row));
            case 1:
                return myTOPM.getTagLength(row);
            case 2:
                return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getMultiMaps(col));
            case 3:
                return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getChromosome(row));
            case 4:
                return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getStrand(row));
            case 5:
                return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getStartPosition(row));
            case 6:
                return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getEndPosition(row));
            case 7:
                return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getDivergence(row));
            default:
                int varIndex = col - 8;
                if (varIndex % 2 == 0) {
                    return AbstractTagsOnPhysicalMap.printWithMissing(myTOPM.getVariantPosOff(row, varIndex / 2));
                } else {
                    byte vd = myTOPM.getVariantDef(row, varIndex / 2);
                    if (vd == TOPMInterface.BYTE_MISSING) {
                        return AbstractTagsOnPhysicalMap.printWithMissing(vd);
                    } else {
                        return NucleotideAlignmentConstants.getHaplotypeNucleotide(vd);
                    }
                }
        }
    }

}
