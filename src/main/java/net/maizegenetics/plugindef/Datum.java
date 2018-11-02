/*
 * Datum.java
 *
 */
package net.maizegenetics.plugindef;

/**
 * This wraps data elements used as input or output to Tassel modules.
 *
 * @author Terry Casstevens
 */
public class Datum {

    private String myName;
    private final Object myData;
    private String myComment;

    /**
     * Creates a new instance of Datum
     *
     * @param name name of this data instance.
     * @param data data element
     * @param comment optional comment
     */
    public Datum(String name, Object data, String comment) {

        if ((name == null) || (name.length() == 0)) {
            throw new IllegalArgumentException("Datum: init: name must be assigned.");
        }

        if (data == null) {
            throw new IllegalArgumentException("Datum: init: data can not be null.");
        }

        myName = name;
        myData = data;
        myComment = comment;

    }

    public String getName() {
        return myName;
    }

    public void setName(String s) {
        myName = s;
    }

    public Object getData() {
        return myData;
    }

    public String getComment() {
        return myComment;
    }

    public Class getDataType() {
        return myData.getClass();
    }

    @Override
    public String toString() {
        return myName;
    }

    @Override
    public boolean equals(Object obj) {

        if (obj == this) {
            return true;
        }

        Datum instance;
        if (obj instanceof Datum) {
            instance = (Datum) obj;
        } else {
            return false;
        }

        return myName.equals(instance.getName()) && (myData == instance.getData()) && (myComment.equals(instance.getComment()));
    }
}
