package net.maizegenetics.util;

/**
 * Based on response in http://stackoverflow.com/questions/2670982/using-pairs-or-2-tuples-in-java
 *
 * @author Eli Rodgers-Melnick
 */
public class Tuple<X, Y> implements Comparable<Tuple<X, Y>> {
    public final X x;
    public final Y y;

    /**
     * Instantiates a tuple object, which just holds 2 values
     *
     * @param x The first object
     * @param y The second object
     */
    public Tuple(X x, Y y) {
        this.x = x;
        this.y = y;
    }

    public X getX() {
        return x;
    }

    public Y getY() {
        return y;
    }

    @Override
    public int hashCode() {
        return (x.hashCode() ^ y.hashCode());
    }

    @Override
    public String toString() {
        return "(" + x.toString() + "," + y.toString() + ")";
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || !(obj instanceof Tuple)) {
            return false;
        }

        return compareTo((Tuple) obj) == 0;
    }

    /**
     * Using this method (i.e. sorting Tuples)
     * seems to imply both x and y are Comparable.
     * If that's not the case, results are unpredictable.
     */
    @Override
    public int compareTo(Tuple<X, Y> o) {
        if (this == o) {
            return 0;
        }
        if (x instanceof Comparable) {
            int i = ((Comparable) x).compareTo(o.x);
            if (i != 0) return i;
        }
        if (y instanceof Comparable) {
            return ((Comparable) y).compareTo(o.y);
        }
        throw new IllegalStateException("Tuple: compareTo: neither x or y is Comparable types");
    }

}
