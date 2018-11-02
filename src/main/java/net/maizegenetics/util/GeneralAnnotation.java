package net.maizegenetics.util;

import com.google.common.collect.SetMultimap;

import java.util.Map;
import java.util.Set;

/**
 * Provide generalized annotations (descriptors).
 *
 */
public interface GeneralAnnotation {

    /**
     * Returns all annotation value for a given annotation key
     *
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public String[] getTextAnnotation(String annoName);

    public Map<String, String> getConcatenatedTextAnnotations();

    /**
     * Returns all annotation value for a given annotation key
     *
     * @param annoName annotation key
     * @return array of annotation values (if not present new double[0])
     */
    public double[] getQuantAnnotation(String annoName);

    /**
     * Returns average annotation for a given annotation key
     *
     * @param annoName annotation key
     * @return average value (if not present - return Double.NaN)
     */
    public double getAverageAnnotation(String annoName);

    /**
     * Returns all keys
     *
     * @return
     */
    public Set<String> getAnnotationKeys();

    public SetMultimap<String, String> getAnnotationAsMap();

    /**
     * Returns whether the entity contains the annotation with the specified
     * value. If either the annotation or the value is missing false is return
     *
     * @param annoName annotation key
     * @param annoValue annotation value;
     * @return
     */
    public boolean isAnnotatedWithValue(String annoName, String annoValue);

    public Map.Entry<String, String>[] getAllAnnotationEntries();

    /**
     * Returns number of annotations including when counts of multiple values
     * for same key.
     *
     * @return number of annotations
     */
    public int numAnnotations();
}
