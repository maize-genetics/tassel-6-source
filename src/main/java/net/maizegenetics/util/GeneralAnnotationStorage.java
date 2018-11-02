/*
 *  GeneralAnnotationStorage
 * 
 *  Created on Feb. 2, 2015
 */
package net.maizegenetics.util;

import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.SetMultimap;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class GeneralAnnotationStorage implements GeneralAnnotation {

    public static final GeneralAnnotationStorage EMPTY_ANNOTATION_STORAGE = new GeneralAnnotationStorage();

    private static final double[] EMPTY_DOUBLE_ARRAY = new double[0];

    private static final int MAX_CACHE_SIZE = 1_000_000;

    private static final Map<Map.Entry<String, String>, Map.Entry<String, String>> CACHE = Collections.synchronizedMap(new LinkedHashMap<Map.Entry<String, String>, Map.Entry<String, String>>((3 * MAX_CACHE_SIZE) / 2) {

        @Override
        protected boolean removeEldestEntry(Map.Entry<Map.Entry<String, String>, Map.Entry<String, String>> eldest) {
            return size() > MAX_CACHE_SIZE;
        }

    });

    private static Map.Entry<String, String> getCanonicalAnnotation(String key, String value) {
        Map.Entry<String, String> temp = new AbstractMap.SimpleImmutableEntry<>(key, value);
        Map.Entry<String, String> entry = CACHE.putIfAbsent(temp, temp);
        return (entry == null) ? temp : entry;
    }

    private final Map.Entry<String, String>[] myAnnotations;

    private GeneralAnnotationStorage(Builder builder) {
        myAnnotations = (Map.Entry<String, String>[]) new Map.Entry<?, ?>[builder.myAnnotations.size()];
        for (int i = 0; i < builder.myAnnotations.size(); i++) {
            myAnnotations[i] = builder.myAnnotations.get(i);
        }
    }

    private GeneralAnnotationStorage() {
        myAnnotations = (Map.Entry<String, String>[]) new Map.Entry<?, ?>[0];
    }

    public static Builder getBuilder() {
        return new Builder();
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        List<String> result = new ArrayList<>(1);
        for (Map.Entry<String, String> me : myAnnotations) {
            if (me.getKey().equals(annoName)) {
                result.add(me.getValue());
            }
        }
        return result.toArray(new String[result.size()]);
    }

    @Override
    public Map<String, String> getConcatenatedTextAnnotations() {
        Map<String, String> result = new HashMap<>();
        for (Map.Entry<String, String> me : myAnnotations) {
            String value = result.get(me.getKey());
            if (value == null) {
                result.put(me.getKey(), me.getValue());
            } else {
                result.put(me.getKey(), value + "," + me.getValue());
            }
        }
        return result;
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        try {
            ArrayList<Double> result = new ArrayList<>(1);
            for (Map.Entry<String, String> me : myAnnotations) {
                if (me.getKey().equals(annoName)) {
                    result.add(Double.parseDouble(me.getValue()));
                }
            }
            if (result.isEmpty()) {
                return EMPTY_DOUBLE_ARRAY;
            }
            double[] d = new double[result.size()];
            for (int i = 0; i < result.size(); i++) {
                d[i] = result.get(i);
            }
            return d;
        } catch (Exception e) {
            return EMPTY_DOUBLE_ARRAY;
        }
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        double[] values = getQuantAnnotation(annoName);
        if (values.length == 0) {
            return Double.NaN;
        } else if (values.length == 1) {
            return values[0];
        } else {
            double result = 0.0;
            for (double current : values) {
                result += current;
            }
            return result / (double) values.length;
        }
    }

    @Override
    public boolean isAnnotatedWithValue(String annoName, String annoValue) {
        for (Map.Entry<String, String> me : myAnnotations) {
            if (me.getKey().equals(annoName) && me.getValue().equals(annoValue)) {
                return true;
            }
        }
        return false;
    }

    @Override
    public Map.Entry<String, String>[] getAllAnnotationEntries() {
        return Arrays.copyOf(myAnnotations, myAnnotations.length);
    }

    @Override
    public Set<String> getAnnotationKeys() {
        Set<String> result = new HashSet<>();
        for (Map.Entry<String, String> me : myAnnotations) {
            result.add(me.getKey());
        }
        return result;
    }

    @Override
    public SetMultimap<String, String> getAnnotationAsMap() {
        ImmutableSetMultimap.Builder<String, String> result = new ImmutableSetMultimap.Builder<String, String>()
                .orderKeysBy(Ordering.natural()).orderValuesBy(Ordering.natural());
        for (Map.Entry<String, String> en : myAnnotations) {
            result.put(en.getKey(), en.getValue());
        }
        return result.build();
    }

    @Override
    public int numAnnotations() {
        return myAnnotations.length;
    }

    public static class Builder {

        private final List<Map.Entry<String, String>> myAnnotations = new ArrayList<>(0);

        private Builder() {
        }

        public Builder addAnnotation(String key, String value) {
            myAnnotations.add(getCanonicalAnnotation(key, value));
            return this;
        }

        public Builder addAnnotation(String key, Number value) {
            myAnnotations.add(getCanonicalAnnotation(key, value.toString()));
            return this;
        }

        public Builder addAnnotations(GeneralAnnotation existing) {
            if (existing == null) {
                return this;
            }
            myAnnotations.addAll(Arrays.asList(existing.getAllAnnotationEntries()));
            return this;
        }

        public GeneralAnnotationStorage build() {
            if (myAnnotations.isEmpty()) {
                return EMPTY_ANNOTATION_STORAGE;
            }
            Collections.sort(myAnnotations, (Map.Entry<String, String> s1, Map.Entry<String, String> s2) -> {
                int keyComp = s1.getKey().compareTo(s2.getKey());
                if (keyComp != 0) {
                    return keyComp;
                }
                return s1.getValue().compareTo(s2.getValue());
            });
            return new GeneralAnnotationStorage(this);
        }
    }
}
