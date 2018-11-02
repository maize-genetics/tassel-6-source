/*
 *  TaxaList
 */
package net.maizegenetics.taxa;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;

/**
 *
 * @author Terry Casstevens
 */
public interface TaxaList extends List<Taxon> {

    /**
     * Returns number of taxa
     *
     * @return number of taxa
     */
    public int numberOfTaxa();

    /**
     * Return taxa name at given index.
     *
     * @param index
     *
     * @return taxa name
     */
    public String taxaName(int index);

    /**
     * Return a list of all matching taxa indices for a given name.
     *
     * @param name name
     *
     * @return Indices for matching taxa (-1 if no match).
     */
    public int indexOf(String name);

    /**
     * Return a list of all matching taxa indices for a given name.
     *
     * @param taxon taxon
     *
     * @return Indices for matching taxa (-1 if no match).
     */
    public int indexOf(Taxon taxon);

    /**
     * Returns TaxaList Collector.
     * 
     * @return collector 
     */
    public static Collector<Taxon, ?, TaxaList> collect() {
        return new TaxaListCollector();
    }

    public static class TaxaListCollector implements Collector<Taxon, TaxaListBuilder, TaxaList> {

        @Override
        public Supplier<TaxaListBuilder> supplier() {
            return TaxaListBuilder::new;
        }

        @Override
        public BiConsumer<TaxaListBuilder, Taxon> accumulator() {
            return TaxaListBuilder::add;
        }

        @Override
        public BinaryOperator<TaxaListBuilder> combiner() {
            return (left, right) -> {
                left.addAll(right);
                return left;
            };
        }

        @Override
        public Function<TaxaListBuilder, TaxaList> finisher() {
            return TaxaListBuilder::build;
        }

        @Override
        public Set<Collector.Characteristics> characteristics() {
            return Collections.EMPTY_SET;
        }

    }
}
