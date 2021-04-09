package net.maizegenetics.plugindef

import kotlin.reflect.KProperty

/**
 * @author Terry Casstevens
 * Created October 20, 2018
 */

class Parameter<T>() {

    operator fun getValue(thisRef: Plugin, property: KProperty<*>): T {
        return thisRef.getParameter(property.name) as T
    }

    operator fun setValue(thisRef: Plugin, property: KProperty<*>, value: T): Plugin {
        thisRef.setParameter(property.name, value)
        return thisRef
    }

}
