# Mapping

This folder contains mappings to make the link between the concrete system and the abstract system, and vice-versa.

## Grids 
```@docs
Dionysos.Mapping.Grid
Dionysos.Mapping.get_pos_by_coord
Dionysos.Mapping.GridFree
Dionysos.Mapping.GridEllipsoidalRectangular
```

## Abstract Mappings 
```@docs
Dionysos.Mapping.AbstractMapping
Dionysos.Mapping.GridMapping
```

## Concrete Mappings 

### Simple Mappings
```@docs
Dionysos.Mapping.ListMapping
```

### Grid Mappings
```@docs
Dionysos.Mapping.ExplicitGridMapping
Dionysos.Mapping.ImplicitGridMapping
Dionysos.Mapping.PeriodicGridMapping
```

### Multi-level Mappings
```@docs
Dionysos.Mapping.AbstractMultiLevelMapping
Dionysos.Mapping.HierarchicalGridMapping
```

## State sets

`MappedStateSet` is the public surface (a state set bundled with its mapping);
the other set types are the building blocks it wraps.

```@docs
Dionysos.Mapping.MappedStateSet
Dionysos.Mapping.AbstractStateSet
Dionysos.Mapping.FullStateSet
Dionysos.Mapping.ExplicitIdSet
Dionysos.Mapping.ImplicitStateSet
```


