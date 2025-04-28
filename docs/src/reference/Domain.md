# Domain

This folder contains structures that are used to encode different kinds of concrete and abstract domains.

## Concrete continuous domains 
```@docs
Dionysos.Domain.ContinuousUnboundedDomain
Dionysos.Domain.ContinuousBoundedDomain
Dionysos.Domain.ContinuousBoundedEllipsoidDomain
```

## Grids 
```@docs
Dionysos.Domain.Grid
Dionysos.Domain.get_pos_by_coord
Dionysos.Domain.GridFree
Dionysos.Domain.GridEllipsoidalRectangular
Dionysos.Domain.DeformedGrid
Dionysos.Domain.get_volume
```

## Abstract domains 
```@docs
Dionysos.Domain.DomainType
Dionysos.Domain.GridDomainType
```

```@docs
Dionysos.Domain.merge_to_hyperrectangles_pos
Dionysos.Domain.merge_to_hyperrectangles_real
```

## Concrete domains 

### Classical domains
```@docs
Dionysos.Domain.DomainList
Dionysos.Domain.CustomList
```

### Periodic domain 
```@docs
Dionysos.Domain.PeriodicDomainList
Dionysos.Domain.is_periodic
Dionysos.Domain.has_same_periodicity
Dionysos.Domain.wrap_pos
Dionysos.Domain.wrap_coord
Dionysos.Domain._make_periodic_index_map
```

### Nested domain 
```@docs
Dionysos.Domain.NestedDomain
```