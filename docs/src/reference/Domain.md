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
Dionysos.Domain.merge_to_hyperrectangles_pos
Dionysos.Domain.merge_to_hyperrectangles_real
Dionysos.Domain.DomainList
Dionysos.Domain.PeriodicDomainList
Dionysos.Domain.has_same_periodicity
Dionysos.Domain.wrap_pos
Dionysos.Domain.wrap_coord
Dionysos.Domain.GeneralDomainList
Dionysos.Domain.RectangularObstacles
Dionysos.Domain.CustomList
```