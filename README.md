# Installation

`m4 libnetcdf-dev libnetcdff-dev libpnetcdf-dev`

# Benchmark

On Surface Pro 6, 8 core, it takes 1 hour to simulate 3 hour.

# List of files needed for a run:

Wildcard: `*.TBL RRTM* namelist.input wrfbdy_* wrfinput_d*`

```
GENPARM.TBL
LANDUSE.TBL
namelist.input
namelist.output
RRTM_DATA
RRTM_DATA_DBL
RRTMG_LW_DATA
RRTMG_LW_DATA_DBL
RRTMG_SW_DATA
RRTMG_SW_DATA_DBL
SOILPARM.TBL
URBPARM.TBL
VEGPARM.TBL
wrfbdy_d01
wrf.exe
wrfinput_d01
wrfinput_d02
```

# Variable naming convention:

- `i`, `j`: horizontal dimension, `k`: vertical dimension,
- `t`: tile (1 tile per processor),
- `m`: memory,
- `d`: domain.

Example:
- `ids`: Start index of `i`-dimension of the whole domain,
- `ite`: End index of `i`-dimension of the current tile.
