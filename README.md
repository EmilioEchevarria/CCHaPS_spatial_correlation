# CCHaPS_spatial_correlation
Spatial correlation analysis of wave variables using CCHaPS output

Written in Fortran to make the analysis tractable: `spatial_correlation.f90` for non-directional variables, and `spatial_correlation_dir.f90` for wave direction. The compilation instructions are written for petrichor. 

The Fortran codes need an input file, `spatial_corr.inp`:

```
output_var=hs
start_year=1981
end_year=2020

the_lon=147.72
the_lat=-43.4

max_dist=1000.0 ! It will select the nodes that are within this distance for calculating correlations
```

It returns a file called `schism_corrs_node_{sch_node}.txt` with the following data on each line:

```
lon  lat  schism_idx  correlation
```