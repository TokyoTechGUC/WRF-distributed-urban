# WRF with distributed urban parameters

This is a modified version of the WRF v3.3.1 used in Tokyo Tech GUC lab.
In this version, urban land use is modeled using 2D maps of urban morphological
parameters.
Hourly varying anthropogenic heat map can also be inputted to the model.

## Compilation

Refer to the documentation of WRF v3.3.1 on how to compile the model.

## How to run a simulation

### WPS

* Download the geogrid files for urban parameters from <https://doi.org/10.6084/m9.figshare.17108981>
* (Optional) Download the anthropogenic heat dataset (AH4GUC) geogrid from <https://urbanclimate.tse.ens.titech.ac.jp/database/AHE/AH4GUC/>.
* Compile the `modify_geo_em.f90` file. NetCDF Fortran library is required.
This script is used to modify `geo_em` files so that most dominant land category and second dominant land category can be simulated.
```bash
gfortran modify_geo_em.f90 -o modify_geo_em.exe $(nf-config --fflags) $(nf-config --flibs)
```
* Edit the GEOGRID.TBL of WPS as follow. Change the `rel_path` according to your directory structure.
* Note that the `rel_path` in `AHE24H` needs to be changed depending on the month of the year.
```
===============================
name = DISP # Displacement height (m)
        priority = 1
        dest_type = continuous
        fill_missing = 0
        masked = water
        interp_option =     30s:average_gcell(4.0)+four_pt+average_4pt
        interp_option = default:average_gcell(4.0)+four_pt+average_4pt
        rel_path = default:urban_params/2010/d
===============================
name = Z0_GRD # Roughness length for momentum (m)
        priority = 1
        dest_type = continuous
        fill_missing = 0
        masked = water
        interp_option =     30s:average_gcell(4.0)+four_pt+average_4pt
        interp_option = default:average_gcell(4.0)+four_pt+average_4pt
        rel_path = default:urban_params/2010/z_0
===============================
name = AVE_BH_GRD # Average building height (m)
        priority = 1
        dest_type = continuous
        fill_missing = 0
        masked = water
        interp_option =     30s:average_gcell(4.0)+four_pt+average_4pt
        interp_option = default:average_gcell(4.0)+four_pt+average_4pt
        rel_path = default:urban_params/2010/H_avg
===============================
name = FAI_GRD # Frontal area index (dimensionless)
        priority = 1
        dest_type = continuous
        fill_missing = -1
        masked = water
        interp_option =     30s:average_gcell(4.0)+four_pt+average_4pt
        interp_option = default:average_gcell(4.0)+four_pt+average_4pt
        rel_path = default:urban_params/2010/lambda_f
===============================
name = PAI_GRD # Plan area index (dimensionless)
        priority = 1
        dest_type = continuous
        masked=water
        fill_missing = -1
        interp_option =     30s:average_gcell(4.0)+four_pt+average_4pt
        interp_option = default:average_gcell(4.0)+four_pt+average_4pt
        rel_path = default:urban_params/2010/lambda_p
===============================
name = AHE24H # Anthropogenic heat (W/m^2)
        priority = 1
        dest_type = continuous
        fill_missing = 0
        interp_option =     30s:average_gcell(4.0)+four_pt+average_4pt
        interp_option = default:average_gcell(4.0)+four_pt+average_4pt
        rel_path = default:AHE/01_24 # 01 means January, 24 means 24 hours
        z_dim_name = utc_hour
==============================
```
* After running `geogrid.exe`, in the directory with `geo_em.d??.nc` files, execute
```
./modify_geo_em.exe
```

### WRF

* In `namelist.input`, modify the `domains` section as follow.
```
&domains
# Your options (e.g., time_step, e_we, e_sn)
offset_zd = 1 # Add zero-plane displacement height to topography.
/
```
* Run WRF as usual.

## References

* [WRF Model]
D.N. Khanh, A.C.G. Varquez and M. Kanda, Impact of urbanization on exposure to
extreme warming in megacities, *Heliyon*, 9, e15511, doi:
<https://doi.org/10.1016/j.heliyon.2023.e15511>
* [AH4GUC anthropogenic heat dataset] Varquez, A.C.G., Kiyomoto, S., Khanh, D.N. et al. Global 1-km present and future hourly anthropogenic heat flux. *Sci Data* 8, 64 (2021). https://doi.org/10.1038/s41597-021-00850-w
