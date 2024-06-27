# SHIP_STATS_v2

Affected files from UWSAM which calculates the environmental accumulation mode aerosol (NAc) on the edges of the domain (controlled by edge_frac) while the LES is running. This differs from the old version which used the background aerosol before the injection began to determine the presence of a ship track. The environmental NAc is updated at each statistics call and used to determine if a column is in the ship track or not and provides more flexibility for environments where the background aerosol may be changing considerably.  

To run conditional ship track statistics add doShipTrackConditionals = .true. to the prm file. This specific version is built for compatibility with specific radiation/microphysics/turbulence routines but could easily be adapted for other packages in the future: MICRO_M2005_PA --> microphysics.f90, SGS_TKE --> (tke_full.f90, diffuse_mom.f90, diffuse_mom3D.f90), RAD_RRTM4PBL --> rad_full.f90. 

This code also includes hyperviscosity in the momentum equations for larger grid spacing runs to potentially limit excessive entrainment. The hyperviscosity functionality can be turned on in params.f90 with doMomentumHyperviscosity = .true.. The timescale over which the hyperviscosity acts is also in params.f90 (tau_MomentumHyperviscosity = 1200.).


