PLUMBER2
=================
Scripts to create the forcing and evaluation dataset for the PLUMBER2 ecosystem 
model intercomparison project.

More details on the data processing steps are available in:
Ukkola et al., A flux tower dataset tailored for land model evaluation, _Earth_ 
_System_ _Science_ _Data_, 14, 449-461, <https://doi.org/10.5194/essd-14-449-2022>

Final PLUMBER2 dataset is available at:
<https://researchdata.edu.au/plumber2-forcing-evaluation-surface-models/1656048>

Also see FluxnetLSM(<https://github.com/aukkola/FluxnetLSM>) R package used for 
pre-processing steps. 

**NB:** there is currently an error in the _Qair_qc_ processing with too many time steps flagged as "observed" at **FLUXNET2015** and **La** **Thuile** sites. As _Qair_ is calculated from _VPD_, _Tair_ and _Psurf_ at these sites, you can use the qc flags for these variables to assess _Qair_ quality 


