# _slPCMlib_ - solid/liquid PCM Modelica library 

The library _slPCMlib_ contains <ins>material data</ins>, <ins>phase transition models</ins> 
and exemplarily <ins>heat transfer models</ins> 
for solid/liqid phase change materials (PCM). 

![Alt text](./Resources/Images/slPCMlib.png?raw=true "Title")


## Description

Model assumptions: 
* Phase transitions are induced by temperature. 
* Phase transitions are pressure independent. 
* Only two phases exist (two-phase model), a solid and a liquid phase. 
* Solid and liquid phases co-exist as homogenous mixture (macroscopic view) over an extended phase transition temperature range (non-isothermal phase transitions). 
* Effective (also apparent) mixture properties are used within the phase transition temperature range. They result from a linear weighting of contributions of solid and liquid phases. The weighting factor is the mass (or volume) liquid phase fraction. 

For a given temperature input signal, liquid phase fractions are computed. The following variables are then determined, which are functions of temperature *T* and liquid mass phase fraction *xi_m*: 

* effective density *rho* 
* liquid volume phase fraction *xi_v*
* effective thermal conductivity *lambda* 
* effective specific heat capacity *cp* and specific enthalpy *h* 
* baseline heat capacity *c_BL*, which describes the mixture heat capacity (without the effect of the phase transition enthalpy) 
* solid *h_S* and liquid *h_L* enthalpies, where the difference *h_L - h_S* is the temperature-dependent phase transition enthalpy

_slPCMlib_ contains <ins>material data</ins> for: 
* generic PCM, for which the user can make copies and adapt the parameters to match certain behavior, 
* specific PCM, which are based on experimental data, e.g. data taken from manufacturer's data sheets. 

Material data means solid and liquid properties mentioned above (which can be functions of temperature), 
and liquid mass phase transition curves for complete melting and solidification. 

_slPCMlib_ also contains the following <ins>phase transition models</ins>: 
* melting curve model - rate-independent, no hysteresis 
* curve track hysteresis model - rate-independent, hysteresis 
* curve switch hysteresis model - rate-independent, hysteresis 
* curve scale hysteresis model - rate-independent, hysteresis 

Phase transition models provide different behavior for temperature cycles with complete and incomplete phase transitions, relevant e.g. for studying cyclic charging/ discharging of thermal storages using PCM. 

_slPCMlib_ also contains exemplarily <ins>heat transfer models</ins> to provide small demonstration projects.  

For more detailed information see the documentation layers of the Modelica packages. 

## References

Barz, T. (2021). Paraffins as phase change material in a compact plate-fin heat exchanger-Part II: Validation of the “curve scale” hysteresis model for incomplete phase transitions. Journal of Energy Storage, 34, 102164. 
![DOI](https://doi.org/10.1016/j.est.2020.102164)

Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase Fraction–Temperature Curves from Heat Capacity Data for Numerical Modeling of Heat Transfer in Commercial Paraffin Waxes. Energies, 13(19), 5149.
![DOI](https://doi.org/10.3390/en13195149)

Barz, T., Emhofer, J., Marx, K., Zsembinszki, G., & Cabeza, L. F. (2019). Phenomenological modelling of phase transitions with hysteresis in solid/liquid PCM. Journal of Building Performance Simulation, 12(6), 770-788.
![DOI](https://doi.org/10.1080/19401493.2019.1657953)

Barz, T., & Sommer, A. (2018). Modeling hysteresis in the phase transition of industrial-grade solid/liquid PCM for thermal energy storages. International Journal of Heat and Mass Transfer, 127, 701-713.
![DOI](https://doi.org/10.1016/j.ijheatmasstransfer.2018.08.032)

## License

The Modelica _slPCMlib_ library is available under a 3-clause BSD-license, see
[slPCMlib license](https://gitlab-intern.ait.ac.at/tes/pcm/slpcmlib/-/blob/master/LICENSE).