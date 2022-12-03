# _slPCMlib_ - solid/liquid PCM Modelica library 

The library _slPCMlib_ contains property models for solid/liqid phase change materials (PCM)
showing a non-isothermal phase transition behavior. 

The library contains generic PCM and specific commercial PCM (media).  
Different phenomenological phase transition models are implemented 
to account for temperature shifts in latent transition changes, 
e.g. due to multi-step transitions and thermal hysteresis. 
The models predict effective properties which are valid over the
PCM functional temperature range where latent heat is absorbed and released. 
Based on the properties and adopting the apparent heat capacity method, heat transfer problems can be solved 
numerically.  

## Description

![Alt text](./Resources/Images/slPCMlib.png?raw=true "Title")

Assumptions for modeling effective PCM properties:
* There are only two phases (two-phase model): a solid and a liquid phase.
* Phase transitions are induced by temperature and are independent of pressure.
* Phase transitions extend over a temperature range (non-isothermal phase transitions) and are continuous.
* Within the phase transition temperature range the solid and liquid phases coexist as a homogenous mixture (macroscopic view). 
  The material is then in a semi-solid or semi-liquid state which produces a mushy zone in the PCM domain.
* Properties of the mushy state are local effective (also apparent) mixture properties, which are defined by a weighting of contributions from
solid and liquid phases. The weighting is based on the phase change progress, i.e. the mass (or volume) phase fraction.  

Temperature is the input to the model using the < inductionAtNode > connector. For a given temperature input *T* the liquid mass phase fraction *xi* is computed, and the following variables are derived: 

* effective density *rho* 
* liquid volume phase fraction *phi*
* effective thermal conductivity *lambda* 
* effective specific heat capacity *cp* and specific enthalpy *h* 
* baseline heat capacity *c_BL*, which describes the mixture heat capacity (without the effect of the phase transition enthalpy) 
* solid *h_S* and liquid *h_L* enthalpies, where the difference *h_L - h_S* is the temperature-dependent phase transition enthalpy

_slPCMlib_ contains <ins>material data</ins> for: 
* generic PCM, for which the user can make copies and adapt the parameters to match certain behavior, 
* specific PCM, which are based on experimental data, e.g. data taken from manufacturer's data sheets. 

Material data means solid and liquid properties mentioned above (which can be functions of temperature), 
and liquid mass phase transition curves for complete melting and solidification. 

_slPCMlib_ also contains the following <ins>phase transition models</ins> which compute the liquid mass phase fraction *xi* : 
* melting curve model - rate-independent, no hysteresis 
* curve track hysteresis model - rate-independent, hysteresis 
* curve switch hysteresis model - rate-independent, hysteresis 
* curve scale hysteresis model - rate-independent, hysteresis 

Phase transition models provide different behavior for temperature cycles with complete and incomplete phase transitions, relevant e.g. for studying cyclic charging/ discharging of thermal storages using PCM. 

_slPCMlib_ also contains exemplarily <ins>heat transfer models</ins> to provide small demonstration projects.  

For more detailed information see the documentation layers of the Modelica packages. 

## References

Barz, T., Bres, A., & Emhofer, J. (2022). slPCMlib: A Modelica Library for the Prediction of Effective Thermal Material Properties of Solid/Liquid Phase Change Materials (PCM). In Proceedings of Asian Modelica Conference 2022 (pp. 63-74). Linkoping University Electronic Press. https://doi.org/10.3384/ecp19363

Barz, T. (2021). Paraffins as phase change material in a compact plate-fin heat exchanger-Part II: Validation of the “curve scale” hysteresis model for incomplete phase transitions. Journal of Energy Storage, 34, 102164. https://doi.org/10.1016/j.est.2020.102164

Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase Fraction–Temperature Curves from Heat Capacity Data for Numerical Modeling of Heat Transfer in Commercial Paraffin Waxes. Energies, 13(19), 5149. https://doi.org/10.3390/en13195149

Barz, T., Emhofer, J., Marx, K., Zsembinszki, G., & Cabeza, L. F. (2019). Phenomenological modelling of phase transitions with hysteresis in solid/liquid PCM. Journal of Building Performance Simulation, 12(6), 770-788. https://doi.org/10.1080/19401493.2019.1657953

Barz, T., & Sommer, A. (2018). Modeling hysteresis in the phase transition of industrial-grade solid/liquid PCM for thermal energy storages. International Journal of Heat and Mass Transfer, 127, 701-713. https://doi.org/10.1016/j.ijheatmasstransfer.2018.08.032

## License

The Modelica _slPCMlib_ library is available under a 3-clause BSD-license, see
[slPCMlib license](https://gitlab-intern.ait.ac.at/tes/pcm/slpcmlib/-/blob/master/LICENSE).
