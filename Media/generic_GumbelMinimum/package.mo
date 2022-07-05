within slPCMlib.Media;
package generic_GumbelMinimum "Generic, Gumbel Minimum distribution, asymmetric peak"

  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "generic_GumbelMinimumHysteresis";

    // --- parameters for phase transition function ---
    constant Real GumbelMinSoli_mu =   273.15 + 28
             "Location parameter";
    constant Real GumbelMinSoli_beta = 0.2
             "Shape parameter, smaller is sharper, beta>0";
  constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
      slPCMlib.BasicUtilities.GumbelMinimum.quantile(
        mu=propData.GumbelMinSoli_mu,
        beta=propData.GumbelMinSoli_beta,
        u=0.001),slPCMlib.BasicUtilities.GumbelMinimum.quantile(
        mu=propData.GumbelMinSoli_mu,
        beta=propData.GumbelMinSoli_beta,
        u=0.999)} "temperature range solidification {startT, endT}";
    // ---
    constant Real GumbelMinMelt_mu =   273.15 + 30
             "Location parameter";
    constant Real GumbelMinMelt_beta = 0.4
             "Shape parameter, smaller is sharper, beta>0";
  constant Modelica.Units.SI.Temperature rangeTmelting[2]={
      slPCMlib.BasicUtilities.GumbelMinimum.quantile(
        mu=propData.GumbelMinMelt_mu,
        beta=propData.GumbelMinMelt_beta,
        u=0.001),slPCMlib.BasicUtilities.GumbelMinimum.quantile(
        mu=propData.GumbelMinMelt_mu,
        beta=propData.GumbelMinMelt_beta,
        u=0.999)} "temperature range melting {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={2000.0,20}
    "solid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={2000.0,10}
    "liquid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=150000.0
    "scalar phase transition enthalpy";

    // --- reference values ---
  constant Modelica.Units.SI.Temperature Tref=273.15 + 26
    "reference Temperature";
  constant Modelica.Units.SI.SpecificEnthalpy href=0.0
    "reference enthalpy at Tref";

   // constant String transFctForHeating="7thSmoothStep";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
     "Returns liquid mass phase fraction for complete melting processes"
  algorithm
    xi  :=slPCMlib.BasicUtilities.GumbelMinimum.cumulative(
        mu=propData.GumbelMinMelt_mu,
        beta=propData.GumbelMinMelt_beta,
        u=T);
    dxi :=slPCMlib.BasicUtilities.GumbelMinimum.density(
        mu=propData.GumbelMinMelt_mu,
        beta=propData.GumbelMinMelt_beta,
        u=T);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
     "Returns liquid mass phase fraction for complete solidification processes"
  algorithm
    xi  :=slPCMlib.BasicUtilities.GumbelMinimum.cumulative(
        mu=propData.GumbelMinSoli_mu,
        beta=propData.GumbelMinSoli_beta,
        u=T);
    dxi :=slPCMlib.BasicUtilities.GumbelMinimum.density(
        mu=propData.GumbelMinSoli_mu,
        beta=propData.GumbelMinSoli_beta,
        u=T);
  end phaseFrac_complSolidification;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.600000000000000e+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
   lambda := 2.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
   lambda := 3.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation (Documentation(info="<html>
      <p>
      In this package, complete phase transitions 
      (complete melting or solidification) 
      are modelled by the 
      Gumbel Minimum distribution 
      (also referred to as Extreme value type I distribution), see
        <blockquote>          
        NIST/SEMATECH (01/2015). e-Handbook of Statistical Methods. 
        <a href>http://www.itl.nist.gov/div898/handbook/</a>
        </blockquote>          
      <p>
      The cumulative distribution function (CDF) is used to model the 
      liquid mass phase fraction <var>xi</var> for complete phase transitions 
      (complete melting or solidification). <br>
      The probability distribution function 
      (PDF) is used to model the derivative w.r.t. temperature <var>d xi / d T</var>. 
      </p>
      <p>
      The Gumbel Minimum distribution 
      is based on the smallest extreme of a distribution. 
      The peak in the PDF is <strong>left-skewed</strong>. 
      The function can be convenient to represent <strong>asymmetric</strong> 
      peaks, see e.g. 
        <blockquote>          
        Barz, T., Emhofer, J., Marx, K., Zsembinszki, G., & Cabeza, L. F. (2019). 
        Phenomenological modelling of phase transitions with hysteresis in 
        solid/liquid PCM. 
        Journal of Building Performance Simulation, 12(6), 770-788.
        <a href>https://doi.org/10.1080/19401493.2019.1657953</a>
        </blockquote> 
      The function has a
      location parameter <var> mu </var> and 
      a shape parameter <var> beta>0 </var>. 
      </p>
      <p>
      The Gumbel CDF reaches 0 and 1 only asymptotically. 
      However, certain phase transition models use the limits of the 
      phase transition temperature range. E.g. the 
      <a href>slPCMlib.Interfaces.phTransModCurveTrackHysteresis</a> 
      switches when the temperature leaves the 
      phase transition temperature range. 
      Therefore, the phase transition temperature range is 
      <strong>approximately</strong> defined by 
      <var>T = invCDF(P=0.001)</var> and <var>T = invCDF(P=0.999)</var>, 
      where <var>invCDF</var> is the the inverse of the CDF and <var>P</var> is the probability. 
      </p></html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end generic_GumbelMinimum;
