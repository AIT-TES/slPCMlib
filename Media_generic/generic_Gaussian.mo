within slPCMlib.Media_generic;
package generic_Gaussian "Generic, Gaussian distribution, symmetric peak"

  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "generic_GaussianHysteresis";

    // --- parameters for phase transition functions ---
    constant Real GaussMelt_mu = 273.15 + 29
             "location parameter";
    constant Real GaussMelt_sigma = 1.0
             "standard deviation of the normal distribution";
  constant Modelica.Units.SI.Temperature rangeTmelting[2]={
      Modelica.Math.Distributions.Normal.quantile(
        mu=propData.GaussMelt_mu,
        sigma=propData.GaussMelt_sigma,
        u=0.001),Modelica.Math.Distributions.Normal.quantile(
        mu=propData.GaussMelt_mu,
        sigma=propData.GaussMelt_sigma,
        u=0.999)} "temperature range melting {startT, endT}";
    // ---
    constant Real GaussSoli_mu = 273.15 + 27
             "location parameter";
    constant Real GaussSoli_sigma = 0.5
             "standard deviation of the normal distribution";
  constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
      Modelica.Math.Distributions.Normal.quantile(
        mu=propData.GaussSoli_mu,
        sigma=propData.GaussSoli_sigma,
        u=0.001),Modelica.Math.Distributions.Normal.quantile(
        mu=propData.GaussSoli_mu,
        sigma=propData.GaussSoli_sigma,
        u=0.999)} "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={2000.0,30}
    "solid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={2000.0,20}
    "liquid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=150000.0
    "scalar phase transition enthalpy";

    // --- reference values ---
  constant Modelica.Units.SI.Temperature Tref=273.15 + 26
    "reference Temperature";
  constant Modelica.Units.SI.SpecificEnthalpy href=0.0
    "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
     "Returns liquid mass phase fraction for complete melting processes"
  algorithm
    dxi :=Modelica.Math.Distributions.Normal.density(
        mu=propData.GaussMelt_mu,
        sigma=propData.GaussMelt_sigma,
        u=T);
     xi :=Modelica.Math.Distributions.Normal.cumulative(
        mu=propData.GaussMelt_mu,
        sigma=propData.GaussMelt_sigma,
        u=T);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
     "Returns liquid mass phase fraction for complete solidification processes"
  algorithm
    dxi :=Modelica.Math.Distributions.Normal.density(
        mu=propData.GaussSoli_mu,
        sigma=propData.GaussSoli_sigma,
        u=T);
     xi :=Modelica.Math.Distributions.Normal.cumulative(
        mu=propData.GaussSoli_mu,
        sigma=propData.GaussSoli_sigma,
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
      Gaussian distribution (also referred to as Normal distribution), see
        <blockquote>          
        <a href>Modelica.Math.Distributions.Normal</a>
        </blockquote>          
      The cumulative distribution function (CDF) is used to model the 
      liquid mass phase fraction <var>xi</var> for complete phase transitions 
      (complete melting or solidification).<br> 
      The probability distribution function 
      (PDF) is used to model the derivative w.r.t. 
      temperature <var>d xi / d T</var>.
      </p>
      <p>
      The peak in the Gaussian PDF is <strong>symmetric</strong>. 
      The function has a
      location parameter <var> mu </var> and 
        a shape parameter <var> sigma>0 </var>. 
      </p>
      <p>
      The Gaussian CDF reaches 0 and 1 only asymptotically. 
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
end generic_Gaussian;
