within slPCMlib.Media_generic;
package generic_uniform "Generic, uniform distribution, linear"

  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

  constant String mediumName = "generic_GaussianHysteresis";

    // --- parameters for phase transition functions ---
  constant Boolean modelForMelting =        true;
  constant Boolean modelForSolidification = true;
  constant Real UniformMelt_min = 273.15 + 26
           "lower limit";                     //62
  constant Real UniformMelt_max = 273.15 + 30
           "upper limit";                     //65
  constant Modelica.Units.SI.Temperature rangeTmelting[2]={UniformMelt_min,
      UniformMelt_max} "temperature range melting {startT, endT}";
  // ---
  constant Real UniformSoli_min = 273.15 + 24
           "lower limit";                     //60.5
  constant Real UniformSoli_max = 273.15 + 29
           "upper limit";                     //62
  constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
      UniformSoli_min,UniformSoli_max}
    "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={2000.0,20}
    "solid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={2000.0,30}
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
    dxi :=Modelica.Math.Distributions.Uniform.density(
        u_min = propData.UniformMelt_min,
        u_max = propData.UniformMelt_max,
        u=T);
     xi :=Modelica.Math.Distributions.Uniform.cumulative(
        u_min = propData.UniformMelt_min,
        u_max = propData.UniformMelt_max,
        u=T);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
     "Returns liquid mass phase fraction for complete solidification processes"
  algorithm
    dxi :=Modelica.Math.Distributions.Uniform.density(
        u_min = propData.UniformSoli_min,
        u_max = propData.UniformSoli_max,
        u=T);
     xi :=Modelica.Math.Distributions.Uniform.cumulative(
        u_min = propData.UniformSoli_min,
        u_max = propData.UniformSoli_max,
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
      Uniform distribution, see
        <blockquote>          
        <a href>Modelica.Math.Distributions.Uniform</a>
        </blockquote>          
      The cumulative distribution function (CDF) is used to model the 
      liquid mass phase fraction <var>xi</var> for complete phase transitions 
      (complete melting or solidification).<br> 
      The probability distribution function 
      (PDF) is used to model the derivative w.r.t. 
      temperature <var>d xi / d T</var>.
      </p>
      <p>
      The function has two parameters: the lower and upper limit. The 
      function returns 0 for inputs less than or equal to the left limit, 
      1 for inputs greater than or equal to the right limit. 
      The limits define the phase transition temperature range 
      for melting rangeTmelting[2] and for solidification rangeTsolidification[2]. 
      </p>
      <p>
      <strong>The uniformn distribution is piecewise-linear. It's first 
      derivative is discontinuous at the lower and upper limits 
      of the transition ranges.</strong>  
        </p></html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end generic_uniform;
