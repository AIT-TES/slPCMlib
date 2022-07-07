within slPCMlib.Media_generic;
package generic_7thOrderSmoothStep "Generic, 7th-order smoothstep, symmetric peak"

  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

  constant String mediumName = "generic_7thOrderSmoothStepHysteresis";

    // --- parameters for phase transition functions ---
  constant Boolean modelForMelting =        true;
  constant Boolean modelForSolidification = true;
  constant Modelica.Units.SI.Temperature rangeTmelting[2]={273.15 + 26,
                                                           273.15 + 30}
                              "temperature range melting {startT, endT}";
  constant Modelica.Units.SI.Temperature rangeTsolidification[2]={273.15 + 24,
                                                                  273.15 + 29}
                              "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={2000.0,20.0}
    "solid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={3000.0,20.0}
    "liquid specific heat capacity, linear coefficients a + b*T";
  constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=150e3
    "scalar phase transition enthalpy";

    // --- reference values ---
  constant Modelica.Units.SI.Temperature Tref=273.15 + 24
    "reference temperature";
  constant Modelica.Units.SI.SpecificEnthalpy href=0.0
    "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  algorithm
    (xi, dxi) := BasicUtilities.smoothStep7thOrder(x = T,
                  startStepX = propData.rangeTmelting[1],
                  endStepX =   propData.rangeTmelting[2]);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
   "Returns liquid mass phase fraction for complete solidification processes"
  algorithm
    (xi, dxi) := BasicUtilities.smoothStep7thOrder(x = T,
                  startStepX = propData.rangeTsolidification[1],
                  endStepX =   propData.rangeTsolidification[2]);
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
      In this package, complete phase transitions (complete melting or solidification) 
      are modelled by the 7th-order smoothstep function, see
        <blockquote>          
        <a href>https://en.wikipedia.org/wiki/Smoothstep</a>. 
        </blockquote>          
      <p>
      The smoothstep function is used to model the 
      liquid mass phase fraction <var>xi</var> for complete phase transitions 
      (complete melting or solidification). <br>
      It's derivative is used to model 
      the derivative of <var>xi</var> w.r.t. temperature: <var>d xi / d T</var>. 
      </p>
      <p>
      The smoothstep function is sigmoid-like and <strong>symmetric</strong>. 
      The 7th-order smoothstep is <strong>3 times continuously differentiable</strong> (C^3 smooth). <br>
        (Note that the effective (or apparent) heat capacity model uses the first 
        derivative of the smoothstep function.) 
      </p>
      <p>
      The function has two parameters: the <strong>left edge</strong> 
      and the <strong>right edge</strong>. 
      The function returns 0 for inputs less than or equal to the left edge, 
      1 for inputs greater than or equal to the right edge. <br>
      The edges define the phase transition temperature range for melting 
      <var> rangeTmelting[2] </var> and 
      for solidification <var> rangeTsolidification[2] </var>. 
      These ranges, 
      specifically the edges, are especially important in the 
      <strong>curve track hysteresis model</strong>: 
      <a href>slPCMlib.Interfaces.phTransModCurveTrackHysteresis</a>. 
      (Note that when using phase transition models which only asymptotically 
      reach 0 and 1, the definition of phase transition temperature ranges is 
      not straightforward. Switches - by the curve track hysteresis model - 
      between the phase transition functions for 
      melting and solification would then produce jump discontinuities in the 
      phase fraction.)  
        </p></html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end generic_7thOrderSmoothStep;
