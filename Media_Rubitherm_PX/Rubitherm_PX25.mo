
// within slPCMlib.Rubitherm_PX;
package Rubitherm_PX25 "Rubitherm GmbH, PX25; data taken from: Rubitherm datasheet; last access: 2020-07-03."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "PX25";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 2.9914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8814999999999998E+02, 2.9914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 6.6145736861479716E+04
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature  Tref = rangeTmelting[1]
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";

  end propData;


  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer len_x =    data_H.len_x;
    constant Real data_x[:] =   data_H.data_x;
    constant Real data_y[:] =   data_H.data_y;
    constant Real m_k[:] =      data_H.m_k;
    constant Real iy_start[:] = data_H.iy_start;
    constant Real iy_scaler =   data_H.iy_scaler;
  algorithm 
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15,
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complMelting;

  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected 
    constant Integer len_x =    data_C.len_x;
    constant Real data_x[:] =   data_C.data_x;
    constant Real data_y[:] =   data_C.data_y;
    constant Real m_k[:] =      data_C.m_k;
    constant Real iy_start[:] = data_C.iy_start;
    constant Real iy_scaler =   data_C.iy_scaler;
  algorithm 
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15,
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {1.5000000000000000E+01, 1.7625000000000000E+01, 1.9875000000000000E+01, 2.2375000000000000E+01, 2.3625000000000000E+01, 2.4625000000000000E+01, 2.6000000000000000E+01};
    constant Real[7] data_y =   {0.0000000000000000E+00, 3.0444594493000000E-02, 7.5811883366999996E-02, 1.4526117675899999E-01, 2.3599448658200001E-01, 2.0505425310100001E-01, 0.0000000000000000E+00};
    constant Real[7] m_k =      {0.0000000000000000E+00, -4.5556537229999999E-03, 3.1018270882000001E-02, 6.5533692652000000E-02, 7.5130670481000000E-02, -1.2638267801799999E-01, 0.0000000000000000E+00};
    constant Real[7] iy_start = {0.0000000000000000E+00, 4.2536823569000001E-02, 1.4697517542600000E-01, 4.0511124734499998E-01, 6.4193682087500004E-01, 8.7904411084599998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9911570270917904E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {1.5000000000000000E+01, 1.7625000000000000E+01, 2.1625000000000000E+01, 2.3625000000000000E+01, 2.4625000000000000E+01, 2.5625000000000000E+01, 2.6000000000000000E+01};
    constant Real[7] data_y =   {0.0000000000000000E+00, 2.7827262318000000E-02, 1.1773081735099999E-01, 2.1767976495999999E-01, 2.3049575595899999E-01, 4.1264868481999999E-02, 0.0000000000000000E+00};
    constant Real[7] m_k =      {0.0000000000000000E+00, -4.3042454870000002E-03, 2.7656114134000001E-02, 7.4887293708000000E-02, -6.2894826385000005E-02, -1.7757016339000001E-01, 0.0000000000000000E+00};
    constant Real[7] iy_start = {0.0000000000000000E+00, 3.9237089677000000E-02, 2.8928309038500000E-01, 6.1093566216600004E-01, 8.4796858417200005E-01, 9.9430860177400004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0062118294654108E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 6.5000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 6.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 6.5000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>PX25</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    other<br>  Data taken from: Rubitherm datasheet - last access 2020-07-03.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-04-20.<br><br>
   See:<br>
    Barz, T., Bres, A., & Emhofer, J. (2022).
    slPCMlib: A Modelica Library for the Prediction of Effective 
    Thermal Material Properties of Solid/Liquid Phase Change  
    Materials (PCM). 
    In Proceedings of Asian Modelica Conference 2022 (pp. 63-74). 
    Linkoping University Electronic Press. 
    <a href>https://doi.org/10.3384/ecp19363</a>.
    </p>
    </blockquote>
    </p></html>",
    revisions="<html>
    <ul>
    <li>file creation date: 2023-04-20 </ul>
    </p></html>"));
end Rubitherm_PX25;