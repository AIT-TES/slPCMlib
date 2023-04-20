
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT22HC "Rubitherm GmbH, RT22HC; data taken from: Rubitherm datasheet; last access: 2020-10-09."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT22HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 2.9814999999999998E+02}
             "temperature range melting {startT, endT}";
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {0.0, 0.0}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.5443129235800411E+05
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
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {1.5000000000000000E+01, 1.7625000000000000E+01, 1.9125000000000000E+01, 2.0375000000000000E+01, 2.1625000000000000E+01, 2.2375000000000000E+01, 2.3125000000000000E+01, 2.5000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 3.3562315386999997E-02, 9.2793335264999999E-02, 1.4765658815400001E-01, 3.2134128824000002E-01, 2.7901425337600000E-01, 1.1472165392899999E-01, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, -3.0372669669999999E-03, 4.5384038138000002E-02, 9.2634268014000001E-02, 1.6163483855999999E-01, -2.3592068740100000E-01, -1.6498345523100000E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 4.5473484958000000E-02, 1.3056038965200001E-01, 2.7367859326999999E-01, 5.5582544532000000E-01, 7.9788489328300005E-01, 9.4119866989400003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9298804514950278E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Integer  len_x =   1;
    constant Real[1] data_x =   {0.0};
    constant Real[1] data_y =   {0.0};
    constant Real[1] m_k =      {0.0};
    constant Real[1] iy_start = {0.0};
    constant Real    iy_scaler = 0.0;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 7.6000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT22HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
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
end Rubitherm_RT22HC;