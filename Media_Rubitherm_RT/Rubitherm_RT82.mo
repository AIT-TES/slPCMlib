
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT82 "Rubitherm GmbH, RT82; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT82";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4114999999999998E+02, 3.5914999999999998E+02}
             "temperature range melting {startT, endT}";
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {0.0, 0.0}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4400000000000000E+05
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
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {6.8000000000000000E+01, 7.0625000000000000E+01, 7.3375000000000000E+01, 7.4625000000000000E+01, 7.7875000000000000E+01, 7.9625000000000000E+01, 8.0625000000000000E+01, 8.2125000000000000E+01, 8.3375000000000000E+01, 8.4625000000000000E+01, 8.6000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.9313336652000000E-02, 5.9659016963000001E-02, 6.7921887441000001E-02, 8.8781562962000002E-02, 8.8796277756999997E-02, 1.0646731939400000E-01, 5.6550235897999998E-02, 2.0425876175999998E-02, 6.8682527219999998E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.7733607340000000E-03, -1.5068907664000000E-02, 2.9934653390999998E-02, 9.6892665820000005E-03, -1.6555086696999999E-02, 3.2437333082000003E-02, -2.7706784944000001E-02, -4.1022352136000001E-02, 4.2874916830000004E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 5.0621402335000003E-02, 1.9744128125200000E-01, 2.7137931237099999E-01, 5.4406313691399999E-01, 7.0627272792300000E-01, 7.9989755299700005E-01, 9.3354584259700002E-01, 9.8343005914799997E-01, 9.9459820378099995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0008095016327239E+00;
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
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.7000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.7000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT82</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-05-15.<br><br>
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
    <li>file creation date: 2023-05-15 </ul>
    </p></html>"));
end Rubitherm_RT82;