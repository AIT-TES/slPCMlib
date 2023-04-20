
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT8 "Rubitherm GmbH, RT8; data taken from: Rubitherm datasheet; last access: 2022-12-02."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT8";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7214999999999998E+02, 2.8814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7014999999999998E+02, 2.8314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.5361354730788307E+05
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
    constant Real[8] data_x =   {-1.0000000000000000E+00, 2.1250000000000000E+00, 6.1250000000000000E+00, 9.1250000000000000E+00, 1.0375000000000000E+01, 1.2625000000000000E+01, 1.3625000000000000E+01, 1.5000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 3.4379633009000002E-02, 1.2927198596500000E-01, 1.2168657420600000E-01, 4.4492596427000000E-02, 9.4567234579999992E-03, 6.5933567609999999E-03, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 1.6828670396999999E-02, 3.3402352215000000E-02, -8.4717436184999995E-02, -1.1030439329000000E-02, -2.5858753351000002E-02, 3.5577838119999999E-03, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 4.0750829149000002E-02, 3.5150629703800002E-01, 8.2499099382499996E-01, 9.2097266433600000E-01, 9.8913887601700001E-01, 9.9481390090900002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0181860469087114E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {-3.0000000000000000E+00, 1.6250000000000000E+00, 3.3750000000000000E+00, 5.3750000000000000E+00, 7.1250000000000000E+00, 8.6250000000000000E+00, 9.3750000000000000E+00, 1.0000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 4.5755191868999998E-02, 9.4197911092000000E-02, 1.2254709247499999E-01, 1.1989051652400000E-01, 1.8827432773000000E-01, 8.9917709496000003E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 8.2543620809999996E-03, -3.1961992149999998E-03, 2.0457179552999999E-02, 1.5673576036000000E-02, 5.3565469070000003E-02, -1.9501725888200000E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 9.1055454433999997E-02, 2.1638219851400001E-01, 4.2515195672900002E-01, 6.3841294126399994E-01, 8.6233446983999995E-01, 9.7825838841699997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9956532832682632E-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT8</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2022-12-02.<br><br>
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
end Rubitherm_RT8;