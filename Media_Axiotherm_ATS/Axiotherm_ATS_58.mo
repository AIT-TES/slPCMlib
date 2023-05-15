
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_58 "Axiotherm GmbH, ATS 58; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 58";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2714999999999998E+02, 3.3414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2514999999999998E+02, 3.3114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9600000000000000E+05
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
    constant Real[11] data_x =   {5.4000000000000000E+01, 5.5625000000000000E+01, 5.6125000000000000E+01, 5.6375000000000000E+01, 5.6625000000000000E+01, 5.6875000000000000E+01, 5.7375000000000000E+01, 5.7625000000000000E+01, 5.8375000000000000E+01, 5.9625000000000000E+01, 6.1000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 5.1278845419999997E-03, 2.6264125547999999E-02, 1.0309374749800000E-01, 4.2913035506000002E-01, 7.7967882024099999E-01, 7.6558180532200004E-01, 4.1026049559099997E-01, 6.5717969227999995E-02, 7.6068799349999997E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -1.1240531450999999E-02, 8.4659611908999996E-02, 1.0113426706830000E+00, 1.3577457356410001E+00, 1.3401131026260000E+00, -1.1414675772880001E+00, -9.6337742246400004E-01, -1.3711658950799999E-01, -7.7048037289999998E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 6.6627493829999997E-03, 1.2532956136000001E-02, 2.3915236037000000E-02, 8.8861711755999995E-02, 2.4057478877499999E-01, 6.8009628332600003E-01, 8.2665142942400005E-01, 9.6689314524500003E-01, 9.9597036377199999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0034399313001168E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {5.2000000000000000E+01, 5.4625000000000000E+01, 5.5125000000000000E+01, 5.6125000000000000E+01, 5.6625000000000000E+01, 5.6875000000000000E+01, 5.7625000000000000E+01, 5.7875000000000000E+01, 5.8000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 6.2846459906000005E-02, 1.0945452261300000E-01, 4.1658070779099998E-01, 5.4441559752699997E-01, 5.1665072517699995E-01, 1.2656057918100000E-01, 2.8008403038999999E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 3.4075511048000003E-02, 1.4299919706899999E-01, 3.3033222723099998E-01, 1.5400634074700001E-01, -3.0600686807999999E-01, -5.6850704428599996E-01, -3.1849402553799999E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 6.2734293048000003E-02, 1.0342038664800000E-01, 3.5009991110700001E-01, 5.9330567624700004E-01, 7.2793808499099999E-01, 9.8070208305100004E-01, 9.9866810588199995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9706149673143596E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3620000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.2800000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.2800000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS 58</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
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
end Axiotherm_ATS_58;