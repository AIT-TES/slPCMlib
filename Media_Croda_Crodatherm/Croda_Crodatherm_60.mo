
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_60 "Croda International Plc, Crodatherm 60; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 60";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2514999999999998E+02, 3.3514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2614999999999998E+02, 3.3414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8138933479725727E+05
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
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {5.2000000000000000E+01, 5.3875000000000000E+01, 5.5875000000000000E+01, 5.8375000000000000E+01, 5.9375000000000000E+01, 5.9625000000000000E+01, 5.9875000000000000E+01, 6.0125000000000000E+01, 6.0625000000000000E+01, 6.0875000000000000E+01, 6.1875000000000000E+01, 6.2000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.9305952070000001E-03, 2.8793211897000001E-02, 1.5288042041600000E-01, 4.4129167645900003E-01, 5.7320968647299997E-01, 5.7293555758000003E-01, 4.4573310380900000E-01, 6.4523329382999997E-02, 1.9277726633000002E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.7847050709999999E-03, 3.6831431440000001E-03, 1.5996592596600001E-01, 4.2360052758799999E-01, 4.2531354321499998E-01, -3.1936530474800001E-01, -7.6246746918700004E-01, -4.5457529229900001E-01, -5.3394378530000002E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.0035891750000000E-03, 3.1718287568999998E-02, 1.7880356014000001E-01, 4.5654577810300001E-01, 5.8455974508300002E-01, 7.3311080844999998E-01, 8.6398951615499997E-01, 9.8629546015900005E-01, 9.9476114109799996E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0095440283550956E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {5.3000000000000000E+01, 5.5125000000000000E+01, 5.6625000000000000E+01, 5.8875000000000000E+01, 5.9375000000000000E+01, 5.9625000000000000E+01, 5.9875000000000000E+01, 6.0625000000000000E+01, 6.0875000000000000E+01, 6.1000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 6.3432678439999998E-03, 3.3872502908000000E-02, 2.6108126975200002E-01, 4.1629617751400000E-01, 5.2278883890700001E-01, 5.3849709318600003E-01, 1.4896172785699999E-01, 3.3418792730000001E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 7.7378183000000001E-03, 8.3389853049999994E-03, 2.4470169604799999E-01, 3.5640300562200000E-01, 3.6216157429099999E-01, -1.5928898240799999E-01, -6.0520183962200003E-01, -3.8341339202199998E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.8757720279999999E-03, 3.4300161753999997E-02, 2.6930640564699998E-01, 4.3840951945399997E-01, 5.5723079734799996E-01, 6.9429813306400001E-01, 9.7647800429700005E-01, 9.9839071129699997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0124889051392110E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.2200000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.2100000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.8999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.2100000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 60</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_60;