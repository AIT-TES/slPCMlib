
// within slPCMlib.Axiotherm_ATP;
package Axiotherm_ATP_28 "Axiotherm GmbH, ATP 28; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 28";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9514999999999998E+02, 3.0314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2900000000000000E+05
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
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {2.2000000000000000E+01, 2.3875000000000000E+01, 2.5625000000000000E+01, 2.6875000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.7875000000000000E+01, 2.8125000000000000E+01, 2.8375000000000000E+01, 2.8625000000000000E+01, 2.8875000000000000E+01, 2.9875000000000000E+01, 3.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 7.8525616080000003E-03, 1.9873480419999999E-02, 1.2887913952800001E-01, 3.7801006773700002E-01, 6.5512109776899996E-01, 7.9392770133500001E-01, 7.1846711025200005E-01, 4.7003867579699998E-01, 1.8446661365899999E-01, 7.1815860294000000E-02, 5.1426498500000003E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 6.1953566799999998E-04, 2.1366675100000002E-03, 2.4554169424100000E-01, 7.9016200967200001E-01, 9.1131856859899996E-01, 1.7441380308400001E-01, -8.9149613054800003E-01, -1.0680109875090000E+00, -8.9837540650600001E-01, -2.1459599032099999E-01, -5.2112990339999998E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 7.2278689259999999E-03, 3.1259222436000003E-02, 9.2942443796000002E-02, 2.0908330056200000E-01, 3.3844554962599999E-01, 5.2464082952900004E-01, 7.2053177755300002E-01, 8.7100524789300005E-01, 9.5247136160300006E-01, 9.8113406957100002E-01, 9.9997447590300004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0066288588739536E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.0000000000000000E+01, 2.4125000000000000E+01, 2.6125000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.7875000000000000E+01, 2.8375000000000000E+01, 2.8625000000000000E+01, 2.8875000000000000E+01, 2.9000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 8.5888347819999995E-03, 5.6764880959000001E-02, 4.0333718064399998E-01, 6.0129332361499999E-01, 6.8919023908400001E-01, 4.4650522901700002E-01, 2.2125760640400000E-01, 5.0446033793999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 3.1835651070000002E-03, 7.6547263429999998E-02, 6.0488047611700002E-01, 6.5716115403900000E-01, 7.0426645184999997E-02, -8.4806664313600000E-01, -8.1820364558699998E-01, -5.8523076741699998E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.3414136496000000E-02, 5.4975902312000002E-02, 2.7729065147499998E-01, 4.0462728796000003E-01, 5.7165657790900004E-01, 8.7962565009000004E-01, 9.6429026989199995E-01, 9.9757040729799995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0162012214951133E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.4400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.6000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 28</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATP_28;