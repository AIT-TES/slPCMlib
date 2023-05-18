
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_minus_21 "Axiotherm GmbH, ATS -21; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5014999999999998E+02, 2.5814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4114999999999998E+02, 2.5714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.7409432142179844E+05
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
    constant Real[12] data_x =   {-2.3000000000000000E+01, -2.1875000000000000E+01, -2.1625000000000000E+01, -2.1125000000000000E+01, -2.0875000000000000E+01, -2.0625000000000000E+01, -2.0375000000000000E+01, -1.9875000000000000E+01, -1.9125000000000000E+01, -1.7875000000000000E+01, -1.5125000000000000E+01, -1.5000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.6219864320000000E-03, 9.7127098820000007E-03, 5.1487218189100004E-01, 7.6873786316899995E-01, 8.0193259614500001E-01, 5.9101568063400001E-01, 2.6643034762300000E-01, 3.6340471114000003E-02, 1.8412213874000001E-02, 4.7273761600000002E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 3.1103104680000001E-03, 7.0683679605999999E-02, 1.1759792039779999E+00, 8.3256525496500000E-01, -7.1924948041299996E-01, -7.3004096780799999E-01, -5.1736484791100001E-01, -7.1554499593000004E-02, -1.3502163653999999E-02, -5.4117154380000001E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 5.8390414599999996E-04, 1.6480262790000001E-03, 1.0968905257100001E-01, 2.7181156854299998E-01, 4.7607988274099999E-01, 6.5012863609899996E-01, 8.5990754014100002E-01, 9.5248222599200005E-01, 9.7912447094900001E-01, 9.9997751667799994E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9927666375219981E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {-3.2000000000000000E+01, -2.8375000000000000E+01, -2.5875000000000000E+01, -2.4375000000000000E+01, -2.3125000000000000E+01, -2.2625000000000000E+01, -2.2375000000000000E+01, -2.1625000000000000E+01, -2.1375000000000000E+01, -2.1125000000000000E+01, -2.0875000000000000E+01, -1.9375000000000000E+01, -1.7625000000000000E+01, -1.6000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 4.6447560927000003E-02, 5.1406397734999998E-02, 8.4890390041999997E-02, 1.9368862261399999E-01, 3.4492104424100001E-01, 4.7735582841599999E-01, 1.3386759786399999E-01, 1.2829409602000000E-02, 1.3085152119999999E-03, 2.3386994399999999E-04, 5.6245300919999996E-03, 5.1289771389999997E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 1.9879817673000000E-02, 9.7326884410000002E-03, 2.3397900520999999E-02, 1.8956773680100000E-01, 3.8986687179300000E-01, 3.9597335687899998E-01, -6.7992407180000003E-01, -1.0310662312100000E-01, -4.7725124920000002E-03, -3.0648688300000002E-04, -9.0548942679999998E-03, 5.8368562039999997E-03, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 6.3118604970000006E-02, 1.9215582375800000E-01, 2.9293680422700002E-01, 4.4712647542299999E-01, 5.7907315202800003E-01, 6.8298134591699999E-01, 9.6576721188600001E-01, 9.8127248855000004E-01, 9.8254168493799998E-01, 9.8271312887399997E-01, 9.8881510516899995E-01, 9.9448699276399999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0112443589550637E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.2340000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1600000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS -21</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-05-18.<br><br>
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
    <li>file creation date: 2023-05-18 </ul>
    </p></html>"));
end Axiotherm_ATS_minus_21;