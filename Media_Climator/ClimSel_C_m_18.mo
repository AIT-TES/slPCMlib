
// within slPCMlib.Media_Climator;
package ClimSel_C_m_18 "Climator Sweden AB, ClimSel C-18; data taken from: Climator Sweden AB datasheet; last access: 2022-10-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C-18";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4914999999999998E+02, 2.5714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4714999999999998E+02, 2.5514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.8500000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.3720909961422652E+05
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
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {-2.4000000000000000E+01, -2.1125000000000000E+01, -1.9375000000000000E+01, -1.8875000000000000E+01, -1.8625000000000000E+01, -1.8125000000000000E+01, -1.7125000000000000E+01, -1.6375000000000000E+01, -1.6125000000000000E+01, -1.6000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 6.8374205020000000E-03, 1.2403491234000001E-02, 3.0564799507000001E-02, 7.2640043747999997E-02, 3.3383332637300001E-01, 5.8580690917599998E-01, 1.3950236538999999E-01, 3.0772974636999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.4180718580000000E-03, -1.0158529393000000E-02, 7.2670855891999994E-02, 3.1663683635300000E-01, 5.7119484606299997E-01, -3.8535077642900001E-01, -6.2060788176199999E-01, -3.4921497652099998E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 8.9460065350000003E-03, 2.8946370259000001E-02, 3.8058562082000001E-02, 4.9811991116000001E-02, 1.4714965680200001E-01, 6.9241033731900004E-01, 9.7843391319999995E-01, 9.9851580333699996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0106173959061922E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-2.6000000000000000E+01, -2.5125000000000000E+01, -2.4125000000000000E+01, -2.2625000000000000E+01, -2.1125000000000000E+01, -1.9875000000000000E+01, -1.9625000000000000E+01, -1.9375000000000000E+01, -1.8875000000000000E+01, -1.8625000000000000E+01, -1.8375000000000000E+01, -1.8125000000000000E+01, -1.8000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.2456204840000001E-03, 4.0762906125999997E-02, 5.2132513563999998E-02, 9.6750079143000003E-02, 2.1759411189700001E-01, 3.1351647618900003E-01, 4.7186581907300001E-01, 4.9811446653300001E-01, 3.5452568957500002E-01, 1.7594593585700000E-01, 4.0147352791999998E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 6.5191402820000003E-03, 5.7931842062000000E-02, 4.0073031347999999E-02, 2.6610434854000001E-02, 2.7236583845500001E-01, 4.6156009793199998E-01, 5.1384024431499997E-01, -4.8479041877899998E-01, -6.7172431368399999E-01, -6.4815149802799998E-01, -4.6602174809699998E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 5.6779183400000003E-04, 1.7826192954000002E-02, 9.1009667895999999E-02, 2.0545134196100001E-01, 3.7028505349500002E-01, 4.3583483106300003E-01, 5.3395437885100006E-01, 7.9784339177100005E-01, 9.0563767762199998E-01, 9.7197194793700004E-01, 9.9809333294000002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0022375097312768E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.1699999999999999E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C-18</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    macroencapsulation<br>  Data taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C_m_18;