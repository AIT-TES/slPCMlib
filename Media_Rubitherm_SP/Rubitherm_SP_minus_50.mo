
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_50 "Rubitherm GmbH, SP-50; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-50";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.1614999999999998E+02, 2.2614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.1414999999999998E+02, 2.2314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8836706993118860E+05
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
    constant Real[8] data_x =   {-5.7000000000000000E+01, -5.4125000000000000E+01, -5.2875000000000000E+01, -5.2625000000000000E+01, -5.1625000000000000E+01, -4.9125000000000000E+01, -4.7625000000000000E+01, -4.7000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 1.3228791352000000E-02, 7.3866917675999996E-02, 1.2373684836000000E-01, 4.0593159064500001E-01, 3.4901876685000000E-02, 1.8394885970000000E-03, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 7.9951688440000005E-03, 1.3396290683400000E-01, 3.4486505582100002E-01, 5.2484900180000002E-02, -6.4342346519999996E-02, -5.3648855010000000E-03, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 1.3602740838000000E-02, 5.1898577938999997E-02, 7.5663852404000004E-02, 3.6686344185300002E-01, 9.8298516685199999E-01, 9.9959702986300003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0069168820392778E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-5.9000000000000000E+01, -5.7125000000000000E+01, -5.5625000000000000E+01, -5.5125000000000000E+01, -5.4625000000000000E+01, -5.3625000000000000E+01, -5.2375000000000000E+01, -5.2125000000000000E+01, -5.1625000000000000E+01, -5.1125000000000000E+01, -5.0000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.5041933882000000E-02, 8.7039125460999997E-02, 7.9826725853000005E-02, 1.4214622716200001E-01, 3.8446630548999999E-01, 1.2458160716500000E-01, 8.5984802550999995E-02, 6.7874898764000005E-02, 8.0380348323999998E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 2.3761562145000002E-02, -2.1532376783000001E-02, 3.9044142509999998E-02, 2.9781080886700001E-01, 3.7393968131000002E-02, -2.9815800146200000E-01, -9.1381617644000004E-02, 2.4601662561999998E-02, -3.4081267078999999E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.6721458331000000E-02, 1.1042952976800000E-01, 1.5138867712400000E-01, 2.0211600120300000E-01, 4.9067932518300000E-01, 8.5704020280899995E-01, 8.8259897686599997E-01, 9.1909731344900003E-01, 9.5786133341699997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0124756441208840E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.1000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP-50</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-06-03.<br><br>
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
end Rubitherm_SP_minus_50;