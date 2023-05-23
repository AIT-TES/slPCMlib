
within slPCMlib.Media_Climator;
package ClimSel_C70 "Climator Sweden AB, ClimSel C70; data taken from: Climator Sweden AB datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C70";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4414999999999998E+02, 3.5114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.3814999999999998E+02, 3.4414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.3500000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.8000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.0028769669063076E+05
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
    constant Real[12] data_x =   {7.1000000000000000E+01, 7.2375000000000000E+01, 7.3125000000000000E+01, 7.3375000000000000E+01, 7.3875000000000000E+01, 7.4375000000000000E+01, 7.4625000000000000E+01, 7.5375000000000000E+01, 7.5625000000000000E+01, 7.6875000000000000E+01, 7.7875000000000000E+01, 7.8000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.4922847730000000E-03, 4.6081877240000003E-03, 2.7968601422999999E-02, 4.0115561950900003E-01, 5.8681812320100002E-01, 4.5208653445500002E-01, 2.6897599137400002E-01, 2.5408482492099999E-01, 2.0717642291999998E-02, 2.2918280299999999E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -1.4007015330999999E-02, 1.5920334150999998E-02, 2.1651454191200001E-01, 8.0881848015699997E-01, -3.9266469649699998E-01, -3.9192593516899998E-01, -9.8633256276999995E-02, -1.1252071945700000E-01, -5.5952358133000001E-02, -2.3750360040000001E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 5.2970365549999997E-03, 7.3075358519999996E-03, 1.0335877097000000E-02, 1.0530875690600000E-01, 3.7742331654500000E-01, 5.0732562246900004E-01, 7.6406110056999998E-01, 8.2953774309600004E-01, 9.9397813620900000E-01, 9.9998876484400001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0003316778541247E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {6.5000000000000000E+01, 6.6375000000000000E+01, 6.7625000000000000E+01, 6.8625000000000000E+01, 6.9375000000000000E+01, 6.9625000000000000E+01, 6.9875000000000000E+01, 7.0375000000000000E+01, 7.0625000000000000E+01, 7.0875000000000000E+01, 7.1000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.4642487907999999E-02, 5.4447628599000000E-02, 1.7149534515900000E-01, 4.0294881363800000E-01, 5.6115005944899998E-01, 6.1719515627499999E-01, 3.8231107856000002E-01, 1.8731054836800001E-01, 4.2449598300000002E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -1.1251663698000000E-02, 9.6247692763999998E-02, 1.3471801904200001E-01, 4.8339467764700000E-01, 5.1611060677599996E-01, -3.2589249894999998E-02, -7.4026101596600002E-01, -7.1442878933400000E-01, -4.9038916825000001E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.8738138402999999E-02, 5.4217032186999997E-02, 1.6412169477100000E-01, 3.6344619382600002E-01, 4.8394059039800003E-01, 6.3428174502099999E-01, 8.9923664994300001E-01, 9.7039482998200000E-01, 9.9798287588199996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0012666645555288E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.7000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.7000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 8.1000000000000005E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.1000000000000005E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C70</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C70;