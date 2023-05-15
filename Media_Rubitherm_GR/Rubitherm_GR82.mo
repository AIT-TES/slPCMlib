
within slPCMlib.Media_Rubitherm_GR;
package Rubitherm_GR82 "Rubitherm GmbH, GR82; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "GR82";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4314999999999998E+02, 3.5714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.4414999999999998E+02, 3.5714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.3585493739764177E+04
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
    constant Real[11] data_x =   {7.0000000000000000E+01, 7.1375000000000000E+01, 7.2375000000000000E+01, 7.3125000000000000E+01, 7.3875000000000000E+01, 7.5875000000000000E+01, 7.8125000000000000E+01, 7.9625000000000000E+01, 8.0625000000000000E+01, 8.2375000000000000E+01, 8.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.2789762333000003E-02, 3.9638974883000000E-02, 6.7491011136000006E-02, 3.1513620653999999E-02, 9.6663400891999998E-02, 1.0396351917799999E-01, 1.0219093122000000E-01, 1.4048953308400000E-01, 6.0880709064999997E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -1.0990086970000000E-02, 5.1126434158000000E-02, -2.5055439245000002E-02, -1.5545076598000000E-02, 2.0901274681000000E-02, 5.3944272649999996E-03, -7.8636041340000001E-03, 5.8769312157999999E-02, -6.2224608404000001E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.4280460378000000E-02, 5.5326108227000001E-02, 9.9081668626000000E-02, 1.3577165522000001E-01, 2.5182851317200000E-01, 4.8413303688699999E-01, 6.4127350714900000E-01, 7.5708955705600001E-01, 9.6421825465300004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0002466655084963E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {7.1000000000000000E+01, 7.2625000000000000E+01, 7.4125000000000000E+01, 7.5125000000000000E+01, 7.5875000000000000E+01, 7.7375000000000000E+01, 8.0625000000000000E+01, 8.2375000000000000E+01, 8.4000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 3.6514663167000000E-02, 3.3826919980000002E-02, 1.1012690297600000E-01, 6.6096461829999995E-02, 1.1307843875100000E-01, 1.1382080820200000E-01, 1.0716546518099999E-01, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, -8.7321106539999995E-03, 3.8118997497999997E-02, -3.2578779160000000E-03, -2.3226109338000001E-02, -1.0920586764000000E-02, 1.1731392000000001E-05, -6.5256395905000000E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 3.1847630069000002E-02, 7.6178286944000001E-02, 1.5221915778199999E-01, 2.1978618367200001E-01, 3.5293852944900000E-01, 7.1495924343799999E-01, 9.2669413090999997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0081655490315686E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>GR82</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: other<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-03.<br><br>
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
end Rubitherm_GR82;