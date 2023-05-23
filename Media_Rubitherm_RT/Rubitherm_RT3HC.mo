
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT3HC "Rubitherm GmbH, RT3HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT3HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6614999999999998E+02, 2.7914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6614999999999998E+02, 2.7714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.6500000000000000E+05
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-7.0000000000000000E+00, -1.8750000000000000E+00, 6.2500000000000000E-01, 1.3750000000000000E+00, 1.6250000000000000E+00, 2.3750000000000000E+00, 3.8750000000000000E+00, 5.3750000000000000E+00, 6.0000000000000000E+00};
    constant Real[9] data_y =   {0.0000000000000000E+00, 2.4093492445000000E-02, 9.2152878514999995E-02, 2.0084139277700000E-01, 2.7615338527900002E-01, 3.5942452611199999E-01, 8.5588209070999996E-02, 1.0043675723000000E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 1.0644343460000001E-03, 6.2348439440999999E-02, 2.5215259729900003E-01, 2.6928240617200000E-01, -1.2225355423500001E-01, -1.0691519691199999E-01, -2.8008969383000001E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 6.0619439859000000E-02, 1.7631748295899999E-01, 2.7934933449900001E-01, 3.4009671880100001E-01, 6.0201850871700002E-01, 9.3963957289900002E-01, 9.9772775751900000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0203620082751961E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {-7.0000000000000000E+00, -5.1250000000000000E+00, -3.6250000000000000E+00, -1.3750000000000000E+00, 2.3750000000000000E+00, 2.6250000000000000E+00, 2.8750000000000000E+00, 3.6250000000000000E+00, 3.8750000000000000E+00, 4.0000000000000000E+00};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.0159544190999999E-02, 1.3909200252000001E-02, 3.2081645237000002E-02, 3.4017269200200001E-01, 3.7359697218999999E-01, 3.5149600934500003E-01, 8.4975704026000007E-02, 1.8777265434999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.0346397718000000E-02, 1.9065467278999999E-02, 5.6359288399999996E-04, 1.3245632352600001E-01, 1.3228313518300000E-01, -2.1899014638500000E-01, -3.7993393273499998E-01, -2.1331906079900001E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.5732286186000002E-02, 3.9443867060999997E-02, 9.8478282437000006E-02, 6.3723208864199998E-01, 7.2568887313999997E-01, 8.1736189185899999E-01, 9.8711433114799996E-01, 9.9911186422900000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9142224333878071E-01;
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
    lambda := 2.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT3HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-09-30.<br><br>
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
end Rubitherm_RT3HC;