
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT5HC "Rubitherm GmbH, RT5HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT5HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7414999999999998E+02, 2.8114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7414999999999998E+02, 2.7914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4100000000000000E+05
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
    constant Real[11] data_x =   {1.0000000000000000E+00, 3.1250000000000000E+00, 4.1250000000000000E+00, 5.1250000000000000E+00, 5.3750000000000000E+00, 5.6250000000000000E+00, 5.8750000000000000E+00, 6.1250000000000000E+00, 6.6250000000000000E+00, 6.8750000000000000E+00, 8.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.2775501183000000E-02, 4.7865682271999999E-02, 4.1420738304799998E-01, 5.8676117736699995E-01, 7.9488917762099998E-01, 7.6335417176200004E-01, 5.1317674039299999E-01, 1.4420638515999999E-02, 4.3537437200000001E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 7.3585690190000001E-03, 9.5368446664999995E-02, 6.1751441253700001E-01, 7.0782872055799995E-01, 6.9836151548699998E-01, -8.2429328322800000E-01, -1.1558656505800000E+00, -1.2715708610300000E-01, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.0863219687000001E-02, 3.3973670932000002E-02, 2.2250976775900000E-01, 3.4783296383899998E-01, 5.2152061438899999E-01, 7.2532518895200004E-01, 8.8748869901799998E-01, 9.9855260976700000E-01, 9.9975378064599996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0053951814138478E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {1.0000000000000000E+00, 2.6250000000000000E+00, 3.8750000000000000E+00, 4.1250000000000000E+00, 4.3750000000000000E+00, 4.6250000000000000E+00, 5.1250000000000000E+00, 5.3750000000000000E+00, 5.6250000000000000E+00, 5.8750000000000000E+00, 6.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 8.6811937619999998E-03, 7.6242049436999998E-02, 1.5101505840700000E-01, 3.0359699309700000E-01, 6.1486094322300000E-01, 8.4331127446800003E-01, 6.4030212697599997E-01, 3.3083308914800003E-01, 7.7108514220999999E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -2.3200930160000000E-03, 1.7607604260900001E-01, 3.5944227591400002E-01, 8.9930760540999999E-01, 1.0838611261640001E+00, -5.2309977066299995E-01, -1.1496412175120001E+00, -1.1226239639019999E+00, -9.0922187880000005E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 7.6079138369999997E-03, 3.7629521515999999E-02, 6.5240963483000000E-02, 1.1956918042400000E-01, 2.3407598483500000E-01, 6.3440756447500002E-01, 8.2421780471100003E-01, 9.4617275110499999E-01, 9.9634350013700002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0058041430994957E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT5HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2022-03-21.<br><br>
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
end Rubitherm_RT5HC;