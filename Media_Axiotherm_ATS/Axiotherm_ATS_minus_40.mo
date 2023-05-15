
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_minus_40 "Axiotherm GmbH, ATS -40; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -40";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.2814999999999998E+02, 2.3914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.2614999999999998E+02, 2.3514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.3659830034074833E+05
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
    constant Real[10] data_x =   {-4.5000000000000000E+01, -4.2375000000000000E+01, -3.9375000000000000E+01, -3.8375000000000000E+01, -3.7625000000000000E+01, -3.7375000000000000E+01, -3.6625000000000000E+01, -3.6375000000000000E+01, -3.5125000000000000E+01, -3.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.2044831941000001E-02, 1.2809182429600000E-01, 1.8978446565099999E-01, 2.8297338644100001E-01, 3.4475414556299999E-01, 2.6601823720099999E-01, 1.7509493172099999E-01, 1.5954250622999998E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -5.5556174499999996E-03, 8.3508979178000001E-02, 4.9636768871999998E-02, 1.9317350367399999E-01, 2.0062447971299999E-01, -3.2623077529000000E-01, -3.0127347944700000E-01, -3.4035567364999998E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.2123556025999998E-02, 1.9052799421500000E-01, 3.5228668024300003E-01, 5.2284033067199998E-01, 6.0126642594000002E-01, 8.5499904785600001E-01, 9.1000747925000003E-01, 9.9461549410899996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9998675119596525E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {-4.7000000000000000E+01, -4.3875000000000000E+01, -4.2375000000000000E+01, -4.1375000000000000E+01, -4.0375000000000000E+01, -3.9375000000000000E+01, -3.9125000000000000E+01, -3.8375000000000000E+01, -3.8125000000000000E+01, -3.8000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.2025659169999998E-02, 5.6480533271000000E-02, 7.3446263780999999E-02, 1.9316398776999999E-01, 4.8565573786799998E-01, 4.8319348353500002E-01, 1.2706122793600000E-01, 2.8345357440999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 8.6066598340000001E-03, 4.0528173920999999E-02, 4.4590727140000003E-03, 2.5877807383500001E-01, 3.0450591853900000E-01, -2.0406934466900001E-01, -5.3911460158500002E-01, -3.2399651401399998E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 4.2733784791999997E-02, 1.0270406818199999E-01, 1.7019595939500001E-01, 2.8152060803800000E-01, 6.1476330699500004E-01, 7.3764930657600003E-01, 9.8048288353000002E-01, 9.9865976304399995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9297820116056568E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.5430000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.4500000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4500000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS -40</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATS_minus_40;