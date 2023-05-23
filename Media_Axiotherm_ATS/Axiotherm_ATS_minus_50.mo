
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_minus_50 "Axiotherm GmbH, ATS -50; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -50";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.1914999999999998E+02, 2.2714999999999998E+02}
             "temperature range melting {startT, endT}";
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {0.0, 0.0}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7724198155780637E+05
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
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {-5.4000000000000000E+01, -5.2875000000000000E+01, -5.1625000000000000E+01, -5.1125000000000000E+01, -5.0875000000000000E+01, -5.0625000000000000E+01, -5.0375000000000000E+01, -5.0125000000000000E+01, -4.9875000000000000E+01, -4.9625000000000000E+01, -4.9375000000000000E+01, -4.9125000000000000E+01, -4.8875000000000000E+01, -4.7875000000000000E+01, -4.6000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 3.2592186440000000E-03, 1.7590683797999999E-02, 1.9387307659999999E-02, 5.5446082297000002E-02, 1.9569670286800001E-01, 7.3423864300300001E-01, 1.1474692295240001E+00, 1.0552936610460000E+00, 5.4268044999099996E-01, 7.9303264008000002E-02, 1.2024547734000000E-02, 2.3073203540000002E-03, 9.6375145129999999E-03, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 8.2261164080000004E-03, -1.3000169993999999E-02, 4.2558279430999997E-02, 1.7216923579000001E-01, 1.2842814236869999E+00, 1.9332512445940000E+00, 1.1887299803619999E+00, -1.9011625157520000E+00, -1.9490455555540001E+00, -9.5163916810100002E-01, -4.3812160798000002E-02, -6.1824865749999999E-03, 8.9267429940000000E-03, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 9.6598724599999996E-04, 1.6765513074000000E-02, 2.4854845917999999E-02, 3.3536433863000002E-02, 5.9144320463999998E-02, 1.7203832348199999E-01, 4.1119760054600002E-01, 7.0271913525300000E-01, 9.0277223624400005E-01, 9.7534603411400000E-01, 9.8203564915599995E-01, 9.8363159907700004E-01, 9.8834625610100002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0002847368721755E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Integer  len_x =   1;
    constant Real[1] data_x =   {0.0};
    constant Real[1] data_y =   {0.0};
    constant Real[1] m_k =      {0.0};
    constant Real[1] iy_start = {0.0};
    constant Real    iy_scaler = 0.0;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3780000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.2400000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATS -50</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
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
end Axiotherm_ATS_minus_50;