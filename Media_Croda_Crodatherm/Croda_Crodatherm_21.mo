within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_21 "Croda International Plc, Crodatherm 21; data taken from: Croda datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 2.9714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8814999999999998E+02, 2.9514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.9000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8667568345569796E+05
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
    constant Real[10] data_x =   {1.5000000000000000E+01, 1.8625000000000000E+01, 1.9125000000000000E+01, 1.9375000000000000E+01, 1.9875000000000000E+01, 2.0125000000000000E+01, 2.0375000000000000E+01, 2.1375000000000000E+01, 2.2625000000000000E+01, 2.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.9393042776000000E-02, 5.5088515254999998E-02, 1.2832964769700000E-01, 5.1686559000700005E-01, 6.3680441831300005E-01, 6.3291322436200004E-01, 1.6510858324800001E-01, 1.0182779104000001E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -3.7446738680000000E-03, 1.3369897946800000E-01, 7.4484991905200004E-01, 7.7600182167700005E-01, 2.9036357783500000E-01, -4.9861755767400001E-01, -3.8114760620200000E-01, -1.4765130035999999E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.9222625051000000E-02, 5.4968413627000001E-02, 7.4698582757000004E-02, 2.3523429173700000E-01, 3.8186818701899999E-01, 5.4457652240900001E-01, 9.3352181655900002E-01, 9.9532893637200004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9928974411366189E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {1.5000000000000000E+01, 1.7625000000000000E+01, 1.8625000000000000E+01, 1.9125000000000000E+01, 1.9375000000000000E+01, 1.9625000000000000E+01, 1.9875000000000000E+01, 2.0125000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.0875000000000000E+01, 2.2000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 3.4189430173000003E-02, 3.5621096641000001E-02, 1.1985557018300000E-01, 3.1494624994699999E-01, 8.8264538611300003E-01, 1.1532800505880001E+00, 8.6759307502799998E-01, 3.0945793646800002E-01, 6.1826593159999997E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.7850889010000001E-02, -1.3089761288000001E-02, 3.2067025716800002E-01, 1.2326948678450000E+00, 1.7941676371799999E+00, -7.8196918594000006E-02, -1.8901456601780000E+00, -1.7448767844149999E+00, -5.2876587374000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 2.8901575620000002E-02, 6.7245692760999995E-02, 9.9184125866999998E-02, 1.4881934783400000E-01, 2.9569789848000000E-01, 5.6012759123699996E-01, 8.2235952424100001E-01, 9.6883794228900000E-01, 9.9950221423700003E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0007081161639779E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.9100000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 1.7999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 1.4999999999999999E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 21</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_21;
