
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_84 "Axiotherm GmbH, ATS 84; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 84";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4814999999999998E+02, 3.6114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.4714999999999998E+02, 3.5814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.0900000000000000E+05
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
    constant Real[15] data_x =   {7.5000000000000000E+01, 7.7375000000000000E+01, 7.8875000000000000E+01, 8.0375000000000000E+01, 8.1625000000000000E+01, 8.3125000000000000E+01, 8.3625000000000000E+01, 8.3875000000000000E+01, 8.4125000000000000E+01, 8.4625000000000000E+01, 8.4875000000000000E+01, 8.5625000000000000E+01, 8.6375000000000000E+01, 8.7625000000000000E+01, 8.8000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.5005167597000000E-02, 2.7593866238999998E-02, 3.3284713291000002E-02, 8.9528749463999996E-02, 3.2035189712800000E-01, 4.7053682623800003E-01, 4.6397037260899998E-01, 3.7352409745600001E-01, 9.5912895932000003E-02, 3.6881637129000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 6.3243272629999996E-03, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, -1.9015732673999999E-02, 2.6428317061999999E-02, 4.9579753799000002E-02, 3.8296248182000003E-02, 2.6292417963600001E-01, 2.8437862569000000E-01, -2.3902500780499999E-01, -5.5435550965199998E-01, -5.2451928789799995E-01, -1.2221315331700000E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, -1.8221445507000000E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 2.6628966649999999E-02, 4.9945357740000002E-02, 9.1065656765000003E-02, 1.6891889313100000E-01, 4.3294199418300000E-01, 6.2927308120500003E-01, 7.4824044577600002E-01, 8.5406072733600003E-01, 9.7023970080099997E-01, 9.8467425939599995E-01, 9.9273735755100001E-01, 9.9273735755100001E-01, 9.9903237436699999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9521416715994548E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {7.4000000000000000E+01, 7.6125000000000000E+01, 7.7875000000000000E+01, 8.0125000000000000E+01, 8.1625000000000000E+01, 8.2375000000000000E+01, 8.2625000000000000E+01, 8.3125000000000000E+01, 8.3625000000000000E+01, 8.3875000000000000E+01, 8.5000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 8.5745534850000006E-03, 2.6396374811000000E-02, 4.7366703151999999E-02, 1.9499551475599999E-01, 4.2565178124399999E-01, 5.8751748958299999E-01, 4.2022998577699999E-01, 3.2881997162999999E-02, 5.6645771000000001E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -3.3441372639999999E-03, 2.2602100366000001E-02, 4.6707785463000001E-02, 1.1715642763299999E-01, 4.7820496735000001E-01, 4.8795967835900000E-01, -8.3534903046999998E-01, -2.9602710267299998E-01, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.0341537272000000E-02, 3.4256199463999999E-02, 1.0687814390200000E-01, 2.7499636765899999E-01, 4.9024607041500001E-01, 6.1650772683999999E-01, 8.9527676905300002E-01, 9.9704991527700004E-01, 9.9968220743100000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9736406879309747E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.7550000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.6500000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATS 84</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATS_84;