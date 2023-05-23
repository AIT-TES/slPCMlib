
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_28 "Rubitherm GmbH, SP-28; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-28";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4114999999999998E+02, 2.4914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4014999999999998E+02, 2.5114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4593246957917645E+05
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
    constant Real[11] data_x =   {-3.2000000000000000E+01, -3.0375000000000000E+01, -2.9875000000000000E+01, -2.9625000000000000E+01, -2.9125000000000000E+01, -2.8625000000000000E+01, -2.7875000000000000E+01, -2.7375000000000000E+01, -2.7125000000000000E+01, -2.5875000000000000E+01, -2.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.0410867530000002E-03, 1.6017176256999999E-02, 5.7100122263000000E-02, 4.1646152469800002E-01, 6.5780981436899999E-01, 3.7949084013200002E-01, 1.2456188220400000E-01, 6.4515691651000004E-02, 6.2516928250000003E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -2.9417689230000002E-03, 5.0562068730999998E-02, 3.5352940851800002E-01, 7.6685217957800000E-01, -5.3295253154000001E-02, -5.2891255029399997E-01, -4.0137349446100001E-01, -1.4407885530600001E-01, -5.4293253069999997E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.1271588159999999E-03, 6.7875179430000000E-03, 1.4370889010000000E-02, 1.2446491864200000E-01, 4.1093753327100002E-01, 8.2339812821199998E-01, 9.4710764987800000E-01, 9.6946614254899999E-01, 9.9571742646700001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0028649146398316E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-3.3000000000000000E+01, -3.1375000000000000E+01, -3.0875000000000000E+01, -3.0625000000000000E+01, -2.9875000000000000E+01, -2.9375000000000000E+01, -2.8875000000000000E+01, -2.8375000000000000E+01, -2.8125000000000000E+01, -2.7375000000000000E+01, -2.4625000000000000E+01, -2.2875000000000000E+01, -2.2000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.0314971009000000E-02, 6.2787356104999997E-02, 1.1968686160800000E-01, 4.7952711720699998E-01, 6.1989870002299996E-01, 3.7044852940900003E-01, 3.5376722606000000E-02, 7.7347871180000003E-03, 5.1204785129999996E-03, 1.9113573246999999E-02, 7.2309594640000002E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 2.9790475715999999E-02, 1.3671987129800001E-01, 3.8194207300100003E-01, 4.6001622806999998E-01, 3.7234383124000001E-02, -7.4717683909199994E-01, -2.3742544409600000E-01, -2.5226431170000001E-02, 1.6446536833000000E-02, -2.9723266546000000E-02, -8.6714452640000005E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 9.9567886469999999E-03, 2.8516475719000001E-02, 5.0062253682999998E-02, 2.7124840923700000E-01, 5.5509329923899997E-01, 8.1918993087400005E-01, 9.1008421758199998E-01, 9.1437067861800003E-01, 9.1723980716800002E-01, 9.7969793669799998E-01, 9.9738805072000003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0006362468668479E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.2000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP-28</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP_minus_28;