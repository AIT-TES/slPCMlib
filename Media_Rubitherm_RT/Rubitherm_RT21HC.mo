within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT21HC "Rubitherm GmbH, RT21HC; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT21HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8514999999999998E+02, 2.9714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8514999999999998E+02, 2.9614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.6600000000000000E+05
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
    constant Real[8] data_x =   {1.2000000000000000E+01, 1.3375000000000000E+01, 1.5625000000000000E+01, 1.7625000000000000E+01, 2.0375000000000000E+01, 2.1375000000000000E+01, 2.2625000000000000E+01, 2.4000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 2.6000966241000002E-02, 3.1466973209000000E-02, 3.8374552392000001E-02, 1.8144523797900000E-01, 2.5776447713000000E-01, 1.8189968543300000E-01, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, -1.3306199038000000E-02, 2.7906012524999999E-02, 2.6358368603000001E-02, 1.1647872033400000E-01, -1.6427528694999999E-02, -1.0171862460400000E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 1.9800644155000002E-02, 6.6659953278999998E-02, 1.3641341520600000E-01, 3.7976403260000002E-01, 6.0846426037300005E-01, 8.9190585680800005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9141604263000827E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {1.2000000000000000E+01, 1.3625000000000000E+01, 1.5875000000000000E+01, 1.7625000000000000E+01, 1.9625000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.1125000000000000E+01, 2.1625000000000000E+01, 2.1875000000000000E+01, 2.3000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.9335111511000001E-02, 2.9462781957000000E-02, 5.8358363295999999E-02, 2.0549759911900001E-01, 3.8277060340100000E-01, 4.8821500169700000E-01, 3.2215254235399998E-01, 2.3416434872999999E-02, 3.8000416500000002E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -5.3956459310000003E-03, 8.5854577980000006E-03, 1.0333249144000000E-02, 1.2297325013000000E-01, 3.2928985060999999E-01, 3.2344820008000003E-01, -6.6281604675899997E-01, -2.0465515870699999E-01, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.7058388471000000E-02, 6.6525452087000003E-02, 1.4365213687299999E-01, 3.7212161356800000E-01, 5.8506447584099996E-01, 6.9500761526800003E-01, 9.2027659362100001E-01, 9.9785734624900002E-01, 9.9978420733499995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0095452609347140E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT21HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT21HC;
