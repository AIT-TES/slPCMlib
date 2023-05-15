
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT18HC "Rubitherm GmbH, RT18HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT18HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8614999999999998E+02, 2.9314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8514999999999998E+02, 2.9214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.3282340012643728E+05
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
    constant Real[12] data_x =   {1.3000000000000000E+01, 1.5625000000000000E+01, 1.6875000000000000E+01, 1.7375000000000000E+01, 1.7625000000000000E+01, 1.7875000000000000E+01, 1.8125000000000000E+01, 1.8375000000000000E+01, 1.8625000000000000E+01, 1.8875000000000000E+01, 1.9875000000000000E+01, 2.0000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.9448945326999999E-02, 1.5081645704999999E-01, 4.3414964352899998E-01, 7.4717518054999998E-01, 8.6513664324999995E-01, 7.2373377804899997E-01, 3.8919691973100001E-01, 9.1244723952000006E-02, 2.5086506498000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.0539198824000000E-02, 2.8339105499000000E-01, 8.6719069089400003E-01, 9.8690855125700006E-01, -4.9461870980000003E-02, -1.2426479763910001E+00, -1.2716350599970001E+00, -8.0541179788399997E-01, -6.9204995616999998E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.9521360472000000E-02, 9.0578638500000003E-02, 2.2497729539700001E-01, 3.7236989398100001E-01, 5.7979994051200001E-01, 7.8511152199000001E-01, 9.2471082584700004E-01, 9.8247515911899996E-01, 9.9320767635200002E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0023838489921419E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {1.2000000000000000E+01, 1.4125000000000000E+01, 1.5625000000000000E+01, 1.6125000000000000E+01, 1.6375000000000000E+01, 1.7125000000000000E+01, 1.7375000000000000E+01, 1.8625000000000000E+01, 1.8875000000000000E+01, 1.9000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 9.0084455370000004E-03, 2.3574224888000001E-02, 6.2824008634000006E-02, 1.3738927522600000E-01, 6.3120546290799995E-01, 6.4294223938799999E-01, 7.2739341543000005E-02, 1.5462587700999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 8.9319464809999999E-03, 3.9255406000000000E-05, 1.4457456209100000E-01, 7.2002980946200001E-01, 3.2560605691799999E-01, -3.3421115433900001E-01, -3.0466424195300001E-01, -1.7165607421400000E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 6.2298101269999996E-03, 3.2415934448999997E-02, 5.1062547391000003E-02, 7.3161027168000006E-02, 3.8083308849699998E-01, 5.4404758540099996E-01, 9.8888992758600003E-01, 9.9925477249899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0031313488064857E+00;
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
    lambda := 7.7000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT18HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT18HC;