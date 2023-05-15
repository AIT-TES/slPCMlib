
within slPCMlib.Media_Rubitherm_PX;
package Rubitherm_PX52 "Rubitherm GmbH, PX52; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "PX52";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1714999999999998E+02, 3.2714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1614999999999998E+02, 3.2714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 7.2492256650577532E+04
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
    constant Real[8] data_x =   {4.4000000000000000E+01, 4.5375000000000000E+01, 4.7875000000000000E+01, 4.9375000000000000E+01, 5.0375000000000000E+01, 5.1625000000000000E+01, 5.2875000000000000E+01, 5.4000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 8.9535834866999994E-02, 1.2950754576000001E-01, 1.3959283712100001E-01, 1.6178825836100000E-01, 1.5085085024799999E-01, 3.0001895824000001E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 1.3905066291000001E-02, 9.3064659910000007E-03, 4.3761373903999999E-02, -2.5681899978000002E-02, 4.1620990710000000E-03, -6.3412880465000002E-02, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 5.8722675901999997E-02, 3.3193302739000002E-01, 5.2518381511699996E-01, 6.7996793160900004E-01, 8.6940891299500000E-01, 9.8992226339000000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.8917818295348392E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {4.3000000000000000E+01, 4.5125000000000000E+01, 4.6375000000000000E+01, 4.7625000000000000E+01, 4.9625000000000000E+01, 5.0625000000000000E+01, 5.1375000000000000E+01, 5.1875000000000000E+01, 5.2875000000000000E+01, 5.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 4.4009664763999999E-02, 6.4522304045999998E-02, 1.0448606508800000E-01, 1.6277355436099999E-01, 1.4840883606399999E-01, 1.6918436919999999E-01, 2.0382849593399999E-01, 5.7743530483000002E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.8428983413000002E-02, 5.1561076132999997E-02, 3.5881595310000001E-03, 4.2449256437000003E-02, -3.5331046039000003E-02, 8.2442162138000000E-02, -3.3517615737000003E-02, -1.1276423672300000E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.9872344450000000E-02, 1.0346561679600000E-01, 2.1547418528100001E-01, 4.7007984521099999E-01, 6.3234376063800002E-01, 7.4605446160599997E-01, 8.4183626758399999E-01, 9.7938810143800004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0011786569875543E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 6.5000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 6.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 6.5000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>PX52</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_PX52;