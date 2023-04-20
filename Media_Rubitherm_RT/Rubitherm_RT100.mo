
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT100 "Rubitherm GmbH, RT100; data taken from: Rubitherm datasheet; last access: 2020-10-09."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT100";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.5814999999999998E+02, 3.9114999999999998E+02}
             "temperature range melting {startT, endT}";
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {0.0, 0.0}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8400000000000000E+05
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
    constant Integer  len_x =    19;
    constant Real[19] data_x =   {8.5000000000000000E+01, 8.7625000000000000E+01, 8.9625000000000000E+01, 9.1875000000000000E+01, 9.3625000000000000E+01, 9.6375000000000000E+01, 9.8625000000000000E+01, 1.0037500000000000E+02, 1.0162500000000000E+02, 1.0337500000000000E+02, 1.0562500000000000E+02, 1.0737500000000000E+02, 1.0837500000000000E+02, 1.1037500000000000E+02, 1.1237500000000000E+02, 1.1337500000000000E+02, 1.1462500000000000E+02, 1.1637500000000000E+02, 1.1800000000000000E+02};
    constant Real[19] data_y =   {0.0000000000000000E+00, 4.1010863819999997E-03, 1.5243502377000000E-02, 3.4102710692000000E-02, 4.0775616113999998E-02, 4.4624273577000002E-02, 5.9668420770000001E-02, 4.3119568039000002E-02, 3.9595388985000000E-02, 4.8736225679000003E-02, 4.3847610348999998E-02, 2.9618273199000001E-02, 4.1536972126000000E-02, 3.6145193713999997E-02, 1.7696180983999998E-02, 2.5461237753000000E-02, 2.1339574383999999E-02, 9.8423283390000008E-03, 0.0000000000000000E+00};
    constant Real[19] m_k =      {0.0000000000000000E+00, -2.1824628880000000E-03, -6.2527824890000003E-03, -8.6032981309999996E-03, 3.1933096018999997E-02, 9.1271002410000003E-03, 2.0958792048000000E-02, 1.4788986113999999E-02, -1.9724366597000002E-02, -1.7818459355000000E-02, 2.2277448950000001E-02, 2.9625926228000001E-02, -3.0431087648000000E-02, 7.6500558460000003E-03, 2.0362748632000001E-02, -2.1869024572999999E-02, 1.4288036811000000E-02, 4.0722243430000001E-03, 0.0000000000000000E+00};
    constant Real[19] iy_start = {0.0000000000000000E+00, 6.6298631590000004E-03, 2.7312433136999999E-02, 8.3767252328999994E-02, 1.3889047827000001E-01, 2.7056821408300002E-01, 3.8280402388399998E-01, 4.7423502666700001E-01, 5.3037479289900002E-01, 6.0710883252299996E-01, 6.9427099205300002E-01, 7.5662159637000004E-01, 7.9716713091299995E-01, 8.6209658798699995E-01, 9.1165536979799999E-01, 9.3673061048200001E-01, 9.6125088848200002E-01, 9.9111507952599998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9909223457417407E-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT100</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-04-20.<br><br>
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
    <li>file creation date: 2023-04-20 </ul>
    </p></html>"));
end Rubitherm_RT100;