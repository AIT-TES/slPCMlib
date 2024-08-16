within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT54HC "Rubitherm GmbH, RT54HC; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT54HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2314999999999998E+02, 3.3214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2314999999999998E+02, 3.2814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.6951398713847113E+05
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
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {5.0000000000000000E+01, 5.1125000000000000E+01, 5.2625000000000000E+01, 5.3125000000000000E+01, 5.3375000000000000E+01, 5.3625000000000000E+01, 5.4125000000000000E+01, 5.4375000000000000E+01, 5.4625000000000000E+01, 5.4875000000000000E+01, 5.5625000000000000E+01, 5.6375000000000000E+01, 5.9000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 6.4641632810000004E-03, 4.1973362939000003E-02, 1.2811356542499999E-01, 2.8646789765699998E-01, 6.7796281884400000E-01, 8.3363052773699997E-01, 4.8459449284400002E-01, 1.2499102251800000E-01, 3.5571035468999999E-02, 2.5370162391000001E-02, 3.7304328837000002E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.3895284603000000E-02, 1.6010581002000000E-02, 3.0861911824200000E-01, 9.9441855206499996E-01, 1.3145552978360000E+00, -1.3629863865920000E+00, -1.4170314055050000E+00, -1.2929825325219999E+00, -1.0414240869999999E-01, 4.7486311267000002E-02, -3.5517212617999998E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.1672316719999998E-03, 3.8043437069999998E-02, 7.4413075901999998E-02, 1.2258959842599999E-01, 2.4129301815599999E-01, 6.7430576656600005E-01, 8.3911124817700000E-01, 9.1454702757100004E-01, 9.2840404162699997E-01, 9.4412515258600005E-01, 9.7147670082799997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9846037179329938E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {5.0000000000000000E+01, 5.1625000000000000E+01, 5.3125000000000000E+01, 5.3375000000000000E+01, 5.3625000000000000E+01, 5.3875000000000000E+01, 5.4375000000000000E+01, 5.4625000000000000E+01, 5.4875000000000000E+01, 5.5000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.5269466506000000E-02, 2.6694476626000002E-01, 4.1060889903600001E-01, 6.2737196843099996E-01, 7.2958760748499996E-01, 4.8007305769699998E-01, 2.3881814854500000E-01, 5.4562100188999998E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.0945965714000000E-02, 4.3999424072900001E-01, 6.6161356540500005E-01, 7.2325525663099999E-01, 1.1459739753600000E-01, -9.0620932108899999E-01, -8.7463718118900002E-01, -6.3390955746699995E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.8232595306000000E-02, 1.5778742088600001E-01, 2.4183367097899999E-01, 3.7204464250500002E-01, 5.4588190847300000E-01, 8.7152561580800003E-01, 9.6176620337399998E-01, 9.9739960658500004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0060606847631548E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.5000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.0000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT54HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_RT54HC;
