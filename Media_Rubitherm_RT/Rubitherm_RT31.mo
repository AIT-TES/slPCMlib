
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT31 "Rubitherm GmbH, RT31; data taken from: Rubitherm datasheet; last access: 2022-10-31."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT31";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9714999999999998E+02, 3.1114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9714999999999998E+02, 3.0814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4920869236091958E+05
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
    constant Real[12] data_x =   {2.4000000000000000E+01, 2.7125000000000000E+01, 2.9125000000000000E+01, 3.0125000000000000E+01, 3.0875000000000000E+01, 3.2375000000000000E+01, 3.3375000000000000E+01, 3.4375000000000000E+01, 3.5125000000000000E+01, 3.6125000000000000E+01, 3.7375000000000000E+01, 3.8000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 2.1173509255999999E-02, 1.2636003421200001E-01, 2.1322143162099999E-01, 1.5128651538900001E-01, 1.5665114827900001E-01, 1.1305184091999999E-01, 6.5945845025000005E-02, 5.7568258198000002E-02, 9.6298375999999998E-03, 8.5832656219999998E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.0439770758999999E-02, 1.1243369516100001E-01, -2.4001399876000001E-02, -5.3288475932000000E-02, 9.9789949860000001E-03, -9.3425999960000006E-02, 1.1269572861999999E-02, -6.9594210808000004E-02, -4.5112625599999998E-03, -1.6792821502999999E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 2.4734076803999999E-02, 1.3894554672100001E-01, 3.2118435661299999E-01, 4.6006958852500002E-01, 6.8046447813599997E-01, 8.2478715443100004E-01, 9.0604223230100001E-01, 9.5644881604800003E-01, 9.8479201945100003E-01, 9.9785165696599998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0059532114838448E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.4000000000000000E+01, 2.5625000000000000E+01, 2.7875000000000000E+01, 2.9375000000000000E+01, 3.0125000000000000E+01, 3.1125000000000000E+01, 3.2125000000000000E+01, 3.2375000000000000E+01, 3.2625000000000000E+01, 3.3625000000000000E+01, 3.3875000000000000E+01, 3.5000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.3500481186000000E-02, 7.6197745383000004E-02, 1.4402439828700001E-01, 2.1589564545699999E-01, 1.5781820575700001E-01, 1.5237830795000001E-01, 1.9405647091400000E-01, 2.6686924419800001E-01, 2.5206172230999999E-02, 7.5586700300000001E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -5.0453110819999999E-03, 4.5286654932999999E-02, 1.1728695916400000E-01, 1.4148927091000000E-02, -6.1451496941000001E-02, 1.1019993732199999E-01, 1.9638997466899999E-01, 2.0631686116600001E-01, -2.0785767776500000E-01, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.2070467591000000E-02, 9.1688438804999994E-02, 2.4324318503199999E-01, 3.8294473650799998E-01, 5.7595930736000001E-01, 7.1664948851499999E-01, 7.5947333755599999E-01, 8.1699491477800001E-01, 9.9741406966199997E-01, 9.9957513823800004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9926282702912039E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.6000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT31</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2022-10-31.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
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
end Rubitherm_RT31;