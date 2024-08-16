within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP24E "Rubitherm GmbH, SP24E; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP24E";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9214999999999998E+02, 2.9914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9114999999999998E+02, 2.9714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8274243091394374E+05
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
    constant Real[13] data_x =   {1.9000000000000000E+01, 2.0125000000000000E+01, 2.2125000000000000E+01, 2.2875000000000000E+01, 2.3125000000000000E+01, 2.3375000000000000E+01, 2.3625000000000000E+01, 2.4125000000000000E+01, 2.4375000000000000E+01, 2.4625000000000000E+01, 2.4875000000000000E+01, 2.5625000000000000E+01, 2.6000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 6.0082734359999997E-03, 3.8600687123999999E-02, 3.1911041727000003E-02, 7.5428390090999997E-02, 2.2301039922900001E-01, 7.2464012696600000E-01, 9.8204156770700002E-01, 4.7484971805999998E-01, 7.1183459244000005E-02, 1.6885350561999999E-02, 6.0427756689999997E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.2699364326000001E-02, -2.4646375775000000E-02, 5.8592111569000002E-02, 2.1054837993400000E-01, 1.1572879317639999E+00, 1.7264359873420001E+00, -1.7061645159890000E+00, -1.8029921538220000E+00, -6.4507270148499996E-01, -3.8089124446000003E-02, -5.3468662159999996E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.0379594820000000E-03, 5.9030948361000000E-02, 8.1545541981999994E-02, 9.4157247722000001E-02, 1.2649453518300000E-01, 2.4185587696600000E-01, 7.3947519477099999E-01, 9.2188431487800004E-01, 9.8403723557400002E-01, 9.9187558608699999E-01, 9.9893084909800001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9886865584828699E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {1.8000000000000000E+01, 2.0625000000000000E+01, 2.1875000000000000E+01, 2.2375000000000000E+01, 2.2625000000000000E+01, 2.2875000000000000E+01, 2.3375000000000000E+01, 2.3625000000000000E+01, 2.3875000000000000E+01, 2.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 4.3967392538999998E-02, 1.8123286369600000E-01, 3.9141432767599998E-01, 5.9445935194800004E-01, 6.8888181314600005E-01, 4.5159121586899997E-01, 2.2443809203000001E-01, 5.1250925760000003E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 4.5132548063000003E-02, 2.3976816822700001E-01, 6.0536378637099997E-01, 6.6846381659300003E-01, 9.9024143415000002E-02, -8.5370190987500005E-01, -8.2386513460900002E-01, -5.9523082498699997E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.1981668206999998E-02, 1.4807990473499999E-01, 2.8443700582499998E-01, 4.0807874033100000E-01, 5.7244083971699999E-01, 8.7923422321900002E-01, 9.6408770974900004E-01, 9.9755731336900000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0059897251221253E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.6000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.0000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.0000000000000000E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP24E</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP24E;
