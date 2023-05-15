
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_7 "Rubitherm GmbH, SP-7; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-7";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6114999999999998E+02, 2.7214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6014999999999998E+02, 2.6914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2354016561450472E+05
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-1.2000000000000000E+01, -9.1250000000000000E+00, -7.3750000000000000E+00, -6.6250000000000000E+00, -6.3750000000000000E+00, -5.6250000000000000E+00, -5.3750000000000000E+00, -2.8750000000000000E+00, -1.0000000000000000E+00};
    constant Real[9] data_y =   {0.0000000000000000E+00, 2.4662135108000000E-02, 7.5509235928000004E-02, 2.2699312896300000E-01, 3.6165700752800001E-01, 4.3764466122599999E-01, 3.2607527111800000E-01, 2.3560382056000001E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 1.0590003565999999E-02, 5.7971603888000003E-02, 4.2120681861199999E-01, 4.7186526989700001E-01, -3.8098612258999998E-01, -3.6289928247499997E-01, -2.1972279498000000E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 2.8625145003999999E-02, 1.0543805836600000E-01, 2.0345135727100000E-01, 2.7798670517599999E-01, 6.2334549226199998E-01, 7.2030055237699997E-01, 9.8408934890999999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0166116992090262E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-1.3000000000000000E+01, -9.8750000000000000E+00, -8.3750000000000000E+00, -6.8750000000000000E+00, -6.6250000000000000E+00, -6.3750000000000000E+00, -5.8750000000000000E+00, -5.3750000000000000E+00, -5.1250000000000000E+00, -4.1250000000000000E+00, -4.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 6.0012801275000000E-02, 9.8481411180000006E-02, 2.8467840139500000E-01, 3.7839428250200002E-01, 5.1347158891400002E-01, 3.9668988338000000E-01, 5.0783536407999999E-02, 1.4873526535999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 2.4711004761000001E-02, 2.3258396837000001E-02, 2.9809290916800002E-01, 4.0524557193600003E-01, 4.1909051811699999E-01, -7.0682519234899999E-01, -3.4844273503000001E-01, -4.0697255111000002E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 7.4522563117999999E-02, 1.9506053791000000E-01, 4.3366017151300001E-01, 5.1695006024599999E-01, 6.2966560985599995E-01, 8.8360128002399996E-01, 9.8922569488300005E-01, 9.9590731114200004E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0117082093924368E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.2500000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1500000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.1500000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP-7</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP_minus_7;