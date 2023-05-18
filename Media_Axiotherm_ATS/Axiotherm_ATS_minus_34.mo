
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_minus_34 "Axiotherm GmbH, ATS -34; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -34";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.3714999999999998E+02, 2.4414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.2614999999999998E+02, 2.3714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9973632927323916E+05
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
    constant Real[11] data_x =   {-3.6000000000000000E+01, -3.4875000000000000E+01, -3.4625000000000000E+01, -3.4125000000000000E+01, -3.3875000000000000E+01, -3.3625000000000000E+01, -3.3375000000000000E+01, -3.2875000000000000E+01, -3.1875000000000000E+01, -3.0125000000000000E+01, -2.9000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 8.2853117230000002E-03, 4.2469681453000001E-02, 8.6992985140900003E-01, 1.0572265788410000E+00, 8.3355272984399997E-01, 3.6379740368899999E-01, 8.3304521487999994E-02, 2.0294394002000001E-02, 5.3521564900000001E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 2.0207535990000001E-02, 3.2101485539000002E-01, 1.6312471633929999E+00, -1.6374534826500001E-01, -1.4862657694950001E+00, -1.1317256879540001E+00, -1.6992360244599999E-01, -3.1370691425999998E-02, -1.1735925760000001E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.5316957270000002E-03, 7.3140336089999999E-03, 2.0831362643800000E-01, 4.5880162621999998E-01, 7.0227484583199995E-01, 8.5024149480900002E-01, 9.4206907484299995E-01, 9.8236179046399996E-01, 9.9711036845500001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0009771550314064E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {-4.7000000000000000E+01, -4.5875000000000000E+01, -4.4625000000000000E+01, -4.3375000000000000E+01, -4.1875000000000000E+01, -4.1625000000000000E+01, -4.1375000000000000E+01, -4.0625000000000000E+01, -4.0375000000000000E+01, -3.9125000000000000E+01, -3.7625000000000000E+01, -3.7375000000000000E+01, -3.6875000000000000E+01, -3.6375000000000000E+01, -3.6125000000000000E+01, -3.6000000000000000E+01};
    constant Real[16] data_y =   {0.0000000000000000E+00, 1.4079888911000000E-02, 6.2082289404999998E-02, 3.5948199966000002E-02, 6.7829546362000004E-02, 9.6491158517999998E-02, 1.4528730059799999E-01, 1.2122622196799999E-01, 6.7414357329999994E-02, 6.4065471306999994E-02, 2.3766268371700000E-01, 3.3764205754299997E-01, 3.3660277619099999E-01, 1.1602997611500000E-01, 2.6343316233000000E-02, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, 3.1786302333999997E-02, -4.4038725542000001E-02, 8.2072130220000004E-03, 7.9104809930000006E-02, 1.4224249827900001E-01, 1.5952510683400001E-01, -1.8708086170400001E-01, -1.5843649403400001E-01, 5.8303741494000000E-02, 2.9762974776800000E-01, 3.2274970624900001E-01, -3.5335422680200002E-01, -4.3825652493900003E-01, -3.0470874598999997E-01, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 4.5946071020000001E-03, 6.2410421555999997E-02, 1.1720015431200000E-01, 1.8212353734500000E-01, 2.0245483971000000E-01, 2.3276612169499999E-01, 3.4964602348800000E-01, 3.7321608776300003E-01, 4.2749008388300003E-01, 6.0999023775200001E-01, 6.8219888952700003E-01, 8.6593053228799999E-01, 9.8154019352599997E-01, 9.9874287559300001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0059400781944892E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.2550000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1800000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATS -34</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATS_minus_34;