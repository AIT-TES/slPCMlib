
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT21 "Rubitherm GmbH, RT21; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8514999999999998E+02, 2.9914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8414999999999998E+02, 2.9614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4100000000000000E+05
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
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {1.2000000000000000E+01, 1.3875000000000000E+01, 1.5375000000000000E+01, 1.7125000000000000E+01, 1.8625000000000000E+01, 2.2125000000000000E+01, 2.3125000000000000E+01, 2.4125000000000000E+01, 2.5625000000000000E+01, 2.6000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.8953763633000000E-02, 1.9817944059000000E-02, 6.1076889402999997E-02, 7.4040760270999997E-02, 1.8773648569500001E-01, 1.7272634744800000E-01, 4.8579984349000002E-02, 5.6466763599999999E-04, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 2.3687065566000001E-02, 4.5165204929999998E-03, 1.6330698469999998E-02, 2.6636722487000001E-02, 2.7097923857000002E-02, -1.1300816999900000E-01, -8.7508927754999993E-02, -2.2704717069999998E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 2.0072472974000000E-02, 5.9983052581000000E-02, 1.2730784545300000E-01, 2.2606372435600000E-01, 6.8071075720800001E-01, 8.7136287182799999E-01, 9.7918147414699996E-01, 9.9992125021400002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9346136270420005E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {1.1000000000000000E+01, 1.4625000000000000E+01, 1.6375000000000000E+01, 1.8625000000000000E+01, 1.9875000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.1125000000000000E+01, 2.1625000000000000E+01, 2.1875000000000000E+01, 2.2875000000000000E+01, 2.3000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 5.9138713194999998E-02, 1.2708144798900001E-01, 8.2432925697999995E-02, 1.1476034521200000E-01, 2.1959281320700000E-01, 3.3485934290399999E-01, 2.8792239953499998E-01, 3.7779174009000001E-02, 1.1082280646000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.2635301734999998E-02, -1.4072832562999999E-02, 5.6278806830000003E-03, 9.9706906095000006E-02, 3.1689138533700001E-01, 3.5115067393299998E-01, -5.0181212959499999E-01, -2.8299886921099998E-01, -3.0314539641000001E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 8.2448928998999998E-02, 2.5485775068099997E-01, 4.8237945547000000E-01, 5.9343847796399996E-01, 6.7254703864999998E-01, 7.4171442374499996E-01, 9.1527853304300000E-01, 9.9218901566999995E-01, 9.9698335732400001E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0005685008552858E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT21</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2022-10-31.<br><br>
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
end Rubitherm_RT21;