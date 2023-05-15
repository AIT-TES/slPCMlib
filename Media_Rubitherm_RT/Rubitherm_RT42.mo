
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT42 "Rubitherm GmbH, RT42; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT42";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0614999999999998E+02, 3.1714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0514999999999998E+02, 3.1714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4000000000000000E+05
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
    constant Real[8] data_x =   {3.3000000000000000E+01, 3.5875000000000000E+01, 3.8125000000000000E+01, 3.9375000000000000E+01, 4.1375000000000000E+01, 4.1875000000000000E+01, 4.2875000000000000E+01, 4.4000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 1.3773842771000000E-02, 5.3742534805999997E-02, 1.7643451924299999E-01, 2.6861712947400002E-01, 2.8444203993500000E-01, 6.2502119689000005E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 1.2438570200000001E-04, 6.1619356728000001E-02, 4.1027787147999997E-02, 8.2575865782999996E-02, -1.1535889525199999E-01, -1.3321028932699999E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 1.9995617245000000E-02, 7.0722219156999999E-02, 2.1935577008900001E-01, 6.5671292241599999E-01, 8.0113376904800004E-01, 9.7859079172499996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0142737256031147E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {3.2000000000000000E+01, 3.4625000000000000E+01, 3.6875000000000000E+01, 3.8625000000000000E+01, 4.0625000000000000E+01, 4.1375000000000000E+01, 4.1625000000000000E+01, 4.2375000000000000E+01, 4.2625000000000000E+01, 4.2875000000000000E+01, 4.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.3236993932000000E-02, 4.0696717815999997E-02, 1.0843419503000000E-01, 1.8644228178800001E-01, 2.8801462780300002E-01, 3.5949053726800001E-01, 1.4211818806500001E-01, 4.0393927045000000E-02, 1.3145174979000000E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -3.8926872010000002E-03, 1.4640720903000001E-02, 6.7933241174000000E-02, 3.9512906559999998E-02, 2.1147083819699999E-01, 2.1342258995499999E-01, -4.6325736544099999E-01, -2.5143437462700002E-01, -3.4293100730999997E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.9700114808000001E-02, 7.2802880220999996E-02, 1.9023601397199999E-01, 4.9600310811400000E-01, 6.6665485964700000E-01, 7.4795967228299998E-01, 9.6880589753400004E-01, 9.9061776102700005E-01, 9.9620510090500003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0046563917117808E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT42</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_RT42;