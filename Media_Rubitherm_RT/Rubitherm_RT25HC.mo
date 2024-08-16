within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT25HC "Rubitherm GmbH, RT25HC; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT25HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 3.0314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8714999999999998E+02, 3.0114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9890373466098134E+05
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
    constant Real[11] data_x =   {1.5000000000000000E+01, 1.7125000000000000E+01, 1.9625000000000000E+01, 2.2375000000000000E+01, 2.4125000000000000E+01, 2.4375000000000000E+01, 2.4625000000000000E+01, 2.5375000000000000E+01, 2.6875000000000000E+01, 2.8125000000000000E+01, 3.0000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 7.3762466820000002E-03, 2.3025303577999999E-02, 1.1570887917000000E-01, 1.7701664369199999E-01, 2.1962604103700001E-01, 2.8245410119600001E-01, 2.8410260319399999E-01, 9.4592305960000001E-03, 5.3077823550000004E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.6296069560000001E-03, -1.6337033880000000E-03, 2.3525705230000001E-03, 1.3243584726499999E-01, 1.9895791228000001E-01, 2.1174051287800000E-01, -2.4735947069500000E-01, -2.2673504726999999E-02, 1.8936142390000000E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 7.2704234410000002E-03, 4.7226927491999998E-02, 2.3668298927500001E-01, 4.6105096966600001E-01, 5.1060096879000005E-01, 5.7369696761599998E-01, 8.0917843247499999E-01, 9.8836440884300003E-01, 9.9443367065699995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0064210399028344E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {1.4000000000000000E+01, 1.7375000000000000E+01, 1.8375000000000000E+01, 1.9375000000000000E+01, 2.3125000000000000E+01, 2.4125000000000000E+01, 2.4375000000000000E+01, 2.4625000000000000E+01, 2.5125000000000000E+01, 2.5625000000000000E+01, 2.5875000000000000E+01, 2.6125000000000000E+01, 2.7625000000000000E+01, 2.8000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 2.7253039481000000E-02, 4.5370042296000002E-02, 8.3265795242999999E-02, 1.3451505739000000E-01, 2.4762909264200000E-01, 3.0489258041299999E-01, 3.8213481498099999E-01, 2.4962282921500001E-01, 1.5715195288000000E-02, 2.4492208750000002E-03, 6.5569506700000003E-04, 5.1552422930000001E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -3.7419850340000001E-03, 6.3751235697999994E-02, -4.1415821062000001E-02, 7.4642544660999996E-02, 1.9069536981800000E-01, 2.3963380091399999E-01, 2.3501969396600000E-01, -5.1760739252499999E-01, -1.1990650344500001E-01, -8.3057949319999998E-03, -8.2267566399999998E-04, -1.5977553719999999E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 4.9371472473999999E-02, 7.9953280163999996E-02, 1.5278435110999999E-01, 4.2418308250100001E-01, 6.0496163878800002E-01, 6.7353584810500000E-01, 7.5914354312099996E-01, 9.3216693745699997E-01, 9.9001682183499995E-01, 9.9170032359500004E-01, 9.9204826548699998E-01, 9.9922330332099996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9656865207932821E-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT25HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_RT25HC;
