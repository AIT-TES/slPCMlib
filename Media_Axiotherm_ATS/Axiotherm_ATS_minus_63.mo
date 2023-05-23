
// within slPCMlib.Axiotherm_ATS;
package Axiotherm_ATS_minus_63 "Axiotherm GmbH, ATS -63; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -63";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.0314999999999998E+02, 2.1914999999999998E+02}
             "temperature range melting {startT, endT}";
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {0.0, 0.0}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4029335036681232E+05
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
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {-7.0000000000000000E+01, -6.7875000000000000E+01, -6.4375000000000000E+01, -6.3375000000000000E+01, -6.2875000000000000E+01, -6.2625000000000000E+01, -6.2375000000000000E+01, -6.1875000000000000E+01, -6.1625000000000000E+01, -6.1375000000000000E+01, -6.1125000000000000E+01, -6.0125000000000000E+01, -5.7625000000000000E+01, -5.5375000000000000E+01, -5.4000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.4790080099999999E-03, 7.5912741675999998E-02, 7.5291418417000003E-02, 1.3830948588000000E-01, 2.4568075978100001E-01, 4.7651566964900000E-01, 5.8282582585499998E-01, 3.9751527775599999E-01, 1.6533560926700000E-01, 7.1119712668999996E-02, 1.0943364038000000E-02, 1.2910015980000000E-02, 2.7945798050000002E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 1.6738168410000000E-03, 3.8856886796000002E-02, -1.7071203789000001E-02, 2.3894243133699999E-01, 6.2430698554200004E-01, 7.6926813018700002E-01, -6.1001425233700002E-01, -8.3975498249799996E-01, -6.7766955138899998E-01, -1.8446541543900000E-01, -5.3008135649999996E-03, -1.3265584920000000E-02, 2.2403599646000000E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 9.5129466099999995E-04, 9.9434308865000001E-02, 1.8052471450400001E-01, 2.2908697518300000E-01, 2.7555290457900000E-01, 3.6599555763600000E-01, 6.6259321298999996E-01, 7.8760838794599997E-01, 8.5783734761200003E-01, 8.8510378694900005E-01, 9.1147409102200005E-01, 9.4578937430800003E-01, 9.7702301483300003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0103117739401264E+00;
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
    rho := 1.4360000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3500000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.3500000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS -63</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATS_minus_63;