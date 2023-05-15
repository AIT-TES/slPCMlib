
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT111HC "Rubitherm GmbH, RT111HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT111HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.7514999999999998E+02, 3.8714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.7514999999999998E+02, 3.8514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8962207544645295E+05
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
    constant Real[10] data_x =   {1.0200000000000000E+02, 1.0362500000000000E+02, 1.0687500000000000E+02, 1.0862500000000000E+02, 1.1062500000000000E+02, 1.1162500000000000E+02, 1.1187500000000000E+02, 1.1262500000000000E+02, 1.1287500000000000E+02, 1.1400000000000000E+02};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.5486616892000000E-02, 2.2853179937000000E-02, 5.2040243283999997E-02, 2.4353700372100001E-01, 3.8089461907200001E-01, 3.5575761908900000E-01, 9.9711084221000001E-02, 5.0169746472000003E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -3.8563098069999999E-03, 8.5882600130000007E-03, 4.6750075940000003E-03, 1.7241616254500000E-01, 6.8862990055999995E-02, -2.1515404326699999E-01, -3.7189056573700002E-01, -1.2453605831400000E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.3484591494000000E-02, 6.5036048442999997E-02, 1.3182964014000001E-01, 3.7244112404000002E-01, 6.9455542448300001E-01, 7.8848627609599997E-01, 9.6733870272300004E-01, 9.8485451018100001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0039553586512093E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.0200000000000000E+02, 1.0362500000000000E+02, 1.0612500000000000E+02, 1.0787500000000000E+02, 1.0962500000000000E+02, 1.1012500000000000E+02, 1.1037500000000000E+02, 1.1062500000000000E+02, 1.1112500000000000E+02, 1.1137500000000000E+02, 1.1162500000000000E+02, 1.1187500000000000E+02, 1.1200000000000000E+02};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.6469813024000000E-02, 2.6411945770000000E-02, 5.7266371280999998E-02, 1.0049209305800000E-01, 1.8023081468400001E-01, 2.8266872065499998E-01, 4.7037293666899999E-01, 5.4570642829299998E-01, 3.9784168412699999E-01, 2.0044270681099999E-01, 4.6103175160000003E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -4.5980728899999998E-03, 1.3453883787000000E-02, 1.4445333208000000E-02, 2.5440870127000002E-02, 2.6532210522600003E-01, 5.3446603816200000E-01, 6.1986497217799996E-01, -4.6187620622199999E-01, -7.3700962562200001E-01, -7.1319343748599995E-01, -5.3824218081499997E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.4615656834000000E-02, 5.9497882315000000E-02, 1.3358937272300000E-01, 2.7090875754100002E-01, 3.3709785297599998E-01, 3.9442979313700000E-01, 4.8956094902199998E-01, 7.7038482393300001E-01, 8.9160351902299995E-01, 9.6741718679800004E-01, 9.9778573694899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0154318016406796E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.0000000000000000E+02;
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
    lambda := 8.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT111HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2022-12-06.<br><br>
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
end Rubitherm_RT111HC;