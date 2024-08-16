within slPCMlib.Media_Rubitherm_PX;
package Rubitherm_PX15 "Rubitherm GmbH, PX15; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "PX15";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7914999999999998E+02, 2.9114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7914999999999998E+02, 2.9114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 5.6000000000000000E+04
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
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {6.0000000000000000E+00, 7.3750000000000000E+00, 8.6250000000000000E+00, 1.2375000000000000E+01, 1.5625000000000000E+01, 1.7125000000000000E+01, 1.8000000000000000E+01};
    constant Real[7] data_y =   {0.0000000000000000E+00, 4.5428283960000003E-02, 7.5340622965999998E-02, 1.0721086282400000E-01, 1.1024670071800000E-01, 5.7711250582000002E-02, 0.0000000000000000E+00};
    constant Real[7] m_k =      {0.0000000000000000E+00, 1.1376448502000000E-02, 4.4498467670999999E-02, -1.8426750900000000E-04, 5.3784761569999999E-03, -8.8398241589000004E-02, 0.0000000000000000E+00};
    constant Real[7] iy_start = {0.0000000000000000E+00, 2.9238215203999999E-02, 9.9919279080000006E-02, 4.9186677709700000E-01, 8.3795562864999995E-01, 9.8052543006600001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9316067079406156E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {6.0000000000000000E+00, 7.3750000000000000E+00, 8.6250000000000000E+00, 1.0625000000000000E+01, 1.2375000000000000E+01, 1.4625000000000000E+01, 1.5625000000000000E+01, 1.6875000000000000E+01, 1.8000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 3.4537498910999997E-02, 4.9809169341999998E-02, 9.2638843200000001E-02, 1.0289062774600000E-01, 1.2036665010300000E-01, 1.7819282575900000E-01, 6.8398373232000001E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, -1.0133116021000000E-02, 4.5268810296000002E-02, 2.2557917545999999E-02, 2.3920357695000002E-02, 1.1015305745000000E-02, 7.5308949141999995E-02, -1.1642270311900001E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 2.5406747047000001E-02, 7.1027635825000004E-02, 2.2143502470900001E-01, 3.9261842854200002E-01, 6.4989271045200003E-01, 7.9418791231999997E-01, 9.7373693386700000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0025935418108689E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 6.5000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 6.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 1.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 1.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>PX15</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: other<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-03.<br><br>
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
end Rubitherm_PX15;
