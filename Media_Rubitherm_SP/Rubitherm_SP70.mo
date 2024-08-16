within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP70 "Rubitherm GmbH, SP70; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP70";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.3414999999999998E+02, 3.5114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.3414999999999998E+02, 3.5114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.2133679025980749E+05
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
    constant Real[12] data_x =   {6.1000000000000000E+01, 6.2875000000000000E+01, 6.4125000000000000E+01, 6.7375000000000000E+01, 6.9625000000000000E+01, 7.0375000000000000E+01, 7.0875000000000000E+01, 7.2125000000000000E+01, 7.3125000000000000E+01, 7.4625000000000000E+01, 7.7875000000000000E+01, 7.8000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.4210397277000001E-02, 3.3856949080000001E-03, 3.1521618461000003E-02, 9.6861341293000003E-02, 1.7036657252500001E-01, 2.5009887974000000E-01, 1.5608450657699999E-01, 5.6051122722000003E-02, 5.7847465833000003E-02, 1.0302161530000001E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 4.2389725660000000E-03, 1.2044268090000000E-03, 4.1064137189999999E-03, 4.1640973961999998E-02, 1.5667397631600000E-01, 9.8272293048999995E-02, -1.3318603195900000E-01, -4.4974457974000001E-02, 3.8175955136999999E-02, -1.1399981393999999E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.2088752853000000E-02, 2.3489347801000000E-02, 7.7697009422000002E-02, 2.0638227866400000E-01, 3.0126644206699998E-01, 4.0767337357400002E-01, 6.9187308156100003E-01, 7.9065850763900003E-01, 8.6054025786099997E-01, 9.9995042079899998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0006946729681980E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {6.1000000000000000E+01, 6.2625000000000000E+01, 6.4375000000000000E+01, 6.5125000000000000E+01, 6.6375000000000000E+01, 6.7125000000000000E+01, 6.7875000000000000E+01, 6.8375000000000000E+01, 6.8625000000000000E+01, 6.9125000000000000E+01, 7.0625000000000000E+01, 7.2625000000000000E+01, 7.4875000000000000E+01, 7.6375000000000000E+01, 7.8000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.0737405405999999E-02, 2.8832756793000000E-02, 2.7593158134000002E-02, 1.5570750239299999E-01, 1.8984404945800001E-01, 2.5661948547199998E-01, 1.6427074267399999E-01, 9.5493563042999996E-02, 4.1365556131000000E-02, 3.0409019648999999E-02, 6.9789856531999997E-02, 2.8099931985000000E-02, 4.3741342632000003E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, -5.2350384280000000E-03, -1.8083636611999999E-02, 3.2863338568000000E-02, 3.9757837652000000E-02, 8.6845423071999994E-02, -7.1068796143000004E-02, -2.5405396266500002E-01, -2.2971961269300001E-01, -4.2805117162000000E-02, -2.9871123150000001E-03, 3.6838418430999997E-02, -1.3664599409999999E-03, -3.8796500056000001E-02, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 9.8608395700000000E-03, 4.7705145471999999E-02, 6.6447674834000006E-02, 1.7993696799400000E-01, 3.0711445543799998E-01, 4.8166996435800002E-01, 5.9053597960299997E-01, 6.2282972394500002E-01, 6.5310353000800003E-01, 6.9939683537800001E-01, 7.8618601903300001E-01, 9.1223435033699996E-01, 9.7303919665100003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9845248635703854E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.5000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.3000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP70</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-11-24.<br><br>
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
end Rubitherm_SP70;
