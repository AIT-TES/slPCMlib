
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT11HC "Rubitherm GmbH, RT11HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT11HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7814999999999998E+02, 2.8714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7814999999999998E+02, 2.8614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9144732148397589E+05
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
    constant Real[11] data_x =   {5.0000000000000000E+00, 6.8750000000000000E+00, 8.8750000000000000E+00, 1.0625000000000000E+01, 1.1375000000000000E+01, 1.1625000000000000E+01, 1.1875000000000000E+01, 1.2125000000000000E+01, 1.2625000000000000E+01, 1.2875000000000000E+01, 1.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.0838443139999999E-02, 4.8991305652000001E-02, 1.9454043995799999E-01, 4.4219867647000000E-01, 5.6445984053800002E-01, 5.6032918990699998E-01, 4.3065862697000001E-01, 7.0534864885999998E-02, 2.3866669524999998E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.6416139550000000E-03, 3.2771516560999997E-02, 2.4374595278399999E-01, 4.0959886796500000E-01, 4.1164257539499999E-01, -3.4845193670399999E-01, -7.2950205114199995E-01, -4.3907314955600002E-01, -6.0818832964000001E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.9034455775999999E-02, 7.8423162653999998E-02, 2.3749848686399999E-01, 4.6825104641300003E-01, 5.9393640955799998E-01, 7.3833725962200003E-01, 8.6405902609300000E-01, 9.8317761151299998E-01, 9.9299707942600002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9891666378907540E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {5.0000000000000000E+00, 7.3750000000000000E+00, 9.8750000000000000E+00, 1.1625000000000000E+01, 1.1875000000000000E+01, 1.2625000000000000E+01, 1.2875000000000000E+01, 1.3000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 3.4542615842000003E-02, 1.2789325898500001E-01, 4.4635221412100001E-01, 4.2648599261699999E-01, 1.0556166228300000E-01, 2.3388876116999999E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, -5.7172034090000002E-03, 1.2491088941200000E-01, 1.9102187999500000E-01, -2.4222942985000001E-01, -4.7600481005500000E-01, -2.6616374317199998E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 4.3606995289999999E-02, 1.7830826506100000E-01, 6.6279279333300001E-01, 7.7389994563300002E-01, 9.8389570367900003E-01, 9.9888730773199996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9771787149267099E-01;
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
    lambda := 7.7000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT11HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-08-10.<br><br>
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
end Rubitherm_RT11HC;