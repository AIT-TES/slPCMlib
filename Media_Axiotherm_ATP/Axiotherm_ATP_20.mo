within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_20 "Axiotherm GmbH, ATP 20; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 20";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8714999999999998E+02, 2.9514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8614999999999998E+02, 2.9414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.5000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.5000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8700000000000000E+05
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
    constant Real[11] data_x =   {1.4000000000000000E+01, 1.6625000000000000E+01, 1.8625000000000000E+01, 1.9125000000000000E+01, 1.9375000000000000E+01, 1.9875000000000000E+01, 2.0125000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.1125000000000000E+01, 2.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 4.0205493090000001E-03, 1.4488358793999999E-02, 4.8177257100000002E-02, 1.2964258981700000E-01, 6.2386168924499996E-01, 7.5576814482300003E-01, 7.0571177780399996E-01, 4.9819540647999999E-01, 2.0409184399499999E-01, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -2.1735052669999999E-03, -7.4844616290000001E-03, 1.2969981756700000E-01, 9.3109824402899999E-01, 9.8615000280099996E-01, 2.1673546530999999E-01, -7.1219271548899998E-01, -6.9447173557999997E-01, -4.4068379717599998E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 6.5860655250000004E-03, 2.7054959034999999E-02, 3.9983150988000003E-02, 5.8205531946999997E-02, 2.4718579832500001E-01, 4.2529729219000001E-01, 6.1457430808799995E-01, 7.6637702708800004E-01, 9.3825422062300001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0093527539797507E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.3000000000000000E+01, 1.5375000000000000E+01, 1.5625000000000000E+01, 1.6375000000000000E+01, 1.7375000000000000E+01, 1.8875000000000000E+01, 1.9375000000000000E+01, 1.9625000000000000E+01, 2.0125000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.0875000000000000E+01, 2.1000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 9.5084260400000000E-04, 0.0000000000000000E+00, 0.0000000000000000E+00, 3.6740259466000000E-02, 1.0965093551200000E-01, 3.2595999195800002E-01, 5.9505175739799998E-01, 7.4909814813900000E-01, 5.5722987870399998E-01, 2.8428671530700000E-01, 6.5820504863999998E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -1.0136474078000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 8.7619399420000004E-03, 1.8611424661100001E-01, 7.7270334857699996E-01, 9.1315392983699994E-01, -5.5147269864500004E-01, -1.0153667195439999E+00, -9.8655024278299996E-01, -7.7220947640100002E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 5.9277670329999999E-03, 5.9942089870000000E-03, 5.9942089870000000E-03, 2.3735850366999999E-02, 1.0071684968600000E-01, 1.9795623203099999E-01, 3.1301053729299999E-01, 6.8167379365199998E-01, 8.4833601495800004E-01, 9.5398093035800002E-01, 9.9687378390199999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0057638034451541E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.9400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.4000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATP 20</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_20;
