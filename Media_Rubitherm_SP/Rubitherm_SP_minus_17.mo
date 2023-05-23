
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_17 "Rubitherm GmbH, SP-17; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-17";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5114999999999998E+02, 2.6414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4914999999999998E+02, 2.5714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.0100000000000000E+05
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
    constant Real[15] data_x =   {-2.2000000000000000E+01, -1.9375000000000000E+01, -1.8875000000000000E+01, -1.8625000000000000E+01, -1.7625000000000000E+01, -1.7375000000000000E+01, -1.6625000000000000E+01, -1.6375000000000000E+01, -1.5625000000000000E+01, -1.5125000000000000E+01, -1.4375000000000000E+01, -1.3125000000000000E+01, -1.1875000000000000E+01, -1.0125000000000000E+01, -9.0000000000000000E+00};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.2349150935999999E-02, 3.3565705602000000E-02, 6.8582267303000002E-02, 2.7217672959799999E-01, 2.2611690645599999E-01, 2.1867981772600001E-01, 2.5822386710599998E-01, 2.2769779184700001E-01, 1.4039406211700001E-01, 1.0887462667300001E-01, 9.6222337319999997E-03, 7.0306564030000003E-03, 1.3955766006000001E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 5.2384924899999998E-03, 7.4181418932000004E-02, 2.9570947046700002E-01, -1.2083031293400000E-01, -1.1697911291099999E-01, 1.0853406974599999E-01, 1.1491505966899999E-01, -1.8192736155299999E-01, -1.0107703892100001E-01, -2.9415442114000000E-02, -1.9291977588000001E-02, -4.3174177300000003E-04, 4.1071574399999998E-04, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 1.3182433859000000E-02, 2.3211305769999999E-02, 3.4810360654999997E-02, 2.3962516542900000E-01, 3.0180791277000002E-01, 4.5782525061800000E-01, 5.1732469845600004E-01, 7.1319553423700000E-01, 8.0341239270900000E-01, 8.9340758942900000E-01, 9.6605195377800002E-01, 9.7399353522099996E-01, 9.9211719965699996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9865258790303291E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {-2.4000000000000000E+01, -2.3125000000000000E+01, -2.2375000000000000E+01, -2.1625000000000000E+01, -1.9875000000000000E+01, -1.7375000000000000E+01, -1.6375000000000000E+01, -1.6000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 4.4142952429999998E-03, 5.3904070604000001E-02, 1.2862126724400000E-01, 1.1367538030000000E-01, 2.6520898513300001E-01, 4.6697196323999997E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 1.2922185685000000E-02, 1.2907039276700000E-01, -4.2432344369999998E-03, 3.7877467957999997E-02, -8.5630043291999997E-02, -2.0054470698800000E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 1.1030909179999999E-03, 1.7473109306000000E-02, 9.1919419055999999E-02, 2.9250645959400001E-01, 8.2864001098700002E-01, 9.9361582794400005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9665620706756752E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.2000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP-17</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP_minus_17;