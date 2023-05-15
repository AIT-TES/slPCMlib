
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_24 "Rubitherm GmbH, SP-24; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-24";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4914999999999998E+02, 2.5614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4514999999999998E+02, 2.5114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.1900000000000000E+05
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-2.4000000000000000E+01, -2.3875000000000000E+01, -2.3625000000000000E+01, -2.2875000000000000E+01, -2.2625000000000000E+01, -2.2375000000000000E+01, -2.1125000000000000E+01, -1.9375000000000000E+01, -1.7000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 3.8320120580999997E-02, 1.7045641657700000E-01, 6.0476212869599999E-01, 5.7951791796899998E-01, 4.5228919477200002E-01, 7.7751272585000003E-02, 2.3982612234999999E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 4.4025541805599999E-01, 6.8272424945300003E-01, 1.4991556741600001E-01, -4.5072006311899998E-01, -4.4534963505500003E-01, -1.0387498446300000E-01, -1.1998133076999999E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 1.8260110270000001E-03, 2.6718192929999999E-02, 3.4313748682500000E-01, 4.9465368057600001E-01, 6.2390261615599996E-01, 9.1138460131800003E-01, 9.7710707839099997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0023344112200103E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-2.8000000000000000E+01, -2.6875000000000000E+01, -2.5375000000000000E+01, -2.4875000000000000E+01, -2.4625000000000000E+01, -2.4375000000000000E+01, -2.4125000000000000E+01, -2.3875000000000000E+01, -2.3625000000000000E+01, -2.3375000000000000E+01, -2.3125000000000000E+01, -2.2000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.7960990159999998E-03, 5.5581516139000003E-02, 1.7891456693600000E-01, 3.8933305339899998E-01, 9.0097960819199996E-01, 1.1050094758480000E+00, 8.4553453448500004E-01, 2.7313695982800001E-01, 1.2018623302000001E-02, 1.4846831999999999E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.1386795361000001E-02, 3.8952741456999999E-02, 4.2659043008399999E-01, 1.1840786012510001E+00, 1.5619057081240000E+00, -1.0296371477900000E-01, -1.7867078800050000E+00, -1.5981609907509999E+00, -7.4945367977999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.4979487010000000E-03, 4.1641865642000002E-02, 9.2227047830000006E-02, 1.5936178360600001E-01, 3.1879946749600002E-01, 5.7840891245100001E-01, 8.3118103527400000E-01, 9.7013444633900003E-01, 9.9786573401199996E-01, 9.9916425528599995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0007309176771411E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.2000000000000000E+03;
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
    lambda := 1.3000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP-24</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP_minus_24;