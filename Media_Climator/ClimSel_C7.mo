
within slPCMlib.Media_Climator;
package ClimSel_C7 "Climator Sweden AB, ClimSel C7; data taken from: Climator Sweden AB datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C7";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7514999999999998E+02, 2.8414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7214999999999998E+02, 2.7814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.8000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 5.5346067535279079E+04
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
    constant Real[10] data_x =   {2.0000000000000000E+00, 3.3750000000000000E+00, 4.6250000000000000E+00, 6.1250000000000000E+00, 7.3750000000000000E+00, 7.6250000000000000E+00, 8.3750000000000000E+00, 8.6250000000000000E+00, 9.8750000000000000E+00, 1.1000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.9291497083000000E-02, 3.0636545449000002E-02, 3.3623136297999998E-02, 3.0045742751700000E-01, 4.4990979852500002E-01, 4.5283755173000001E-01, 3.0420987014500001E-01, 4.0391029342999998E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -7.1758359060000000E-03, 2.9928144181000000E-02, 3.5157629840000003E-02, 4.7586917115400001E-01, 5.1677701668700005E-01, -5.1051499285400004E-01, -4.6822444954600001E-01, -7.4739337244999995E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.4494157144000000E-02, 4.1052425787000003E-02, 8.8596930835999996E-02, 2.4107215115700001E-01, 3.3530962201900000E-01, 7.2469912146500004E-01, 8.1977020602999995E-01, 9.8505892088199998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0069951851390007E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-1.0000000000000000E+00, 8.7500000000000000E-01, 1.8750000000000000E+00, 2.3750000000000000E+00, 2.6250000000000000E+00, 3.1250000000000000E+00, 3.3750000000000000E+00, 3.6250000000000000E+00, 3.8750000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 4.5623041996000001E-02, 1.5778715978600000E-01, 4.1989257430900001E-01, 7.2055035272900003E-01, 6.9198718449800001E-01, 3.8014580010499999E-01, 9.7991858690999994E-02, 2.9615175123000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 4.4870453136000001E-02, 2.4951399974999999E-01, 8.2377903305699995E-01, 9.5134510462499999E-01, -1.1781740324709999E+00, -1.1950876699730000E+00, -8.0992552410700003E-01, -8.0138616383999994E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.9677333294000002E-02, 1.1447559258900000E-01, 2.4716135265399999E-01, 3.8929835357800002E-01, 7.8748699149500001E-01, 9.2182424660100004E-01, 9.7968556006600005E-01, 9.9185653398700002E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0017340193360045E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 7.8000000000000003E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C7</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C7;