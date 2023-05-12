
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT44HC_DSC_1_Kmin "Rubitherm GmbH, RT44HC; data taken from: AIT Austrian Institute of Technology GmbH; last access: 2020-05-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT44HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0314999999999998E+02, 3.1889999999999998E+02}
             "temperature range melting {startT, endT}";
    // -> These are just dummy variables - there is no data for cooling! <-
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {0.0, 0.0}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4078322410616104E+05
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
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {3.0000000000000000E+01, 3.6625000000000000E+01, 4.0125000000000000E+01, 4.1375000000000000E+01, 4.2125000000000000E+01, 4.3125000000000000E+01, 4.4125000000000000E+01, 4.4375000000000000E+01, 4.4625000000000000E+01, 4.4875000000000000E+01, 4.5125000000000000E+01, 4.5375000000000000E+01, 4.5750000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 3.7412744840000000E-03, 2.6704161782999999E-02, 1.4164157015500001E-01, 2.5892743334500001E-01, 2.0980775312700001E-01, 3.1915305086599999E-01, 3.1918481489599998E-01, 2.6933308689500002E-01, 1.1467994493399999E-01, 1.1747829422000000E-02, 9.1382293199999997E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.5072701329999999E-03, 2.1880604859000001E-02, 1.6348829008800000E-01, 1.2750756524300000E-01, 5.7164106410999997E-02, 1.0582622869800000E-01, -4.2626705314999999E-02, -4.0521903832599998E-01, -5.0849677637500001E-01, -1.0308825865900000E-01, -3.5578936280000000E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 6.8928855460000001E-03, 3.9435205484000001E-02, 1.2637464678700000E-01, 2.7855794686699997E-01, 5.1923557663200004E-01, 7.8014655283400003E-01, 8.6086225088699997E-01, 9.3645622744599999E-01, 9.8508629835499995E-01, 9.9880380602600005E-01, 9.9987011044399998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0018652227827409E+00;
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
    rho := 8.0000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT44HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: AIT Austrian Institute of Technology GmbH - last access 2020-05-14.<br><br>
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
end Rubitherm_RT44HC_DSC_1_Kmin;