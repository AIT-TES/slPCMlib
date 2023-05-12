
// within slPCMlib.Axiotherm_ATS;
package Axiotherm_ATS_minus_24 "Axiotherm GmbH, ATS -24; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -24";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4714999999999998E+02, 2.5814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4414999999999998E+02, 2.5614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.0000000000000000E+05
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
    constant Real[10] data_x =   {-2.6000000000000000E+01, -2.4875000000000000E+01, -2.4625000000000000E+01, -2.4125000000000000E+01, -2.3625000000000000E+01, -2.3375000000000000E+01, -2.2875000000000000E+01, -2.0875000000000000E+01, -1.7125000000000000E+01, -1.5000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 0.0000000000000000E+00, 1.3333333333000000E-02, 4.7517167569300001E-01, 7.0715095976499998E-01, 4.8875055820899999E-01, 2.3079386379100000E-01, 4.5480319941999997E-02, 7.0500681810000002E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 0.0000000000000000E+00, 1.1625806790500000E-01, 1.0702004605030000E+00, -6.6339222087700001E-01, -6.5357262132700000E-01, -3.3663152611699998E-01, -1.9211329425000000E-02, -2.3054942999999999E-04, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 0.0000000000000000E+00, 1.0487250830000000E-03, 1.0210335051500000E-01, 4.2991488745099998E-01, 5.7760087273299998E-01, 7.4885412898899995E-01, 9.1732465515999995E-01, 9.9268279157899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.8828559168279906E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-2.9000000000000000E+01, -2.6875000000000000E+01, -2.6125000000000000E+01, -2.5625000000000000E+01, -2.5375000000000000E+01, -2.4875000000000000E+01, -2.4625000000000000E+01, -2.4375000000000000E+01, -2.4125000000000000E+01, -2.3125000000000000E+01, -1.9625000000000000E+01, -1.7000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.8751945587000002E-02, 9.1305037330000005E-02, 3.1025072125999997E-01, 6.3937676431699997E-01, 7.8070544004300002E-01, 5.0476879139200004E-01, 1.8149288270899999E-01, 6.8103221884999995E-02, 1.1229643213000001E-02, 3.2629086430000002E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 3.9093026209000001E-02, 1.4616083618500000E-01, 8.7970581803799996E-01, 1.1005281972509999E+00, -9.6624399382600001E-01, -1.2004979352080001E+00, -9.6524757086599999E-01, -1.9115923660600001E-01, -9.2957339629999997E-03, -2.3376682999999999E-04, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 3.6950353041999999E-02, 8.4276421113999994E-02, 1.6906695744200001E-01, 2.8618350437899998E-01, 6.8278273987799998E-01, 8.4408553188400004E-01, 9.2832880359299996E-01, 9.5539566396300002E-01, 9.7981573288900004E-01, 9.9586707894000004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9628447580912294E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.1910000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1200000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.1200000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS -24</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
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
end Axiotherm_ATS_minus_24;