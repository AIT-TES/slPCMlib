
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP50 "Rubitherm GmbH, SP50; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP50";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1914999999999998E+02, 3.2714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1614999999999998E+02, 3.2214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0630096070419959E+05
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
    constant Real[11] data_x =   {4.6000000000000000E+01, 4.8625000000000000E+01, 4.9125000000000000E+01, 4.9375000000000000E+01, 4.9875000000000000E+01, 5.0125000000000000E+01, 5.0375000000000000E+01, 5.0625000000000000E+01, 5.1875000000000000E+01, 5.3125000000000000E+01, 5.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.5027983852000000E-02, 4.4648646593000003E-02, 1.1516903378399999E-01, 5.3591151281299998E-01, 6.6282285540399999E-01, 6.4551608430499996E-01, 4.9172471619500002E-01, 4.9571956958999998E-02, 1.5923075746000000E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -1.1167668184000000E-02, 1.1543138842800001E-01, 7.4522760416800005E-01, 8.3906462777799995E-01, 2.8225724133600000E-01, -5.4673540280599997E-01, -5.4344256329600005E-01, -8.4128095693000005E-02, -2.3492869507999999E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.6250285176000001E-02, 3.8585235646999998E-02, 5.5354682345000003E-02, 2.1686743565600000E-01, 3.7027180565000001E-01, 5.3885995559699995E-01, 6.8161444536500004E-01, 9.6132633766599995E-01, 9.9450883267099999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0043376164630935E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {4.3000000000000000E+01, 4.4375000000000000E+01, 4.6625000000000000E+01, 4.7375000000000000E+01, 4.7625000000000000E+01, 4.8125000000000000E+01, 4.8375000000000000E+01, 4.8625000000000000E+01, 4.8875000000000000E+01, 4.9000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 6.3040649943000002E-02, 1.3607155186100001E-01, 3.4592061283600001E-01, 5.1280282811099998E-01, 5.3337800206300001E-01, 3.7817055217500001E-01, 1.8722582807400001E-01, 4.2666029041000003E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.4170250509000001E-02, 1.0060215209100000E-01, 4.9372232040300001E-01, 5.4354913218400003E-01, -5.2967038606000005E-01, -7.1932651547199999E-01, -6.9395651882099996E-01, -4.9480958029000000E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 4.1225574584000002E-02, 2.2930021054999999E-01, 3.9208443924500003E-01, 4.9947189699400002E-01, 7.8418858094199995E-01, 8.9944895843799999E-01, 9.7019331347899995E-01, 9.9797186709700003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0028627162406392E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP50</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP50;