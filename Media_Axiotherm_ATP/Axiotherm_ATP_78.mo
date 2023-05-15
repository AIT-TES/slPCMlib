
within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_78 "Axiotherm GmbH, ATP 78; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 78";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.4114999999999998E+02, 3.5314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.4114999999999998E+02, 3.5314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9900000000000000E+05
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
    constant Real[11] data_x =   {6.8000000000000000E+01, 7.0625000000000000E+01, 7.4625000000000000E+01, 7.6625000000000000E+01, 7.7375000000000000E+01, 7.7625000000000000E+01, 7.8125000000000000E+01, 7.8375000000000000E+01, 7.8625000000000000E+01, 7.9375000000000000E+01, 8.0000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.0028276298999999E-02, 3.1069674104000002E-02, 9.9928273757000002E-02, 2.6461654513799998E-01, 4.2454339611000003E-01, 5.0836344580100001E-01, 3.9752716675600003E-01, 2.3222294369700000E-01, 3.3548414393999998E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -2.9878918540000000E-03, -8.1782535919999993E-03, 3.6600567948999997E-02, 4.6817467482500003E-01, 5.3492323579099998E-01, -3.0594819897699999E-01, -5.7607104394700004E-01, -5.1218442379700002E-01, -1.0218996193100000E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.8376254976000000E-02, 1.3894778788900000E-01, 2.5656736577500000E-01, 3.7459490940700002E-01, 4.6153642601000000E-01, 7.1562516447600000E-01, 8.3179723516899995E-01, 9.1122859324799999E-01, 9.9274716751900005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0133357590081709E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {6.8000000000000000E+01, 7.0875000000000000E+01, 7.3875000000000000E+01, 7.5875000000000000E+01, 7.6875000000000000E+01, 7.7375000000000000E+01, 7.7625000000000000E+01, 7.7875000000000000E+01, 7.8125000000000000E+01, 7.8625000000000000E+01, 7.8875000000000000E+01, 7.9875000000000000E+01, 8.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.0270767508000000E-02, 2.4188899257000001E-02, 8.0689613431000007E-02, 1.5268781184800001E-01, 3.8706843085100001E-01, 6.7892004722499999E-01, 7.4054244574600003E-01, 5.3018112453099997E-01, 1.8245688239999999E-02, 2.8139837100000001E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 3.5363044700000002E-04, 9.2148205919999997E-03, 3.4641086090000002E-02, 1.9732027265600000E-01, 7.5117121054299996E-01, 8.5805448338900003E-01, -5.0170444391500002E-01, -1.1331827051990000E+00, -1.3599865722800000E-01, -8.3804121159999998E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.4804228415999999E-02, 6.0727516756000001E-02, 1.5901331534499999E-01, 2.6415955153900000E-01, 3.8997000563200002E-01, 5.2525329258700004E-01, 7.1337167231999998E-01, 8.7866737991800004E-01, 9.9727132624100001E-01, 9.9927753671899999E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0195295370364665E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8900000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 78</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATP_78;