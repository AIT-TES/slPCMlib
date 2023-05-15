
within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_2 "Axiotherm GmbH, ATP 2; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 2";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6714999999999998E+02, 2.7814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6714999999999998E+02, 2.7614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8246363715007584E+05
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
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {-6.0000000000000000E+00, -4.3750000000000000E+00, -2.3750000000000000E+00, -1.2500000000000000E-01, 1.1250000000000000E+00, 1.3750000000000000E+00, 1.6250000000000000E+00, 2.1250000000000000E+00, 2.3750000000000000E+00, 2.6250000000000000E+00, 2.8750000000000000E+00, 3.3750000000000000E+00, 3.6250000000000000E+00, 4.1250000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[16] data_y =   {0.0000000000000000E+00, 6.6357869970000002E-03, 1.6133053583999999E-02, 3.8966876333999999E-02, 1.4630574507100000E-01, 2.3911093530800001E-01, 4.0826057722300002E-01, 4.7319383837700002E-01, 3.3678234536899998E-01, 1.6011491877299999E-01, 9.4084649174000001E-02, 9.9650830620000005E-02, 1.7630796397699999E-01, 2.1602703242499999E-01, 1.8739038590999998E-02, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, -8.8794314210000006E-03, 6.5139532540000001E-03, 3.9237880193000001E-02, 2.4200750226600001E-01, 4.8125393394799998E-01, 5.6010679205500002E-01, -4.3861463749399998E-01, -6.2513682865499998E-01, -5.1050352828199996E-01, -1.1205226103300001E-01, 2.0572316474100000E-01, 2.4457336207300001E-01, -1.6692934783800001E-01, -2.1949258712600000E-01, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 7.3546854659999996E-03, 2.5014419162000000E-02, 7.3256590681999995E-02, 1.6276129512400000E-01, 2.0975090011799999E-01, 2.9036216788000002E-01, 5.3183357832500000E-01, 6.3417969463099999E-01, 6.9577160811299998E-01, 7.2550836824400000E-01, 7.6737412259399995E-01, 8.0170944260900001E-01, 9.0849932972000003E-01, 9.9911350228200002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0012485372180440E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-6.0000000000000000E+00, -4.3750000000000000E+00, -2.3750000000000000E+00, -6.2500000000000000E-01, 8.7500000000000000E-01, 1.3750000000000000E+00, 1.6250000000000000E+00, 2.1250000000000000E+00, 2.3750000000000000E+00, 2.6250000000000000E+00, 2.8750000000000000E+00, 3.0000000000000000E+00};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.0401151181000001E-02, 2.1111743578000001E-02, 4.1510433147000002E-02, 1.3249816409400000E-01, 3.3941459637100002E-01, 5.5713943533800003E-01, 6.3835652414900002E-01, 4.6393398944899999E-01, 2.3328509953200000E-01, 5.3601824254000001E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -4.0570919869999996E-03, -6.8652098619999997E-03, -9.0956387959999999E-03, 2.1628823376299999E-01, 6.4039030177199996E-01, 7.2635190518699999E-01, -5.5101715884199998E-01, -8.6181670718299996E-01, -8.3353642789299998E-01, -6.2531498298400001E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 9.3871148699999993E-03, 4.1986796697999998E-02, 9.7607627793000001E-02, 1.8626456793800000E-01, 2.9591434189799998E-01, 4.0805443854700002E-01, 7.3505239115900001E-01, 8.7510508609500004E-01, 9.6251437802799999E-01, 9.9745231718100003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0046456927746406E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.4399999999999997E-01;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.6000000000000001E-01;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 2</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_2;