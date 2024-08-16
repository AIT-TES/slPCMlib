within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_23 "Axiotherm GmbH, ATP 23; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 23";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8814999999999998E+02, 2.9914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8714999999999998E+02, 2.9714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9700000000000000E+05
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
    constant Real[10] data_x =   {1.5000000000000000E+01, 1.8125000000000000E+01, 2.0625000000000000E+01, 2.1625000000000000E+01, 2.2125000000000000E+01, 2.2375000000000000E+01, 2.3625000000000000E+01, 2.4625000000000000E+01, 2.5125000000000000E+01, 2.6000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 5.1315788150000001E-03, 3.6816523682000002E-02, 5.0738953869999999E-02, 9.3483392183999994E-02, 1.5862107620999999E-01, 5.3533482190899995E-01, 4.7305927732999999E-02, 2.3498966670000001E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 2.7306167450000000E-03, 3.2967314002000002E-02, 3.2934169129999998E-03, 1.5372530103400001E-01, 5.1121683304900001E-01, -2.6908185542500002E-01, -2.8126542142200001E-01, -4.2413846950000002E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 5.8760262349999998E-03, 4.3069990916000001E-02, 8.9959868521000005E-02, 1.2333652725600000E-01, 1.5339751343499999E-01, 6.9612103632400002E-01, 9.9249761710100004E-01, 9.9923205899400003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0138228385290353E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.4000000000000000E+01, 1.5625000000000000E+01, 1.7625000000000000E+01, 1.9625000000000000E+01, 2.0375000000000000E+01, 2.0625000000000000E+01, 2.0875000000000000E+01, 2.1125000000000000E+01, 2.1625000000000000E+01, 2.1875000000000000E+01, 2.2375000000000000E+01, 2.3625000000000000E+01, 2.4000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 9.2007465739999995E-03, 1.8472352184000000E-02, 9.0448828982999996E-02, 4.2049172726099998E-01, 7.9052300376100004E-01, 9.0200799156500000E-01, 6.5916671556899997E-01, 1.8990633518999999E-02, 2.1966686800000001E-04, 2.1966686800000001E-04, 5.8342907749999997E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -3.5408481310000000E-03, -6.4200675900000000E-03, 8.6228865143000002E-02, 9.4702528455400004E-01, 1.1298369243419999E+00, -4.5414570699099999E-01, -1.3950677266930001E+00, -1.8206741433699999E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, -1.6385241602999999E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 8.3455974309999993E-03, 3.7293457166999999E-02, 1.1619024434300000E-01, 2.6910721901500001E-01, 4.2118690297200001E-01, 6.4342162616300003E-01, 8.4567004102200005E-01, 9.9152579390399997E-01, 9.9299479992799999E-01, 9.9310584176000005E-01, 9.9908816207700002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0110020931747230E+00;
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
    lambda := 2.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 23</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_23;
