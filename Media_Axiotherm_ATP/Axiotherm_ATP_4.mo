
within slPCMlib.Media_Axiotherm_ATP;
package Axiotherm_ATP_4 "Axiotherm GmbH, ATP 4; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 4";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6914999999999998E+02, 2.8014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6814999999999998E+02, 2.7814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2500000000000000E+05
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
    constant Real[13] data_x =   {-4.0000000000000000E+00, -1.6250000000000000E+00, 8.7500000000000000E-01, 3.6250000000000000E+00, 4.3750000000000000E+00, 4.6250000000000000E+00, 4.8750000000000000E+00, 5.1250000000000000E+00, 5.3750000000000000E+00, 5.6250000000000000E+00, 5.8750000000000000E+00, 6.6250000000000000E+00, 7.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 6.5285921510000004E-03, 2.2572500727000001E-02, 8.9607126046999996E-02, 3.5144917882200000E-01, 6.0771916914799995E-01, 7.3460404580500005E-01, 6.6241803453100001E-01, 4.3081430607200000E-01, 1.6684381486100000E-01, 6.4091472009000006E-02, 2.6003803530000001E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -7.5851082739999999E-03, 4.6448578960000004E-03, 9.4197944325000005E-02, 7.1865347805900004E-01, 8.3619521875199998E-01, 1.5080049008800001E-01, -8.3723607166900005E-01, -9.9078494832099995E-01, -8.3648425257500003E-01, -1.9300133621500001E-01, -1.0388679249999999E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.1325287634000000E-02, 4.1350937887999997E-02, 1.3922294573100000E-01, 2.7543416435899998E-01, 3.9479377250100001E-01, 5.6626277522699997E-01, 7.4615073944800003E-01, 8.8369181730500002E-01, 9.5764237165400001E-01, 9.8317401470800003E-01, 9.9963393866799999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0006351338712201E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-5.0000000000000000E+00, -3.7500000000000000E-01, 1.1250000000000000E+00, 2.6250000000000000E+00, 3.1250000000000000E+00, 3.3750000000000000E+00, 3.6250000000000000E+00, 4.1250000000000000E+00, 4.3750000000000000E+00, 4.6250000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[12] data_y =   {0.0000000000000000E+00, 6.5586133100000001E-03, 4.0892234557000001E-02, 4.9448181257999997E-02, 1.2591199515599999E-01, 2.5239045110800001E-01, 5.3791735882100000E-01, 7.6794790784099998E-01, 5.8857128832200001E-01, 3.0585687570900000E-01, 7.1500826993999997E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -4.3645661570000001E-03, 2.4866948151999999E-02, 1.3764727653999999E-02, 2.6960261098999999E-01, 8.1416690192499996E-01, 1.0029925081810001E+00, -4.3464701628199998E-01, -1.0503596663860000E+00, -1.0281879612600000E+00, -8.4507995048200002E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 2.3082606790000001E-02, 5.3367939443000001E-02, 1.2361805310300000E-01, 1.6235595589700000E-01, 2.0707045180399999E-01, 3.0545404762000000E-01, 6.6397964888700001E-01, 8.3777345297200001E-01, 9.5012221014700005E-01, 9.9661163612699999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0059157120639206E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.3300000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.5000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.5000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 4</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_4;