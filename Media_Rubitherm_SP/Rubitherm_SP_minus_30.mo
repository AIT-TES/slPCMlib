within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_30 "Rubitherm GmbH, SP-30; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-30";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.3914999999999998E+02, 2.4714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.3614999999999998E+02, 2.4614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2961894212444883E+05
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
    constant Real[13] data_x =   {-3.4000000000000000E+01, -3.2125000000000000E+01, -3.0625000000000000E+01, -3.0125000000000000E+01, -2.9875000000000000E+01, -2.9625000000000000E+01, -2.9375000000000000E+01, -2.9125000000000000E+01, -2.8625000000000000E+01, -2.8375000000000000E+01, -2.7875000000000000E+01, -2.7125000000000000E+01, -2.6000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.9512584410000000E-03, 2.4640126644000001E-02, 2.8212791790000000E-02, 6.2206059346000001E-02, 1.6759205369700000E-01, 4.8714920438900000E-01, 7.6939387421500005E-01, 6.5010300992600001E-01, 3.2315992100599999E-01, 1.0009414359200000E-01, 5.3154594414999999E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.6604219040000001E-03, -7.5311399459999998E-03, 5.0865050771000002E-02, 1.6620622811200000E-01, 9.8417618908899995E-01, 1.2035062518199999E+00, 1.0204883991320000E+00, -1.0830043629749999E+00, -8.7128829406599995E-01, -1.5330468413900000E-01, -5.7803167545999998E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.2861052250000000E-03, 2.4759608424000001E-02, 3.6786512131000003E-02, 4.7515127951999997E-02, 7.2041343324000001E-02, 1.5294522220300000E-01, 3.1136493158799999E-01, 7.1106762747600005E-01, 8.3192690515900003E-01, 9.2301161154400002E-01, 9.7613692557099996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0025224895297680E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-3.7000000000000000E+01, -3.4375000000000000E+01, -3.2375000000000000E+01, -3.1625000000000000E+01, -3.1375000000000000E+01, -3.0375000000000000E+01, -2.9625000000000000E+01, -2.9125000000000000E+01, -2.8875000000000000E+01, -2.7375000000000000E+01, -2.7000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.2145052895999998E-02, 1.2976582069900000E-01, 2.7850840553000000E-01, 3.7055112348800001E-01, 3.8336087821499998E-01, 7.4388583773000005E-02, 2.6780353649999999E-03, 9.1441908499999996E-04, 4.3819690180000002E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -3.5141148360000002E-03, 1.0401637937600000E-01, 3.0782668600000002E-01, 3.2598083318900001E-01, -3.5208248268800002E-01, -3.6315050590499998E-01, -8.4433211090000002E-03, -9.9488734400000000E-04, -1.3795097293000000E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.1186374823000001E-02, 1.4763881722600000E-01, 2.9166428502300001E-01, 3.7297102582300001E-01, 8.0787036068100004E-01, 9.8061643243899999E-01, 9.9253275726199996E-01, 9.9294438123600004E-01, 9.9933785261899999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0033176151127667E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.1000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.2000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP-30</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-06-03.<br><br>
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
end Rubitherm_SP_minus_30;
