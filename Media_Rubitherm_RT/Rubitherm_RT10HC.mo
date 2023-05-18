
within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT10HC "Rubitherm GmbH, RT10HC; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT10HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7914999999999998E+02, 2.8614999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7814999999999998E+02, 2.8314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7561812317773499E+05
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
    constant Real[11] data_x =   {6.0000000000000000E+00, 7.6250000000000000E+00, 8.1250000000000000E+00, 8.3750000000000000E+00, 9.1250000000000000E+00, 9.3750000000000000E+00, 1.0125000000000000E+01, 1.0375000000000000E+01, 1.0875000000000000E+01, 1.1875000000000000E+01, 1.3000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.6497660713999999E-02, 5.5821535396000002E-02, 1.3245254759899999E-01, 6.7010948899199996E-01, 6.7924701829599998E-01, 2.6348730889900002E-01, 1.3145595465200000E-01, 1.3501657005999999E-02, 1.6918622122000002E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 7.1940056350000003E-03, 1.4201688426199999E-01, 7.3843111431700004E-01, 3.5186740069400002E-01, -3.7829448673499999E-01, -5.9659170094600000E-01, -4.9191255742000001E-01, -3.4764907847000003E-02, 1.1744375878000000E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.1815862193000000E-02, 2.7079835100000001E-02, 4.7498385887999997E-02, 3.6643272418900003E-01, 5.3882597282199995E-01, 9.0241690135099994E-01, 9.5121717417399998E-01, 9.7792039347600002E-01, 9.8924955190099995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9954055701916911E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {5.0000000000000000E+00, 6.8750000000000000E+00, 8.1250000000000000E+00, 8.3750000000000000E+00, 8.6250000000000000E+00, 9.1250000000000000E+00, 9.3750000000000000E+00, 9.6250000000000000E+00, 9.8750000000000000E+00, 1.0000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 9.0955950450000000E-03, 2.3471566783600001E-01, 3.9622184779800002E-01, 6.5278548860200003E-01, 7.5044997517400003E-01, 5.4585568609799995E-01, 2.7462472998400000E-01, 6.3118711785999995E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.3766132622000000E-02, 4.7440567586200000E-01, 7.7548613751899997E-01, 8.6622918009799998E-01, -6.4442186218800002E-01, -1.0132409706620000E+00, -9.8012506348700001E-01, -7.3648230359599998E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 4.5057356260000004E-03, 9.7148449585000005E-02, 1.7464809955800001E-01, 3.0564043658700002E-01, 6.8891320227899999E-01, 8.5329781009700001E-01, 9.5595107501300003E-01, 9.9700629340600000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0025949593350489E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.7000000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT10HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-06-10.<br><br>
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
end Rubitherm_RT10HC;