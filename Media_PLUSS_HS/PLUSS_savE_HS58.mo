
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS58 "Pluss Advanced Technologies Pvt Ltd, HS58; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS58";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2614999999999998E+02, 3.3514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2114999999999998E+02, 3.3014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.4000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.1500000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4294479887492157E+05
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
    constant Real[11] data_x =   {5.3000000000000000E+01, 5.5125000000000000E+01, 5.6625000000000000E+01, 5.7125000000000000E+01, 5.7375000000000000E+01, 5.7875000000000000E+01, 5.8625000000000000E+01, 5.9375000000000000E+01, 5.9875000000000000E+01, 6.1125000000000000E+01, 6.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.7890271316000000E-02, 4.9854704844000002E-02, 9.8222253608999999E-02, 1.6581807992299999E-01, 4.2127860407099998E-01, 4.7402690165200001E-01, 1.8015034910100000E-01, 4.6489328883000003E-02, 2.1908727658000001E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.5589558450000001E-02, 1.4959077534000000E-02, 1.6495138861000000E-01, 4.7460534961799999E-01, 4.7191401682900003E-01, -2.9957881672800002E-01, -4.2510770459699998E-01, -8.6147777650999999E-02, -2.9835297988000001E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.3248200970000000E-02, 6.4586573738999994E-02, 9.8754797520999998E-02, 1.3040066802200001E-01, 2.7841711153100002E-01, 6.5332491208700005E-01, 9.0655493191799996E-01, 9.5655387758400001E-01, 9.9225642887599996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0080787382329168E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {4.8000000000000000E+01, 4.9625000000000000E+01, 5.1625000000000000E+01, 5.3625000000000000E+01, 5.4625000000000000E+01, 5.5375000000000000E+01, 5.5625000000000000E+01, 5.6125000000000000E+01, 5.6625000000000000E+01, 5.6875000000000000E+01, 5.7000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 4.3904390378000001E-02, 9.0456158043000001E-02, 1.0870033504300000E-01, 1.1803055949100000E-01, 2.3151281748600000E-01, 3.4074613359700001E-01, 3.5191234540900002E-01, 1.2315665416200000E-01, 2.8048535547999999E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 2.3797812764999999E-02, -1.8892805540000001E-03, 2.4830041290000000E-02, 3.1038844020000000E-03, 3.1227830690000002E-01, 3.4812424367599998E-01, -3.5281484790899997E-01, -4.5786148131400001E-01, -3.2514259167100001E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.0589156685999999E-02, 1.7423335225400000E-01, 3.6544353120400003E-01, 4.8120074494199999E-01, 5.9837533072100002E-01, 6.7008105904500004E-01, 8.5879617451699997E-01, 9.8036231779899996E-01, 9.9866361888800004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0050466555242761E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3200000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 6.0999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.1000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS58</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: multiple options available<br>  The data is taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS58;