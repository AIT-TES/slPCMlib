
// within slPCMlib.Croda_Crodatherm;
package Croda_Crodatherm_32 "Croda International Plc, Crodatherm 32; data taken from: Croda datasheet; last access: 2023-02-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 32";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9514999999999998E+02, 3.0714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7640157456873904E+05
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
    constant Real[10] data_x =   {2.2000000000000000E+01, 2.4125000000000000E+01, 2.5875000000000000E+01, 2.7625000000000000E+01, 2.9375000000000000E+01, 3.0875000000000000E+01, 3.1625000000000000E+01, 3.2625000000000000E+01, 3.3125000000000000E+01, 3.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 2.2170956243000001E-02, 8.6343232588999994E-02, 8.2624594176000002E-02, 9.6911267509999996E-02, 1.8186228384700001E-01, 2.5942175480700003E-01, 7.3378559120999998E-02, 1.7623465724999999E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 3.0425731612999999E-02, 9.5092856360000002E-03, 8.5350543590000001E-03, 5.9840037040000002E-02, 6.9981090905000001E-02, 1.0470097378600000E-01, -2.5767168694499998E-01, -4.7667972251000001E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.2166183283000001E-02, 1.1294127289800000E-01, 2.6175608906199999E-01, 4.0645595401399998E-01, 6.1464098739600004E-01, 7.7929088624800003E-01, 9.7684368870100002E-01, 9.9530837145499995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0048572434632574E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.0000000000000000E+01, 2.2125000000000000E+01, 2.4875000000000000E+01, 2.8125000000000000E+01, 3.0375000000000000E+01, 3.0625000000000000E+01, 3.0875000000000000E+01, 3.1625000000000000E+01, 3.1875000000000000E+01, 3.2000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 4.8262520589999998E-03, 5.6780079141999999E-02, 8.4525520054000006E-02, 2.7201161049799999E-01, 3.5863680352200000E-01, 3.8117358159699999E-01, 1.1027140195500000E-01, 2.4860865180000000E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 5.0628029359999997E-03, 2.7244508306000000E-02, 3.7571664230000002E-03, 2.6921926538699997E-01, 2.8147923385800000E-01, -6.8716725805000003E-02, -4.3358441254500002E-01, -2.8617577670299998E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.2063088660000002E-03, 7.3575096037999999E-02, 3.2259358277700001E-01, 6.1023114266500000E-01, 6.8859651450699999E-01, 7.8241568321400001E-01, 9.8278332071800001E-01, 9.9882484634799995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9489856286613787E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.1600000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.3600000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.2000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.3600000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 32</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
       material class: paraffin-based;  encapsulation:    none<br>  Data taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_32;