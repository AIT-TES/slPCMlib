
// within slPCMlib.Croda_Crodatherm;
package Croda_Crodatherm_37 "Croda International Plc, Crodatherm 37; data taken from: Croda datasheet; last access: 2023-02-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 37";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0214999999999998E+02, 3.1314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0214999999999998E+02, 3.1014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9044172893814975E+05
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
    constant Real[13] data_x =   {2.9000000000000000E+01, 3.1875000000000000E+01, 3.3875000000000000E+01, 3.4875000000000000E+01, 3.5125000000000000E+01, 3.5375000000000000E+01, 3.5625000000000000E+01, 3.5875000000000000E+01, 3.6375000000000000E+01, 3.6625000000000000E+01, 3.7375000000000000E+01, 3.8625000000000000E+01, 4.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.8558059859999999E-03, 3.1469059840999999E-02, 4.6456112867999999E-02, 8.9947697090000001E-02, 2.0570863748900001E-01, 5.0655940057500004E-01, 7.5305330679799998E-01, 6.3243906568099995E-01, 3.3404247314399999E-01, 4.8434393728000000E-02, 5.9760216000000000E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 2.7505969710000001E-03, 2.3538128576000000E-02, 7.6793229106999997E-02, 2.1416396849700001E-01, 8.8558888869100005E-01, 1.0974381669650000E+00, 8.3401754392799998E-01, -9.9944578111500004E-01, -8.4011992395500001E-01, -1.1313820900500000E-01, -3.3232200249999999E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.2153128750000001E-03, 2.9669355634999998E-02, 6.4267556136000006E-02, 8.0637354123999999E-02, 1.1416866264400000E-01, 2.0228608113999999E-01, 3.6144795382500000E-01, 7.4683735268700002E-01, 8.6707328683200002E-01, 9.7665776810200000E-01, 9.9640742930500004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0021300481421873E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {2.9000000000000000E+01, 3.2625000000000000E+01, 3.4625000000000000E+01, 3.5125000000000000E+01, 3.5375000000000000E+01, 3.5625000000000000E+01, 3.6125000000000000E+01, 3.6375000000000000E+01, 3.6625000000000000E+01, 3.6875000000000000E+01, 3.7000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 8.5834945349999997E-03, 8.2405178470000007E-02, 1.9248571286899999E-01, 3.3183930825699998E-01, 5.9481956796799995E-01, 7.3713420851699996E-01, 5.4623558937200001E-01, 2.7801675429200001E-01, 6.4288595791999997E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 6.2361341129999999E-03, 6.8380876660999998E-02, 3.5237766545100002E-01, 7.5295498633900004E-01, 8.8755415434700002E-01, -5.5831070332199995E-01, -9.9830650302299995E-01, -9.6914664688300001E-01, -7.5353091485400003E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 8.7831530290000005E-03, 7.9495370171000002E-02, 1.4269335881799999E-01, 2.0654355590900000E-01, 3.2239321574300001E-01, 6.8776939669799997E-01, 8.5149747890700000E-01, 9.5501904364199997E-01, 9.9694417485300002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0062392894073970E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.5700000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.1900000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.3999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.1900000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 37</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
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
end Croda_Crodatherm_37;