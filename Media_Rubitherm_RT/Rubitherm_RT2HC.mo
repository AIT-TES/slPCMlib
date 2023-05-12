
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT2HC "Rubitherm GmbH, RT2HC; data taken from: Rubitherm datasheet; last access: 2020-09-30."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT2HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6614999999999998E+02, 2.7814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6614999999999998E+02, 2.7614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7066833224395153E+05
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
    constant Real[11] data_x =   {-7.0000000000000000E+00, -4.1250000000000000E+00, 1.2500000000000000E-01, 1.6250000000000000E+00, 2.3750000000000000E+00, 2.6250000000000000E+00, 2.8750000000000000E+00, 3.6250000000000000E+00, 3.8750000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.4756871680000002E-02, 1.0624829598400000E-01, 1.9904045603900000E-01, 2.7020569742299999E-01, 3.0557999340300002E-01, 2.9191500689200001E-01, 6.7203790018000004E-02, 2.9764152709000000E-02, 3.0691769500000002E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 3.3280218650000000E-03, 5.1045991840999998E-02, 6.1162516110000002E-02, 1.1830595377800000E-01, 1.1673429186100000E-01, -1.7172347702800000E-01, -3.1746815987499999E-01, -8.1956896545999997E-02, -3.1703918320000000E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 4.8401960803000001E-02, 2.7970725212899999E-01, 5.1026039979100002E-01, 6.8620749677399995E-01, 7.5929314635300005E-01, 8.3665121007700005E-01, 9.8032328507400002E-01, 9.9138478516899997E-01, 9.9998471481399998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0153408037493750E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-7.0000000000000000E+00, -4.3750000000000000E+00, -1.8750000000000000E+00, 3.7500000000000000E-01, 8.7500000000000000E-01, 1.3750000000000000E+00, 1.6250000000000000E+00, 2.1250000000000000E+00, 2.3750000000000000E+00, 2.6250000000000000E+00, 2.8750000000000000E+00, 3.0000000000000000E+00};
    constant Real[12] data_y =   {0.0000000000000000E+00, 3.2890756755999999E-02, 5.9549842114000003E-02, 9.9293689685000006E-02, 1.1217643958699999E-01, 2.4165927291700001E-01, 4.0738703665900000E-01, 4.7822165575600001E-01, 3.4966731678700003E-01, 1.7649436065999999E-01, 4.0634125558000000E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.4617826719999999E-03, 1.8457985670999998E-02, -5.8261600970000004E-03, 9.8360610739000001E-02, 4.6634341821600001E-01, 5.4646516044899995E-01, -3.9718288291800002E-01, -6.4612439293599999E-01, -6.2556827136500004E-01, -4.7472922614599999E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 4.2759057256000002E-02, 1.5053978850100000E-01, 3.4139996612200002E-01, 3.9261112602699999E-01, 4.7422317086100002E-01, 5.5575528143399999E-01, 7.9926171349599995E-01, 9.0510714579200002E-01, 9.7143626946799999E-01, 9.9805901572599998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0101423273621888E+00;
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
    lambda := 7.7000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT2HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-09-30.<br><br>
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
end Rubitherm_RT2HC;