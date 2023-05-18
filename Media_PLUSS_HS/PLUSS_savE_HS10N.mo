
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS10N "Pluss Advanced Technologies Pvt Ltd, HS10N; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS10N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.5814999999999998E+02, 2.6914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5714999999999998E+02, 2.6514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.9000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.0231630885743711E+05
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
    constant Real[13] data_x =   {-1.5000000000000000E+01, -1.2375000000000000E+01, -1.1875000000000000E+01, -1.1625000000000000E+01, -1.1375000000000000E+01, -1.1125000000000000E+01, -1.0625000000000000E+01, -1.0375000000000000E+01, -9.8750000000000000E+00, -8.8750000000000000E+00, -7.6250000000000000E+00, -6.3750000000000000E+00, -4.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 4.4341300740000002E-03, 2.3170789750000000E-02, 9.1629218814000005E-02, 3.8384866496100001E-01, 6.9697679916699995E-01, 6.7541036857199999E-01, 3.5514762121600002E-01, 1.2009493786800000E-01, 6.0803795329000002E-02, 2.1055104916999998E-02, 6.1763127709999997E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -1.0743205848000001E-02, 7.4961214146999994E-02, 9.1429830706899995E-01, 1.2143249199960000E+00, 1.1987879971380000E+00, -1.0216963403649999E+00, -8.3975206326499996E-01, -1.7765215820100000E-01, -2.2646684845999999E-02, -4.6535322103000001E-02, 5.4881699810000002E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.2056874340999999E-02, 1.7201666630000001E-02, 2.7236824848999999E-02, 8.5437791180000000E-02, 2.2139010568699999E-01, 6.1295958434499997E-01, 7.4155836487899995E-01, 8.4717203557999998E-01, 9.2514487148299995E-01, 9.7972559677000004E-01, 9.9002956406200004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0056826948730080E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-1.6000000000000000E+01, -1.4375000000000000E+01, -1.3625000000000000E+01, -1.3375000000000000E+01, -1.2625000000000000E+01, -1.2375000000000000E+01, -1.2125000000000000E+01, -1.1875000000000000E+01, -1.1625000000000000E+01, -1.1375000000000000E+01, -1.1125000000000000E+01, -1.0125000000000000E+01, -8.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.2163334075999999E-02, 1.9168825930000000E-03, 1.9378562097000000E-02, 1.9378562097000000E-02, 6.9385106827900001E-01, 1.3081961783730001E+00, 1.2907351144690000E+00, 5.0746013136699997E-01, 3.8407113854000002E-02, 5.9666241570000000E-03, 2.5657095039999999E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.4258573469000000E-02, -1.8457819806999998E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 2.5737074099610000E+00, 2.2734019451909999E+00, -2.1156437071790002E+00, -2.3326000579619999E+00, -2.9319383893200002E-01, -1.7493838890000001E-03, 3.7185539969999998E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 6.7354215949999999E-03, 1.3539327031000000E-02, 1.6101449075000000E-02, 3.0614558956000001E-02, 1.0625506873200000E-01, 3.5771447461099998E-01, 7.0494257647799996E-01, 9.3052348186300005E-01, 9.8805248396500001E-01, 9.9207549238199999E-01, 9.9588054131000003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9856806083973493E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.0570000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1250000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 4.2500000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS10N</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS10N;