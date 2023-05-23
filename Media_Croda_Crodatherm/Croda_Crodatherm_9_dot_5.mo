
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_9_dot_5 "Croda International Plc, Crodatherm 9.5; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 9.5";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7714999999999998E+02, 2.8414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7614999999999998E+02, 2.8314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.2000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.1000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9228338115194495E+05
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
    constant Real[11] data_x =   {4.0000000000000000E+00, 6.8750000000000000E+00, 7.8750000000000000E+00, 8.3750000000000000E+00, 8.6250000000000000E+00, 9.1250000000000000E+00, 9.3750000000000000E+00, 9.6250000000000000E+00, 9.8750000000000000E+00, 1.0625000000000000E+01, 1.1000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 2.1744502448000001E-02, 8.5000057276999996E-02, 3.1828905086699999E-01, 6.5856453195800002E-01, 8.2161386757199995E-01, 5.4457046363100003E-01, 2.0752646445600001E-01, 7.8433796058999999E-02, 3.1419755390000000E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 2.1881931236000001E-02, 1.6549697601900001E-01, 9.2406348624699997E-01, 1.1476718317679999E+00, -9.3533677723099995E-01, -1.2340437845429999E+00, -9.9116206520700001E-01, -2.3801975628600000E-01, -1.2621768506999999E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.6277268471000000E-02, 5.7916633145999999E-02, 1.4341799579000000E-01, 2.6504651670700002E-01, 6.8083375894700005E-01, 8.5414067962399998E-01, 9.4741420596899995E-01, 9.7941724181699996E-01, 9.9955628669700003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0056758566973205E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {3.0000000000000000E+00, 5.1250000000000000E+00, 6.8750000000000000E+00, 8.1250000000000000E+00, 8.3750000000000000E+00, 8.6250000000000000E+00, 8.8750000000000000E+00, 9.3750000000000000E+00, 9.6250000000000000E+00, 9.8750000000000000E+00, 1.0000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 6.3041928520000002E-03, 4.1620528431000003E-02, 2.6707373308900001E-01, 4.0286500796600000E-01, 6.0928687471400000E-01, 7.0431866349899996E-01, 4.6048847481100003E-01, 2.2870895360300000E-01, 5.2208100408000002E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 7.3878370649999997E-03, 4.2047023246000001E-02, 4.1344530327899998E-01, 6.2633056707000001E-01, 6.8505246684599996E-01, 9.4631188071000003E-02, -8.7143235722900003E-01, -8.4091517229699997E-01, -6.0619448098299999E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.9521790380000000E-03, 3.7328403048000003E-02, 1.8315892678300000E-01, 2.6651023854299999E-01, 3.9381965485999998E-01, 5.6254839817100000E-01, 8.7658067014100005E-01, 9.6331831159000003E-01, 9.9750482323599998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0086858958712355E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.6300000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.5800000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.3999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4999999999999999E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 9.5</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_9_dot_5;