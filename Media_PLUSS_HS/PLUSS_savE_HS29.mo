
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS29 "Pluss Advanced Technologies Pvt Ltd, HS29; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS29";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9714999999999998E+02, 3.0814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9614999999999998E+02, 3.0314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.5100000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.6200000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9331117724030360E+05
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {2.4000000000000000E+01, 2.5625000000000000E+01, 2.7125000000000000E+01, 2.7875000000000000E+01, 2.9375000000000000E+01, 3.0625000000000000E+01, 3.3375000000000000E+01, 3.4875000000000000E+01, 3.5000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 5.8381403005999997E-02, 1.4724120247600000E-01, 2.9826970592199997E-01, 1.8177434577599999E-01, 6.0457180404999998E-02, 1.2906343385000000E-02, 7.0699799000000001E-05, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 3.5581970373000003E-02, 1.3389457683100001E-01, 1.5918154292100001E-01, -1.7350899497399999E-01, -3.2798098908999998E-02, -2.3980743244000000E-02, -7.5161993999999996E-04, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 3.9706097697000002E-02, 1.7583603430799999E-01, 3.4214072162999998E-01, 7.6563147292900002E-01, 8.9904412334499995E-01, 9.9460550340499998E-01, 9.9999655115300001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0025525680059000E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {2.3000000000000000E+01, 2.4375000000000000E+01, 2.6625000000000000E+01, 2.7125000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.8375000000000000E+01, 2.9625000000000000E+01, 3.0000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 4.9673914240999997E-02, 9.7021127010000000E-02, 1.5017833773100001E-01, 2.1547497662800000E-01, 3.2838143355700000E-01, 4.4402742384100002E-01, 5.1225991318000000E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 1.2805384203000000E-02, 1.6564447286000001E-02, 1.7773060382600001E-01, 3.6512517746500001E-01, 4.0246293268000000E-01, -2.2590503022600000E-01, -2.1476607151300001E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 3.2107005674000003E-02, 1.9541932072500001E-01, 2.5381373404200003E-01, 2.9850778082800000E-01, 3.6623989251099998E-01, 6.8508682893700001E-01, 9.9291771686800001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9918167735503061E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.6810000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5300000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 4.7799999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 3.8200000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS29</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
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
end PLUSS_savE_HS29;