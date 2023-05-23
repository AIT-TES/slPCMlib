
// within slPCMlib.Axiotherm_ATS;
package Axiotherm_ATS_44 "Axiotherm GmbH, ATS 44; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 44";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1514999999999998E+02, 3.2214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1114999999999998E+02, 3.1814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9371000000000003E+05
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
    constant Real[9] data_x =   {4.2000000000000000E+01, 4.4625000000000000E+01, 4.5125000000000000E+01, 4.5375000000000000E+01, 4.5875000000000000E+01, 4.6125000000000000E+01, 4.6375000000000000E+01, 4.7875000000000000E+01, 4.9000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 1.4775345547000000E-02, 4.7123350165000001E-02, 1.1341889062799999E-01, 4.7729210749899997E-01, 6.0274359217200002E-01, 6.1934929985400000E-01, 7.4485334562999994E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 3.2059737500000001E-04, 1.1932639838600000E-01, 6.2664213872100005E-01, 7.2708356325400003E-01, 3.5064444681000001E-01, -3.6808861003799997E-01, -1.5862283245000000E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 1.9305962962999999E-02, 3.2367254663000000E-02, 4.9881137674000001E-02, 1.9620468391099999E-01, 3.3386437666300001E-01, 4.9116309493799998E-01, 9.7470411196100004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0050714377745864E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {3.8000000000000000E+01, 4.0125000000000000E+01, 4.1625000000000000E+01, 4.2625000000000000E+01, 4.3125000000000000E+01, 4.3375000000000000E+01, 4.3625000000000000E+01, 4.4125000000000000E+01, 4.4375000000000000E+01, 4.4625000000000000E+01, 4.4875000000000000E+01, 4.5000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.1491014755999999E-02, 7.7974253239999999E-02, 1.1888814153400000E-01, 2.0240534292400000E-01, 3.1872570696100000E-01, 5.3867768740099997E-01, 6.3380142615000001E-01, 4.6369141433599997E-01, 2.3413147200400000E-01, 5.3913839604000001E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.0114184619999999E-02, 7.9176666093999995E-02, 1.8621237844999999E-02, 2.9232833043899997E-01, 6.2143443927099995E-01, 7.2708513688400001E-01, -5.2442019828800002E-01, -8.5640028266699997E-01, -8.2924099856000000E-01, -6.2996545221500000E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 8.4589077240000005E-03, 6.2967527172999999E-02, 1.6713081941100000E-01, 2.4224651979299999E-01, 3.0609417849999998E-01, 4.1342600950899999E-01, 7.3473430349799995E-01, 8.7457065139199996E-01, 9.6223423295800004E-01, 9.9743375644099996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0066276173649826E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4890000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS 44</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATS_44;