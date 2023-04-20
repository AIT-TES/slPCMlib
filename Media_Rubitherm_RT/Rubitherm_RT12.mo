
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT12 "Rubitherm GmbH, RT12; data taken from: Rubitherm datasheet; last access: 2020-08-30."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT12";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7614999999999998E+02, 2.9114999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7314999999999998E+02, 2.8814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4228400413534552E+05
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
    constant Real[9] data_x =   {3.0000000000000000E+00, 4.3750000000000000E+00, 5.6250000000000000E+00, 7.8750000000000000E+00, 1.0625000000000000E+01, 1.2875000000000000E+01, 1.4375000000000000E+01, 1.6625000000000000E+01, 1.8000000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 6.6775323760000005E-02, 7.7280599086999996E-02, 1.0574489659700000E-01, 1.1415862610100000E-01, 7.6797680384000006E-02, 2.2947056552000002E-02, 7.3047384219999999E-03, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, -6.4546775760000001E-03, 3.0904526977999999E-02, 1.7818702344000001E-02, -1.7604058094000001E-02, -5.0469905230999999E-02, 1.4801361600000000E-03, 3.3600190890000002E-03, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 4.6889281211999999E-02, 1.3199495327499999E-01, 3.4325836099300000E-01, 6.6770238993900000E-01, 8.9621951971600000E-01, 9.6123793425600002E-01, 9.9445283796899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9923917344995683E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {0.0000000000000000E+00, 2.1250000000000000E+00, 4.6250000000000000E+00, 8.3750000000000000E+00, 9.6250000000000000E+00, 1.0625000000000000E+01, 1.1625000000000000E+01, 1.2875000000000000E+01, 1.4625000000000000E+01, 1.5000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 4.2374966898000001E-02, 7.2889295827000003E-02, 8.6593929700000002E-02, 1.1577997234500000E-01, 1.0662189154000000E-01, 1.3000260852800000E-01, 4.4723157559999997E-02, 8.0903149099999998E-04, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.7208604225999999E-02, 2.5410599552000000E-02, -3.3275418099999997E-04, 4.1700856733000002E-02, -3.2089366585000001E-02, 4.7204034380000003E-02, -6.8247675674999994E-02, -3.2854599080000000E-03, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.8347146089999999E-02, 1.7742796608100000E-01, 5.0491368374000001E-01, 6.2529444425299996E-01, 7.4203381048600003E-01, 8.5315690901999996E-01, 9.7674669690399996E-01, 9.9988739718800002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9479548039609611E-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT12</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-08-30.<br><br>
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
end Rubitherm_RT12;