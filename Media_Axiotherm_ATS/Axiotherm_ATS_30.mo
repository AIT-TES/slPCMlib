within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_30 "Axiotherm GmbH, ATS 30; data taken from: Axiotherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS 30";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9714999999999998E+02, 3.0914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9414999999999998E+02, 3.0414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7200000000000000E+05
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
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {2.4000000000000000E+01, 2.6375000000000000E+01, 2.8625000000000000E+01, 2.9375000000000000E+01, 3.0125000000000000E+01, 3.0625000000000000E+01, 3.1125000000000000E+01, 3.1375000000000000E+01, 3.1625000000000000E+01, 3.2125000000000000E+01, 3.2375000000000000E+01, 3.2625000000000000E+01, 3.2875000000000000E+01, 3.3875000000000000E+01, 3.6000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 6.6065773240000002E-03, 3.3282500239999999E-02, 7.8499805382000007E-02, 1.4670070359199999E-01, 1.0598711000600000E-01, 1.5308198285600000E-01, 2.5771029346000002E-01, 4.8953080802299997E-01, 5.8501389165800000E-01, 3.9424974187099998E-01, 1.6060892876399999E-01, 6.8211585524000001E-02, 1.2289985053000000E-02, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 6.0613310280000004E-03, 1.8893113864000000E-02, 1.1393548739700000E-01, -1.1172732802999999E-02, -7.2723801322000001E-02, 2.2409040067800001E-01, 6.1805797461199996E-01, 7.6317690637000002E-01, -6.4316292481100001E-01, -8.5182828308299996E-01, -6.8984633089799996E-01, -1.7619804702700001E-01, -8.6802273310000004E-03, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 5.0508520419999996E-03, 4.4944562783999997E-02, 8.2817848994999999E-02, 1.7412097550299999E-01, 2.3928069057699999E-01, 2.9850552874100000E-01, 3.4834220259300003E-01, 4.4200555315500001E-01, 7.4320135281300004E-01, 8.6804774775100002E-01, 9.3731130330800005E-01, 9.6352238886099995E-01, 9.9010110875799995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0109449346873389E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.1000000000000000E+01, 2.4375000000000000E+01, 2.7125000000000000E+01, 2.8125000000000000E+01, 2.8375000000000000E+01, 2.8625000000000000E+01, 2.9125000000000000E+01, 2.9375000000000000E+01, 2.9625000000000000E+01, 2.9875000000000000E+01, 3.0875000000000000E+01, 3.1000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 1.0803414371000001E-02, 6.8695253684999999E-02, 2.0435898453599999E-01, 3.3838955711899998E-01, 5.8630784769099997E-01, 6.5983865065699998E-01, 4.4652163102999998E-01, 1.8933332261800001E-01, 7.9385669490999997E-02, 7.0606932900000002E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -1.5140923850000001E-03, 5.4806763644000001E-02, 3.4489456075000002E-01, 6.9429402533600004E-01, 8.1313762298699999E-01, -7.3397999296799998E-01, -9.4375060895399998E-01, -7.8690817316499995E-01, -2.2689543377900001E-01, -7.2396504129999999E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.9646360471000001E-02, 9.3382065112000007E-02, 2.0561177607000000E-01, 2.7156302475700000E-01, 3.8640492546100003E-01, 7.2979552190499997E-01, 8.6902998876199999E-01, 9.4760855320000004E-01, 9.7824799834700005E-01, 9.9996533541800003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9890145424661236E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.3830000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.3000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.9999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS 30</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATS_30;
