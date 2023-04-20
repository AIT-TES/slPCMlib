
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT_minus_9HC "Rubitherm GmbH, RT-9HC; data taken from: Rubitherm datasheet; last access: 2020-09-29."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT-9HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6014999999999998E+02, 2.6514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5914999999999998E+02, 2.6514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.2100000000000000E+05
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
    constant Real[10] data_x =   {-1.3000000000000000E+01, -1.1375000000000000E+01, -1.0875000000000000E+01, -1.0625000000000000E+01, -1.0125000000000000E+01, -9.6250000000000000E+00, -8.8750000000000000E+00, -8.3750000000000000E+00, -8.1250000000000000E+00, -8.0000000000000000E+00};
    constant Real[10] data_y =   {0.0000000000000000E+00, 7.2362444369999997E-03, 2.7266418644999998E-02, 8.1026703828999999E-02, 4.6774564806500002E-01, 6.8794811936199995E-01, 3.5438508406199998E-01, 9.4016588797000006E-02, 2.0184941242000001E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, -5.3052855050000000E-03, 7.7973621098000004E-02, 4.7495327045800001E-01, 7.9702378317199996E-01, -1.4272157251199999E-01, -5.6913097758400000E-01, -3.9791458851500000E-01, -2.2534950782700000E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 7.0600743850000001E-03, 1.3963657294999999E-02, 2.5454156966000001E-02, 1.5618161245400000E-01, 4.6526037305899998E-01, 8.7689210188699995E-01, 9.8562860764700000E-01, 9.9903005336600004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0018712731078587E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {-1.4000000000000000E+01, -1.2875000000000000E+01, -1.1375000000000000E+01, -1.0875000000000000E+01, -1.0625000000000000E+01, -1.0375000000000000E+01, -1.0125000000000000E+01, -9.8750000000000000E+00, -9.6250000000000000E+00, -9.3750000000000000E+00, -8.3750000000000000E+00, -8.0000000000000000E+00};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.7873713070000004E-03, 2.9683662537000000E-02, 8.6428950177999994E-02, 1.9068912296200000E-01, 4.4697120420399999E-01, 6.7277696080399996E-01, 7.6011280467300002E-01, 6.6716509547000002E-01, 4.3842954702199999E-01, 3.2886685002000003E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.0392891519000001E-02, 8.9416797959999995E-03, 2.0430887619900001E-01, 8.2586867888799997E-01, 9.6133211871400004E-01, 8.2313085593099999E-01, -1.4749856980000000E-02, -7.7020413019400003E-01, -7.2195791251899999E-01, -1.3399014488100000E-01, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.5986936310000001E-03, 2.7755526746000000E-02, 5.2743579988999999E-02, 8.4183857761000000E-02, 1.6328098006899999E-01, 3.0413869079400002E-01, 4.8783478171200001E-01, 6.7039872359700003E-01, 8.0851286403300004E-01, 9.9539840981899996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0012040221691234E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT-9HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-09-29.<br><br>
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
end Rubitherm_RT_minus_9HC;