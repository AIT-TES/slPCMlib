
// within slPCMlib.Rubitherm_GR;
package Rubitherm_GR42 "Rubitherm GmbH, GR42; data taken from: Rubitherm datasheet; last access: 2020-07-03."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "GR42";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0914999999999998E+02, 3.1714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0814999999999998E+02, 3.1714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 4.2000000000000000E+04
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
    constant Real[10] data_x =   {3.6000000000000000E+01, 3.7625000000000000E+01, 3.8625000000000000E+01, 3.9875000000000000E+01, 4.1375000000000000E+01, 4.1625000000000000E+01, 4.2375000000000000E+01, 4.2625000000000000E+01, 4.3125000000000000E+01, 4.4000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.2998142425000002E-02, 9.9830062979999995E-02, 2.0923509193299999E-01, 2.8950843772700002E-01, 3.1686996158699998E-01, 1.7609547742100001E-01, 9.6228944255000004E-02, 2.5930692208000001E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.3101862474000000E-02, 1.2442475007100000E-01, 5.0312931307000000E-02, 8.9083213122000002E-02, 8.7708577760000003E-02, -3.0920183252000000E-01, -2.9432615304699999E-01, -6.7339722953999995E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 2.4005603531000001E-02, 8.1328349809000000E-02, 2.8480268432000000E-01, 6.5278203455899997E-01, 7.2883266613999997E-01, 9.3296063331500001E-01, 9.6703400232600001E-01, 9.9292884780699997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0032474620785703E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {3.5000000000000000E+01, 3.7125000000000000E+01, 3.8375000000000000E+01, 3.9375000000000000E+01, 4.1375000000000000E+01, 4.1625000000000000E+01, 4.1875000000000000E+01, 4.2625000000000000E+01, 4.2875000000000000E+01, 4.3625000000000000E+01, 4.4000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 4.9274778016999997E-02, 8.9354080372000003E-02, 1.5787378972800001E-01, 2.9606643837800001E-01, 3.7151991480100000E-01, 3.7042230291299999E-01, 6.3551690458999999E-02, 1.9165237325000001E-02, 1.4581174970000001E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.7036755030999998E-02, 8.7127470717000000E-02, 5.0683685479999997E-03, 2.2863216320800001E-01, 2.3225553271600000E-01, -2.0133485202699999E-01, -3.8604696312900000E-01, -6.8868326142999994E-02, -3.5733114529999999E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 4.6093167524999999E-02, 1.2386235529600000E-01, 2.5473955761599998E-01, 6.3539465475500001E-01, 7.1909588737300001E-01, 8.1440645802500000E-01, 9.8636349657699995E-01, 9.9507944452700003E-01, 9.9976772342099995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0032579655739240E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.0000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.0000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>GR42</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    other<br>  Data taken from: Rubitherm datasheet - last access 2020-07-03.<br><br>
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
end Rubitherm_GR42;