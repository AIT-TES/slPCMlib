
// within slPCMlib.Rubitherm_SP;
package Rubitherm_SP58 "Rubitherm GmbH, SP58; data taken from: Rubitherm datasheet; last access: 2020-07-12."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP58";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2214999999999998E+02, 3.3514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2114999999999998E+02, 3.3014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.1844434476354884E+05
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
    constant Real[11] data_x =   {4.9000000000000000E+01, 5.0625000000000000E+01, 5.3625000000000000E+01, 5.5875000000000000E+01, 5.7375000000000000E+01, 5.8625000000000000E+01, 5.9125000000000000E+01, 5.9625000000000000E+01, 6.0125000000000000E+01, 6.1625000000000000E+01, 6.2000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.0604854553000000E-02, 2.6995057558000000E-02, 1.1533825790700000E-01, 2.5958830273400002E-01, 2.6857082744100003E-01, 1.7596682447699999E-01, 4.1349907557999999E-02, 8.5566724239999996E-03, 1.1050314039999999E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -8.7411027300000003E-03, 5.6283146139999999E-03, 8.0059072612000001E-02, -1.0338094500000001E-03, 1.2839996665999999E-02, -2.7653056431799999E-01, -1.9200186127800001E-01, -1.5831899523000001E-02, -4.8107433099999998E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.0609452576000000E-02, 5.6533137827000000E-02, 1.8610657328399999E-01, 4.8446115107200000E-01, 8.1491919779400002E-01, 9.3285484796200002E-01, 9.8576970073500003E-01, 9.9463421620299997E-01, 9.9984818787600005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0065949494014310E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {4.8000000000000000E+01, 5.0875000000000000E+01, 5.2875000000000000E+01, 5.4125000000000000E+01, 5.4375000000000000E+01, 5.4625000000000000E+01, 5.5375000000000000E+01, 5.5875000000000000E+01, 5.6875000000000000E+01, 5.7000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.2868615002000001E-02, 3.5080257294000002E-02, 1.6379839170999999E-01, 2.5826097338900000E-01, 4.1647191153300001E-01, 5.1516405230300000E-01, 2.7205516496499998E-01, 8.3951719159999997E-03, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 6.5550233999999994E-05, 2.4875258359000001E-02, 2.6549081763200000E-01, 4.9964746570899998E-01, 5.5827251636899999E-01, -4.7526873855899998E-01, -4.3208039366200002E-01, -9.2158001466999998E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 1.8589289456000000E-02, 5.8560272378000003E-02, 1.5221347297700000E-01, 2.0413061533099999E-01, 2.8878534352200003E-01, 6.8952372384799998E-01, 8.8687051421700003E-01, 9.9959232079799998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0073593981286764E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.2000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.2000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP58</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP58;