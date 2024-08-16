within slPCMlib.Media_Knauf_SmartBoard;
package SmartBoard_26 "Knauf Gips KG, SmartBoard 26; data taken from: DBU-Abschlussbericht-AZ-23836.pdf."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SmartBoard 26";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9414999999999998E+02, 3.0289999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0264999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.2000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.2000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.8103435654108758E+04
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
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {2.1000000000000000E+01, 2.4875000000000000E+01, 2.6375000000000000E+01, 2.7375000000000000E+01, 2.8125000000000000E+01, 2.8375000000000000E+01, 2.9375000000000000E+01, 2.9750000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 3.9524494463000000E-02, 1.0281523645500000E-01, 2.4022445100400000E-01, 4.4356390479600000E-01, 4.7611561813300002E-01, 5.6959529910000001E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 1.9764773778000001E-02, 7.6598784335999995E-02, 1.7530288815200001E-01, 2.7412477235700000E-01, 3.9902207219000001E-02, -4.1577317347100001E-01, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 5.2395786946999999E-02, 1.4951137468600001E-01, 3.1453428882099999E-01, 5.6898773381000001E-01, 6.8639730445100000E-01, 9.9413095902100002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0105846331077069E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {2.0000000000000000E+01, 2.4125000000000000E+01, 2.6375000000000000E+01, 2.7125000000000000E+01, 2.7875000000000000E+01, 2.8375000000000000E+01, 2.9375000000000000E+01, 2.9500000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 3.6361261967000000E-02, 1.6546862689700001E-01, 4.1897048432799999E-01, 3.9885247544500002E-01, 1.7009147588800000E-01, 1.0442696199000000E-02, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 2.1911077574000001E-02, 1.2753961943799999E-01, 4.0240379798400000E-01, -4.6209066945799998E-01, -3.1795292637599998E-01, -9.4118475100000001E-02, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 4.4311964735999999E-02, 2.2841259688599999E-01, 4.3650610071800000E-01, 7.8676465249399996E-01, 9.2722155601599998E-01, 9.9946522212500000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0087894701505709E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 7.6700000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.6700000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 1.7999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 1.7999999999999999E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SmartBoard 26</strong>  from manufacturer: <strong>Knauf Gips KG</strong>.<br>
  Basic characteristics are the material class: paraffin-based composite, and encapsulation: microencapsulated<br>  The data is taken from: DBU-Abschlussbericht-AZ-23836.pdf - last access 2010-10-01.<br><br>
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
end SmartBoard_26;
