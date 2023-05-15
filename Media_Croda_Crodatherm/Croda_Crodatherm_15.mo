
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_15 "Croda International Plc, Crodatherm 15; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 15";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7614999999999998E+02, 2.9014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7614999999999998E+02, 2.8714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.9000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7579231659097533E+05
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
    constant Real[11] data_x =   {3.0000000000000000E+00, 4.6250000000000000E+00, 5.8750000000000000E+00, 8.8750000000000000E+00, 1.1125000000000000E+01, 1.2125000000000000E+01, 1.3125000000000000E+01, 1.3875000000000000E+01, 1.4875000000000000E+01, 1.6125000000000000E+01, 1.7000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.0760318410000000E-03, 7.4903060380000002E-03, 2.0048088395000000E-02, 5.9674702927000001E-02, 1.3593036339600001E-01, 3.7105614748299998E-01, 3.8200443424199998E-01, 6.6809916405000003E-02, 1.8973883010000000E-03, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -6.2318991209999999E-03, 8.0285911239999994E-03, 1.3002032283000000E-02, 4.0703089384000003E-02, 1.4717716555300001E-01, 1.8285434389800001E-01, -2.3783134136500000E-01, -1.4582821530199999E-01, -5.3373494650000002E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 3.8412594020000002E-03, 8.5523786909999996E-03, 4.5844861972000002E-02, 1.2325497110100000E-01, 2.1151012793300000E-01, 4.6013007591499999E-01, 7.5995585866000004E-01, 9.7505212002999997E-01, 9.9951414007499995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9241494463877411E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {3.0000000000000000E+00, 7.8750000000000000E+00, 9.6250000000000000E+00, 1.0625000000000000E+01, 1.1125000000000000E+01, 1.1375000000000000E+01, 1.1625000000000000E+01, 1.1875000000000000E+01, 1.2125000000000000E+01, 1.2375000000000000E+01, 1.2625000000000000E+01, 1.3125000000000000E+01, 1.4000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.6976500685999998E-02, 4.2286522599000002E-02, 1.8117090993099999E-01, 3.8761096345899998E-01, 5.3406851808000000E-01, 7.2402674877399997E-01, 6.8312405922599995E-01, 4.4889845706400000E-01, 1.3262876466099999E-01, 5.4727632999999999E-03, 6.2827826700000004E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 6.0784195110000000E-03, 1.3549278935999999E-02, 3.2205759189799998E-01, 5.0957707653700002E-01, 6.0605025112800004E-01, 5.9021120282799999E-01, -7.3414435803300004E-01, -1.1027965969819999E+00, -9.6137952701600005E-01, -3.1597818095999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.9717079880000000E-02, 8.0303920143999993E-02, 1.6742290475200000E-01, 3.0747898835900001E-01, 4.2365237094600000E-01, 5.8300756774300000E-01, 7.6813510070100000E-01, 9.1339086131100000E-01, 9.8626476231799998E-01, 9.9884356346799996E-01, 9.9972161551600003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0127795673308342E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.9600000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.5900000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.8999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 8.5900000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 15</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
  <br><br>
  The package contains phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  <p>
   Code export from <strong><u>slPCMlib database</u></strong> on 2023-05-15.<br><br>
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
    <li>file creation date: 2023-05-15 </ul>
    </p></html>"));
end Croda_Crodatherm_15;