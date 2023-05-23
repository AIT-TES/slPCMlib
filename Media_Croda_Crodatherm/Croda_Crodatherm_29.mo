
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_29 "Croda International Plc, Crodatherm 29; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 29";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9514999999999998E+02, 3.0414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.3000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9626543688347557E+05
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
    constant Real[11] data_x =   {2.2000000000000000E+01, 2.4875000000000000E+01, 2.6625000000000000E+01, 2.7125000000000000E+01, 2.7375000000000000E+01, 2.7875000000000000E+01, 2.8375000000000000E+01, 2.9375000000000000E+01, 2.9875000000000000E+01, 3.0875000000000000E+01, 3.1000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 5.5584935650000002E-03, 1.5317296344000000E-02, 4.5908811405999998E-02, 1.0836182652900000E-01, 4.5465516346000001E-01, 6.3957937902999995E-01, 2.0148850734599999E-01, 3.6569108401999999E-02, 3.3454414200000002E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 4.8599537159999999E-03, -9.6893404699999999E-04, 1.1356508755000000E-01, 5.4225275330800005E-01, 6.9783573875799998E-01, -8.0171229534999999E-02, -5.5000216485800002E-01, -1.0384485132800000E-01, -3.4338493500000000E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 4.6351176949999996E-03, 2.4356378293000001E-02, 3.7255433402000000E-02, 5.4278344024000003E-02, 1.9156408799999999E-01, 4.8085246480499999E-01, 9.3977951709600005E-01, 9.8991600716600003E-01, 9.9998358930700004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9834785488300670E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.0000000000000000E+01, 2.3625000000000000E+01, 2.5625000000000000E+01, 2.6125000000000000E+01, 2.6375000000000000E+01, 2.6625000000000000E+01, 2.7125000000000000E+01, 2.7375000000000000E+01, 2.7625000000000000E+01, 2.7875000000000000E+01, 2.8625000000000000E+01, 3.0000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 2.9985493550000000E-03, 2.2440489691000001E-02, 8.4433422738999994E-02, 2.3209832109100001E-01, 6.7792517661200002E-01, 9.7295895560400003E-01, 5.9613109662899999E-01, 1.6758918897600000E-01, 4.7876215063999999E-02, 3.0481911200000000E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -7.6289301999999998E-04, -3.7038648240000000E-03, 2.3422808692899999E-01, 1.1458316226759999E+00, 1.5790225197970000E+00, -1.3397697510340001E+00, -1.6157375460990000E+00, -1.2512890854609999E+00, -1.5872834908899999E-01, -5.5562073260000001E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 6.2714332790000001E-03, 3.2695662689000003E-02, 5.4461233600999998E-02, 8.9286180049000005E-02, 2.0080345229900001E-01, 6.7441988143499998E-01, 8.7202986315800002E-01, 9.6561396474899996E-01, 9.8686063265299995E-01, 9.9877953587699997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0001841993523972E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.1700000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.5100000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.2000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4999999999999999E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 29</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
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
end Croda_Crodatherm_29;