
// within slPCMlib.PLUSS_HS;
package PLUSS_savE_HS24 "Pluss Advanced Technolgies Pvt Ltd, HS24; data taken from: PLUSS datasheet; last access: 2022-02-13."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS24";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9414999999999998E+02, 3.0514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0700000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.4200000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.9165966741499995E+05
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
    constant Real[10] data_x =   {2.1000000000000000E+01, 2.4625000000000000E+01, 2.5375000000000000E+01, 2.5875000000000000E+01, 2.6375000000000000E+01, 2.6875000000000000E+01, 2.8125000000000000E+01, 3.0125000000000000E+01, 3.1125000000000000E+01, 3.2000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 4.2860868305000001E-02, 1.0223385587700000E-01, 2.0667156100799999E-01, 2.3264419085499999E-01, 1.8889784781700000E-01, 2.0282478780700000E-01, 6.9072872566000004E-02, 3.6183549681000002E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 1.1042356789999999E-02, 1.8950672681700001E-01, 1.6788539908799999E-01, -7.8343324017999996E-02, -4.2701968463999998E-02, -1.8563336937000002E-02, -3.4399926848000001E-02, -5.4802384500000002E-02, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 6.5960160057000006E-02, 1.1226264369600000E-01, 1.9037380373700000E-01, 3.0597534522999997E-01, 4.1120348055700001E-01, 6.5423854921799995E-01, 9.3296501844599999E-01, 9.8759723386400000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0055919077528848E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.0000000000000000E+01, 2.1375000000000000E+01, 2.4875000000000000E+01, 2.5625000000000000E+01, 2.6375000000000000E+01, 2.6625000000000000E+01, 2.6875000000000000E+01, 2.7125000000000000E+01, 2.7625000000000000E+01, 2.7875000000000000E+01, 2.8375000000000000E+01, 2.9000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 5.4634358111000002E-02, 1.0778338806600000E-01, 1.1897276322000000E-01, 3.2768911786600002E-01, 5.5427175920200000E-01, 6.0035025005300002E-01, 4.3040565973799999E-01, 1.3516602159000000E-02, 1.1762405130000000E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, -5.0170206499999996E-03, 2.2445947561999999E-02, 1.3319352493000000E-02, 5.8639056028600001E-01, 6.6849475904099998E-01, -4.0145363931799999E-01, -9.1841427722299995E-01, -1.1224578111200000E-01, -5.0551063019999998E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 3.8064114737999998E-02, 2.9233984484600001E-01, 3.7716066986899999E-01, 5.1674208636700003E-01, 6.2573647775399999E-01, 7.7451335944400002E-01, 9.0506446955099995E-01, 9.9854392981200002E-01, 9.9981266925099999E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9250490737423913E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.6210000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5100000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.0500000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.5100000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS24</strong>  from manufacturer: <strong>Pluss Advanced Technolgies Pvt Ltd</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    multiple options available<br>  Data taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS24;