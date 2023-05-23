
// within slPCMlib.Rubitherm_RT;
package Rubitherm_RT35HC "Rubitherm GmbH, RT35HC; data taken from: Rubitherm datasheet; last access: 2020-10-09."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT35HC";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0214999999999998E+02, 3.1214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0114999999999998E+02, 3.1014999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.1547052462262398E+05
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
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.9000000000000000E+01, 3.0875000000000000E+01, 3.3625000000000000E+01, 3.4125000000000000E+01, 3.4375000000000000E+01, 3.4625000000000000E+01, 3.5125000000000000E+01, 3.5375000000000000E+01, 3.5875000000000000E+01, 3.6875000000000000E+01, 3.8125000000000000E+01, 3.9000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 8.5990677689999993E-03, 5.8664880296999997E-02, 1.2887392311099999E-01, 2.1508194567500000E-01, 3.7369194780400000E-01, 5.5810456663300001E-01, 5.2358236868600005E-01, 3.0677999948099999E-01, 4.5543229647999998E-02, 2.7058323910000002E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.6647167500000002E-04, 4.6653446451000002E-02, 2.2300136774900001E-01, 5.0947629664000005E-01, 5.7382353071000003E-01, 5.1547785766999997E-02, -4.2137455796300000E-01, -3.8999389853700001E-01, -9.2169118905000000E-02, -6.6978311910000001E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 8.0741543650000008E-03, 7.2046427018000000E-02, 1.1574756227800000E-01, 1.5772095219000001E-01, 2.3181390848299999E-01, 4.7841072615000002E-01, 6.1764703396200005E-01, 8.2693214996499997E-01, 9.7999244894500004E-01, 9.9923494965000004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0113478481285403E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {2.8000000000000000E+01, 3.0125000000000000E+01, 3.1625000000000000E+01, 3.2875000000000000E+01, 3.4375000000000000E+01, 3.4625000000000000E+01, 3.4875000000000000E+01, 3.5125000000000000E+01, 3.5625000000000000E+01, 3.5875000000000000E+01, 3.6875000000000000E+01, 3.7000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 4.7713430079999999E-03, 1.3166066707000001E-02, 9.3312958557000003E-02, 4.4825031706000001E-01, 5.3792118969100000E-01, 5.3324626059299995E-01, 4.3590924170399997E-01, 1.2420481717300000E-01, 5.4195911526999999E-02, 5.3697604100000004E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 4.8448336830000001E-03, -6.1860905159999997E-03, 1.5331410561500000E-01, 3.2224115205699999E-01, 3.2189919208200002E-01, -2.4884562757500001E-01, -6.1526281022200002E-01, -6.1139001999099996E-01, -1.5081460721800000E-01, -5.5349425900000004E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 3.2549672520000000E-03, 1.8817135387000000E-02, 6.4718650483999995E-02, 4.4020200980899998E-01, 5.6379936942900000E-01, 7.0102781842899997E-01, 8.2440423985599998E-01, 9.6472005666700000E-01, 9.8467364676400004E-01, 9.9997357665700004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0026294456836473E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT35HC</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT35HC;