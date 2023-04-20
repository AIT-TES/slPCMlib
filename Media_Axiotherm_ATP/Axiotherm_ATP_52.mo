
// within slPCMlib.Axiotherm_ATP;
package Axiotherm_ATP_52 "Axiotherm GmbH, ATP 52; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 52";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1714999999999998E+02, 3.2714999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1614999999999998E+02, 3.2614999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0000000000000000E+05
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
    constant Real[12] data_x =   {4.4000000000000000E+01, 4.9125000000000000E+01, 5.0625000000000000E+01, 5.1375000000000000E+01, 5.1625000000000000E+01, 5.1875000000000000E+01, 5.2125000000000000E+01, 5.2375000000000000E+01, 5.2625000000000000E+01, 5.2875000000000000E+01, 5.3875000000000000E+01, 5.4000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 3.7011306331999998E-02, 1.2263743454000001E-01, 3.7784793225000002E-01, 5.8660472041800005E-01, 6.7260610928900000E-01, 5.9313182020999999E-01, 3.8765735008399999E-01, 1.5813916475900000E-01, 6.3887235737000003E-02, 5.1238146200000001E-04, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.1343316362000000E-02, 1.2657196278800001E-01, 6.0706056786399998E-01, 6.7486539656900002E-01, 1.7371147643000000E-02, -7.5795693139499998E-01, -8.6827984213300002E-01, -7.5464929081900001E-01, -1.8686762403700000E-01, -5.2276371379999997E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 4.8647379498999999E-02, 1.4973872283799999E-01, 3.1668995994400001E-01, 4.3819770893499999E-01, 6.0076861943399995E-01, 7.6478464780199995E-01, 8.8929443533999997E-01, 9.5766105119599998E-01, 9.8272621556299999E-01, 9.9997450934999998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0108509178289009E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {4.3000000000000000E+01, 4.7625000000000000E+01, 4.9375000000000000E+01, 4.9875000000000000E+01, 5.0125000000000000E+01, 5.0375000000000000E+01, 5.0625000000000000E+01, 5.0875000000000000E+01, 5.1125000000000000E+01, 5.1375000000000000E+01, 5.1625000000000000E+01, 5.1875000000000000E+01, 5.2875000000000000E+01, 5.3000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 1.9121799163999999E-02, 3.5629263033000000E-02, 3.8671186065000002E-02, 9.0188953254999996E-02, 2.6055436189400000E-01, 8.1146929268099999E-01, 1.1367568858849999E+00, 8.7345698667600002E-01, 3.4663366644100002E-01, 1.6347535637000001E-02, 2.1582414609999999E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 9.0341505159999991E-03, -2.4845220750999999E-02, 7.2990969756000001E-02, 2.5054295596600001E-01, 1.2066366506190001E+00, 1.8492727821449999E+00, 2.1517492242000000E-01, -1.8434261336430000E+00, -1.7498720084179999E+00, -1.3338294089300001E-01, -6.2817041820000001E-03, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 2.8366937855000000E-02, 8.5426546704000003E-02, 1.0211140343800000E-01, 1.1743005076900000E-01, 1.5664112151199999E-01, 2.8846633440899999E-01, 5.4276118555300001E-01, 8.0710458547999997E-01, 9.6048921215600003E-01, 9.9777336735699995E-01, 9.9943938180900005E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0089496097545447E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.8900000000000000E+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATP 52</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
       material class: unknown;  encapsulation:    macroencapsulation<br>  Data taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATP_52;