
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS30N "Pluss Advanced Technolgies Pvt Ltd, HS30N; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS30N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.3814999999999998E+02, 2.4914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.3714999999999998E+02, 2.4514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.1000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.7000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7323360421788169E+05
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
    constant Real[12] data_x =   {-3.5000000000000000E+01, -3.3375000000000000E+01, -3.1875000000000000E+01, -3.0875000000000000E+01, -3.0625000000000000E+01, -3.0375000000000000E+01, -2.9875000000000000E+01, -2.9625000000000000E+01, -2.9375000000000000E+01, -2.8375000000000000E+01, -2.5375000000000000E+01, -2.4000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 3.7419815751000000E-02, 1.8560093638000000E-02, 1.2169528179100000E-01, 2.2492136061500001E-01, 4.2151096347599998E-01, 5.8569117177300001E-01, 4.7938310791299998E-01, 2.9192840392199998E-01, 7.9148305634999999E-02, 1.4146744879000001E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 3.3730834110000003E-02, 1.1089160149000001E-02, 2.5915567950200002E-01, 5.8883099275200002E-01, 6.8935183226399999E-01, -2.1579400263699999E-01, -6.2479976543100002E-01, -5.4487752935599998E-01, -4.8740907586999997E-02, -1.3068989110000000E-03, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 2.2968245378000001E-02, 6.9172712573999998E-02, 1.1860061318800000E-01, 1.6018743216199999E-01, 2.4042316096499999E-01, 5.1092997664799999E-01, 6.4611907432600002E-01, 7.4206322183999995E-01, 8.8617645349300000E-01, 9.9048532603399997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9944238423216236E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-3.6000000000000000E+01, -3.3375000000000000E+01, -3.2125000000000000E+01, -3.1625000000000000E+01, -3.1375000000000000E+01, -3.1125000000000000E+01, -3.0875000000000000E+01, -3.0375000000000000E+01, -3.0125000000000000E+01, -2.9625000000000000E+01, -2.8875000000000000E+01, -2.8125000000000000E+01, -2.8000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 4.2721588118000003E-02, 2.1291369177800001E-01, 3.9241554834400000E-01, 5.1605475965799996E-01, 5.2535672019000001E-01, 4.1598381366400000E-01, 7.6599076712999994E-02, 3.1748780919999997E-02, 3.7760513482999997E-02, 1.6530213240700001E-01, 1.6443767445000000E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.6196404328999998E-02, 2.8064944710799999E-01, 4.1275037646500001E-01, 4.1764551498399999E-01, -2.6099504198599999E-01, -6.7833637509400002E-01, -5.1502844961200001E-01, -5.5805593850000001E-02, 1.9226989338799999E-01, -5.9487165240999999E-02, -1.9605131543099999E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 4.7046647971000000E-02, 1.7312122518299999E-01, 3.2257452630700001E-01, 4.3677496880599997E-01, 5.7127170882000000E-01, 6.9181714476500000E-01, 8.1226426097399995E-01, 8.2348148757799999E-01, 8.3576231543900004E-01, 9.2422873324499999E-01, 9.9922300051799995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0058762455882106E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4600000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.4250000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.0000000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4250000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS30N</strong>  from manufacturer: <strong>Pluss Advanced Technolgies Pvt Ltd</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: multiple options available<br>  The data is taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS30N;