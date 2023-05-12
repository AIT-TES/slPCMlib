
// within slPCMlib.Axiotherm_ATP;
package Axiotherm_ATP_18 "Axiotherm GmbH, ATP 18; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 18";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8614999999999998E+02, 2.9314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8414999999999998E+02, 2.9414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.3700000000000000E+05
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
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {1.3000000000000000E+01, 1.5125000000000000E+01, 1.6375000000000000E+01, 1.6875000000000000E+01, 1.7125000000000000E+01, 1.7375000000000000E+01, 1.7625000000000000E+01, 1.7875000000000000E+01, 1.8125000000000000E+01, 1.8375000000000000E+01, 1.8625000000000000E+01, 1.8875000000000000E+01, 2.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 4.0403288440000002E-03, 1.7231036838000001E-02, 1.9843062277000001E-02, 5.6598649883999999E-02, 1.9839049596600000E-01, 7.5225300666899997E-01, 1.1651232441999999E+00, 1.0811793815660000E+00, 5.1199289077099996E-01, 6.8161205969999999E-02, 1.4981489863000000E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 5.6165215970000000E-03, -1.1728527424000001E-02, 4.3709026689000001E-02, 1.7464664109000000E-01, 1.2708639930020000E+00, 1.9597574315460000E+00, 1.2040644987110001E+00, -1.8856785542450001E+00, -1.9989614181840001E+00, -6.2739933915699997E-01, -3.4595973725000000E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.1929149599999999E-03, 1.7842863021000001E-02, 2.6006975344999999E-02, 3.4935490937999997E-02, 6.1262636191000003E-02, 1.7722289176900000E-01, 4.2234817133899999E-01, 7.2107757314900001E-01, 9.2205821032199997E-01, 9.8784112205800001E-01, 9.9519194375599995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0062286457863401E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {1.1000000000000000E+01, 1.3375000000000000E+01, 1.3625000000000000E+01, 1.4375000000000000E+01, 1.6125000000000000E+01, 1.6875000000000000E+01, 1.7125000000000000E+01, 1.7375000000000000E+01, 1.7625000000000000E+01, 1.7875000000000000E+01, 1.8125000000000000E+01, 1.8625000000000000E+01, 1.8875000000000000E+01, 1.9375000000000000E+01, 2.0875000000000000E+01, 2.1000000000000000E+01};
    constant Real[16] data_y =   {0.0000000000000000E+00, 7.6629087499999996E-04, 0.0000000000000000E+00, 0.0000000000000000E+00, 1.7179004851000001E-02, 2.9816346620999998E-02, 8.2102377628000003E-02, 2.6030499703900001E-01, 8.6147035888500001E-01, 1.2169485505080000E+00, 9.7812075266700005E-01, 1.6877637131000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 7.4742157400000005E-04, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, -8.1690571679999997E-03, 0.0000000000000000E+00, 0.0000000000000000E+00, -7.5531552899999995E-04, 7.1425170752000000E-02, 2.5109676616600002E-01, 1.2860880542369999E+00, 1.9960515450450000E+00, 4.4602107677600000E-01, -1.9827078705240000E+00, -1.6720050018900001E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, -7.6963429410000003E-03, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 4.7960370689999997E-03, 4.8497939170000001E-03, 4.8497939170000001E-03, 2.0222217537000001E-02, 3.4600472069999999E-02, 4.7781451689999997E-02, 8.5555546437000002E-02, 2.2340720713600001E-01, 4.9388728341900001E-01, 7.8371149134899998E-01, 9.9668892363299999E-01, 9.9793983824599997E-01, 9.9793983824599997E-01, 9.9996295065800000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0097233149222435E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 8.4400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.6000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.6000000000000000E+02;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATP 18</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_18;