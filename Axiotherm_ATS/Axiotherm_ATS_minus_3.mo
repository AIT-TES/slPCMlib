
// within slPCMlib.Axiotherm_ATS;
package Axiotherm_ATS_minus_3 "Axiotherm GmbH, ATS -3; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -3";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6614999999999998E+02, 2.7314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6214999999999998E+02, 2.7114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.7756900306783005E+05
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
    constant Real[11] data_x =   {-7.0000000000000000E+00, -4.3750000000000000E+00, -3.8750000000000000E+00, -3.6250000000000000E+00, -3.3750000000000000E+00, -3.1250000000000000E+00, -2.6250000000000000E+00, -2.3750000000000000E+00, -2.1250000000000000E+00, -1.6250000000000000E+00, 0.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.0072765009999999E-03, 2.3790534400999999E-02, 1.1072370337900000E-01, 5.3615997488500000E-01, 9.7086346357999997E-01, 7.3775780888499998E-01, 2.5905639232500000E-01, 9.5754864078000004E-02, 2.3532882073999999E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -7.4284887850000004E-03, 8.2099827525999994E-02, 1.3286844405470000E+00, 1.7201696564419999E+00, 1.7090563327770001E+00, -1.6506660045970001E+00, -1.1747558128219999E+00, -2.6731824807299998E-01, -6.3369889339999997E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 8.2819106870000003E-03, 1.3156972729000000E-02, 2.3565698974000002E-02, 1.0305212086000000E-01, 2.9307760013000000E-01, 7.9442114791100005E-01, 9.1757446936300002E-01, 9.5753392187800002E-01, 9.8212446223500005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0084361223752527E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {-1.1000000000000000E+01, -9.8750000000000000E+00, -8.3750000000000000E+00, -5.6250000000000000E+00, -5.1250000000000000E+00, -4.8750000000000000E+00, -4.6250000000000000E+00, -4.3750000000000000E+00, -4.1250000000000000E+00, -3.8750000000000000E+00, -3.6250000000000000E+00, -3.3750000000000000E+00, -3.1250000000000000E+00, -2.1250000000000000E+00, -2.0000000000000000E+00};
    constant Real[15] data_y =   {0.0000000000000000E+00, 3.8155658570000000E-03, 1.5753424296999999E-02, 2.3298442538000001E-02, 4.6018958190999999E-02, 1.0581442335500001E-01, 2.7338025420599998E-01, 7.4070921071100004E-01, 1.0311698975500001E+00, 9.0126103215999998E-01, 4.6703757657799999E-01, 8.5057757306000001E-02, 2.1304922510000001E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, 8.9037073369999994E-03, -1.0411524605000000E-02, 8.3124499849999995E-03, 1.0153376636500000E-01, 2.8979525263700001E-01, 1.1242506274240001E+00, 1.5748223973919999E+00, 5.4378220093399998E-01, -1.6023355733460001E+00, -1.6347416036119999E+00, -8.7904678626900001E-01, -5.7580137725000001E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 1.2157090130000001E-03, 1.9643142698000000E-02, 6.1835010890000000E-02, 7.7330799073000001E-02, 9.5456414017000002E-02, 1.3881334440099999E-01, 2.6410547710499999E-01, 4.9256070679199998E-01, 7.4707515599399998E-01, 9.1948903117699998E-01, 9.8502411388699995E-01, 9.9410458588899997E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0070544661143290E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.1700000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.1000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ATS -3</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATS_minus_3;