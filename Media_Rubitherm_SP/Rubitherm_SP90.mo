
// within slPCMlib.Rubitherm_SP;
package Rubitherm_SP90 "Rubitherm GmbH, SP90; data taken from: Rubitherm datasheet; last access: 2020-11-24."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP90";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.5214999999999998E+02, 3.6514999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.5214999999999998E+02, 3.6314999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.2716665211715113E+05
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
    constant Real[12] data_x =   {7.9000000000000000E+01, 8.0875000000000000E+01, 8.2625000000000000E+01, 8.6875000000000000E+01, 8.8625000000000000E+01, 8.9375000000000000E+01, 8.9625000000000000E+01, 8.9875000000000000E+01, 9.0125000000000000E+01, 9.0625000000000000E+01, 9.0875000000000000E+01, 9.2000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 6.4933552380000001E-03, 2.5672854554999999E-02, 6.0277313954000003E-02, 1.7099044582199999E-01, 3.6202353184500002E-01, 4.8858246503500002E-01, 4.9052401662000000E-01, 3.7632210653999998E-01, 3.7296650404999999E-02, 1.0834334076999999E-02, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 1.8209466290000000E-03, -6.3311246849999999E-03, 2.0124952912000001E-02, 1.2794802829700000E-01, 3.7516575430799998E-01, 3.7975203116200001E-01, -2.7136361266400000E-01, -6.8159828551099999E-01, -2.3544998994300001E-01, -2.3064304465000001E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 5.6272409850000003E-03, 3.6251521286000003E-02, 1.8095608301900001E-01, 3.5810240005900001E-01, 5.4887595271199996E-01, 6.5657884835500002E-01, 7.8401612981799995E-01, 8.9596479735300005E-01, 9.9131507663399998E-01, 9.9628998946799996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0131797586239597E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {7.9000000000000000E+01, 8.0625000000000000E+01, 8.2625000000000000E+01, 8.5125000000000000E+01, 8.6375000000000000E+01, 8.6875000000000000E+01, 8.7375000000000000E+01, 8.7625000000000000E+01, 8.8125000000000000E+01, 8.8375000000000000E+01, 8.8625000000000000E+01, 8.9375000000000000E+01, 9.0000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.3235387754000001E-02, 1.6425752650000000E-02, 3.7176613571000001E-02, 1.0604331535200000E-01, 1.2591225489300001E-01, 2.7046199707800000E-01, 4.4984085629699999E-01, 5.4756545274799995E-01, 4.2318542440899998E-01, 2.3920169632100000E-01, 3.2008604635999999E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -6.7124455940000003E-03, 1.3313783871999999E-02, 3.9303238749999997E-02, 2.0042928924000000E-02, 1.1598048073099999E-01, 5.1615065857700004E-01, 5.9965040910900003E-01, -3.4496819551000002E-01, -6.4050232482000002E-01, -5.6239526668399997E-01, -9.8429086444000002E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.0432700205000001E-02, 5.3542979017000003E-02, 1.0721164340600001E-01, 1.9957940407999999E-01, 2.5578100661900000E-01, 3.4688037194200000E-01, 4.3682165886600000E-01, 7.0686860390799999E-01, 8.3021568845399996E-01, 9.1291836685500005E-01, 9.9317571251100001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0037758214502601E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.7000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.6500000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.6500000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP90</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    multiple options available<br>  Data taken from: Rubitherm datasheet - last access 2020-11-24.<br><br>
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
end Rubitherm_SP90;