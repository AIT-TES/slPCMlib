
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP15_gel "Rubitherm GmbH, SP15_gel; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP15_gel";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8314999999999998E+02, 2.9814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8014999999999998E+02, 2.9714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.4256237306750476E+05
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
    constant Integer  len_x =    17;
    constant Real[17] data_x =   {1.0000000000000000E+01, 1.1625000000000000E+01, 1.3625000000000000E+01, 1.4625000000000000E+01, 1.5125000000000000E+01, 1.5375000000000000E+01, 1.5625000000000000E+01, 1.5875000000000000E+01, 1.6375000000000000E+01, 1.6625000000000000E+01, 1.6875000000000000E+01, 1.7375000000000000E+01, 1.8375000000000000E+01, 1.9375000000000000E+01, 2.0375000000000000E+01, 2.1625000000000000E+01, 2.5000000000000000E+01};
    constant Real[17] data_y =   {0.0000000000000000E+00, 1.4970506587000000E-02, 4.4394386918000003E-02, 1.6674473812000001E-02, 3.5808103902000003E-02, 1.0496246219700001E-01, 3.4845137509099999E-01, 5.6323088629200002E-01, 4.2200508440399997E-01, 1.6725416005700000E-01, 7.3831724890999997E-02, 3.5079592207999999E-02, 4.8824856213999998E-02, 3.2387636442000003E-02, 4.6714647549000003E-02, 5.4714110320000002E-02, 0.0000000000000000E+00};
    constant Real[17] m_k =      {0.0000000000000000E+00, -3.6639572959999999E-03, 4.3254410733000001E-02, -5.3910073719000001E-02, 9.2573917360999999E-02, 7.2184682126999999E-01, 9.1785862733900003E-01, 7.6968063234899997E-01, -8.7606278531199999E-01, -6.5025998100100002E-01, -1.6380442022300001E-01, 3.9144761043000000E-02, -4.0866040748000002E-02, 2.7617535996999999E-02, -1.9555047079999999E-02, 3.5689961888999998E-02, 0.0000000000000000E+00};
    constant Real[17] iy_start = {0.0000000000000000E+00, 1.2928564290000000E-02, 5.6514989315999999E-02, 9.5023642386000001E-02, 1.0506052584200000E-01, 1.1933386086700000E-01, 1.7481275396300000E-01, 2.8918003944800003E-01, 5.6888328111200004E-01, 6.4113419676500005E-01, 6.6864855668099998E-01, 6.9157515700399996E-01, 7.4004037537000000E-01, 7.7482870472700005E-01, 8.1817265617599999E-01, 8.7419359877199998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9682078721056588E-01;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    16;
    constant Real[16] data_x =   {7.0000000000000000E+00, 8.8750000000000000E+00, 1.0375000000000000E+01, 1.1375000000000000E+01, 1.3625000000000000E+01, 1.4375000000000000E+01, 1.4625000000000000E+01, 1.5125000000000000E+01, 1.5625000000000000E+01, 1.5875000000000000E+01, 1.6625000000000000E+01, 1.7375000000000000E+01, 1.8625000000000000E+01, 2.0875000000000000E+01, 2.2375000000000000E+01, 2.4000000000000000E+01};
    constant Real[16] data_y =   {0.0000000000000000E+00, 2.1835558299999999E-02, 1.7032709965999999E-02, 3.6292905317000003E-02, 8.5021923175000000E-02, 2.5357527478899999E-01, 3.8724293549900002E-01, 3.4856364339599999E-01, 5.6783922644999998E-02, 1.9199706909999999E-02, 2.9318109492000002E-02, 5.1335494382000003E-02, 3.3295015293000003E-02, 7.1889732572000001E-02, 4.9873693060000002E-02, 0.0000000000000000E+00};
    constant Real[16] m_k =      {0.0000000000000000E+00, 1.0525474773000001E-02, 3.1130334156000000E-02, -1.9178995436000001E-02, 8.7772417074999998E-02, 3.8513331282300001E-01, 4.2098754096600000E-01, -5.8032131266099996E-01, -4.5624982527500002E-01, -4.6740014210000000E-02, 6.5975529076999995E-02, -3.7750285744000002E-02, 1.5891917802999998E-02, 1.9791024795999999E-02, -5.9498183473000003E-02, 0.0000000000000000E+00};
    constant Real[16] iy_start = {0.0000000000000000E+00, 1.7164007374000000E-02, 4.2127187024000000E-02, 7.2586360845999995E-02, 1.6277267187600000E-01, 2.7435683679599998E-01, 3.5324652485000002E-01, 5.5542966834200003E-01, 6.5291409182399995E-01, 6.6018463886300005E-01, 6.7292954960200002E-01, 7.0758613928799996E-01, 7.5290622310199995E-01, 8.6809624732099999E-01, 9.7292242867000001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.8716335374285369E-01;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.9999999999999998E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP15_gel</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: macroencapsulation<br>  The data is taken from: Rubitherm datasheet - last access 2020-07-12.<br><br>
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
end Rubitherm_SP15_gel;