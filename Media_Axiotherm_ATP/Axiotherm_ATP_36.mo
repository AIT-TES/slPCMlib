
// within slPCMlib.Axiotherm_ATP;
package Axiotherm_ATP_36 "Axiotherm GmbH, ATP 36; data taken from: Axiotherm datasheet; last access: 2023-03-28."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATP 36";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.0414999999999998E+02, 3.1314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.0314999999999998E+02, 3.1214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0300000000000000E+05
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
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {3.1000000000000000E+01, 3.2125000000000000E+01, 3.3125000000000000E+01, 3.3375000000000000E+01, 3.4375000000000000E+01, 3.4625000000000000E+01, 3.5125000000000000E+01, 3.5375000000000000E+01, 3.5625000000000000E+01, 3.6125000000000000E+01, 3.6375000000000000E+01, 3.6625000000000000E+01, 3.6875000000000000E+01, 3.7625000000000000E+01, 4.0000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 6.7916420849999997E-03, 0.0000000000000000E+00, 0.0000000000000000E+00, 3.3513539140999998E-02, 1.7333005068999999E-02, 5.3579074268000002E-02, 1.8370464850599999E-01, 6.9116378045899995E-01, 1.0434752382750001E+00, 5.7644290878700000E-01, 1.0943527406799999E-01, 2.1443523291000001E-02, 3.3986984190000000E-03, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, -1.9135062772999999E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, -3.4655981671000001E-02, -3.8845097752999998E-02, 1.5673453143200000E-01, 1.2635965255449999E+00, 1.8482733345639999E+00, -1.8680408748830000E+00, -1.8680778257580000E+00, -1.3132232888169999E+00, -7.6084737476000003E-02, 4.8827161220000000E-03, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 5.8636030730000004E-03, 7.6725956329999996E-03, 7.6725956329999996E-03, 2.7401997397000000E-02, 3.3807109911000000E-02, 4.7519375996999999E-02, 7.1517882130999993E-02, 1.7828926397200001E-01, 6.9157408239600004E-01, 8.9493641182899997E-01, 9.7813823027900004E-01, 9.8809737198000003E-01, 9.9364163911799996E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0043082066593070E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {3.0000000000000000E+01, 3.2125000000000000E+01, 3.4125000000000000E+01, 3.4875000000000000E+01, 3.5375000000000000E+01, 3.5625000000000000E+01, 3.5875000000000000E+01, 3.6125000000000000E+01, 3.6625000000000000E+01, 3.6875000000000000E+01, 3.7375000000000000E+01, 3.8875000000000000E+01, 3.9000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 4.7321998749999998E-03, 3.0783071879000000E-02, 9.1378968071999994E-02, 3.9154665672500000E-01, 8.6951534974199995E-01, 1.0529545544479999E+00, 7.9172922301899995E-01, 1.9417475728000001E-02, 0.0000000000000000E+00, 0.0000000000000000E+00, 8.5989763600000004E-04, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 3.8921708530000000E-03, 2.8004916325000000E-02, 1.8892205719800001E-01, 1.1452432597640001E+00, 1.4681251951609999E+00, -2.9321987181499998E-01, -1.6574632146410000E+00, -1.8451086316300000E-01, 0.0000000000000000E+00, 0.0000000000000000E+00, -8.8545304710000005E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 3.6023985510000001E-03, 3.1381352559000000E-02, 7.0068692606999994E-02, 1.7198199503600001E-01, 3.2964290921799999E-01, 5.8186063268499999E-01, 8.2215757067899997E-01, 9.9614463605100001E-01, 9.9762690168000001E-01, 9.9762690168000001E-01, 9.9995732289899997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0109639521945237E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATP 36</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
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
end Axiotherm_ATP_36;