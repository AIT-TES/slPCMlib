
within slPCMlib.Media_Axiotherm_ATS;
package Axiotherm_ATS_minus_6 "Axiotherm GmbH, ATS -6; data taken from: Axiotherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ATS -6";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6314999999999998E+02, 2.7014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.5914999999999998E+02, 2.6814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 3.0756552360903344E+05
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
    constant Real[11] data_x =   {-1.0000000000000000E+01, -7.3750000000000000E+00, -6.8750000000000000E+00, -6.6250000000000000E+00, -6.3750000000000000E+00, -6.1250000000000000E+00, -5.6250000000000000E+00, -5.3750000000000000E+00, -5.1250000000000000E+00, -4.6250000000000000E+00, -3.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.5637983720000002E-03, 2.3559359148999998E-02, 1.0917122926300001E-01, 5.3096274998000004E-01, 9.6445792889799997E-01, 7.4140795276000004E-01, 2.6496806357000002E-01, 9.9265544240999995E-02, 2.3748595353000002E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -1.2257648245999999E-02, 8.0243024767000001E-02, 1.3100547511609999E+00, 1.7105290677350000E+00, 1.6969689316879999E+00, -1.6363184602900001E+00, -1.1702060792180000E+00, -2.7683767689700001E-01, -1.0733803480000000E-02, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.1786968029000000E-02, 1.6670036935000000E-02, 2.6917741977999999E-02, 1.0532036121700000E-01, 2.9345037492499998E-01, 7.9236181833099995E-01, 9.1647784228200002E-01, 9.5760148550400004E-01, 9.8296377096499998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0060524811235363E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-1.4000000000000000E+01, -1.1125000000000000E+01, -8.3750000000000000E+00, -7.8750000000000000E+00, -7.6250000000000000E+00, -7.3750000000000000E+00, -7.1250000000000000E+00, -6.8750000000000000E+00, -6.6250000000000000E+00, -6.3750000000000000E+00, -6.1250000000000000E+00, -5.1250000000000000E+00, -5.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.5485380534000000E-02, 2.7304511633999998E-02, 1.0784287133900000E-01, 2.7494082503299999E-01, 7.3582871007999995E-01, 1.0202473120850000E+00, 8.9126092546199998E-01, 4.6329392912900003E-01, 8.5764500896000004E-02, 2.1797227294000000E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -1.0494778380000000E-03, 1.2669930103000000E-02, 2.9231540770900000E-01, 1.1115185382260000E+00, 1.5503090573260001E+00, 5.2391230133699995E-01, -1.5818251107040000E+00, -1.6133941758879999E+00, -8.8192060904000003E-01, -5.8339860955000003E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.3144284976999999E-02, 7.3686256429999999E-02, 1.0184322909900000E-01, 1.4573011967300001E-01, 2.7066091448899998E-01, 4.9709302017099999E-01, 7.4875138469799996E-01, 9.1942366010499998E-01, 9.8470077049600002E-01, 9.9392070756799999E-01, 1.0000000000000000E+00, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0070124490665773E+00;
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
  This package contains solid and liquid properties for the PCM:  <strong>ATS -6</strong>  from manufacturer: <strong>Axiotherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: macroencapsulation<br>  The data is taken from: Axiotherm datasheet - last access 2023-03-28.<br><br>
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
end Axiotherm_ATS_minus_6;