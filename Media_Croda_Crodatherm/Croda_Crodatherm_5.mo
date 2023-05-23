
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_5 "Croda International Plc, Crodatherm 5; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 5";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7314999999999998E+02, 2.8214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7114999999999998E+02, 2.7814999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.8000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.7000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.0394451724012539E+05
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
    constant Real[11] data_x =   {0.0000000000000000E+00, 2.6250000000000000E+00, 3.1250000000000000E+00, 3.3750000000000000E+00, 3.8750000000000000E+00, 4.1250000000000000E+00, 4.6250000000000000E+00, 5.3750000000000000E+00, 5.8750000000000000E+00, 7.1250000000000000E+00, 9.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.8878354005999999E-02, 5.1532356111000002E-02, 1.1483555258800000E-01, 4.4451097716100002E-01, 5.6898322103200005E-01, 5.4659787106299995E-01, 2.0547520895499999E-01, 5.0844390679999997E-02, 1.1432347339000001E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, -4.1607998999999997E-04, 1.2070395892800000E-01, 5.5868094216200004E-01, 6.6052057839199996E-01, 3.8914463932499999E-01, -3.6564234034100002E-01, -4.8800065333900000E-01, -1.0896027788500000E-01, -7.5057679350000002E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 2.5086402580999999E-02, 4.0207724058000001E-02, 5.8774124451999997E-02, 1.9687247071500000E-01, 3.2532926898400000E-01, 6.2076943812800001E-01, 9.0933346767699996E-01, 9.6567309676500002E-01, 9.9145741486700001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0027838143111678E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-2.0000000000000000E+00, 1.2500000000000000E-01, 2.6250000000000000E+00, 3.1250000000000000E+00, 3.3750000000000000E+00, 3.6250000000000000E+00, 4.1250000000000000E+00, 4.3750000000000000E+00, 4.6250000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 3.4161522679999998E-03, 7.7492575041000003E-02, 1.7980209288800000E-01, 3.1794449000500002E-01, 5.8925767234899995E-01, 7.5140216503199997E-01, 5.6068825123599997E-01, 2.8660378155799998E-01, 6.6424205595999999E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 4.7220073499999998E-03, 4.9623042637000001E-02, 3.3657763677800001E-01, 7.7259061492699999E-01, 9.2165105870999997E-01, -5.4011076746400000E-01, -1.0192449954489999E+00, -9.9102948611399999E-01, -7.7988844903499999E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.8628684939999999E-03, 8.0036994064000000E-02, 1.3870074950200001E-01, 1.9897522033500001E-01, 3.1221356469400002E-01, 6.7982656122100005E-01, 8.4724174777299999E-01, 9.5358330471200003E-01, 9.9684685797299999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0054555957278977E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.2400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.7000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.3000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4999999999999999E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 5</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
  Basic characteristics are the material class: paraffin-based, and encapsulation: none<br>  The data is taken from: Croda datasheet - last access 2023-02-28.<br><br>
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
end Croda_Crodatherm_5;