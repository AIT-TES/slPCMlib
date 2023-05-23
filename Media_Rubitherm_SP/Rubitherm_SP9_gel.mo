
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP9_gel "Rubitherm GmbH, SP9_gel; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP9_gel";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.7314999999999998E+02, 2.9214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7314999999999998E+02, 2.9214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.3600000000000000E+05
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
    constant Real[13] data_x =   {0.0000000000000000E+00, 3.6250000000000000E+00, 5.6250000000000000E+00, 7.6250000000000000E+00, 1.0375000000000000E+01, 1.0625000000000000E+01, 1.0875000000000000E+01, 1.1625000000000000E+01, 1.1875000000000000E+01, 1.3125000000000000E+01, 1.4875000000000000E+01, 1.6875000000000000E+01, 1.9000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 3.0272280632000000E-02, 4.0757446847999999E-02, 3.9530342865999997E-02, 2.3785554978000001E-01, 3.0648462780899999E-01, 3.1678285773900000E-01, 7.3477965706000001E-02, 3.5098019345000002E-02, 3.1912485885000001E-02, 2.7149197859000000E-02, 5.5388095780999999E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -7.3101328200000004E-04, 2.3065434635000000E-02, 1.2341203822000000E-02, 2.2363514506400001E-01, 2.2858177856199999E-01, -1.0486095061800001E-01, -3.5756536028100000E-01, -7.2344037801000002E-02, 1.6307124669000000E-02, -3.1612942000000001E-04, -1.7782823994999999E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 5.5915331967000000E-02, 1.1929210598800000E-01, 2.0352571648100001E-01, 4.5287056360500000E-01, 5.2118828445499998E-01, 6.0118581686499994E-01, 7.6007912291699997E-01, 7.7221907386300004E-01, 8.0269176188199998E-01, 8.5886057048099995E-01, 9.4761107122900001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0044248266681155E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    18;
    constant Real[18] data_x =   {0.0000000000000000E+00, 2.6250000000000000E+00, 4.8750000000000000E+00, 6.3750000000000000E+00, 7.3750000000000000E+00, 7.8750000000000000E+00, 8.3750000000000000E+00, 8.6250000000000000E+00, 9.1250000000000000E+00, 9.6250000000000000E+00, 9.8750000000000000E+00, 1.0875000000000000E+01, 1.2625000000000000E+01, 1.5375000000000000E+01, 1.6375000000000000E+01, 1.7125000000000000E+01, 1.8125000000000000E+01, 1.9000000000000000E+01};
    constant Real[18] data_y =   {0.0000000000000000E+00, 3.2692935537999998E-02, 3.9508411435000002E-02, 4.5563123584999997E-02, 7.7421754125000006E-02, 8.9686238582999997E-02, 2.0612503506200000E-01, 3.5945843872099997E-01, 3.4830067415499999E-01, 4.6660107803999998E-02, 1.3104947980000000E-02, 3.1502650567000003E-02, 2.5963464671999999E-02, 3.6838772196999998E-02, 3.9414116694000002E-02, 5.9525734917000001E-02, 9.8310967509999993E-03, 0.0000000000000000E+00};
    constant Real[18] m_k =      {0.0000000000000000E+00, -4.1770253820000002E-03, 4.9544988950000002E-03, 4.2506450695000003E-02, 5.3880373949999998E-03, 8.6303649982999997E-02, 4.0393993183100002E-01, 4.7777595932500000E-01, -5.9819119417400002E-01, -4.2896221397500001E-01, -3.6202489371999998E-02, 4.9747465232000000E-02, -8.6317731040000002E-03, -2.0573293430000000E-02, 4.3706930937000003E-02, -3.7236694787999999E-02, -2.3360207256000000E-02, 0.0000000000000000E+00};
    constant Real[18] iy_start = {0.0000000000000000E+00, 4.5039266662999997E-02, 1.2195448720000000E-01, 1.7838046843200001E-01, 2.4258302851200000E-01, 2.8243648958099998E-01, 3.4937249581199997E-01, 4.1926881434000002E-01, 6.1744212660200004E-01, 7.1209196749299997E-01, 7.1748479828900003E-01, 7.3253629301300005E-01, 7.9733139652899998E-01, 8.9065329051300002E-01, 9.2322868076099995E-01, 9.6388278516299997E-01, 9.9720599615200001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9406865235005637E-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP9_gel</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP9_gel;