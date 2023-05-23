
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP_minus_21 "Rubitherm GmbH, SP-21; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP-21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4714999999999998E+02, 2.5814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4614999999999998E+02, 2.5214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.7800000000000000E+05
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
    constant Real[17] data_x =   {-2.6000000000000000E+01, -2.4375000000000000E+01, -2.2875000000000000E+01, -2.2125000000000000E+01, -2.1875000000000000E+01, -2.1625000000000000E+01, -2.1375000000000000E+01, -2.1125000000000000E+01, -2.0875000000000000E+01, -2.0625000000000000E+01, -2.0375000000000000E+01, -2.0125000000000000E+01, -1.9625000000000000E+01, -1.9125000000000000E+01, -1.8375000000000000E+01, -1.7625000000000000E+01, -1.5000000000000000E+01};
    constant Real[17] data_y =   {0.0000000000000000E+00, 4.0122033767999998E-02, 1.9565742140000000E-03, 2.4965489190000001E-03, 1.2511082466000000E-02, 7.0464136120000007E-02, 4.0546330613900000E-01, 7.7469555642599996E-01, 8.7427281411799995E-01, 6.2182731612300002E-01, 2.2599269773399999E-01, 9.4751077584999996E-02, 5.5640378753000000E-02, 1.1854484439300000E-01, 3.4839288918999997E-02, 9.0102649720000001E-03, 0.0000000000000000E+00};
    constant Real[17] m_k =      {0.0000000000000000E+00, -2.0707092097000002E-02, -3.8736115539999998E-03, 7.8579320149999993E-03, 4.5686102632000000E-02, 8.4556963343799996E-01, 1.4093884686350000E+00, 1.3813330758889999E+00, -5.2096027333600003E-01, -1.3372224310220000E+00, -9.3796781808700003E-01, -2.0596683166200000E-01, 1.1165609876099999E-01, 3.8979709760000002E-02, -1.2460886014899999E-01, 5.3059257289999997E-03, 0.0000000000000000E+00};
    constant Real[17] iy_start = {0.0000000000000000E+00, 3.7428654733000000E-02, 6.6039915822000000E-02, 6.7168145949999994E-02, 6.8859407841000003E-02, 7.5110824379000005E-02, 1.3208052047499999E-01, 2.8083092495000000E-01, 4.9844621487099999E-01, 6.9111469428600003E-01, 7.9577574817399999E-01, 8.3232265031200003E-01, 8.6353088525599997E-01, 9.0892219561999998E-01, 9.7458818002500003E-01, 9.8501804325900000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0073437688181723E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {-2.7000000000000000E+01, -2.5125000000000000E+01, -2.4375000000000000E+01, -2.3625000000000000E+01, -2.3375000000000000E+01, -2.2875000000000000E+01, -2.2625000000000000E+01, -2.2375000000000000E+01, -2.1875000000000000E+01, -2.1625000000000000E+01, -2.1375000000000000E+01, -2.1125000000000000E+01, -2.1000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.4708527489000001E-02, 7.9887319016999997E-02, 9.6363790370999999E-02, 7.5905273081999999E-02, 1.3037959691500001E-01, 2.4673993812000000E-01, 5.2258355880899998E-01, 7.4247155761299999E-01, 5.6841278878699997E-01, 2.9517729687400002E-01, 6.8978563569000007E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 3.5076926118999999E-02, 9.2828250003999996E-02, -4.9313567351000002E-02, -5.2376831858000003E-02, 2.3441337560200001E-01, 7.7643423436799996E-01, 9.6433725509500001E-01, -4.2503250837800000E-01, -1.0150986390570000E+00, -9.9337963344799995E-01, -8.1504928326799997E-01, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 1.2931381113000000E-02, 4.9571209404999997E-02, 1.2257429938099999E-01, 1.4419675622899999E-01, 1.8994736094600001E-01, 2.3441413411100001E-01, 3.2992278146100001E-01, 6.7629907987399995E-01, 8.4379736692999996E-01, 9.5199765145199999E-01, 9.9673911241199997E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0033814760348174E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.3000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.2000000000000000E+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>SP-21</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP_minus_21;