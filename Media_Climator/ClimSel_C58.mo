
// within slPCMlib.Media_Climator;
package ClimSel_C58 "Climator Sweden AB, ClimSel C58; data taken from: Climator Sweden AB datasheet; last access: 2022-10-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C58";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.2614999999999998E+02, 3.3414999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.2414999999999998E+02, 3.2914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {3.2000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.2500000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.1007872024036245E+05
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
    constant Real[12] data_x =   {5.3000000000000000E+01, 5.4875000000000000E+01, 5.6625000000000000E+01, 5.7375000000000000E+01, 5.7625000000000000E+01, 5.7875000000000000E+01, 5.8125000000000000E+01, 5.8375000000000000E+01, 5.8625000000000000E+01, 5.8875000000000000E+01, 5.9875000000000000E+01, 6.1000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 2.3535997850000000E-02, 1.4187980821599999E-01, 3.7438622930600002E-01, 5.8703235478200000E-01, 6.7424064684700002E-01, 5.9166170153599995E-01, 3.8128893781899997E-01, 1.4997320109800000E-01, 6.0624504147999997E-02, 5.7747994770000000E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.4919669411000000E-02, 6.9157430827999994E-02, 6.0202645295799995E-01, 6.7860724042499998E-01, 1.2417984069000000E-02, -7.8536187672500002E-01, -8.8117631546499997E-01, -7.6852471065899997E-01, -1.6645633449700001E-01, -1.0647002787999999E-02, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 1.4859356507000000E-02, 1.4916739806500001E-01, 3.1887440099100001E-01, 4.3942392054700002E-01, 6.0159002287100005E-01, 7.6502829830300001E-01, 8.8793227607799996E-01, 9.5417702880999999E-01, 9.7751524428799996E-01, 9.9786091945400002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0064373409260492E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {5.1000000000000000E+01, 5.2125000000000000E+01, 5.3625000000000000E+01, 5.4375000000000000E+01, 5.4625000000000000E+01, 5.4875000000000000E+01, 5.5375000000000000E+01, 5.5625000000000000E+01, 5.5875000000000000E+01, 5.6000000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 1.1273801003000000E-02, 1.5930345105600000E-01, 4.1047207683300002E-01, 5.9488384393799998E-01, 6.7040758874399997E-01, 4.2646984319600001E-01, 2.1035967172100001E-01, 4.7844473307000000E-02, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 2.7570365871999999E-02, 1.2878661069399999E-01, 5.5494805308999995E-01, 6.0262802707200003E-01, 2.5753984837000001E-02, -8.1617800961599996E-01, -7.8731972174700005E-01, -5.5409463272399995E-01, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.4433153340000001E-03, 1.1270328001900000E-01, 3.0693511214899999E-01, 4.3270744580499998E-01, 5.9432468875300004E-01, 8.8690121775300002E-01, 9.6657707144299998E-01, 9.9772484520299998E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0027999802170684E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.4000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.6999999999999995E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C58</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    macroencapsulation<br>  Data taken from: Climator Sweden AB datasheet - last access 2022-10-14.<br><br>
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
end ClimSel_C58;