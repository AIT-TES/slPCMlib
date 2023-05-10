
// within slPCMlib.Media_Climator;
package ClimSel_C10 "Climator Sweden AB, ClimSel C10; data taken from: Climator Sweden AB datasheet; last access: 2022-10-14."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "ClimSel C10";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.8014999999999998E+02, 2.8814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.7614999999999998E+02, 2.8214999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.5000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 7.7700000000000000E+04
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
    constant Real[13] data_x =   {7.0000000000000000E+00, 9.1250000000000000E+00, 9.8750000000000000E+00, 1.0125000000000000E+01, 1.0375000000000000E+01, 1.0625000000000000E+01, 1.1125000000000000E+01, 1.1375000000000000E+01, 1.1625000000000000E+01, 1.1875000000000000E+01, 1.2875000000000000E+01, 1.4125000000000000E+01, 1.5000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.5753051366000000E-02, 5.6291148558999998E-02, 1.1425968803200000E-01, 2.5837504036400000E-01, 6.1313959637799997E-01, 8.2191192425000004E-01, 5.3634209678300004E-01, 1.8896852990800000E-01, 7.0537354928999998E-02, 2.1692592430000001E-02, 3.5300473058999998E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, 1.4044256563000001E-02, 1.1156227007699999E-01, 2.8225468645700003E-01, 9.4087287015099996E-01, 1.2227967792570000E+00, -9.5667370981199995E-01, -1.2719644917910000E+00, -9.8625518607600005E-01, -1.9367071014900000E-01, 1.1003686020000000E-02, -4.2411622319000003E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.2149189108000001E-02, 4.8429389200999998E-02, 6.8925342338999995E-02, 1.1221403308800000E-01, 2.2003284099000001E-01, 6.2550945596700003E-01, 7.9748816500399999E-01, 8.8695254320000005E-01, 9.1535436117000002E-01, 9.4450718448200000E-01, 9.8722076607900000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0032365203654969E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    13;
    constant Real[13] data_x =   {3.0000000000000000E+00, 4.3750000000000000E+00, 4.8750000000000000E+00, 5.1250000000000000E+00, 5.3750000000000000E+00, 5.6250000000000000E+00, 5.8750000000000000E+00, 6.1250000000000000E+00, 6.3750000000000000E+00, 6.6250000000000000E+00, 7.3750000000000000E+00, 8.3750000000000000E+00, 9.0000000000000000E+00};
    constant Real[13] data_y =   {0.0000000000000000E+00, 2.9283791678000000E-02, 3.1156102251000000E-02, 6.5115631271999999E-02, 1.6790622356400001E-01, 4.7055108647200000E-01, 7.4663239602300002E-01, 8.3840364662699995E-01, 6.9333683070600005E-01, 3.9444045760300001E-01, 6.2134968394999997E-02, 2.4331487650000000E-03, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -1.1096656922999999E-02, 5.1768252394000000E-02, 1.6706777423600000E-01, 9.8797403556499996E-01, 1.1553026945510001E+00, 1.0299759599070000E+00, -1.4928913323500001E-01, -9.8264234608000001E-01, -8.5546447972700002E-01, -1.6162241871399999E-01, -7.8827411690000002E-03, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 2.1930020333999999E-02, 3.5761283228999997E-02, 4.7220394181999997E-02, 7.2128353805999995E-02, 1.5124118680699999E-01, 3.0438482972800002E-01, 5.0911483227099996E-01, 7.0536226811900005E-01, 8.4097574995400004E-01, 9.7997898494000002E-01, 9.9949510996000002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0022445201295596E+00;
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
    lambda := 8.2999999999999996E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.4000000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>ClimSel C10</strong>  from manufacturer: <strong>Climator Sweden AB</strong>.<br>
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
end ClimSel_C10;