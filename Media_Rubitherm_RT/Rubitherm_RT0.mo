within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT0 "Rubitherm GmbH, RT0; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT0";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6414999999999998E+02, 2.7814999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6414999999999998E+02, 2.7414999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8379605832213687E+05
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
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {-9.0000000000000000E+00, -6.1250000000000000E+00, -3.3750000000000000E+00, -1.6250000000000000E+00, -8.7500000000000000E-01, -1.2500000000000000E-01, 3.7500000000000000E-01, 1.1250000000000000E+00, 1.8750000000000000E+00, 2.3750000000000000E+00, 2.8750000000000000E+00, 3.1250000000000000E+00, 4.8750000000000000E+00, 5.0000000000000000E+00};
    constant Real[14] data_y =   {0.0000000000000000E+00, 1.5916425509000000E-02, 3.4741471879999999E-02, 1.2470160227500000E-01, 1.8042937108700000E-01, 1.3869025900300000E-01, 1.6743389252499999E-01, 2.7269322368600002E-01, 2.3911433185500000E-01, 7.4161959412999995E-02, 3.5983665869999998E-03, 1.2328110529999999E-03, 5.7487734300000004E-04, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -7.9417791600000000E-04, 1.6714576703999999E-02, 1.0508629409700000E-01, -1.0756729729000001E-02, -2.3443594060000000E-03, 1.4686919626700001E-01, 8.5727258087999997E-02, -3.1333607895400001E-01, -3.2973252080600002E-01, -1.1280972909999999E-02, -1.9025031810000000E-03, -6.8685442150000003E-03, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 2.3783547727999999E-02, 8.3296440334999999E-02, 2.0203654095800000E-01, 3.2371547751899998E-01, 4.4480687814800002E-01, 5.1934709463800000E-01, 6.8981713575400005E-01, 9.0365778660899998E-01, 9.8351599414299995E-01, 9.9651662396800000E-01, 9.9708012515900002E-01, 9.9997260273800004E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0152241552172538E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    11;
    constant Real[11] data_x =   {-9.0000000000000000E+00, -7.1250000000000000E+00, -4.8750000000000000E+00, -3.8750000000000000E+00, -1.3750000000000000E+00, -6.2500000000000000E-01, -3.7500000000000000E-01, -1.2500000000000000E-01, 6.2500000000000000E-01, 8.7500000000000000E-01, 1.0000000000000000E+00};
    constant Real[11] data_y =   {0.0000000000000000E+00, 1.7468771978000001E-02, 8.5002249547000003E-02, 3.7030189871999998E-02, 1.6835286464200000E-01, 2.9189763644000000E-01, 3.7396569220300002E-01, 3.9022760316600003E-01, 1.0997699133500000E-01, 2.4723502279000002E-02, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 1.3451312537000000E-02, -5.0444102933000001E-02, -3.2928751770000000E-03, 6.6013289476999998E-02, 2.5164026122799998E-01, 2.6159114519100002E-01, -9.6837353035000004E-02, -4.4052155208999999E-01, -2.8404340630000002E-01, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.2464324567000000E-02, 1.5502223369699999E-01, 2.1223847331699999E-01, 4.3337001140999998E-01, 5.9763386444500000E-01, 6.8100333915000000E-01, 7.7861488446000005E-01, 9.8276310852700000E-01, 9.9882196734700002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0022647713425652E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.8000000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.7000000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 2.0000000000000001E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 2.0000000000000001E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT0</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-09-30.<br><br>
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
end Rubitherm_RT0;
