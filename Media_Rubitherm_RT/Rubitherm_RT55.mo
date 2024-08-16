within slPCMlib.Media_Rubitherm_RT;
package Rubitherm_RT55 "Rubitherm GmbH, RT55; data taken from: Rubitherm datasheet."
  extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT55";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1914999999999998E+02, 3.3214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1914999999999998E+02, 3.3114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.2942853374215616E+05
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
    constant Real[11] data_x =   {4.6000000000000000E+01, 4.7125000000000000E+01, 4.8875000000000000E+01, 5.1125000000000000E+01, 5.3375000000000000E+01, 5.4375000000000000E+01, 5.5125000000000000E+01, 5.6125000000000000E+01, 5.7375000000000000E+01, 5.8625000000000000E+01, 5.9000000000000000E+01};
    constant Real[11] data_y =   {0.0000000000000000E+00, 4.2148719249999996E-03, 2.5608494432999999E-02, 6.9223633142000002E-02, 1.6574514612800001E-01, 1.8086723803599999E-01, 2.0352251369200000E-01, 1.1980373229000001E-01, 3.2792927200000002E-02, 5.8809248900000002E-04, 0.0000000000000000E+00};
    constant Real[11] m_k =      {0.0000000000000000E+00, 9.2668454159999993E-03, 5.8777529730000000E-03, 3.5777084107000000E-02, -7.7283883430000000E-03, 5.9144131092999999E-02, -4.4158118679999997E-02, -6.9433480236999995E-02, -7.6316251747000005E-02, -2.3631892329999999E-03, 0.0000000000000000E+00};
    constant Real[11] iy_start = {0.0000000000000000E+00, 1.4030480269999999E-03, 2.8548090519000001E-02, 1.2326482763000000E-01, 4.0789496222900001E-01, 5.7677737988099997E-01, 7.2678636612799996E-01, 8.9167755174899999E-01, 9.8860607922999999E-01, 9.9991686067100005E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0068497685044213E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {4.6000000000000000E+01, 4.8625000000000000E+01, 5.2625000000000000E+01, 5.5125000000000000E+01, 5.5875000000000000E+01, 5.6625000000000000E+01, 5.7125000000000000E+01, 5.8000000000000000E+01};
    constant Real[8] data_y =   {0.0000000000000000E+00, 2.5119151045000000E-02, 1.4688390918999999E-01, 1.9978946452200000E-01, 2.0381297590799999E-01, 3.5151770449000003E-02, 5.0662695150000001E-03, 0.0000000000000000E+00};
    constant Real[8] m_k =      {0.0000000000000000E+00, 5.5196369800000001E-03, 5.7293365153000000E-02, 4.0163821327999998E-02, -1.6980808391600000E-01, -1.7545693185100000E-01, -1.5782676636000002E-02, 0.0000000000000000E+00};
    constant Real[8] iy_start = {0.0000000000000000E+00, 2.9620361969999998E-02, 3.0294270688899999E-01, 7.4254879636799997E-01, 9.0277364143399996E-01, 9.9211020279499995E-01, 9.9879774196700000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 9.9399166799052152E-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT55</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
  Basic characteristics are the material class: unknown, and encapsulation: multiple options available<br>  The data is taken from: Rubitherm datasheet - last access 2020-10-09.<br><br>
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
end Rubitherm_RT55;
