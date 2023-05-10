
// within slPCMlib.PLUSS_HS;
package PLUSS_savE_HS7N "Pluss Advanced Technolgies Pvt Ltd, HS7N; data taken from: PLUSS datasheet; last access: 2022-02-13."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS7N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.6114999999999998E+02, 2.7214999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.6014999999999998E+02, 2.6914999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.5000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.6732080512716109E+05
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
    constant Real[14] data_x =   {-1.2000000000000000E+01, -9.3750000000000000E+00, -7.3750000000000000E+00, -6.8750000000000000E+00, -6.6250000000000000E+00, -6.3750000000000000E+00, -6.1250000000000000E+00, -5.8750000000000000E+00, -5.6250000000000000E+00, -5.3750000000000000E+00, -4.6250000000000000E+00, -3.8750000000000000E+00, -2.1250000000000000E+00, -1.0000000000000000E+00};
    constant Real[14] data_y =   {0.0000000000000000E+00, 1.4922426934000000E-02, 4.3550766941999997E-02, 9.9627670532000001E-02, 1.9865158684799999E-01, 4.2742861604999999E-01, 6.1993009974199997E-01, 6.8892047085800001E-01, 6.0311215626600001E-01, 4.0141851697399999E-01, 6.1014008496000000E-02, 2.2983301359999998E-03, 7.3001569879999999E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 2.6809138929999999E-03, -1.0542691070000000E-03, 2.0600562414699999E-01, 7.1421380288699998E-01, 8.4076891721600000E-01, 6.7621133030500002E-01, -4.3497871153000002E-02, -6.9535845000200003E-01, -6.6219765247999995E-01, -2.3960185673400000E-01, -2.4865358419999999E-03, -8.8719640900000001E-03, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 1.8222317477999999E-02, 7.8523196537000003E-02, 1.1031119220900000E-01, 1.4528711748600001E-01, 2.2364509260200000E-01, 3.5670764681799999E-01, 5.2569520318100005E-01, 6.9220318079599996E-01, 8.1882016970000004E-01, 9.7392177303200000E-01, 9.8667231446799997E-01, 9.9679844337099999E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0097562158235354E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {-1.3000000000000000E+01, -1.1125000000000000E+01, -9.1250000000000000E+00, -7.8750000000000000E+00, -7.1250000000000000E+00, -6.6250000000000000E+00, -6.1250000000000000E+00, -5.1250000000000000E+00, -4.0000000000000000E+00};
    constant Real[9] data_y =   {0.0000000000000000E+00, 2.0390878134999998E-02, 4.3808415747000000E-02, 2.2602718809799999E-01, 3.7833916608399998E-01, 3.4020740640800001E-01, 2.2796361043600000E-01, 1.1667510834100001E-01, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 3.8999634259999998E-03, 4.3181745890999997E-02, 2.2462517869199999E-01, 6.8364738878000000E-02, -2.2116788311800001E-01, -1.6330191251699999E-01, -1.3487882840700000E-01, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 1.8047443202999999E-02, 6.9361970403000003E-02, 2.1497731188700001E-01, 4.4989694805000002E-01, 6.3632541337299997E-01, 7.7773903615600004E-01, 9.4838536894500003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0040927366648043E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.0260000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1200000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.7600000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.1200000000000000E+03;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS7N</strong>  from manufacturer: <strong>Pluss Advanced Technolgies Pvt Ltd</strong>.<br>
       material class: salt hydrate-based;  encapsulation:    multiple options available<br>  Data taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS7N;