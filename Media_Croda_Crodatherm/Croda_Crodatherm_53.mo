
within slPCMlib.Media_Croda_Crodatherm;
package Croda_Crodatherm_53 "Croda International Plc, Crodatherm 53; data taken from: Croda datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "Crodatherm 53";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {3.1814999999999998E+02, 3.3014999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.1714999999999998E+02, 3.2714999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.9000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.2000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.8527040821919346E+05
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
    constant Real[12] data_x =   {4.5000000000000000E+01, 4.7875000000000000E+01, 4.9625000000000000E+01, 5.1125000000000000E+01, 5.1375000000000000E+01, 5.1625000000000000E+01, 5.2375000000000000E+01, 5.3625000000000000E+01, 5.4125000000000000E+01, 5.5625000000000000E+01, 5.6875000000000000E+01, 5.7000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 6.4872159200000000E-03, 2.9205371378000000E-02, 1.3595090341400001E-01, 2.0331617946200001E-01, 3.0938664826399997E-01, 4.6755811344300002E-01, 1.1547324068000001E-01, 4.0601708737000003E-02, 1.9403670009999999E-03, 1.7549406000000001E-05, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 5.5293062260000003E-03, 7.5851621810000000E-03, 1.9933605933500001E-01, 3.6508970782400002E-01, 3.9177467364399998E-01, -7.9366208339000002E-02, -2.9836772528700001E-01, -7.4919766440999999E-02, -4.4180467259999998E-03, -1.8605863500000000E-04, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 5.5741613050000003E-03, 3.6599914171000003E-02, 1.2542831667400001E-01, 1.6740555805800000E-01, 2.3201963008700000E-01, 5.4871906965600004E-01, 9.4571658336499997E-01, 9.8043760833899996E-01, 9.9931948184399999E-01, 9.9999913653700001E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0104021035505359E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {4.4000000000000000E+01, 4.6125000000000000E+01, 4.9125000000000000E+01, 5.0625000000000000E+01, 5.1375000000000000E+01, 5.1625000000000000E+01, 5.1875000000000000E+01, 5.2125000000000000E+01, 5.2375000000000000E+01, 5.2625000000000000E+01, 5.3125000000000000E+01, 5.4000000000000000E+01};
    constant Real[12] data_y =   {0.0000000000000000E+00, 2.9242716719999999E-03, 3.9060089430999999E-02, 2.3006329704500000E-01, 5.2468007592300003E-01, 6.6372551409400005E-01, 6.1034244659000003E-01, 4.0889371867199997E-01, 1.2718372733100000E-01, 8.8098906919999997E-03, 1.6221408240000000E-03, 0.0000000000000000E+00};
    constant Real[12] m_k =      {0.0000000000000000E+00, 2.6711292010000000E-03, 3.3300662048000002E-02, 3.0847140224899999E-01, 4.6256799240599999E-01, 4.4383752249500003E-01, -6.3147319782800004E-01, -9.5844574924100001E-01, -8.5293560208399999E-01, -4.7497410015999997E-02, 0.0000000000000000E+00, 0.0000000000000000E+00};
    constant Real[12] iy_start = {0.0000000000000000E+00, 2.1180539660000000E-03, 4.2430142281999997E-02, 1.9383381049100001E-01, 4.7176067293500001E-01, 6.2155226360399995E-01, 7.8767935927300003E-01, 9.1777990044299995E-01, 9.8475123158900002E-01, 9.9765392765799998E-01, 9.9928485478999995E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0076915628705319E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 9.0400000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 8.2900000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 2.8000000000000003E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.6000000000000000E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Crodatherm 53</strong>  from manufacturer: <strong>Croda International Plc</strong>.<br>
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
end Croda_Crodatherm_53;