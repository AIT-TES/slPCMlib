
within slPCMlib.Media_Knauf_SmartBoard;
package SmartBoard_21 "Knauf Gips KG, SmartBoard 21; data taken from: DBU-Abschlussbericht-AZ-23836.pdf."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SmartBoard 21";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9314999999999998E+02, 3.0139999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.9314999999999998E+02, 3.0039999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.2000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.2000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.5905813368276304E+04
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
    constant Integer  len_x =    9;
    constant Real[9] data_x =   {2.0000000000000000E+01, 2.3875000000000000E+01, 2.4875000000000000E+01, 2.5375000000000000E+01, 2.6125000000000000E+01, 2.6375000000000000E+01, 2.6875000000000000E+01, 2.7625000000000000E+01, 2.8250000000000000E+01};
    constant Real[9] data_y =   {0.0000000000000000E+00, 3.5495152774999997E-02, 9.7485655835000004E-02, 1.7815058052200000E-01, 5.3699678911600002E-01, 6.1530212280300001E-01, 4.4913942753000002E-01, 1.6062598016000001E-02, 0.0000000000000000E+00};
    constant Real[9] m_k =      {0.0000000000000000E+00, 2.0900746365999999E-02, 1.3929809375100000E-01, 1.9692141481200001E-01, 5.5176984768399995E-01, 1.7377224372600000E-01, -7.2302176744299995E-01, -7.0878688539000007E-02, 0.0000000000000000E+00};
    constant Real[9] iy_start = {0.0000000000000000E+00, 4.2794530236000002E-02, 9.9652085370999996E-02, 1.6763998351199999E-01, 4.2022445761900001E-01, 5.6683289261400005E-01, 8.5280137609399997E-01, 9.9727649772200000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0041253994008581E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {2.0000000000000000E+01, 2.3125000000000000E+01, 2.3875000000000000E+01, 2.4625000000000000E+01, 2.5375000000000000E+01, 2.5625000000000000E+01, 2.6375000000000000E+01, 2.6625000000000000E+01, 2.6875000000000000E+01, 2.7250000000000000E+01};
    constant Real[10] data_y =   {0.0000000000000000E+00, 3.1555221236999997E-02, 1.7871430272899999E-01, 4.2748980250800001E-01, 4.9089503724299999E-01, 4.3837743024900000E-01, 6.0417624481000000E-02, 6.4919349899999999E-03, 1.9710545100000000E-04, 0.0000000000000000E+00};
    constant Real[10] m_k =      {0.0000000000000000E+00, 2.2872079550000000E-02, 3.8231510490400000E-01, 1.9029717288199999E-01, -4.8858455188999998E-02, -4.8166105008200000E-01, -4.1684469009100000E-01, -4.1592185629000002E-02, -8.8725551200000000E-04, 0.0000000000000000E+00};
    constant Real[10] iy_start = {0.0000000000000000E+00, 3.0873562148000001E-02, 9.3243158689000002E-02, 3.3097098408800002E-01, 6.8868298655200000E-01, 8.0779792664899996E-01, 9.9289822395000005E-01, 9.9934545916799999E-01, 9.9997328286400000E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0059258701712850E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 7.6700000000000000E+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.6700000000000000E+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.7999999999999999E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.7999999999999999E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SmartBoard 21</strong>  from manufacturer: <strong>Knauf Gips KG</strong>.<br>
  Basic characteristics are the material class: paraffin-based composite, and encapsulation: microencapsulated<br>  The data is taken from: DBU-Abschlussbericht-AZ-23836.pdf - last access 2010-10-01.<br><br>
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
end SmartBoard_21;