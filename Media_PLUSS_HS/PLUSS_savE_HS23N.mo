
within slPCMlib.Media_PLUSS_HS;
package PLUSS_savE_HS23N "Pluss Advanced Technologies Pvt Ltd, HS23N; data taken from: PLUSS datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "HS23N";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.4814999999999998E+02, 2.5914999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.4714999999999998E+02, 2.5514999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.5800000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {3.4000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.4814184369412460E+05
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
    constant Real[13] data_x =   {-2.5000000000000000E+01, -2.2375000000000000E+01, -2.1875000000000000E+01, -2.1625000000000000E+01, -2.1375000000000000E+01, -2.1125000000000000E+01, -2.0625000000000000E+01, -2.0375000000000000E+01, -2.0125000000000000E+01, -1.9625000000000000E+01, -1.6875000000000000E+01, -1.5625000000000000E+01, -1.4000000000000000E+01};
    constant Real[13] data_y =   {0.0000000000000000E+00, 1.6383583660000000E-03, 1.5064789065000000E-02, 7.6579435540999996E-02, 4.0322122678799999E-01, 7.6288848190100000E-01, 6.5860304994800001E-01, 2.7550096880800001E-01, 1.2391549778900000E-01, 4.6254639372000000E-02, 3.8074503354999997E-02, 2.5754600611000000E-02, 0.0000000000000000E+00};
    constant Real[13] m_k =      {0.0000000000000000E+00, -4.5269433950000000E-03, 5.3158499040000001E-02, 8.8461786746500004E-01, 1.3752472774450000E+00, 1.3505403654340000E+00, -1.2638682362880000E+00, -9.4150123737699998E-01, -2.8853032121200001E-01, 7.4128363360000003E-03, 9.2685490280000008E-03, -4.1303099583999997E-02, 0.0000000000000000E+00};
    constant Real[13] iy_start = {0.0000000000000000E+00, 4.7541473890000000E-03, 7.7308755290000000E-03, 1.4862405726000000E-02, 7.2334668315000003E-02, 2.1836056102400001E-01, 6.2857530935600003E-01, 7.4376462331500004E-01, 7.9033336472899995E-01, 8.2674370300599997E-01, 9.4163181969900001E-01, 9.8815238899500002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0009150394011734E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    14;
    constant Real[14] data_x =   {-2.6000000000000000E+01, -2.4375000000000000E+01, -2.3625000000000000E+01, -2.3375000000000000E+01, -2.3125000000000000E+01, -2.2875000000000000E+01, -2.2375000000000000E+01, -2.2125000000000000E+01, -2.1375000000000000E+01, -2.0625000000000000E+01, -1.9625000000000000E+01, -1.9125000000000000E+01, -1.8125000000000000E+01, -1.8000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 9.6982058809999994E-02, 3.2960103423300002E-01, 5.5515741833499999E-01, 6.4708672659800004E-01, 5.5153408345199995E-01, 9.8627216988000005E-02, 3.3762145068000002E-02, 2.8184265499000001E-02, 4.5967605268999998E-02, 6.2224114306999999E-02, 1.1205980841500000E-01, 8.4416750220000000E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, 7.2892453597000001E-02, 6.2244595025799998E-01, 7.1691530749099996E-01, -1.0248222397000000E-02, -9.0486224952799998E-01, -8.8734923046100000E-01, -8.7944132098999997E-02, 4.8695074844000000E-02, -1.6586269696000001E-02, 9.2683988733000003E-02, 1.9702561701000001E-02, -9.8124693535999993E-02, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 6.2856643636999998E-02, 1.9727639061900001E-01, 3.0755260504699999E-01, 4.6186312346399999E-01, 6.1659352339100004E-01, 7.7902445548199994E-01, 7.9142906674000002E-01, 8.0828051080999996E-01, 8.3919614761700001E-01, 8.8425702155300001E-01, 9.2941947740800002E-01, 9.9959953200700002E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0015752169148342E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.0780000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.1550000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 4.9760000000000000E+00;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 7.0199999999999996E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>HS23N</strong>  from manufacturer: <strong>Pluss Advanced Technologies Pvt Ltd</strong>.<br>
  Basic characteristics are the material class: salt hydrate-based, and encapsulation: multiple options available<br>  The data is taken from: PLUSS datasheet - last access 2022-02-13.<br><br>
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
end PLUSS_savE_HS23N;