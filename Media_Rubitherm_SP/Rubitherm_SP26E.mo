
within slPCMlib.Media_Rubitherm_SP;
package Rubitherm_SP26E "Rubitherm GmbH, SP26E; data taken from: Rubitherm datasheet."
  extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP26E";
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] = {2.9114999999999998E+02, 3.0314999999999998E+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.8914999999999998E+02, 3.0114999999999998E+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.0000000000000000E+03, 0.0000000000000000E+00}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.7500000000000000E+05
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
    constant Real[14] data_x =   {1.8000000000000000E+01, 1.9375000000000000E+01, 2.1875000000000000E+01, 2.3625000000000000E+01, 2.4625000000000000E+01, 2.5125000000000000E+01, 2.5375000000000000E+01, 2.5625000000000000E+01, 2.6125000000000000E+01, 2.6375000000000000E+01, 2.6625000000000000E+01, 2.6875000000000000E+01, 2.7875000000000000E+01, 3.0000000000000000E+01};
    constant Real[14] data_y =   {0.0000000000000000E+00, 1.0244429522000000E-02, 1.4511103316000000E-02, 8.9664905958999994E-02, 6.2123171876000000E-02, 1.1424382499000001E-01, 2.3336285278499999E-01, 5.3420332572499996E-01, 7.0774801649999997E-01, 4.6772460149000000E-01, 1.7179881218099999E-01, 6.4600205640000005E-02, 5.5642737180000004E-03, 0.0000000000000000E+00};
    constant Real[14] m_k =      {0.0000000000000000E+00, -1.0858215235000000E-02, 7.8617077279999996E-03, 6.3650567934999999E-02, -6.1942019369999997E-02, 2.2539564524299999E-01, 7.9012431494699997E-01, 1.0286229817170001E+00, -7.9456589113400000E-01, -1.0783214418280000E+00, -8.4431417677499998E-01, -1.8699820051300001E-01, -4.9804360520000002E-03, 0.0000000000000000E+00};
    constant Real[14] iy_start = {0.0000000000000000E+00, 8.8223637490000006E-03, 3.0182874022000001E-02, 1.0770172254500000E-01, 1.9473842351400000E-01, 2.3314252032999999E-01, 2.7396944361199999E-01, 3.6941501807499999E-01, 7.2061615158500003E-01, 8.7019089887599999E-01, 9.4952930521500001E-01, 9.7586035623599998E-01, 9.9593047160000003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0078347894912689E+00;
  end data_H;

  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    15;
    constant Real[15] data_x =   {1.6000000000000000E+01, 1.7375000000000000E+01, 1.9875000000000000E+01, 2.2875000000000000E+01, 2.3625000000000000E+01, 2.4375000000000000E+01, 2.4625000000000000E+01, 2.4875000000000000E+01, 2.5125000000000000E+01, 2.5375000000000000E+01, 2.5625000000000000E+01, 2.5875000000000000E+01, 2.6375000000000000E+01, 2.7875000000000000E+01, 2.8000000000000000E+01};
    constant Real[15] data_y =   {0.0000000000000000E+00, 1.0731890765000000E-02, 1.1087070358000000E-02, 8.5618349129999993E-02, 1.1289996892900001E-01, 3.8149461141899998E-01, 6.8587417830300001E-01, 7.7318854247100000E-01, 5.6637579334599997E-01, 2.3962172684399999E-01, 2.1772428343999999E-02, 2.3239571979999998E-03, 0.0000000000000000E+00, 1.1255014150000000E-03, 0.0000000000000000E+00};
    constant Real[15] m_k =      {0.0000000000000000E+00, -1.1566784774000000E-02, 5.1830253700000002E-04, 4.6021091574000000E-02, 3.3990181831999998E-02, 7.8184588554900003E-01, 9.2652220217199999E-01, -4.0852920518000002E-01, -1.1773770791439999E+00, -1.1428886874380000E+00, -1.8562653055299999E-01, -1.0088887556000001E-02, 1.8130196690000001E-03, -9.0029353110000002E-03, 0.0000000000000000E+00};
    constant Real[15] iy_start = {0.0000000000000000E+00, 9.2604104309999997E-03, 3.0376300496000000E-02, 1.4202912448100000E-01, 2.1752549355699999E-01, 3.6884594542200000E-01, 5.0237673966700003E-01, 6.9294491215400000E-01, 8.6551043052300003E-01, 9.6673486507499995E-01, 9.9460355291600000E-01, 9.9671499207900005E-01, 9.9705019189800004E-01, 9.9994099730700003E-01, 1.0000000000000000E+00};
    constant Real    iy_scaler = 1.0065066238267109E+00;
  end data_C;

  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.6000000000000000E+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.5000000000000000E+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 5.0000000000000000E-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 5.0000000000000000E-01;
  end conductivity_liquid;


annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>SP26E</strong>  from manufacturer: <strong>Rubitherm GmbH</strong>.<br>
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
end Rubitherm_SP26E;