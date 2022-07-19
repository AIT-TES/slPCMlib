within slPCMlib.Media_Rubitherm_RT;
package RT42 "Rubitherm RT42; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT42";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.091500000000000e+02, 3.171500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.081500000000000e+02, 3.171500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.297905932028317e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 3.091500000000000e+02
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer pieces =   data_H.pieces;
    constant Integer order[:] = data_H.order;
    constant Real breaks[:] =   data_H.breaks;
    constant Real coefs[:,:] =  data_H.coefs;
  algorithm
    (xi, dxi) := BasicUtilities.quartQuintSplineEval(T-273.15,
                 pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer pieces =   data_C.pieces;
    constant Integer order[:] = data_C.order;
    constant Real breaks[:] =   data_C.breaks;
    constant Real coefs[:,:] =  data_C.coefs;
  algorithm
    (xi, dxi) := BasicUtilities.quartQuintSplineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 5, 5, 5, 5, 5, 6, 6, 1};
    constant Real[11] breaks = {-6.400000000000000e+01, 3.600000000000000e+01, 3.700000000000000e+01, 3.800000000000000e+01, 3.900000000000000e+01, 4.000000000000000e+01, 4.100000000000000e+01, 4.200000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 1.440000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.355916461745797e-13, 8.705795673669145e-16, -9.273025621602680e-19, 0.000000000000000e+00, 1.232017896876534e-03, -4.775029035763142e-04},  {7.545149934366814e-04, 3.295072063062393e-03, 5.157635415121497e-03, 2.770120897239057e-03, -1.002454069262043e-03, 1.325017830068663e-03, 0.000000000000000e+00},  {1.229990712966625e-02, 2.453597845829518e-02, 2.070345199195289e-02, 1.201048292087750e-02, 5.622635081081270e-03, -1.910632190504796e-03, 0.000000000000000e+00},  {7.326182339136829e-02, 1.149117105765009e-01, 7.136438933602421e-02, 1.539470134015463e-02, -3.930525871442709e-03, 9.359962483359539e-04, 0.000000000000000e+00},  {2.719380950209414e-01, 2.927824710248037e-01, 1.033253006111906e-01, 9.032560337743328e-03, 7.494553702370607e-04, -6.982621924245278e-04, 0.000000000000000e+00},  {6.771296201724916e-01, 5.260372637791566e-01, 1.279370919215972e-01, 5.047759894446295e-03, -2.741855591885578e-03, 6.851886780538671e-04, 0.000000000000000e+00},  {1.334095068853860e+00, 7.895132483283730e-01, 1.334811248341609e-01, 9.322243074426720e-04, 6.840877983837571e-04, -2.919681104060542e-02, 1.215558426035860e-02},  {2.241664527341973e+00, 9.889579724721040e-01, 3.074797804611892e-02, -4.518784969790461e-02, 3.703379650073559e-02, -1.607068229121709e-02, 2.887974330356656e-03},  {3.240033716702167e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 6, 5, 5, 5, 5, 5, 6, 5, 1};
    constant Real[12] breaks = {-6.500000000000000e+01, 3.500000000000000e+01, 3.600000000000000e+01, 3.700000000000000e+01, 3.800000000000000e+01, 3.900000000000000e+01, 4.000000000000000e+01, 4.100000000000000e+01, 4.200000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 1.440000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.129391969962751e-12, -7.696746953434361e-15, 1.108449787174081e-18, 5.782411586589357e-19, 2.153461634099179e-03, -8.879739805445195e-04},  {1.265487652417572e-03, 5.439464286098824e-03, 8.215006632816402e-03, 3.775136730101402e-03, -2.552301537671898e-03, 1.919136736971039e-03, -4.572505042650397e-04},  {1.760467999646830e-02, 2.983786225061091e-02, 1.655921740282399e-02, 3.612287863823399e-03, 1.846245832076978e-04, 4.731972798304760e-04, 0.000000000000000e+00},  {6.827186937676477e-02, 7.689764537971691e-02, 3.323580129184516e-02, 9.082758994958950e-03, 2.550610982360078e-03, -9.163193165360722e-04, 0.000000000000000e+00},  {1.891223667091098e-01, 1.762383722951227e-01, 6.662455100552230e-02, 1.012200975903854e-02, -2.030985600320283e-03, 3.335757539572740e-04, 0.000000000000000e+00},  {4.404098899224303e-01, 3.333974399517480e-01, 8.814042422028869e-02, 5.333824897330149e-03, -3.631068305339133e-04, 2.436138916388434e-03, 0.000000000000000e+00},  {8.693546110776517e-01, 5.364080303441661e-01, 1.263246470929603e-01, 2.824278673907884e-02, 1.181758775140826e-02, -1.248150554029884e-02, 0.000000000000000e+00},  {1.559666157464966e+00, 8.586485080514048e-01, 1.571434784156575e-01, -4.930191765827657e-02, -5.058993995008596e-02, 4.982908359452495e-02, -1.287480229063413e-02},  {2.512520567627556e+00, 9.945665563369692e-01, 1.086688732605357e-02, -1.086688732605359e-02, 5.433443663026796e-03, -1.086688732605359e-03, 0.000000000000000e+00},  {3.511433878894947e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.600000000000000e+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 2.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 2.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT42</strong>.<br><br>
  Information taken from: data_sheet - last access 01.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  These functions are modelled by piece-wise splines using <strong>variable order quartic and quintic</strong> method,
  see also 
  <blockquote>
  <p>
  Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase
  Fraction–Temperature Curves from Heat Capacity Data for Numerical
  Modeling of Heat Transfer in Commercial Paraffin Waxes.
  Energies, 13(19), 5149.
  <a href>doi.org/10.3390/en13195149</a>.
  </p>
  </blockquote>
  </p></html>",
  revisions="<html>
  <ul>
  <li>file creation date: 19-Jul-2022  </ul>
  </html>"));
end RT42;
