within slPCMlib.Media_Rubitherm_RT;
package RT25HC "Rubitherm RT25HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT25HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.901500000000000e+02, 3.001500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.891500000000000e+02, 2.991500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.795288869777121e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+25
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
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 5, 6, 1};
    constant Real[13] breaks = {-8.300000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 1.270000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.509605903512589e-12, -5.045224517099588e-15, -1.121095649494404e-20, 1.445602896647339e-19, 4.935495879238497e-04, -1.430032612364554e-04},  {3.505463251727433e-04, 1.609728370665352e-03, 2.790446960686569e-03, 2.075430654509389e-03, 3.226990210724172e-04, -2.552843245505609e-04, 0.000000000000000e+00},  {6.893567007555910e-03, 1.343128871716477e-02, 8.400089805143835e-03, 8.133834932934482e-04, -9.537226016803872e-04, 6.129047195537346e-04, 0.000000000000000e+00},  {2.919751114103131e-02, 3.192125199836533e-02, 1.124695187047915e-02, 3.127540282109244e-03, 2.110800996088286e-03, 3.684776490356964e-05, 0.000000000000000e+00},  {7.764090405297690e-02, 7.242521939459531e-02, 3.366285634237254e-02, 1.193922191549807e-02, 2.295039820606134e-03, -4.984595549911602e-03, 1.524083276319423e-03},  {1.945027292524568e-01, 1.689702790165458e-01, 5.626605465817867e-02, 1.755091225195063e-03, 2.333112158394752e-04, -6.993680840507049e-05, 0.000000000000000e+00},  {4.216575285598108e-01, 2.873512228299239e-01, 6.223182754475035e-02, 1.988968004502259e-03, -1.163728261858772e-04, 2.305488349488382e-03, 0.000000000000000e+00},  {7.754186624622897e-01, 4.288437323755471e-01, 9.055537809602542e-02, 2.457836019464257e-02, 1.141106892125603e-02, -6.890062257909263e-03, 0.000000000000000e+00},  {1.323917139791851e+00, 6.948835335472355e-01, 1.638562496283975e-01, 1.322013300574064e-03, -2.303924236829029e-02, 4.323494710226515e-03, 0.000000000000000e+00},  {2.165263188609995e+00, 9.560225767839270e-01, 7.282178242264377e-02, -4.760000907032192e-02, -1.421768817157711e-03, 1.541741777482275e-02, -5.044354670463734e-03},  {3.155458833033445e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[13] breaks = {-8.400000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 1.260000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 3.001593945259695e-14, 4.600302141616118e-16, -8.166602907170884e-19, 0.000000000000000e+00, 1.159888751162218e-03, -5.517705462378433e-04},  {6.081182049548496e-04, 2.488820478412940e-03, 3.322329318054955e-03, 5.634765868653103e-04, -2.477114437756561e-03, 1.655922322283309e-03, 0.000000000000000e+00},  {6.161552472814803e-03, 9.195062735506427e-03, 6.709295674944571e-03, 7.214242058672155e-03, 5.802497173659983e-03, -2.620779814223965e-03, 0.000000000000000e+00},  {3.246187030137398e-02, 5.436246988492144e-02, 3.695920675068112e-02, 4.216432611072439e-03, -7.301401897459841e-03, 1.904747853081163e-03, 0.000000000000000e+00},  {1.226033255036703e-01, 1.212483128950683e-01, 2.484757172995105e-02, -5.941696447955285e-03, 2.222337367945977e-03, -1.463836045386099e-04, 0.000000000000000e+00},  {2.648334674441417e-01, 1.612757984601916e-01, 1.889267054837490e-02, 1.483816978442522e-03, 1.490419345252927e-03, 4.851909737533689e-05, 0.000000000000000e+00},  {4.480246918737790e-01, 2.097168633601536e-01, 3.277182852897334e-02, 7.930685333207598e-03, 1.733014832129612e-03, 1.024033588690638e-03, -4.528127915809126e-04},  {7.007483047253528e-01, 3.083879269402172e-01, 7.041011753456662e-02, 1.604682471701417e-02, 6.099090186910983e-05, 2.354165911091217e-03, 0.000000000000000e+00},  {1.098008330730111e+00, 5.093634293233461e-01, 1.424581962077363e-01, 3.983244743540278e-02, 1.183182045732519e-02, -1.585059585374245e-02, 0.000000000000000e+00},  {1.785643628300179e+00, 8.818514666056567e-01, 1.744405027204719e-01, -7.134622927272095e-02, -6.742115881138706e-02, 7.534079583092593e-02, -2.061885468954950e-02},  {2.757890150683576e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.700000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT25HC</strong>.<br><br>
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
  <li>file creation date: 07-Jul-2022  </ul>
  </html>"));
end RT25HC;
