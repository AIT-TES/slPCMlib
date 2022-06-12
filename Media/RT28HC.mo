within slPCMlib.Media;
package RT28HC "Rubitherm RT28HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT28HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.961500000000000e+02, 3.031500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.961500000000000e+02, 3.021500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 2.221900304037665e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.SIunits.Temp_K            Tref = 273.15+25
             "reference temperature";
    constant Modelica.SIunits.SpecificEnthalpy  href = 0.0
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
    (xi, dxi) := BasicUtilities.splineEval(T-273.15,
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
    (xi, dxi) := BasicUtilities.splineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 6, 5, 5, 6, 5, 6, 5, 1};
    constant Real[10] breaks = {-7.700000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01, 1.300000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 5.613211432370316e-15, -1.804507938597140e-18, -1.711461988436696e-18, 5.782411586589357e-19, 1.803715262963865e-03, -7.525081622380055e-04},  {1.051207100731470e-03, 4.503527341392023e-03, 6.749530196068566e-03, 2.986989384878540e-03, -2.269046118750757e-03, 6.898223170455629e-04, 0.000000000000000e+00},  {1.371203022136540e-02, 2.133648299838932e-02, 8.994444808655270e-03, 8.090280803311392e-04, 1.180065466477057e-03, -5.141711106298510e-05, 0.000000000000000e+00},  {4.598063446415521e-02, 4.621563316728100e-02, 1.798775073788118e-02, 5.015118835609507e-03, 9.229799111621315e-04, 1.442923721019392e-02, -4.092899998242897e-03},  {1.264584543280401e-01, 1.485171968559642e-01, 1.214698588399782e-01, 7.114141061733928e-02, 1.167566598848827e-02, -1.283468277687080e-02, 0.000000000000000e+00},  {4.664279038529393e-01, 5.874103964577365e-01, 2.766012588542177e-01, -1.050275319741563e-02, -5.249774789586573e-02, 1.042412741554585e-02, 2.340133845233861e-03},  {1.280203319332392e+00, 9.652751031396894e-01, 6.944979372074281e-02, -6.944979372074284e-02, 3.472489686037143e-02, -6.944979372074285e-03, 0.000000000000000e+00},  {2.273258339960378e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  7;
    constant Integer[7] order =  {1, 6, 5, 6, 5, 6, 1};
    constant Real[8] breaks = {-7.700000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 1.290000000000000e+02};
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.531143249103176e-15, -1.846954464302997e-18, 1.193892646759935e-18, 2.891205793294678e-19, 1.259130858783489e-03, -5.180079941413973e-04},  {7.411228646436224e-04, 3.187606329067825e-03, 4.821188675713934e-03, 2.231148705006946e-03, -1.474465618203514e-03, 5.457702945229771e-04, 0.000000000000000e+00},  {3.812347939916962e-02, 4.572486927133115e-02, 2.648252963070943e-02, 1.226623554029791e-02, 3.983237327026256e-03, -1.498661180253427e-02, 1.108341345335742e-02},  {1.226771528193575e-01, 1.429890061692127e-01, 1.035657439887793e-01, 1.000013358902086e-01, 9.530138011471621e-02, -5.974965656620094e-02, 0.000000000000000e+00},  {5.047849624160734e-01, 7.325817394453907e-01, 3.778814666856932e-01, -1.162897093129359e-01, -2.034469027162885e-01, 1.976444349669114e-01, -5.231835147455125e-02},  {1.440837640010293e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT28HC</strong>.
  It also contains the phase transition functions for complete melting
  (and solidification), which are modelled by piece-wise splines,
  see
  <blockquote>
  <p>
  Barz, T., Krämer, J., & Emhofer, J. (2020). Identification of Phase
  Fraction–Temperature Curves from Heat Capacity Data for Numerical
  Modeling of Heat Transfer in Commercial Paraffin Waxes.
  Energies, 13(19), 5149.
  <a href>doi.org/10.3390/en13195149</a>.
  </p>
  </blockquote>
  <p>
  </p></html>",
  revisions="<html>
  <ul>
  <li>01-Jun-2022  </ul>
  </html>"));
end RT28HC;
