within slPCMlib.Media_Rubitherm_RT;
package RT10 "Rubitherm RT10; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT10";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.751500000000000e+02, 2.861500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.751500000000000e+02, 2.841500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.281475284698701e+05
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
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[14] breaks = {-9.800000000000000e+01, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.130000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 7.809405953495779e-15, 7.617593221210582e-18, -2.262331057677846e-19, -1.156482317317871e-18, 4.091923620565962e-03, -1.556700101917943e-03},  {2.535223518655835e-03, 1.111941749132584e-02, 1.756873467689048e-02, 9.785234167300757e-03, -2.890883425939334e-03, 6.685894766854380e-04, 5.784647536078798e-05},  {3.884416238027979e-02, 6.773908187884728e-02, 3.713272852042294e-02, 6.064524737613560e-03, 1.319761087899675e-03, 1.296099649700389e-04, 0.000000000000000e+00},  {1.512298685700333e-01, 1.661252073089879e-01, 6.454096891036207e-02, 1.263966873891265e-02, 1.967810912749870e-03, -1.835022065510217e-03, 0.000000000000000e+00},  {3.946685023755355e-01, 3.318222846698919e-01, 9.591661994849707e-02, 2.160691734809957e-03, -7.207299414801215e-03, 2.133043590035307e-03, 0.000000000000000e+00},  {8.194938429039685e-01, 5.119736200622793e-01, 8.048533456447271e-02, -5.338070024041832e-03, 3.457918535375319e-03, -8.484304051325250e-04, 0.000000000000000e+00},  {1.409224215636922e+00, 6.665196012349495e-01, 7.673433165327388e-02, 9.300066134192739e-06, -7.842334902873056e-04, -6.884449315962432e-04, 0.000000000000000e+00},  {2.151014770169395e+00, 8.134370061207601e-01, 6.517238159399019e-02, -1.001208321097746e-02, -4.226458148268521e-03, 2.037051103119413e-03, 0.000000000000000e+00},  {3.017422667628019e+00, 9.070249425983283e-01, 3.014789410264080e-02, -6.547404772857413e-03, 5.958797367328544e-03, -1.995915251123099e-03, 0.000000000000000e+00},  {3.952010981672336e+00, 9.615341296987358e-01, 2.629931147680884e-02, -2.671367814774232e-03, -4.020778888286953e-03, 1.266398974580186e-03, 0.000000000000000e+00},  {4.934418675119399e+00, 9.963675285277839e-01, 6.824524448566290e-03, -6.090493622120178e-03, 2.311215984613978e-03, -2.182470105512859e-05, -1.468061652892223e-04},  {5.933662819591898e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 6, 1};
    constant Real[12] breaks = {-9.800000000000000e+01, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.110000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 2.758757886337934e-15, -2.039835191391896e-18, -1.633019570838444e-18, 5.782411586589357e-19, 1.933648699960442e-03, -5.500095603714484e-04},  {1.383639139591749e-03, 6.368186137577340e-03, 1.108634359403269e-02, 8.336295792175451e-03, 1.418100094230483e-03, -2.500697534335132e-04, 0.000000000000000e+00},  {2.834249500417420e-02, 5.797181231193798e-02, 4.210313400160679e-02, 1.150799863476225e-02, 1.677513270629167e-04, -4.562156635276227e-04, 0.000000000000000e+00},  {1.396369756160165e-01, 1.750920032100561e-01, 7.307148123299481e-02, 7.616847307737686e-03, -2.113326990575197e-03, -1.415918690204291e-04, 0.000000000000000e+00},  {3.931623885072094e-01, 3.349242402918600e-01, 8.182614252255239e-02, -2.252379344767384e-03, -2.821286335677342e-03, 1.469891548879946e-03, 0.000000000000000e+00},  {8.063089971900570e-01, 4.878836997043435e-01, 7.284020196298566e-02, 1.161390801322712e-03, 4.528171408722389e-03, -1.746116811859302e-03, 0.000000000000000e+00},  {1.370976344255572e+00, 6.464303776099094e-01, 8.603223470069510e-02, 1.812908317619261e-03, -4.202412650574119e-03, -3.382279833273593e-03, 2.382398889042997e-03},  {2.100049571288991e+00, 8.045069155297188e-01, 6.816966875301721e-02, -1.171562836553194e-03, 1.462217151870288e-02, -8.432862191055506e-03, 0.000000000000000e+00},  {2.977743902062821e+00, 9.536559396456638e-01, 6.805938744501983e-02, -2.701149867229675e-02, -2.754213943657465e-02, 3.013716115094874e-02, -8.209577754544600e-03},  {3.966833174441038e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT10</strong>.<br><br>
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
end RT10;
