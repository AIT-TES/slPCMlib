within slPCMlib.Media;
package RT50 "Rubitherm RT50, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT50";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.161500000000000e+02, 3.251500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {3.161500000000000e+02, 3.251500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.245030499214933e+05
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
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
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
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-5.700000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 4.600000000000000e+01, 4.700000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 5.000000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 1.520000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 3.801388074488828e-14, 2.835159746775900e-16, -4.117988252282101e-19, 0.000000000000000e+00, 1.283535914855655e-03, -3.179861722597900e-04},  {9.655497426341624e-04, 4.509762540773469e-03, 8.065566564660101e-03, 6.475635703360749e-03, 1.647886990381426e-03, -8.953426994264155e-04, 0.000000000000000e+00},  {2.076905884238349e-02, 4.218263724452865e-02, 2.842636862276645e-02, 4.113756670622306e-03, -2.828826506750651e-03, 1.717530751176897e-03, 0.000000000000000e+00},  {9.438052562472715e-02, 1.086489922308264e-01, 4.096998710589855e-02, 9.973758155388677e-03, 5.758827249133835e-03, -2.742673389225778e-03, 0.000000000000000e+00},  {2.569894169767489e-01, 2.298321829592075e-01, 7.801749117460992e-02, 5.582333259666246e-03, -7.954539696995050e-03, 7.399813812277354e-03, -1.982094177649436e-03},  {5.678846043078656e-01, 3.959025102950011e-01, 9.130397822967080e-02, 8.120429041470844e-03, -6.868833003498309e-04, 1.368074685677885e-04, 0.000000000000000e+00},  {1.062661446042226e+00, 6.008082580201727e-01, 1.129120402376621e-01, 6.740970525749406e-03, -2.845957510888448e-06, -2.247060883690399e-03, 0.000000000000000e+00},  {1.780872807984609e+00, 8.356085618242404e-01, 1.106472672329409e-01, -1.574102214119814e-02, -1.123815037596288e-02, 3.542776242998940e-03, 0.000000000000000e+00},  {2.703692240767628e+00, 9.824413095776574e-01, 3.142306098355849e-02, -2.526586121506028e-02, 6.475730839031816e-03, 2.399173693292630e-03, -1.231439953699664e-03},  {3.699934214692409e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-5.700000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 4.600000000000000e+01, 4.700000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 5.000000000000000e+01, 5.100000000000000e+01, 5.200000000000000e+01, 1.520000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.473917928848470e-13, 6.621302804326147e-16, 1.334511447493068e-18, 0.000000000000000e+00, 2.401632274534671e-03, -1.000314117963546e-03},  {1.401318156719180e-03, 6.006276665063048e-03, 9.011610975894275e-03, 4.010040386075778e-03, -2.996550396779842e-03, 1.054634770462324e-03, 0.000000000000000e+00},  {1.848733055743476e-02, 2.934659204026179e-02, 1.360877745806576e-02, 2.570186503579651e-03, 2.276623455531780e-03, -3.210436296908661e-04, 0.000000000000000e+00},  {6.596846638518288e-02, 7.177598214086632e-02, 3.176864140508701e-02, 8.466244028798108e-03, 6.714053070774489e-04, -6.568752950034046e-04, 0.000000000000000e+00},  {1.779938639720084e-01, 1.601132417907245e-01, 5.462705238391197e-02, 4.583112307073849e-03, -2.612971167939574e-03, 5.314336968071577e-03, -1.878597167273063e-03},  {3.981400390865776e-01, 2.879649006447320e-01, 7.766297446911591e-02, 9.702653970570067e-03, -4.220243836677630e-03, 2.847636509780093e-03, 0.000000000000000e+00},  {7.720979608440981e-01, 4.697560186968709e-01, 1.099258384585613e-01, 2.129804372166047e-02, 1.001793871222283e-02, -7.329624911324233e-03, 0.000000000000000e+00},  {1.375766175522090e+00, 7.569254570713703e-01, 1.606313527836379e-01, -1.192645054269052e-02, -2.663018584439833e-02, 7.974336354521676e-03, 0.000000000000000e+00},  {2.262740685344530e+00, 9.757597494050054e-01, 4.481424963439055e-02, -3.870383037506708e-02, 1.324149592821005e-02, 1.017952369952084e-03, -1.222083851864698e-03},  {3.257648218455156e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT50</strong>.
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
end RT50;
