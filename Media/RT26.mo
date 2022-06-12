within slPCMlib.Media;
package RT26 "Rubitherm RT26, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT26";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.931500000000000e+02, 3.021500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.931500000000000e+02, 3.001500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.314215812126732e+05
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
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 6, 5, 5, 5, 6, 5, 1};
    constant Real[12] breaks = {-8.000000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 1.290000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 7.545361039186382e-15, 1.770246788197564e-16, 5.770294779206345e-19, -2.891205793294678e-19, 1.369503115559058e-03, -6.548002588672594e-04},  {7.147028566995211e-04, 2.918714024599369e-03, 3.873027272581856e-03, 5.990259782453862e-04, -2.974488305213604e-03, 1.233474954759951e-03, 0.000000000000000e+00},  {6.364456781672480e-03, 6.731268057443984e-03, 1.579249236358872e-04, 1.035822304990479e-03, 3.192886468586150e-03, -6.745763490858778e-04, 0.000000000000000e+00},  {1.680778218724310e-02, 1.955324894860543e-02, 1.567694715926551e-02, 7.061604688476303e-03, -1.799952768432381e-04, 3.668660115383658e-03, -1.662868164856364e-03},  {6.092537965727440e-02, 7.973806781297070e-02, 4.752536824462605e-02, 9.770861437812649e-03, -6.779717172770411e-03, 7.155484987330248e-03, 0.000000000000000e+00},  {1.983354449672436e-01, 2.127599448612343e-01, 1.077144993947441e-01, 5.420684262003349e-02, 2.899770776388083e-02, -1.864345501455196e-02, 0.000000000000000e+00},  {5.833709845925843e-01, 6.135830274935897e-01, 2.578867236926100e-01, -1.623687646996280e-02, -6.421956730887897e-02, 2.285859176234379e-02, 0.000000000000000e+00},  {1.397242883762286e+00, 9.380605350451504e-01, 5.244460805288630e-02, -4.452922808204081e-02, 5.007339150283997e-02, -3.265397101643723e-02, 7.943365987874911e-03},  {2.368581585252560e+00, 9.940459737612255e-01, 1.190805247755507e-02, -1.190805247755500e-02, 5.954026238777503e-03, -1.190805247755501e-03, 0.000000000000000e+00},  {3.367390780004807e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 6, 5, 5, 5, 5, 6, 6, 1};
    constant Real[10] breaks = {-8.000000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 1.270000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.174118090876502e-12, -1.020564022758184e-14, -2.548662049770398e-18, 0.000000000000000e+00, 2.146920402893336e-03, -8.827375738336793e-04},  {1.264182827875331e-03, 5.438176570246651e-03, 8.228140421417753e-03, 3.814452552259772e-03, -2.506461593038510e-03, 7.875436792489751e-04, 0.000000000000000e+00},  {1.702603445800997e-02, 2.724968709394307e-02, 1.250816531245568e-02, 1.664042972595483e-03, 1.431256803206365e-03, -5.110407448380292e-04, 0.000000000000000e+00},  {5.936814589537254e-02, 6.042797012528257e-02, 2.097742760110008e-02, 2.278662737040648e-03, -1.123946920983781e-03, 2.083544345955428e-03, 0.000000000000000e+00},  {1.440118037837675e-01, 1.151407475844164e-01, 4.190517774587334e-02, 1.861831851265979e-02, 9.293774808793358e-03, -5.304491855245208e-03, 0.000000000000000e+00},  {3.236653305802651e-01, 2.654586985731106e-01, 1.004778635841610e-01, 2.748499195381417e-03, -1.722868446743268e-02, 1.413432418520379e-01, -6.881950779079793e-02},  {7.476454415267255e-01, 6.995443479725791e-01, 3.864910560241133e-01, -2.912397596992874e-02, -3.428050920692121e-01, 2.829812664463483e-01, -7.147341601083529e-02},  {1.673259627919790e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.500000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT26</strong>.
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
end RT26;
