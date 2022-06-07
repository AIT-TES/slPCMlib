within slPCMlib.Media;
package RT0 "Rubitherm RT0, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT0";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.671500000000000e+02, 2.761500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.661500000000000e+02, 2.741500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.536745580039079e+05
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
    constant Integer[11] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[12] breaks = {-1.060000000000000e+02, -6.000000000000000e+00, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 1.030000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.443037381270635e-13, -1.486784290574459e-15, 1.838578744097348e-18, -2.891205793294678e-19, 8.456124221093658e-04, -3.453521416137809e-04},  {5.002602803497960e-04, 2.155949260715228e-03, 3.275842096885445e-03, 1.549081388818040e-03, -9.522200136598851e-04, 5.658558785840719e-04, -8.257782544456180e-05},  {7.012191066248135e-03, 1.187981000655564e-02, 6.629657585552568e-03, 1.747203611127983e-03, 6.383919975920473e-04, 9.105427331534389e-05, 0.000000000000000e+00},  {2.799830854039171e-02, 3.338957536799292e-02, 1.661216313764227e-02, 5.211314334649610e-03, 1.093663364168767e-03, 7.674742184578017e-04, 0.000000000000000e+00},  {8.507249896330307e-02, 9.045986919617105e-02, 4.648282851118152e-02, 1.726070997590269e-02, 4.931034456457775e-03, -3.429964756663433e-03, 0.000000000000000e+00},  {2.407769763463527e-01, 2.377819701888040e-01, 9.355151761100242e-02, 2.685200235099472e-03, -1.221878932685939e-02, 5.402813128801622e-03, 0.000000000000000e+00},  {5.679796881832008e-01, 4.110795144526838e-01, 8.232251364316079e-02, 7.838174215678133e-03, 1.479527631714872e-02, -4.816697537230210e-03, 0.000000000000000e+00},  {1.079198469274642e+00, 6.343366819685348e-01, 1.464417188207859e-01, 1.885230411197093e-02, -9.288211369002326e-03, -3.684344650920772e-03, 0.000000000000000e+00},  {1.865856618156011e+00, 9.282024632154169e-01, 1.104259164334772e-01, -5.514398787324582e-02, -2.770993462360619e-02, 3.871114406085871e-02, -1.105638571204582e-02},  {2.849285833656865e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 5, 6, 5, 5, 5, 5, 6, 1};
    constant Real[11] breaks = {-1.070000000000000e+02, -7.000000000000000e+00, -6.000000000000000e+00, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 1.010000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.147125376522839e-13, -2.059702921661063e-15, 3.472085124598320e-18, -1.156482317317871e-18, 6.098257784813401e-03, -2.349627463488703e-03},  {3.748630321107928e-03, 1.639352414295474e-02, 2.573816589580177e-02, 1.399002857835994e-02, -4.753123028263544e-03, 5.339396201062829e-06, 0.000000000000000e+00},  {5.512256530616191e-02, 9.085414653760478e-02, 3.924290742331112e-02, -4.969069572683606e-03, -4.726426047258230e-03, 1.623250259123931e-03, 1.310866531474540e-04},  {1.772784605594074e-01, 1.444298196916139e-01, 1.417594481016173e-02, -5.020538107528135e-03, 5.356125045573237e-03, -3.970715172074520e-04, 0.000000000000000e+00},  {3.358227404820207e-01, 1.771592375855820e-01, 2.728036558894196e-02, 1.243324690269029e-02, 3.370767459535977e-03, -1.822988176420006e-03, 0.000000000000000e+00},  {5.542433698423510e-01, 2.733878384275991e-01, 6.657482929002882e-02, 7.686434976634137e-03, -5.744173422564054e-03, 7.057937697163228e-03, 0.000000000000000e+00},  {9.032062368112123e-01, 4.419097967331677e-01, 1.257484706561796e-01, 5.528911825801019e-02, 2.954551506325208e-02, -2.528803995465033e-02, 0.000000000000000e+00},  {1.530411097567171e+00, 8.510159532994733e-01, 2.160085162632210e-01, -7.940922103548474e-02, -9.689468470999955e-02, 1.013385140786451e-01, -2.731985904554838e-02},  {2.495150316417478e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT0</strong>.
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
end RT0;
