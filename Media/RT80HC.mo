within slPCMlib.Media;
package RT80HC "Rubitherm RT80HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT80HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.471500000000000e+02, 3.541500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {3.461500000000000e+02, 3.521500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.753526313167934e+05
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
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 6, 5, 6, 5, 5, 6, 5, 1};
    constant Real[10] breaks = {-2.600000000000000e+01, 7.400000000000000e+01, 7.500000000000000e+01, 7.600000000000000e+01, 7.700000000000000e+01, 7.800000000000000e+01, 7.900000000000000e+01, 8.000000000000000e+01, 8.100000000000000e+01, 1.810000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.745723405613770e-13, 1.393770227173784e-15, 2.145013632201093e-19, -2.312964634635743e-18, 6.351107318643658e-03, -2.666922511583314e-03},  {3.684184807236308e-03, 1.575400152389592e-02, 2.350723551268825e-02, 1.017262295477029e-02, -8.248301080531428e-03, 2.513064898794242e-03, 0.000000000000000e+00},  {4.738280861685358e-02, 7.285846158540611e-02, 2.966594688175277e-02, 2.310067620587003e-03, 4.317023413439781e-03, -4.842366305654340e-03, 3.362851386966313e-03},  {1.550547931993512e-01, 1.523539286579923e-01, 6.451739797210404e-02, 3.841152595712899e-02, 3.054796268966278e-02, -1.336002463409647e-02, 0.000000000000000e+00},  {4.275255838421428e-01, 4.520150300614141e-01, 2.294395056405003e-01, 2.700313037481545e-02, -3.625216048081955e-02, 5.146840382100800e-03, 0.000000000000000e+00},  {1.104877929820154e+00, 8.726289924540186e-01, 1.444043377010368e-01, -6.653710772745473e-02, -1.051795857031555e-02, 2.269537958456591e-02, -6.485254650839403e-03},  {2.061066318611165e+00, 9.943198804100523e-01, 1.136023917984570e-02, -1.136023917984590e-02, 5.680119589922951e-03, -1.136023917984590e-03, 0.000000000000000e+00},  {3.059930294693155e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 5, 5, 6, 5, 6, 1};
    constant Real[9] breaks = {-2.700000000000000e+01, 7.300000000000000e+01, 7.400000000000000e+01, 7.500000000000000e+01, 7.600000000000000e+01, 7.700000000000000e+01, 7.800000000000000e+01, 7.900000000000000e+01, 1.790000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -8.881802547356324e-12, -2.848854291911812e-14, 1.065324299747723e-18, 1.156482317317871e-18, 3.335915895561276e-03, -1.304631941168768e-03},  {2.031283945482219e-03, 8.851787822325492e-03, 1.378967983805427e-02, 7.266520132237404e-03, -2.889899639725140e-03, 4.359474665644642e-04, 0.000000000000000e+00},  {2.948531956493871e-02, 4.885084666895526e-02, 2.260931706205992e-02, 6.639623898148462e-05, -7.101623069028194e-04, 2.451746449340269e-03, 0.000000000000000e+00},  {1.027534636773728e-01, 1.036867525291871e-01, 4.306499643099040e-02, 2.174321150477296e-02, 1.154856993979853e-02, -2.609452145523439e-02, 1.328819246017028e-02},  {2.699906650870577e-01, 2.504972071475397e-01, 1.159637229343043e-01, 7.275612591502884e-02, 8.039884956618083e-02, -4.961837779707749e-02, 0.000000000000000e+00},  {7.399881928530339e-01, 7.741965400384712e-01, 3.204414201056944e-01, -1.018322537910227e-01, -1.676930394192066e-01, 1.647041076726721e-01, -4.372183326294360e-02},  {1.686083134196698e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 9.000000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 8.000000000000000e+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 1.400000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 1.400000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>RT80HC</strong>.
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
end RT80HC;
