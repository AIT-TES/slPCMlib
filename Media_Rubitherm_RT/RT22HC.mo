within slPCMlib.Media_Rubitherm_RT;
package RT22HC "Rubitherm RT22HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT22HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.881500000000000e+02, 2.981500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.881500000000000e+02, 2.971500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.483780742635985e+05
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
    constant Integer[12] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 6, 5, 1};
    constant Real[13] breaks = {-8.500000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 1.250000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 8.524271564785069e-14, 5.170470398933253e-16, -4.579742095884125e-19, 0.000000000000000e+00, 4.319125968407920e-04, -6.256892959395687e-05},  {3.693436673325945e-04, 1.784149406725838e-03, 3.380592024499079e-03, 3.067747376528783e-03, 1.221029040294607e-03, -6.312792797619157e-04, 0.000000000000000e+00},  {9.191582235618986e-03, 1.947629534770183e-02, 1.359721559823404e-02, 1.639070740088052e-03, -1.935367358514972e-03, 1.024364830025732e-03, 0.000000000000000e+00},  {4.299316139315367e-02, 4.896829348049861e-02, 1.714587196766567e-02, 4.141249606285487e-03, 3.186456791613689e-03, -7.616530380319266e-04, 0.000000000000000e+00},  {1.156733802011852e-01, 1.046213482110227e-01, 4.107183115588524e-02, 9.270546392420977e-03, -6.218083985459437e-04, -1.740137145479899e-03, 9.448447759202651e-04},  {2.692200051924086e-01, 2.090577990339341e-01, 6.192392012587714e-02, 8.278836861843513e-03, 4.850177512858538e-03, -2.378703665465077e-04, 0.000000000000000e+00},  {5.530928683603754e-01, 3.759535080898656e-01, 1.134827921230935e-01, 2.530084324781258e-02, 3.660825680125999e-03, -4.891315285376591e-03, 0.000000000000000e+00},  {1.066599522215896e+00, 6.690083483729277e-01, 1.624371230935202e-01, -8.969006885449330e-03, -2.079575074675695e-02, 6.218913108456780e-03, 0.000000000000000e+00},  {1.874499149158595e+00, 9.148871364588783e-01, 7.294472904119828e-02, -2.996287878790934e-02, 1.029881479552695e-02, -5.633392475858955e-03, 1.616756790638535e-03},  {2.838650314981068e+00, 9.936167957241937e-01, 1.276640855162042e-02, -1.276640855162040e-02, 6.383204275810200e-03, -1.276640855162040e-03, 0.000000000000000e+00},  {3.837373674125910e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[12] breaks = {-8.500000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.700000000000000e+01, 1.800000000000000e+01, 1.900000000000000e+01, 2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 1.240000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 3.649041819748621e-13, 3.119155242828259e-15, -1.509318440554790e-19, 1.445602896647339e-19, 6.335987574666615e-04, -1.933114983939279e-04},  {4.402872594407569e-04, 2.008124797342884e-03, 3.436315098760833e-03, 2.469757606788057e-03, 2.683213114243890e-04, 1.375381530176118e-04, 0.000000000000000e+00},  {8.760344226774529e-03, 1.805100382599750e-02, 1.383089731784731e-02, 4.918424382661732e-03, 9.560120765124482e-04, -6.447036972699762e-04, 0.000000000000000e+00},  {4.587197813252356e-02, 6.106860142937407e-02, 2.787520595220740e-02, 2.295435716011766e-03, -2.267506409837433e-03, 1.407506845714349e-03, 0.000000000000000e+00},  {1.362512216659937e-01, 1.216728290710553e-01, 3.523154309836168e-02, 7.300478533805521e-03, 4.770027818734312e-03, -2.231114498921882e-03, 0.000000000000000e+00},  {3.029949856890286e-01, 2.219618896495508e-01, 6.344200062296555e-02, 4.069444819523945e-03, -6.385544675875096e-03, 7.846845992507106e-03, -2.043162219763665e-03},  {5.918864598779373e-01, 3.624873032944904e-01, 8.515809365490273e-02, 1.613248164582131e-02, 2.201251990205458e-03, -1.541194717162644e-04, 0.000000000000000e+00},  {1.057711470991641e+00, 5.892353461439163e-01, 1.452218558164361e-01, 2.339629488948050e-02, 1.430654631624136e-03, -8.970123101007045e-03, 0.000000000000000e+00},  {1.808025499372091e+00, 9.107399454666527e-01, 1.342934372645516e-01, -6.058231759409338e-02, -4.341996087341108e-02, 5.291066397695690e-02, -1.474222393409156e-02},  {2.787225043678656e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 7.600000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.000000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT22HC</strong>.<br><br>
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
end RT22HC;
