within slPCMlib.Media_Rubitherm_RT;
package RT31 "Rubitherm RT31, data taken from data_sheet"
   extends slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT31";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        2.961500000000000e+02,3.091500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        2.961500000000000e+02,3.071500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.314279314663342e+05
      "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature Tref=273.15 + 25
      "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy href=0.0
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
    (xi, dxi) :=BasicUtilities.quartQuintSplineEval(
        T - 273.15,
        pieces,
        order,
        breaks,
        coefs[:, :]);
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
    (xi, dxi) :=BasicUtilities.quartQuintSplineEval(
        T - 273.15,
        pieces,
        order,
        breaks,
        coefs[:, :]);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  15;
    constant Integer[15] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[16] breaks = {-7.700000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 3.500000000000000e+01, 3.600000000000000e+01, 1.360000000000000e+02};
    constant Real[15,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 6.553354657205833e-14, 2.493136194144738e-17, 1.519645716719446e-18, -1.445602896647339e-19, 2.128774398602879e-03, -9.371691867991833e-04},  {1.191605211869255e-03, 5.020856872295501e-03, 7.230206184041063e-03, 2.544360250045116e-03, -3.413665808973358e-03, 9.049479569923476e-04, 0.000000000000000e+00},  {1.347831066626992e-02, 1.798442653958119e-02, 3.430771650259737e-03, -2.060823415924841e-03, 1.111073975988380e-03, 3.322973950208230e-04, 0.000000000000000e+00},  {3.427605681119521e-02, 2.476928247136820e-02, 7.237719208623723e-03, 5.706446438236906e-03, 2.772560951092495e-03, -1.471943422829922e-03, 0.000000000000000e+00},  {7.329012245768661e-02, 6.009458689356712e-02, 2.627299000159018e-02, 2.077256014307661e-03, -4.587156163057118e-03, 2.131060332219163e-03, 0.000000000000000e+00},  {1.592788595363136e-01, 1.111790119485373e-01, 2.629242438836209e-02, 5.039234684270815e-03, 6.068145498038697e-03, -1.345048135213742e-03, 0.000000000000000e+00},  {3.065126279203088e-01, 1.964289060941473e-01, 6.436852007726930e-02, 1.586133532428818e-02, -6.570951780300120e-04, -3.153950123957094e-03, 9.727287214178167e-04},  {5.803330728354442e-01, 3.601881932180794e-01, 9.106138456365007e-02, 1.148027800953535e-03, -1.835914976548227e-03, 7.640662482082493e-04, 0.000000000000000e+00},  {1.031658829689787e+00, 5.422317170830833e-01, 9.113064058930380e-02, 1.445030376843120e-03, 1.984416264493019e-03, -1.617065679803630e-03, 0.000000000000000e+00},  {1.666833568323707e+00, 7.286804260512144e-01, 9.120157250875498e-02, -6.787961363221105e-03, -6.100912134525131e-03, 1.898166568052661e-03, 0.000000000000000e+00},  {2.475724859953983e+00, 8.758068712811647e-01, 5.321388129246747e-02, -1.220994422079503e-02, 3.389920705738170e-03, -6.505414835032659e-04, 0.000000000000000e+00},  {3.395275047529055e+00, 9.559117766091263e-01, 3.041815802947875e-02, -5.155676232875004e-03, 1.372132882218411e-04, -2.366710588779653e-03, 1.033332165582740e-03},  {4.375253140799810e+00, 9.961963571719418e-01, 7.607285656129363e-03, -7.607285656129363e-03, 3.803642828064682e-03, -7.607285656129364e-04, 0.000000000000000e+00},  {5.374492412234204e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 5, 5, 5, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[14] breaks = {-7.700000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 1.340000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 6.670490802923986e-14, 6.565783520778227e-16, -7.642499098570628e-20, -1.807003620809174e-20, -6.445201727569402e-05, 4.292555080145325e-05},  {-2.152646640687938e-05, -6.470678150217405e-05, -6.369107344894920e-07, 2.139908432721268e-04, 3.216231756433288e-04, 6.427097159275156e-04, 0.000000000000000e+00},  {1.091453576199428e-03, 5.076033209057201e-03, 8.998171832217031e-03, 7.927580705120599e-03, 3.535171755280907e-03, -1.768576859409423e-03, 0.000000000000000e+00},  {2.485983421846574e-02, 5.215292171291359e-02, 3.630617588516988e-02, 4.382499132149996e-03, -5.307712541766210e-03, 1.942956154028865e-03, 0.000000000000000e+00},  {1.143366745609619e-01, 1.263967014828016e-01, 3.703695957131144e-02, 2.581210505373812e-03, 4.407068228378117e-03, -1.089026864109300e-03, 0.000000000000000e+00},  {2.836695874847175e-01, 2.203973907345180e-01, 6.033273181660864e-02, 9.319214777793280e-03, -1.038066092168383e-03, -2.422723952198790e-04, 0.000000000000000e+00},  {5.724385863262491e-01, 3.636568723563433e-01, 7.963925564477942e-02, 2.744226456920957e-03, -2.249428068267778e-03, 6.857893042737777e-04, 0.000000000000000e+00},  {1.016915302020299e+00, 5.255992972649781e-01, 8.123325964867353e-02, 6.044072265876314e-04, 1.179518453101111e-03, 9.744346821983579e-04, -8.332250981545278e-04},  {1.625672994197683e+00, 6.944699348765604e-01, 8.736956239670879e-02, -1.597674102114903e-03, -6.446684608225017e-03, 2.263296250443235e-03, 0.000000000000000e+00},  {2.401731429011055e+00, 8.499457801829647e-01, 6.652939494544646e-02, -4.751450030582626e-03, 4.869796643991156e-03, -3.615174205990336e-03, 0.000000000000000e+00},  {3.314709776546884e+00, 9.701535355280995e-01, 4.534208265774192e-02, -2.142400551452136e-02, -1.320607438596052e-02, 1.699206116312483e-02, -4.783615428644240e-03},  {4.307783760566725e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT31</strong>.
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
end RT31;
