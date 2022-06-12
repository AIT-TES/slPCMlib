within slPCMlib.Media;
package RT47 "Rubitherm RT47, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT47";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.121500000000000e+02, 3.231500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {3.121500000000000e+02, 3.221500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.027538889307089e+05
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
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[14] breaks = {-6.100000000000000e+01, 3.900000000000000e+01, 4.000000000000000e+01, 4.100000000000000e+01, 4.200000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 4.600000000000000e+01, 4.700000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 5.000000000000000e+01, 1.500000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 2.831955329168070e-15, 4.685296205997820e-17, -4.918413272072547e-19, 0.000000000000000e+00, 8.758552506314323e-04, -2.563702078362071e-04},  {6.194850427981035e-04, 2.841055006143232e-03, 4.912999388771267e-03, 3.631148349590179e-03, 5.337231356140540e-04, 8.876576093621810e-05, 0.000000000000000e+00},  {1.262717668385305e-02, 2.613922017959149e-02, 1.989644086058827e-02, 6.653698501408575e-03, 9.775519402951446e-04, -5.353423334537842e-04, 0.000000000000000e+00},  {6.575874583228275e-02, 8.712669349889957e-02, 4.036942467204693e-02, 5.210482928051310e-03, -1.699159726973776e-03, 5.825633289140015e-04, 0.000000000000000e+00},  {1.973487505332208e-01, 1.796131693638243e-01, 5.163154838349824e-02, 4.239477309296221e-03, 1.213656917596231e-03, -3.431612750048158e-04, 0.000000000000000e+00},  {4.337034412324311e-01, 2.987335193540726e-01, 6.820030906691617e-02, 5.662492229632988e-03, -5.021494574278477e-04, -1.752605248328604e-04, 0.000000000000000e+00},  {8.056223519007922e-01, 4.492367137229357e-01, 8.042228376291956e-02, 1.901289151592993e-03, -1.378452081592150e-03, 8.446288831996823e-04, -2.110872868536522e-04},  {1.336437728052994e+00, 6.132279610720619e-01, 8.313541825733768e-02, 6.120239201481718e-04, -3.216169683985218e-04, -2.029210089373730e-04, 0.000000000000000e+00},  {2.032888593325207e+00, 7.790337964289015e-01, 8.101257811801735e-02, -2.703654042819645e-03, -1.336222013085387e-03, -5.791606525659566e-04, 0.000000000000000e+00},  {2.888315931163654e+00, 9.247072992213006e-01, 5.909267738538643e-02, -1.384014862082076e-02, -4.232025275915170e-03, 1.940143999840137e-03, 0.000000000000000e+00},  {3.855983877873446e+00, 9.941448270251538e-01, 1.158151986583454e-02, -1.136680972608007e-02, 5.468694723285514e-03, -9.649128608043916e-04, -4.294202795090379e-05},  {4.854804254872884e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 5, 5, 6, 5, 6, 1};
    constant Real[13] breaks = {-6.100000000000000e+01, 3.900000000000000e+01, 4.000000000000000e+01, 4.100000000000000e+01, 4.200000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 4.600000000000000e+01, 4.700000000000000e+01, 4.800000000000000e+01, 4.900000000000000e+01, 1.490000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 2.618134848531676e-14, 2.412113999755596e-16, -1.211737292493671e-18, 2.891205793294678e-19, 8.021885816639241e-04, -2.497692780769274e-04},  {5.524193036134183e-04, 2.512327239891887e-03, 4.275346645485633e-03, 3.026500255100692e-03, 2.644037371657083e-04, 2.273322617890866e-04, 0.000000000000000e+00},  {1.085832944304643e-02, 2.233679755377816e-02, 1.721459245167287e-02, 6.357437821654389e-03, 1.401065046111142e-03, -5.699409444675944e-04, 0.000000000000000e+00},  {5.759828137179539e-02, 7.859285138415617e-02, 3.899388674862660e-02, 6.262288561423010e-03, -1.448639676226830e-03, 3.515978512900692e-04, 0.000000000000000e+00},  {1.803502662410644e-01, 1.713309211172331e-01, 5.260489288843544e-02, 3.983708369416380e-03, 3.093495802235155e-04, 4.037993258156908e-06, 0.000000000000000e+00},  {4.085831761896309e-01, 2.897494202895329e-01, 6.645249541060720e-02, 5.261486622892010e-03, 3.295395465143000e-04, -3.639828028036128e-04, 0.000000000000000e+00},  {7.700121352563737e-01, 4.379371151514531e-01, 8.057436453033281e-02, 2.939816780913084e-03, -1.490374467503764e-03, 5.913246523688127e-04, 0.000000000000000e+00},  {1.290564381903938e+00, 6.049004199466909e-01, 8.636481459173763e-02, 2.891565434586161e-03, 1.466248794340299e-03, 2.922457992480470e-03, -1.862062956106386e-03},  {1.987247825707667e+00, 7.956096528370413e-01, 1.051306392447467e-01, 7.398814146243333e-04, -1.185240558485314e-02, 1.520705566417773e-03, 0.000000000000000e+00},  {2.878396299185644e+00, 9.682844810631487e-01, 5.144290564367921e-02, -3.146268526061051e-02, -4.248877752764280e-03, 1.283790778039458e-02, -3.996044076613906e-03},  {3.871253986582878e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT47</strong>.
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
end RT47;
