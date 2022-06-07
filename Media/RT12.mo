within slPCMlib.Media;
package RT12 "Rubitherm RT12, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT12";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.771500000000000e+02, 2.891500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.761500000000000e+02, 2.871500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.102483641764689e+05
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
    constant Integer pieces =  14;
    constant Integer[14] order =  {1, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 1};
    constant Real[15] breaks = {-9.600000000000000e+01, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.500000000000000e+01, 1.600000000000000e+01, 1.160000000000000e+02};
    constant Real[14,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.745879713567152e-13, -1.245315399497431e-15, 8.023357469919646e-18, 0.000000000000000e+00, 4.616584304644744e-03, -1.561838884948557e-03},  {3.054745419420362e-03, 1.371188821331673e-02, 2.273825977221815e-02, 1.492906534747632e-02, -3.446617510046263e-04, -6.048109942027332e-04, 0.000000000000000e+00},  {5.348448600722419e-02, 9.957290182515263e-02, 5.940937536659202e-02, 7.502308401430487e-03, -3.368716722018292e-03, 7.372352962149856e-04, 0.000000000000000e+00},  {2.173375901745960e-01, 2.311098873556099e-01, 6.907635320092351e-02, 1.399794475507176e-03, 3.174597590566363e-04, -9.196564368565102e-05, 0.000000000000000e+00},  {5.191491193220077e-01, 3.742719880017791e-01, 7.426083874492834e-02, 1.749977074877209e-03, -1.423684593716189e-04, -3.709327194404918e-04, 0.000000000000000e+00},  {9.689186219647803e-01, 5.256194592815779e-01, 7.494723201892534e-02, -2.528823957014181e-03, -1.997032056574078e-03, 6.587126597015912e-04, 0.000000000000000e+00},  {1.565618169911397e+00, 6.632328865205985e-01, 6.196569440545424e-02, -3.929825586294578e-03, 1.296531241933879e-03, -4.538323358038228e-04, 0.000000000000000e+00},  {2.287729624157285e+00, 7.782917618613363e-01, 5.341708174013554e-02, -3.282023976597292e-03, -9.726304370852352e-04, 2.504106298380233e-04, 0.000000000000000e+00},  {3.115434223974913e+00, 8.726413848126671e-01, 4.023933348621250e-02, -4.668439426558000e-03, 2.794227121048812e-04, -9.725274777071377e-05, 0.000000000000000e+00},  {4.023828672811569e+00, 9.397461606150185e-01, 2.693802400146080e-02, -4.523276055845613e-03, -2.068410267486876e-04, 1.387421492085183e-04, 0.000000000000000e+00},  {4.985921482494662e+00, 9.799187270894576e-01, 1.351457116551705e-02, -3.963218670755181e-03, 4.868697192939040e-04, 1.914700691589952e-04, -1.292249262118224e-04},  {5.975940676941123e+00, 9.971876930739256e-01, 4.522460267427545e-03, -2.685537626226061e-03, -4.941538280884558e-04, 1.200984350338583e-03, -3.673845282402974e-04},  {6.975304738650260e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  13;
    constant Integer[13] order =  {1, 6, 5, 5, 5, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[14] breaks = {-9.700000000000000e+01, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.140000000000000e+02};
    constant Real[13,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.895356452973497e-16, -4.118918141340779e-19, 1.372972713780260e-19, 0.000000000000000e+00, 2.202887474206507e-03, -8.899037508484347e-04},  {1.312983723358261e-03, 5.675014865942241e-03, 8.680318479338547e-03, 4.230799725096372e-03, -2.334118891693988e-03, 1.429074247688963e-03, 0.000000000000000e+00},  {1.899407214973040e-02, 3.353694667157645e-02, 2.165874678135336e-02, 9.185066635210041e-03, 4.811252346750824e-03, -1.748262768028254e-03, 0.000000000000000e+00},  {8.643782181659282e-02, 1.149133356867760e-01, 6.059883308720589e-02, 1.094744834193080e-02, -3.930061493390445e-03, 3.472636079936632e-04, 0.000000000000000e+00},  {2.693146410471087e-01, 2.549694189533852e-01, 7.333344523259225e-02, -1.300161551694351e-03, -2.193743453422129e-03, 4.337947696516171e-04, 2.488801731046236e-04},  {5.948062751707259e-01, 3.926231058366806e-01, 6.434165015006195e-02, -7.595842067742261e-04, 3.708432991405310e-03, -1.157355228653175e-03, 0.000000000000000e+00},  {1.053562524713446e+00, 5.280746093388388e-01, 7.273994319163937e-02, 2.500595472315258e-03, -2.078343151860567e-03, 3.072252949515734e-05, 0.000000000000000e+00},  {1.654830052093874e+00, 6.728965221790906e-01, 6.807889599237332e-02, -5.505551840175437e-03, -1.924730504384781e-03, 1.910835589567896e-03, 0.000000000000000e+00},  {2.390286023510345e+00, 7.943929145736189e-01, 5.912221334121728e-02, 5.903882037964399e-03, 7.629447443454698e-03, -5.509758635722727e-03, 0.000000000000000e+00},  {3.251824722270878e+00, 9.333179839651493e-01, 6.751295775861141e-02, -1.867591454544407e-02, -1.991934573515894e-02, 1.827496091754598e-02, -4.546144587890436e-03},  {4.227789220043693e+00, 9.967367099657850e-01, 6.526580068428793e-03, -6.526580068428793e-03, 3.263290034214397e-03, -6.526580068428794e-04, 0.000000000000000e+00},  {5.227136562036850e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT12</strong>.
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
end RT12;
