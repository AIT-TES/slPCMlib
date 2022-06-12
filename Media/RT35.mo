within slPCMlib.Media;
package RT35 "Rubitherm RT35, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT35";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {3.001500000000000e+02, 3.101500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.991500000000000e+02, 3.091500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.239109491869517e+05
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
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 6, 5, 5, 5, 6, 1};
    constant Real[13] breaks = {-7.300000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 3.500000000000000e+01, 3.600000000000000e+01, 3.700000000000000e+01, 1.370000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 8.375196103604028e-14, 6.759044425499380e-16, -3.656701316533700e-19, 1.445602896647339e-19, 6.098776540219499e-04, -1.363743070935728e-04},  {4.735033470128048e-04, 2.231142427627160e-03, 4.053161933816532e-03, 3.371290398348045e-03, 1.003773663706158e-03, 1.842998980454681e-05, 0.000000000000000e+00},  {1.115130176031525e-02, 2.455858209415326e-02, 2.037397500914309e-02, 7.570684951218144e-03, 1.095923612728892e-03, -4.959866731486979e-05, 0.000000000000000e+00},  {6.470086876024375e-02, 9.215428808042263e-02, 4.916558486602209e-02, 1.145839272898501e-02, 8.479302761545431e-04, -6.107804580829336e-04, 0.000000000000000e+00},  {2.177162842537451e-01, 2.251984548135920e-01, 8.252054012907478e-02, 8.742309252773849e-03, -2.205972014260124e-03, 5.891302566598043e-05, 0.000000000000000e+00},  {5.320305294605916e-01, 4.079371399013489e-01, 9.610076605849537e-02, 5.075514523931544e-04, -1.911406885930222e-03, 1.277922370166366e-03, -3.555019048233244e-04},  {1.035587000452242e+00, 5.982722992536870e-01, 9.360167422940725e-02, -1.468890486130563e-03, -8.543236074482585e-04, -1.129765879876585e-04, 0.000000000000000e+00},  {1.725024783253769e+00, 7.770867988843726e-01, 8.293929524644938e-02, -6.015950795800182e-03, -1.419206547386551e-03, -5.379440490378245e-05, 0.000000000000000e+00},  {2.577561925636501e+00, 9.189717387757956e-01, 5.583825952569162e-02, -1.223072103438421e-02, -1.688178571905463e-03, 7.048768854602922e-04, 0.000000000000000e+00},  {3.539157901217159e+00, 9.907277648636960e-01, 1.606579384570906e-02, -1.193466646740314e-02, 1.836205855395998e-03, 2.111435255904144e-03, -8.262254756611145e-04},  {4.537138209094800e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 5, 5, 5, 6, 6, 1};
    constant Real[13] breaks = {-7.400000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01, 2.900000000000000e+01, 3.000000000000000e+01, 3.100000000000000e+01, 3.200000000000000e+01, 3.300000000000000e+01, 3.400000000000000e+01, 3.500000000000000e+01, 3.600000000000000e+01, 1.360000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 8.786929941290932e-16, 2.772255214336691e-18, 2.323972458723076e-19, 2.891205793294678e-19, 8.015067885048010e-04, -2.814520657776388e-04},  {5.200547227280441e-04, 2.318821547859274e-03, 3.793286898383431e-03, 2.386026569495233e-03, -2.142470441405780e-04, 1.745914108747436e-04, 0.000000000000000e+00},  {8.978534105200149e-03, 1.707944393092223e-02, 1.141179845077310e-02, 3.274952501680356e-03, 6.587100102331398e-04, 1.480420010905230e-04, 0.000000000000000e+00},  {4.155148099989951e-02, 5.310294838389594e-02, 2.666933602811824e-02, 7.390212553518146e-03, 1.398920015685755e-03, -3.812490118185106e-04, 0.000000000000000e+00},  {1.297316489692991e-01, 1.323016931043340e-01, 5.342100366460210e-02, 9.173402498076060e-03, -5.073250434067980e-04, -1.345045437314563e-04, 0.000000000000000e+00},  {3.239859186491730e-01, 2.639620850354852e-01, 7.655221546107494e-02, 5.799056887134306e-03, -1.179847762064079e-03, 1.565337306433971e-04, 0.000000000000000e+00},  {6.692759620014469e-01, 4.305269642240038e-01, 8.843563685652736e-02, 2.645003145311961e-03, -3.971791088470938e-04, -1.189793589339082e-04, 0.000000000000000e+00},  {1.190367407759509e+00, 6.131496341429401e-01, 9.279777805004161e-02, -1.335068794154950e-04, -9.920759035166346e-04, -6.482226784602742e-05, 0.000000000000000e+00},  {1.895124414901712e+00, 7.940522546514828e-01, 8.579657931223504e-02, -4.750033171942316e-03, -1.316187242746772e-03, 6.219695783723251e-03, -3.898364489908069e-03},  {2.771228359744556e+00, 9.538388567683246e-01, 6.737084682853892e-02, -2.578511410385827e-02, -2.869317567275155e-02, 3.069007476935873e-02, -8.317146544936139e-03},  {3.760332701789232e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.600000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>RT35</strong>.
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
end RT35;
