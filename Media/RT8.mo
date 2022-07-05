within slPCMlib.Media;
package RT8 "Rubitherm RT8, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT8";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        2.741500000000000e+02,2.841500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        2.731500000000000e+02,2.831500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.091421721879785e+05
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
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 5, 5, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[13] breaks = {-9.900000000000000e+01, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.110000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.835712667563832e-14, 4.036670989619963e-16, -5.849286898929360e-19, 5.782411586589357e-19, 2.432659550145455e-03, -1.001930680124871e-03},  {1.430728870039345e-03, 6.151713669997660e-03, 9.297635299581896e-03, 4.287981898957130e-03, -2.865662451145791e-03, 9.454882282625354e-04, 0.000000000000000e+00},  {1.924788551569278e-02, 3.087572130276229e-02, 1.442248857220389e-02, 2.280214376999317e-03, 1.861778690166886e-03, -3.356199910735921e-04, 0.000000000000000e+00},  {6.835246846675158e-02, 7.233035638346781e-02, 2.907760393346724e-02, 6.371129226930937e-03, 1.836787347989253e-04, 4.454688185951891e-04, 0.000000000000000e+00},  {1.767607055640117e-01, 1.525610109633664e-01, 5.374775220900550e-02, 1.156053235207853e-02, 2.411022827774871e-03, -1.389373176665475e-03, 0.000000000000000e+00},  {3.956516507395715e-01, 2.974353378653943e-01, 8.900175446523576e-02, 7.310891896523270e-03, -4.535843055552501e-03, 1.401522162428473e-03, 0.000000000000000e+00},  {7.862653140736007e-01, 4.862357610753648e-01, 9.773459344577523e-02, 3.182741298597989e-03, 2.471767756589864e-03, 1.406889192805683e-03, -1.333315501679656e-03},  {1.375963751341055e+00, 7.001747958430165e-01, 1.161825832839703e-01, 4.723942194211390e-04, -1.049351880457657e-02, 2.017170895429473e-03, 0.000000000000000e+00},  {2.184317176778315e+00, 9.020689243280554e-01, 7.481036206906889e-02, -2.132997204459042e-02, -4.076643274292058e-04, 8.824716129389317e-04, 0.000000000000000e+00},  {3.140341298416359e+00, 9.904814330873881e-01, 1.719917610011148e-02, -1.413591322491792e-02, 4.004693737265453e-03, 1.037018977663014e-03, -6.126525750387015e-04},  {4.138315054518829e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  12;
    constant Integer[12] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[13] breaks = {-1.000000000000000e+02, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.641473909444156e-14, 2.509022158693751e-16, -4.786554619532341e-18, 1.156482317317871e-18, 3.845154238590491e-03, -1.680681911499732e-03},  {2.164472327107420e-03, 9.141679723979187e-03, 1.324131371340929e-02, 4.837904155910250e-03, -5.984457479543532e-03, 4.142941282860520e-03, -1.012310420219809e-03},  {2.653154330350333e-02, 4.084103359334667e-02, 1.809303782918703e-02, 2.083278661945125e-03, -4.544073685380766e-04, 5.876767305937309e-04, 0.000000000000000e+00},  {8.768216275003782e-02, 8.439769941637401e-02, 2.749319690973129e-02, 6.142416493730127e-03, 2.483976284430577e-03, -1.029984275580370e-03, 0.000000000000000e+00},  {2.071694675787235e-01, 1.625973264768472e-01, 5.052446134170143e-02, 5.778478875648734e-03, -2.665945093471273e-03, 9.110779822445350e-04, 0.000000000000000e+00},  {4.243148671616941e-01, 2.748732953245343e-01, 6.097500723026535e-02, 4.225478324208991e-03, 1.889444817751402e-03, -8.805507946402629e-04, 0.000000000000000e+00},  {7.653975420638139e-01, 4.126547700554951e-01, 7.618260316299807e-02, 2.977749648811968e-03, -2.513309155449913e-03, 4.122730640180192e-04, 0.000000000000000e+00},  {1.255111628839687e+00, 5.659613540262174e-01, 7.415872781691470e-02, -2.952756332807486e-03, -4.519438353598170e-04, 2.258092278143227e-03, 0.000000000000000e+00},  {1.894085102792795e+00, 7.149032267109059e-01, 8.516971858776567e-02, 1.782039110718552e-02, 1.083851755535632e-02, -1.020874022632444e-02, 0.000000000000000e+00},  {2.712608216527685e+00, 9.310142062977943e-01, 1.015745949782157e-01, -4.091294093463347e-02, -4.020518357626587e-02, 4.443802914140274e-02, -1.213233080871652e-02},  {3.696384591625481e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT8</strong>.
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
end RT8;
