within slPCMlib.Media;
package RT5HC "Rubitherm RT5HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT5HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        2.741500000000000e+02,2.801500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        2.741500000000000e+02,2.791500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.930316874611757e+05
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
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 6, 5, 5, 6, 5, 6, 1};
    constant Real[9] breaks = {-9.900000000000000e+01, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 1.070000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.040210351808191e-12, 4.907662778061449e-15, -5.230132868042571e-19, -7.228014483236696e-20, 3.135633689174258e-04, 9.382357660847924e-05},  {4.073869465710224e-04, 2.130758305296816e-03, 4.542987338306393e-03, 5.012105221343842e-03, 2.975170493714317e-03, -1.721650750205493e-03, 0.000000000000000e+00},  {1.334675755502690e-02, 2.954547686971473e-02, 2.021381846256863e-02, -3.037203058538207e-04, -5.633083257313148e-03, 6.083578398421081e-03, 0.000000000000000e+00},  {6.325282772256438e-02, 7.694751184026662e-02, 4.633994198533967e-02, 3.799973064910438e-02, 2.478480873479226e-02, 2.017770627127395e-02, -1.411411579988374e-02},  {2.553884114034575e-01, 3.989696592540999e-01, 2.991133120558879e-01, 5.663371230333811e-02, -8.603839690709410e-02, 1.881780474571047e-02, 0.000000000000000e+00},  {9.428845028553997e-01, 9.170328563762042e-01, 1.409621149804430e-01, -9.934182786793358e-02, 8.050626821458249e-03, 2.336204690321348e-02, -8.324057422501710e-03},  {1.924626262646283e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  7;
    constant Integer[7] order =  {1, 6, 5, 6, 5, 6, 1};
    constant Real[8] breaks = {-9.900000000000000e+01, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 1.060000000000000e+02};
    constant Real[7,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.895899314635641e-12, -9.090063069783639e-15, 1.626598647730058e-18, -1.156482317317871e-18, 4.539157485430136e-03, -1.972389351394549e-03},  {2.566768132130598e-03, 1.086145131669953e-02, 1.580573458337322e-02, 5.943787826410374e-03, -6.890052843767560e-03, 2.916331826097704e-03, 0.000000000000000e+00},  {3.120402084094387e-02, 4.732573171807752e-02, 2.146009926097593e-02, 7.546894712317163e-03, 7.691606286720956e-03, -2.417833233198080e-02, 1.539020426549288e-02},  {1.064402247525475e-01, 1.151026034564021e-01, 7.932016178083579e-02, 1.043340818492504e-01, 1.176530086092101e-01, -6.927708806740693e-02, 0.000000000000000e+00},  {4.535729923808389e-01, 7.109717666675283e-01, 4.054695783097874e-01, -1.178247643879784e-01, -2.287324317278245e-01, 2.183333746986532e-01, -5.752896278436274e-02},  {1.384261553156642e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT5HC</strong>.
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
end RT5HC;
