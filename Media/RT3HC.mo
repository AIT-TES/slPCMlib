within slPCMlib.Media;
package RT3HC "Rubitherm RT3HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT3HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]={
        2.691500000000000e+02,2.791500000000000e+02}
      "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]={
        2.691500000000000e+02,2.771500000000000e+02}
      "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef={
        2.000000000000000e+03,0.0}
      "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef={
        2.000000000000000e+03,0.0}
      "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy phTrEnth=1.445107429767713e+05
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
    constant Integer[12] order =  {1, 6, 6, 5, 5, 5, 5, 5, 5, 5, 6, 1};
    constant Real[13] breaks = {-1.040000000000000e+02, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 1.060000000000000e+02};
    constant Real[12,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -4.787978431753636e-14, -4.141541947960278e-16, 2.975266023813900e-18, -5.782411586589357e-19, 1.659944632637339e-03, -7.584728906590663e-04},  {9.014717419299810e-04, 3.748885819181120e-03, 5.222352966486960e-03, 1.429988513192066e-03, -3.077370196699303e-03, 2.319211096621740e-03, -5.951641240429145e-04},  {9.949375816669650e-03, 1.419914724377401e-02, 5.312746431440918e-03, 4.093362117539564e-04, -4.087765742343250e-04, 4.793721589038088e-04, 0.000000000000000e+00},  {2.994120128830802e-02, 2.681440323949634e-02, 8.881817210334892e-03, 3.567951503854743e-03, 1.988084220284719e-03, -5.039375948972187e-04, 0.000000000000000e+00},  {7.068951986738149e-02, 6.071454107839670e-02, 2.647480109463536e-02, 6.480912436021430e-03, -5.316037542013747e-04, 2.946145062702801e-03, 0.000000000000000e+00},  {1.667743157849364e-01, 1.457111908724604e-01, 7.218936650451956e-02, 3.381594804624394e-02, 1.419912155931263e-02, -7.362493266119716e-03, 0.000000000000000e+00},  {4.253274495013532e-01, 4.115217879269283e-01, 1.852070073379304e-01, 1.698750162229729e-02, -2.261334477128595e-02, 3.931071898307392e-03, 0.000000000000000e+00},  {1.020361473515531e+00, 7.621002878760776e-01, 1.398001625601805e-01, -3.415515847977258e-02, -2.957985279748990e-03, 2.487219685434370e-03, 0.000000000000000e+00},  {1.887635999877702e+00, 9.398392948653168e-01, 4.445897229671270e-02, -2.111490274442484e-02, 9.478113147422860e-03, -2.421912695443262e-03, 0.000000000000000e+00},  {2.857875564747286e+00, 9.912154203379465e-01, 1.376381599354277e-02, -7.421577109166031e-03, -2.631450329793450e-03, 4.331633396584570e-03, -1.268447776875293e-03},  {3.855864959259525e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[11] breaks = {-1.040000000000000e+02, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, -1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 1.040000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.585463599063179e-13, -1.178441953197812e-15, 1.922961145830211e-18, 0.000000000000000e+00, 1.328294887183507e-03, -5.121219948661642e-04},  {8.161728921576201e-04, 3.568742466558090e-03, 5.601118948841425e-03, 3.040508974511790e-03, -1.040355487074928e-03, 9.607031560732252e-05, 0.000000000000000e+00},  {1.208225811060132e-02, 2.021143691750852e-02, 9.441216106000421e-03, -1.602098177146947e-04, -5.600039090383150e-04, 7.705813431150362e-04, 0.000000000000000e+00},  {4.178527875047228e-02, 4.022613075578489e-02, 1.330637662977679e-02, 5.305587977282403e-03, 3.292902806536865e-03, -9.279885233835600e-04, 0.000000000000000e+00},  {1.029882883964697e-01, 9.128731655642559e-02, 3.970067216700967e-02, 9.197313969594254e-03, -1.347039810380935e-03, 7.537108073729103e-03, -2.796041525951636e-03},  {2.465676178268957e-01, 2.138017347706247e-01, 9.264083306152318e-02, 2.325940494632880e-02, -5.602122331009974e-03, 3.631022574303562e-03, 0.000000000000000e+00},  {5.742984908486660e-01, 4.646082392801462e-01, 1.651165396574854e-01, 3.716114136532452e-02, 1.255299054050783e-02, -1.651768509426819e-02, 0.000000000000000e+00},  {1.237219716597862e+00, 8.739482793818570e-01, 1.867410560538246e-01, -7.780374741532609e-02, -7.003543493083314e-02, 7.936947216926435e-02, -2.178746172769924e-02},  {2.207651880128950e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT3HC</strong>.
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
end RT3HC;
