within slPCMlib.Media_Rubitherm_RT;
package RT_minus_4 "Rubitherm RT_minus_4; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT_minus_4";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.631500000000000e+02, 2.711500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.631500000000000e+02, 2.711500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.501818667202909e+05
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
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[11] breaks = {-1.100000000000000e+02, -1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, -7.000000000000000e+00, -6.000000000000000e+00, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, 9.800000000000000e+01};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.188875207306421e-14, -1.457483020492170e-17, 2.545312100338156e-18, 0.000000000000000e+00, 2.333357041125656e-03, -8.894310103048620e-04},  {1.443926030798893e-03, 6.330199143776314e-03, 9.992105256683616e-03, 5.544950205159310e-03, -1.674679948944653e-03, 2.077697218127126e-03, -7.128844969555208e-04},  {2.300131340864508e-02, 4.236171958574014e-02, 2.666258090543207e-02, 5.365512651541541e-03, -1.979461312641836e-03, 5.811007729066834e-04, 0.000000000000000e+00},  {9.599276601162371e-02, 1.067710779652001e-01, 3.669335871327252e-02, 3.258675130041033e-03, 9.260425518915816e-04, 4.702906650666134e-04, 0.000000000000000e+00},  {2.441122110370955e-01, 1.959894443147712e-01, 5.672854606541124e-02, 1.166575198827349e-02, 3.277495877224648e-03, 2.535691183611608e-03, 0.000000000000000e+00},  {5.143091404663878e-01, 3.702322318373799e-01, 1.367476891296957e-01, 5.013264733328816e-02, 1.595595179528269e-02, -1.659610332401558e-02, 0.000000000000000e+00},  {1.070781557238019e+00, 7.749688426575595e-01, 2.169203086611005e-01, -5.200457872573692e-02, -6.702456482479523e-02, 5.723954618627395e-02, -1.381277912101941e-02},  {1.987068332071401e+00, 9.880185207087174e-01, 2.396295858256662e-02, -2.396295858256663e-02, 1.198147929128331e-02, -2.396295858256663e-03, 0.000000000000000e+00},  {2.984672036213145e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[11] breaks = {-1.100000000000000e+02, -1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, -7.000000000000000e+00, -6.000000000000000e+00, -5.000000000000000e+00, -4.000000000000000e+00, -3.000000000000000e+00, -2.000000000000000e+00, 9.800000000000000e+01};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -2.188875207306421e-14, -1.457483020492170e-17, 2.545312100338156e-18, 0.000000000000000e+00, 2.333357041125656e-03, -8.894310103048620e-04},  {1.443926030798893e-03, 6.330199143776314e-03, 9.992105256683616e-03, 5.544950205159310e-03, -1.674679948944653e-03, 2.077697218127126e-03, -7.128844969555208e-04},  {2.300131340864508e-02, 4.236171958574014e-02, 2.666258090543207e-02, 5.365512651541541e-03, -1.979461312641836e-03, 5.811007729066834e-04, 0.000000000000000e+00},  {9.599276601162371e-02, 1.067710779652001e-01, 3.669335871327252e-02, 3.258675130041033e-03, 9.260425518915816e-04, 4.702906650666134e-04, 0.000000000000000e+00},  {2.441122110370955e-01, 1.959894443147712e-01, 5.672854606541124e-02, 1.166575198827349e-02, 3.277495877224648e-03, 2.535691183611608e-03, 0.000000000000000e+00},  {5.143091404663878e-01, 3.702322318373799e-01, 1.367476891296957e-01, 5.013264733328816e-02, 1.595595179528269e-02, -1.659610332401558e-02, 0.000000000000000e+00},  {1.070781557238019e+00, 7.749688426575595e-01, 2.169203086611005e-01, -5.200457872573692e-02, -6.702456482479523e-02, 5.723954618627395e-02, -1.381277912101941e-02},  {1.987068332071401e+00, 9.880185207087174e-01, 2.396295858256662e-02, -2.396295858256663e-02, 1.198147929128331e-02, -2.396295858256663e-03, 0.000000000000000e+00},  {2.984672036213145e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT_minus_4</strong>.<br><br>
  Information taken from: data_sheet - last access 01.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  false</li>
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
end RT_minus_4;
