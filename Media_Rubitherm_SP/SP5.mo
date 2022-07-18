within slPCMlib.Media_Rubitherm_SP;
package SP5 "Rubitherm SP5; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP5";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.731500000000000e+02+0.000000000000000e+00, 2.731500000000000e+02+1.100000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.731500000000000e+02+ERROR
                                                                                            -1.000000000000000e+00, 2.731500000000000e+02+8.000000000000000e+00}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.700000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+0.000000000000000e+00
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";

  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer len_x =    data_H.len_x;
    constant Real data_x[:] =   data_H.data_x;
    constant Real data_y[:] =   data_H.data_y;
    constant Real m_k[:] =      data_H.m_k;
    constant Real iy_start[:] = data_H.iy_start;
    constant Real iy_scaler =   data_H.iy_scaler;
  algorithm
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15,
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer len_x =    data_C.len_x;
    constant Real data_x[:] =   data_C.data_x;
    constant Real data_y[:] =   data_C.data_y;
    constant Real m_k[:] =      data_C.m_k;
    constant Real iy_start[:] = data_C.iy_start;
    constant Real iy_scaler =   data_C.iy_scaler;
  algorithm
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15,
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    12;
    constant Real[12] data_x =   {0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01};
    constant Real[12] data_y =   {0.000000000000000e+00, 4.668949331757315e-02, 4.814276584174094e-02, 8.076760228282799e-02, 2.544764100653031e-01, 2.157078331792126e-01, 2.068588616972009e-01, 8.879230314426721e-02, 1.374025937784933e-02, 2.973110764930954e-02, 1.509336344472946e-02, 0.000000000000000e+00};
    constant Real[12] m_k =      {4.668949331757315e-02, 3.558520160684479e-03, 2.388988251908355e-03, 9.784534894618843e-02, 0.000000000000000e+00, -9.325403077355016e-03, -2.485509050576917e-02, -9.655930115967576e-02, 0.000000000000000e+00, 0.000000000000000e+00, -1.486555382465477e-02, -1.509336344472946e-02};
    constant Real[12] iy_start = {0.000000000000000e+00, 2.680100752112997e-02, 7.407122400865337e-02, 1.302823045847730e-01, 3.051577312181072e-01, 5.398188004942434e-01, 7.513074242723368e-01, 9.043205585344883e-01, 9.473188515838652e-01, 9.689432005926619e-01, 9.924730869916505e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.948778006123833e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    10;
    constant Real[10] data_x =   {-1.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00, 2.000000000000000e+00, 3.000000000000000e+00, 4.000000000000000e+00, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00};
    constant Real[10] data_y =   {0.000000000000000e+00, 4.300234878828089e-02, 5.560213798666030e-02, 1.547484452913692e-01, 4.409182590340598e-01, 1.910470471628772e-02, 9.175664694926443e-03, 2.568337453611176e-01, 2.061469412712218e-02, 0.000000000000000e+00};
    constant Real[10] m_k =      {4.300234878828089e-02, 1.683872458838713e-02, 3.384153579299901e-02, 1.926580605236998e-01, 0.000000000000000e+00, -2.978712006408383e-02, 0.000000000000000e+00, 0.000000000000000e+00, -6.106230774549692e-02, -9.802300668081847e-03};
    constant Real[10] iy_start = {0.000000000000000e+00, 2.357772528425130e-02, 7.125327683905511e-02, 1.627910563871068e-01, 4.753040685885306e-01, 7.067792300724177e-01, 7.183860801352814e-01, 8.508080770692923e-01, 9.939907631632776e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.956188911220224e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.400000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.350000000000000e+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 6.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 6.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP5</strong>.<br><br>
  Information taken from: data sheet - last access 02.12.2019.<br><br>
  It also contains the phase transition functions for
  <ul>
  <li>complete melting       :  true</li>
  <li>complete solidification:  true</li>
  </ul></p><p>
  These functions are modelled by piece-wise splines using <strong>pchip</strong> method,
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
  <li>file creation date: 18-Jul-2022  </ul>
  </html>"));
end SP5;
