within slPCMlib.Media_Rubitherm_SP;
package SP26E "Rubitherm SP26E; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SP26E";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.731500000000000e+02+2.100000000000000e+01, 2.731500000000000e+02+2.800000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.731500000000000e+02+2.000000000000000e+01, 2.731500000000000e+02+2.600000000000000e+01}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.800000000000000e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+2.100000000000000e+01
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
    constant Integer  len_x =    8;
    constant Real[8] data_x =   {2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01, 2.700000000000000e+01, 2.800000000000000e+01};
    constant Real[8] data_y =   {0.000000000000000e+00, 6.353657044215330e-03, 4.451596826734281e-02, 8.927440677221198e-02, 1.216099109567450e-01, 6.687970844266601e-01, 6.944897253282195e-02, 0.000000000000000e+00};
    constant Real[8] m_k =      {5.232060790212295e-03, 1.832883412569078e-02, 4.146037486399833e-02, 1.279208881393659e-02, 9.615937781406425e-02, 0.000000000000000e+00, -2.039939921574014e-01, -4.236613313093360e-02};
    constant Real[8] iy_start = {0.000000000000000e+00, 2.077191530378272e-03, 2.549150254826848e-02, 9.450198272073121e-02, 1.926077294719790e-01, 5.942314616243278e-01, 9.788284791653192e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.960491549307208e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    7;
    constant Real[7] data_x =   {2.000000000000000e+01, 2.100000000000000e+01, 2.200000000000000e+01, 2.300000000000000e+01, 2.400000000000000e+01, 2.500000000000000e+01, 2.600000000000000e+01};
    constant Real[7] data_y =   {0.000000000000000e+00, 1.355154071296485e-02, 4.074644426674210e-02, 8.845292531131370e-02, 2.313434395067131e-01, 6.259056502022700e-01, 0.000000000000000e+00};
    constant Real[7] m_k =      {1.355154071296485e-02, 2.037322213337105e-02, 3.745069229917443e-02, 9.529849761998552e-02, 2.687263624454782e-01, 0.000000000000000e+00, -6.259056502022700e-01};
    constant Real[7] iy_start = {0.000000000000000e+00, 5.893256469187733e-03, 3.031759962328510e-02, 8.707228467915318e-02, 2.251597276539317e-01, 6.533601494718448e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.494078597477356e-01;
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 1.500000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 1.400000000000000e+03;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm
    lambda := 5.000000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm
    lambda := 5.000000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------

annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP26E</strong>.<br><br>
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
end SP26E;
