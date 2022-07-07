within slPCMlib.Media_Rubitherm_SP;
package SP_minus_21 "Rubitherm SP-21; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP_minus_21";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.501500000000000e+02, 2.571500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.471500000000000e+02, 2.521500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.850000000000000e+05
             "scalar phase transition enthalpy";
      
    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 273.15+2.500000000000000e+01
             "reference temperature";
    constant Modelica.Units.SI.SpecificEnthalpy  href = 0.0
             "reference enthalpy at Tref";
      
  end propData;
  // ----------------------------------
  redeclare function extends phaseFrac_complMelting
    "Returns liquid mass phase fraction for complete melting processes"
  protected
    constant Integer len_x    = data_H.len_x;
    constant Real data_x[:]   = data_H.data_x;
    constant Real data_y[:]   = data_H.data_y;
    constant Real m_k[:]      = data_H.m_k;
    constant Real iy_start[:] = data_H.iy_start;
    constant Real iy_scaler   = data_H.iy_scaler;
  algorithm 
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15, 
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complMelting;
  // ----------------------------------
  redeclare function extends phaseFrac_complSolidification
    "Returns liquid mass phase fraction for complete solidification processes"
  protected
    constant Integer len_x    = data_C.len_x;
    constant Real data_x[:]   = data_C.data_x;
    constant Real data_y[:]   = data_C.data_y;
    constant Real m_k[:]      = data_C.m_k;
    constant Real iy_start[:] = data_C.iy_start;
    constant Real iy_scaler   = data_C.iy_scaler;
  algorithm
    (xi, dxi) := slPCMlib.BasicUtilities.cubicHermiteSplineEval(T-273.15, 
                 len_x, data_x, data_y, m_k, iy_start, iy_scaler);
  end phaseFrac_complSolidification;
  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 8; 
    constant Real[8] data_x   = {-2.300000000000000e+01, -2.200000000000000e+01, -2.100000000000000e+01, -2.000000000000000e+01, -1.900000000000000e+01, -1.800000000000000e+01, -1.700000000000000e+01, -1.600000000000000e+01}; 
    constant Real[8] data_y   = {0.000000000000000e+00, 1.984126984126984e-02, 7.341269841269841e-01, 1.150793650793651e-01, 1.071428571428571e-01, 1.587301587301587e-02, 7.936507936507936e-03, 0.000000000000000e+00}; 
    constant Real[8] m_k      = {3.212812977083360e-03, 5.943704007604216e-02, 0.000000000000000e+00, -2.351695777616457e-02, -3.721037622810851e-03, -2.351049123611203e-02, -3.761678597777924e-03, -7.936507936507936e-03}; 
    constant Real[8] iy_start = {0.000000000000000e+00, 5.230423023654640e-03, 3.868131043666715e-01, 8.129800694298381e-01, 9.223399131954029e-01, 9.854383454485866e-01, 9.956878500488212e-01, 9.999999999999999e-01}; 
    constant Real    iy_scaler = 9.990717523680953e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 6; 
    constant Real[6] data_x   = {-2.600000000000000e+01, -2.500000000000000e+01, -2.400000000000000e+01, -2.300000000000000e+01, -2.200000000000000e+01, -2.100000000000000e+01}; 
    constant Real[6] data_y   = {0.000000000000000e+00, 1.908396946564886e-02, 9.160305343511450e-02, 1.335877862595420e-01, 7.557251908396947e-01, 0.000000000000000e+00}; 
    constant Real[6] m_k      = {1.908396946564886e-02, 4.580152671755725e-02, 2.140048948378750e-02, 1.241228390059675e-01, 0.000000000000000e+00, -7.557251908396947e-01}; 
    constant Real[6] iy_start = {0.000000000000000e+00, 6.871825515386915e-03, 6.076876957133184e-02, 1.584941189663688e-01, 5.858978189423365e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.393486704511503e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.100000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.200000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-21</strong>.<br><br>
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
  <li>file creation date: 07-Jul-2022  </ul>
  </html>"));
  end SP_minus_21;
