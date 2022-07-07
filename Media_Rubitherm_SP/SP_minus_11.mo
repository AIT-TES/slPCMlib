within slPCMlib.Media_Rubitherm_SP;
package SP_minus_11 "Rubitherm SP-11; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP_minus_11";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.591500000000000e+02, 2.671500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.571500000000000e+02, 2.641500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.198513964714595e+05
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
    constant Integer  len_x    = 9; 
    constant Real[9] data_x   = {-1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.100000000000000e+01, -1.000000000000000e+01, -9.000000000000000e+00, -8.000000000000000e+00, -7.000000000000000e+00, -6.000000000000000e+00}; 
    constant Real[9] data_y   = {0.000000000000000e+00, 4.517703596818905e-03, 5.866891782127964e-02, 8.304010865341245e-01, 6.017236434223561e-02, 2.321521284723033e-02, 1.386516730185005e-02, 9.159547556460975e-03, 0.000000000000000e+00}; 
    constant Real[9] m_k      = {2.062948832758760e-03, 5.266973729725725e-03, 1.623682389064479e-01, 0.000000000000000e+00, -1.106894553835870e-01, -6.350109172907781e-03, -7.027832645384679e-03, -6.932583650925025e-03, -9.159547556460975e-03}; 
    constant Real[9] iy_start = {0.000000000000000e+00, 1.989988670273312e-03, 2.047424066032776e-02, 4.781119429584329e-01, 9.321981246340578e-01, 9.651661357624135e-01, 9.837454272903196e-01, 9.952390983297918e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.990656657644471e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 8; 
    constant Real[8] data_x   = {-1.600000000000000e+01, -1.500000000000000e+01, -1.400000000000000e+01, -1.300000000000000e+01, -1.200000000000000e+01, -1.100000000000000e+01, -1.000000000000000e+01, -9.000000000000000e+00}; 
    constant Real[8] data_y   = {0.000000000000000e+00, 1.060068461751247e-02, 2.644552428660119e-02, 3.690603857286015e-02, 1.427172701102885e-01, 7.780018967760808e-01, 5.328585636652013e-03, 0.000000000000000e+00}; 
    constant Real[8] m_k      = {1.060068461751247e-02, 1.322276214330059e-02, 6.924760199103562e-03, 2.613165003068152e-02, 3.163562664811295e-01, 0.000000000000000e+00, -1.598425734709627e-02, -2.189544381966993e-04}; 
    constant Real[8] iy_start = {0.000000000000000e+00, 5.077258006688492e-03, 2.410803710736032e-02, 5.415615190475861e-02, 1.197233037990121e-01, 6.060074578635025e-01, 9.986506990032293e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.990991756282421e-01; 
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP-11</strong>.<br><br>
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
  end SP_minus_11;
