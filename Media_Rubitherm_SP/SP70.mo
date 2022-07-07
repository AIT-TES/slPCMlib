within slPCMlib.Media_Rubitherm_SP;
package SP70 "Rubitherm SP70; data taken from: data sheet; last access: 02.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SP70";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.371500000000000e+02, 3.501500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.361500000000000e+02, 3.481500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 1.500000000000000e+05
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
    constant Integer  len_x    = 14; 
    constant Real[14] data_x   = {6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 6.800000000000000e+01, 6.900000000000000e+01, 7.000000000000000e+01, 7.100000000000000e+01, 7.200000000000000e+01, 7.300000000000000e+01, 7.400000000000000e+01, 7.500000000000000e+01, 7.600000000000000e+01, 7.700000000000000e+01}; 
    constant Real[14] data_y   = {0.000000000000000e+00, 9.042173908512759e-03, 1.816610872601731e-02, 2.737254374887167e-02, 3.666222495828009e-02, 7.316242680627691e-02, 1.371191915526174e-01, 2.649494818480094e-01, 1.859656359632346e-01, 7.009890877395542e-02, 4.360623357234820e-02, 7.112705044534286e-02, 6.272801969655643e-02, 0.000000000000000e+00}; 
    constant Real[14] m_k      = {9.042173908512759e-03, 9.083054363008655e-03, 9.165184920179455e-03, 9.248058116131389e-03, 2.289494152870262e-02, 5.022848329716864e-02, 9.589352752086623e-02, 0.000000000000000e+00, -9.742528653702698e-02, -7.117970119544320e-02, 0.000000000000000e+00, 0.000000000000000e+00, -2.519709224635927e-02, -6.272801969655643e-02}; 
    constant Real[14] iy_start = {0.000000000000000e+00, 4.490821323396127e-03, 1.800727853167224e-02, 4.063436944482385e-02, 7.133092205370578e-02, 1.236525251565759e-01, 2.243854455507878e-01, 4.321681922148459e-01, 6.643558455955614e-01, 7.894527998004640e-01, 8.400709894348281e-01, 8.970965700303570e-01, 9.657134754065492e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.940547084268252e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 13; 
    constant Real[13] data_x   = {6.300000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 6.800000000000000e+01, 6.900000000000000e+01, 7.000000000000000e+01, 7.100000000000000e+01, 7.200000000000000e+01, 7.300000000000000e+01, 7.400000000000000e+01, 7.500000000000000e+01}; 
    constant Real[13] data_y   = {0.000000000000000e+00, 2.639392868018648e-02, 2.534896947342418e-02, 1.695119888758530e-01, 2.683765725793594e-01, 3.237361261056857e-01, 4.697983696837036e-02, 5.528970248567900e-03, 5.310073378948366e-03, 3.149377151206269e-02, 6.983779799360550e-02, 2.748196418418981e-02, 0.000000000000000e+00}; 
    constant Real[13] m_k      = {2.639392868018648e-02, 0.000000000000000e+00, 0.000000000000000e+00, 1.215138015529676e-01, 7.711206861491635e-02, 0.000000000000000e+00, -1.232998994432190e-01, -6.566906088586022e-04, 0.000000000000000e+00, 3.226386230732857e-02, 0.000000000000000e+00, -3.491889899680275e-02, -2.748196418418981e-02}; 
    constant Real[13] iy_start = {0.000000000000000e+00, 1.532764252788027e-02, 4.108345681375403e-02, 1.279975712325542e-01, 3.496468668960558e-01, 6.507772496033328e-01, 8.455358234224549e-01, 8.614982936002940e-01, 8.668391127466345e-01, 8.824821481250154e-01, 9.355981155558889e-01, 9.869374090382248e-01, 9.999999999999997e-01}; 
    constant Real    iy_scaler = 9.955304091940006e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 1.500000000000000e+03;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 1.300000000000000e+03;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm SP70</strong>.<br><br>
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
  end SP70;
