within slPCMlib.Media_Knauf_SmartBoard;
package SmartBoard_26 "Knauf Gips KG Knauf SmartBoard 26; data taken from: DBU-Abschlussbericht-AZ-23836.pdf & Buildings Library, https://simulationresearch.lbl.gov; last access: 13.07.2022."
   extends  slPCMlib.Interfaces.partialPCM;
      
  // ----------------------------------
  redeclare replaceable record propData "PCM record"
      
    constant String mediumName = "SmartBoard_26";
      
    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting        = true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {2.951500000000000e+02, 3.031500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {2.931500000000000e+02, 3.026500000000000e+02}
             "temperature range solidification {startT, endT}";
      
    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.200000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.200000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.893480000000000e+04
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
    constant Integer  len_x    = 17; 
    constant Real[17] data_x   = {2.200000000000000e+01, 2.250000000000000e+01, 2.300000000000000e+01, 2.350000000000000e+01, 2.400000000000000e+01, 2.450000000000000e+01, 2.500000000000000e+01, 2.550000000000000e+01, 2.600000000000000e+01, 2.650000000000000e+01, 2.700000000000000e+01, 2.750000000000000e+01, 2.800000000000000e+01, 2.850000000000000e+01, 2.900000000000000e+01, 2.950000000000000e+01, 3.000000000000000e+01}; 
    constant Real[17] data_y   = {0.000000000000000e+00, 3.069415253380637e-03, 7.527587486567503e-03, 1.167059293333463e-02, 1.800601429200370e-02, 2.611582877139498e-02, 3.654092290476624e-02, 5.269330000320109e-02, 7.654621830223128e-02, 1.129953517804195e-01, 1.819727394381412e-01, 2.717291260748688e-01, 4.291832362123292e-01, 4.946238901337682e-01, 2.620019493424283e-01, 1.532382707148857e-02, 0.000000000000000e+00}; 
    constant Real[17] m_k      = {6.138830506761275e-03, 7.527587486567503e-03, 8.601177679953994e-03, 1.047842680543620e-02, 1.444523583806035e-02, 1.853490861276254e-02, 2.657747123180611e-02, 4.000529539746504e-02, 6.030205177721838e-02, 1.054265211359099e-01, 1.587337742944493e-01, 2.472104967741880e-01, 2.228947640588994e-01, 0.000000000000000e+00, -4.793000630622797e-01, -9.132031133225117e-02, -1.068218509428702e-02}; 
    constant Real[17] iy_start = {0.000000000000000e+00, 7.381626961359974e-04, 3.364126681652719e-03, 8.122894777505309e-03, 1.545683463905845e-02, 2.639825959960090e-02, 4.188946542626510e-02, 6.391055779275528e-02, 9.578641777082195e-02, 1.422154466919247e-01, 2.148214578297902e-01, 3.263645702306575e-01, 5.020376760662065e-01, 7.375505659456925e-01, 9.366226813960126e-01, 9.978497577216463e-01, 9.999999999999999e-01}; 
    constant Real    iy_scaler = 9.996496849384645e-01; 
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;    
    constant Integer  len_x    = 20; 
    constant Real[20] data_x   = {2.000000000000000e+01, 2.050000000000000e+01, 2.100000000000000e+01, 2.150000000000000e+01, 2.200000000000000e+01, 2.250000000000000e+01, 2.300000000000000e+01, 2.350000000000000e+01, 2.400000000000000e+01, 2.450000000000000e+01, 2.500000000000000e+01, 2.550000000000000e+01, 2.600000000000000e+01, 2.650000000000000e+01, 2.700000000000000e+01, 2.750000000000000e+01, 2.800000000000000e+01, 2.850000000000000e+01, 2.900000000000000e+01, 2.950000000000000e+01}; 
    constant Real[20] data_y   = {0.000000000000000e+00, 1.317353145197244e-03, 3.883961516791748e-03, 6.504699670153018e-03, 9.127299408814713e-03, 1.233704347511479e-02, 1.695757427126962e-02, 2.533944757282883e-02, 3.376979291289783e-02, 4.649727061622753e-02, 6.172976110134171e-02, 8.436078983202527e-02, 1.247870665861878e-01, 1.848297421683142e-01, 3.674174957226387e-01, 4.986079712485443e-01, 3.399799150677990e-01, 1.334823284274044e-01, 4.907048725860752e-02, 0.000000000000000e+00}; 
    constant Real[20] m_k      = {2.634706290394488e-03, 3.883961516791748e-03, 5.187346524955774e-03, 5.243337892022965e-03, 5.832343804961775e-03, 7.830274862454904e-03, 1.300240409771404e-02, 1.681221864162822e-02, 2.115782304339870e-02, 2.795996818844387e-02, 3.786351921579774e-02, 6.305730548484610e-02, 1.004689523362889e-01, 2.426304291364508e-01, 3.137782290802301e-01, 0.000000000000000e+00, -3.651256428211399e-01, -2.909094278091915e-01, -1.334823284274044e-01, -9.814097451721504e-02}; 
    constant Real[20] iy_start = {0.000000000000000e+00, 3.026766680971235e-04, 1.573184057440995e-03, 4.163744003979547e-03, 8.051310896358880e-03, 1.336461780112385e-02, 2.056540154771904e-02, 3.103829809574134e-02, 4.569430460391651e-02, 6.557761425260140e-02, 9.237179400002560e-02, 1.282941424543548e-01, 1.796937840437301e-01, 2.539803247008229e-01, 3.902737412196005e-01, 6.128498573256587e-01, 8.296484459348422e-01, 9.462230879766711e-01, 9.884928156488472e-01, 1.000000000000000e+00}; 
    constant Real    iy_scaler = 9.979049052861844e-01; 
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm 
    rho := 7.670000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm 
    rho := 7.670000000000000e+02;
  end density_liquid;
  // ----------------------------------
  redeclare function extends conductivity_solid "Returns solid thermal conductivity"
  algorithm 
    lambda := 1.800000000000000e-01;
  end conductivity_solid;
  // ----------------------------------
  redeclare function extends conductivity_liquid "Returns liquid thermal conductivity"
  algorithm 
    lambda := 1.800000000000000e-01;
  end conductivity_liquid;
  // ----------------------------------
      
annotation(Documentation(
  info="<html>
  <p>
  This package contains solid and liquid properties for the PCM:  <strong>Knauf Gips KG Knauf SmartBoard 26</strong>.<br><br>
  Information taken from: DBU-Abschlussbericht-AZ-23836.pdf & Buildings Library, https://simulationresearch.lbl.gov - last access 13.07.2022.<br><br>
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
  <li>file creation date: 13-Jul-2022  </ul>
  </html>"));
  end SmartBoard_26;
