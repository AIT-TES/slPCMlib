within slPCMlib.Media_Rubitherm_RT;
package RT44HC "Rubitherm RT44HC; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT44HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.121500000000000e+02, 3.181500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.111500000000000e+02, 3.181500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 2.197454929466062e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature            Tref = 3.121500000000000e+02
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
    constant Integer pieces =  8;
    constant Integer[8] order =  {1, 5, 6, 5, 6, 5, 6, 1};
    constant Real[9] breaks = {-6.100000000000000e+01, 3.900000000000000e+01, 4.000000000000000e+01, 4.100000000000000e+01, 4.200000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 1.450000000000000e+02};
    constant Real[8,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -6.532119647411234e-15, -5.945168257897701e-18, -2.644206516638919e-18, -2.891205793294678e-19, 8.195023607598293e-04, 0.000000000000000e+00},  {8.195023607532883e-04, 4.097511803785665e-03, 8.195023607598276e-03, 8.195023607598216e-03, 4.097511803799146e-03, 1.287625446787967e-02, -6.244511842118311e-03},  {3.203631580929595e-02, 8.837687834365940e-02, 9.246003230020979e-02, 2.845737865922522e-02, -2.518889348857720e-02, 8.870949209189554e-03, 0.000000000000000e+00},  {2.250126608330027e-01, 3.022682510133954e-01, 1.154082994383178e-01, 1.641129679681200e-02, 1.916585255737057e-02, 9.586772563539118e-03, -9.127258145353468e-03},  {6.787258750570843e-01, 6.521524644555254e-01, 2.385961586280664e-01, 6.397269754616117e-03, -6.980915680523586e-02, 2.133783886480758e-02, 0.000000000000000e+00},  {1.527400449954864e+00, 9.759891580786230e-01, 5.231141570857534e-02, -5.946096881825156e-02, 3.688003751880201e-02, -1.166573936956614e-02, 1.429910621935245e-03},  {2.522884263694982e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  9;
    constant Integer[9] order =  {1, 5, 6, 5, 6, 5, 6, 5, 1};
    constant Real[10] breaks = {-6.200000000000000e+01, 3.800000000000000e+01, 3.900000000000000e+01, 4.000000000000000e+01, 4.100000000000000e+01, 4.200000000000000e+01, 4.300000000000000e+01, 4.400000000000000e+01, 4.500000000000000e+01, 1.450000000000000e+02};
    constant Real[9,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -6.177192886879588e-14, -5.714274072162773e-17, 2.020406255786045e-17, -2.891205793294678e-19, 1.388232103908729e-03, 0.000000000000000e+00},  {1.388232103846919e-03, 6.941160519452244e-03, 1.388232103908729e-02, 1.388232103908733e-02, 6.941160519543642e-03, 2.681633435048170e-02, -1.454102731687836e-02},  {5.531050225462077e-02, 1.509529156441285e-01, 1.472241810252527e-01, 1.898976028451173e-02, -7.709257748122324e-02, 2.658799103400826e-02, 0.000000000000000e+00},  {3.219727727612987e-01, 3.269402037933992e-01, 7.517907331531140e-03, -2.350063930029859e-02, 5.584737768881807e-02, -4.604333401557558e-02, 1.350191596241308e-02},  {6.562362042215861e-01, 3.456584370074427e-01, 1.419565484398431e-02, 9.493850547479633e-03, 2.815944704713646e-02, -2.634252200491522e-03, 0.000000000000000e+00},  {1.051109341467138e+00, 5.019978255239106e-01, 1.852913667643267e-01, 9.578911673111025e-02, 1.498818604467885e-02, -9.805448850602133e-02, 3.550743074242486e-02},  {1.786628778767567e+00, 9.426727953491633e-01, 1.146544093018904e-01, -1.146544093018904e-01, 5.732720465094519e-02, -1.146544093018904e-02, 0.000000000000000e+00},  {2.775163337837486e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.000000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.000000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT44HC</strong>.<br><br>
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
  <li>file creation date: 19-Jul-2022  </ul>
  </html>"));
end RT44HC;
