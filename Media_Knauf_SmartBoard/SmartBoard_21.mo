within slPCMlib.Media_Knauf_SmartBoard;
package SmartBoard_21 "Knauf Gips KG Knauf SmartBoard 21; data taken from: DBU-Abschlussbericht-AZ-23836.pdf; last access: 13.07.2022."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "SmartBoard_21";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.Units.SI.Temperature rangeTmelting[2]=
    {2.731500000000000e+02+1.900000000000000e+01, 2.731500000000000e+02+2.800000000000000e+01}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2]=
    {2.731500000000000e+02+1.600000000000000e+01, 2.731500000000000e+02+2.700000000000000e+01}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {1.200000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {1.200000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 28696
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.Units.SI.Temperature       Tref = 273.15 + 17.0
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
    constant Integer  len_x =    19;
    constant Real[19] data_x =   {1.900000000000000e+01, 1.950000000000000e+01, 2.000000000000000e+01, 2.050000000000000e+01, 2.100000000000000e+01, 2.150000000000000e+01, 2.200000000000000e+01, 2.250000000000000e+01, 2.300000000000000e+01, 2.350000000000000e+01, 2.400000000000000e+01, 2.450000000000000e+01, 2.500000000000000e+01, 2.550000000000000e+01, 2.600000000000000e+01, 2.650000000000000e+01, 2.700000000000000e+01, 2.750000000000000e+01, 2.800000000000000e+01};
    constant Real[19] data_y =   {0.000000000000000e+00, 5.465061430767304e-04, 2.933057449351532e-03, 4.432183541059063e-03, 7.509162376516408e-03, 9.695162854904836e-03, 1.293629506992100e-02, 1.653483897247048e-02, 2.063972033583406e-02, 3.017736513435309e-02, 4.011253025533176e-02, 6.001591143955424e-02, 1.160904994390771e-01, 2.045504481203158e-01, 4.716995759049793e-01, 6.121300457233732e-01, 3.585771915493322e-01, 3.141950569042121e-02, 0.000000000000000e+00};
    constant Real[19] m_k =      {1.093012286153461e-03, 2.933057449351532e-03, 3.885677397982333e-03, 4.576104927164876e-03, 5.262979313845773e-03, 5.427132693404590e-03, 6.839676117565639e-03, 7.703425265913061e-03, 1.364252616188262e-02, 1.947280991949770e-02, 2.983854630520115e-02, 7.597796918374530e-02, 1.445345366807616e-01, 3.556090764659022e-01, 4.075795976030573e-01, 0.000000000000000e+00, -5.807105400329520e-01, -1.856872788559082e-01, -3.254084560952551e-02};
    constant Real[19] iy_start = {0.000000000000000e+00, 9.822343584666052e-05, 9.476728712406471e-04, 2.773319971253579e-03, 5.742266211009954e-03, 1.003691837094851e-02, 1.566141375002042e-02, 2.300605605636930e-02, 3.216954370974744e-02, 4.474354016547029e-02, 6.208790787987803e-02, 8.614192554754693e-02, 1.287104384586429e-01, 2.044202388390700e-01, 3.722824037342108e-01, 6.515353766995126e-01, 9.061319250990250e-01, 9.953389402996023e-01, 9.999999999999999e-01};
    constant Real    iy_scaler = 9.992997852714598e-01;
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer  len_x =    23;
    constant Real[23] data_x =   {1.600000000000000e+01, 1.650000000000000e+01, 1.700000000000000e+01, 1.750000000000000e+01, 1.800000000000000e+01, 1.850000000000000e+01, 1.900000000000000e+01, 1.950000000000000e+01, 2.000000000000000e+01, 2.050000000000000e+01, 2.100000000000000e+01, 2.150000000000000e+01, 2.200000000000000e+01, 2.250000000000000e+01, 2.300000000000000e+01, 2.350000000000000e+01, 2.400000000000000e+01, 2.450000000000000e+01, 2.500000000000000e+01, 2.550000000000000e+01, 2.600000000000000e+01, 2.650000000000000e+01, 2.700000000000000e+01};
    constant Real[23] data_y =   {0.000000000000000e+00, 6.127548215583842e-04, 3.373741688658988e-03, 3.663154801552382e-03, 3.672732438269915e-03, 4.179981374891020e-03, 5.190660776712669e-03, 7.291788462862735e-03, 9.262672841823860e-03, 1.094131505552181e-02, 1.388440817360499e-02, 1.726322054813351e-02, 2.149534094693003e-02, 3.078932399406117e-02, 3.572072831988318e-02, 7.050150349602093e-02, 2.226779427317542e-01, 3.918978134087402e-01, 4.644011220087688e-01, 4.562544370814194e-01, 2.195775730336422e-01, 7.347783995192542e-03, 0.000000000000000e+00};
    constant Real[23] m_k =      {1.225509643116768e-03, 3.373741688658988e-03, 1.728196857534788e-03, 1.789790140546889e-05, 5.460756017832844e-05, 1.517928338442753e-03, 3.111807087971714e-03, 4.072012065111191e-03, 3.649526592659075e-03, 4.621735331781125e-03, 6.321905492611701e-03, 7.610932773325040e-03, 1.352610344592766e-02, 9.978077451446250e-03, 2.785521353110468e-02, 1.869572144118710e-01, 3.213963099127193e-01, 2.417231792770146e-01, 0.000000000000000e+00, -4.888010956409594e-02, -4.489066530862269e-01, -4.398829867850505e-02, -2.943984784421927e-03};
    constant Real[23] iy_start = {0.000000000000000e+00, 1.084244528751315e-04, 1.139241223061662e-03, 2.933940678255164e-03, 4.766988476932510e-03, 6.699513212943943e-03, 9.008767351432545e-03, 1.210910608169725e-02, 1.625616295665892e-02, 2.128646862744147e-02, 2.745694322805039e-02, 3.521632165919660e-02, 4.478189840188898e-02, 5.792584010386174e-02, 7.417950098176207e-02, 9.741841528147888e-02, 1.679063394528415e-01, 3.231966461664485e-01, 5.422882483137992e-01, 7.734503939256399e-01, 9.507268837499399e-01, 9.990182291549967e-01, 1.000000000000000e+00};
    constant Real    iy_scaler = 9.999131430775384e-01;
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
  This package contains solid and liquid properties for the PCM:  <strong>Knauf Gips KG Knauf SmartBoard 21</strong>.<br><br>
  Information taken from: DBU-Abschlussbericht-AZ-23836.pdf - last access 13.07.2022.<br><br>
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
end SmartBoard_21;
