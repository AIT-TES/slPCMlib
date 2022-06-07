within slPCMlib.Media;
package RT11HC "Rubitherm RT11HC, data taken from data_sheet"
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT11HC";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = true;
    constant Modelica.SIunits.Temp_K rangeTmelting[2] =  {2.781500000000000e+02, 2.871500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.SIunits.Temp_K rangeTsolidification[2] = {2.781500000000000e+02, 2.861500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.SIunits.SpecificEnthalpy   phTrEnth = 1.833430289519252e+05
             "scalar phase transition enthalpy";

    // --- reference values ---
    constant Modelica.SIunits.Temp_K            Tref = 273.15+25
             "reference temperature";
    constant Modelica.SIunits.SpecificEnthalpy  href = 0.0
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
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
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
    (Xi, dXi) := BasicUtilities.splineEval(T-273.15,
                     pieces, order, breaks, coefs[:,:]);
  end phaseFrac_complSolidification;

  // ----------------------------------
  package data_H "spline interpolation data for heating"
    extends Modelica.Icons.Package;
    constant Integer pieces =  11;
    constant Integer[11] order =  {1, 6, 5, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[12] breaks = {-9.500000000000000e+01, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.400000000000000e+01, 1.140000000000000e+02};
    constant Real[11,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 7.721545428012638e-16, -1.203294092601427e-18, -4.224831238404343e-18, 5.782411586589357e-19, 2.173162078411849e-03, -9.039144247716040e-04},  {1.269247653641012e-03, 5.442323843430779e-03, 8.172904412544419e-03, 3.653332288686406e-03, -2.692905979514815e-03, 7.892651527386838e-04, 0.000000000000000e+00},  {1.663416737152648e-02, 2.592283138021264e-02, 1.086811692890158e-02, 7.743598980139888e-04, 1.253419784178604e-03, -1.686608085814124e-03, 6.474261978685898e-04},  {5.441371347488776e-02, 5.044734082691206e-02, 1.355702743790278e-02, 1.870482133958963e-03, 2.531772323136832e-03, -7.227585966942320e-04, 0.000000000000000e+00},  {1.220975776001042e-01, 8.968613841367278e-02, 2.713152181165834e-02, 4.769985459563971e-03, -1.082020660334328e-03, 1.661326273693404e-03, 0.000000000000000e+00},  {2.442645288983583e-01, 1.622376871428303e-01, 5.156261696527830e-02, 1.705516555516070e-02, 7.224610708132693e-03, 1.427538157875385e-03, 0.000000000000000e+00},  {4.837721474276357e-01, 3.525645513607957e-01, 1.603511594583104e-01, 6.022898996644532e-02, 1.436230149750962e-02, -1.966730652616339e-02, 0.000000000000000e+00},  {1.051611843184533e+00, 8.130665135359836e-01, 2.305388730810701e-01, -7.899486930515008e-02, -8.397423113330731e-02, 8.619912269985126e-02, -2.282284395784063e-02},  {1.995624408105140e+00, 9.953212770016573e-01, 9.357445996679252e-03, -9.357445996679250e-03, 4.678722998339625e-03, -9.357445996679251e-04, 0.000000000000000e+00},  {2.994688663505469e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  10;
    constant Integer[10] order =  {1, 6, 5, 5, 5, 6, 5, 5, 6, 1};
    constant Real[11] breaks = {-9.500000000000000e+01, 5.000000000000000e+00, 6.000000000000000e+00, 7.000000000000000e+00, 8.000000000000000e+00, 9.000000000000000e+00, 1.000000000000000e+01, 1.100000000000000e+01, 1.200000000000000e+01, 1.300000000000000e+01, 1.130000000000000e+02};
    constant Real[10,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, 1.831721285800849e-15, -9.818292607668563e-18, -1.135388605705447e-18, 5.782411586589357e-19, 1.907950460385240e-03, -7.021216944860938e-04},  {1.205828765900967e-03, 5.327022135064702e-03, 8.547679186560695e-03, 5.037070714130522e-03, -9.920731153652089e-04, -7.319402757694609e-05, 0.000000000000000e+00},  {1.905233365871473e-02, 3.319933005117764e-02, 1.697451236099184e-02, 3.368379769002258e-04, -1.358043253249939e-03, 6.960819111115242e-04, 0.000000000000000e+00},  {6.890105270564602e-02, 6.620710524633891e-02, 1.679758588330856e-02, 1.865484075015710e-03, 2.122366302307681e-03, -7.132135904479766e-04, 0.000000000000000e+00},  {1.551803806221689e-01, 1.103221264950529e-01, 2.799610001772170e-02, 3.222813379766656e-03, -1.443701649932202e-03, 1.001030936472486e-02, -3.532423059828399e-03},  {3.017556051696744e-01, 1.990649685345994e-01, 7.611907800725175e-02, 2.690263923071849e-02, -4.378500723733874e-03, 3.579223016097065e-03, 0.000000000000000e+00},  {6.030430132346070e-01, 4.323931544272257e-01, 1.663482215179724e-01, 4.518086649675365e-02, 1.351761435675145e-02, -1.817970357970465e-02, 0.000000000000000e+00},  {1.242303166453606e+00, 8.638041364815686e-01, 2.011994713516975e-01, -8.254571187328702e-02, -7.738090354177178e-02, 8.666843639540352e-02, -2.373075189568305e-02},  {2.210317843371533e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
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
  This package contains solid and liquid properties for the PCM:  <strong>RT11HC</strong>.
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
end RT11HC;
