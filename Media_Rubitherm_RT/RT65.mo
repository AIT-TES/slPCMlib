within slPCMlib.Media_Rubitherm_RT;
package RT65 "Rubitherm RT65; data taken from: data_sheet; last access: 01.12.2019."
   extends  slPCMlib.Interfaces.partialPCM;

  // ----------------------------------
  redeclare replaceable record propData "PCM record"

    constant String mediumName = "RT65";

    // --- parameters for phase transition functions ---
    constant Boolean modelForMelting =        true;
    constant Boolean modelForSolidification = false;
    constant Modelica.Units.SI.Temperature rangeTmelting[2] =  {3.281500000000000e+02, 3.401500000000000e+02}
             "temperature range melting {startT, endT}";
    constant Modelica.Units.SI.Temperature rangeTsolidification[2] = {3.281500000000000e+02, 3.401500000000000e+02}
             "temperature range solidification {startT, endT}";

    // --- parameters for heat capacity and enthalpy ---
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpS_linCoef = {2.000000000000000e+03, 0.0}
             "solid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificHeatCapacity[2] cpL_linCoef = {2.000000000000000e+03, 0.0}
             "liquid specific heat capacity, linear coefficients a + b*T";
    constant Modelica.Units.SI.SpecificEnthalpy   phTrEnth = 6.673358460067782e+04
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
    constant Integer pieces =  14;
    constant Integer[14] order =  {1, 6, 5, 5, 5, 5, 5, 5, 5, 5, 6, 5, 6, 1};
    constant Real[15] breaks = {-4.500000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 6.300000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 1.670000000000000e+02};
    constant Real[14,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.091188406467113e-15, -7.358135874889196e-17, -1.944315952336992e-18, 0.000000000000000e+00, 3.514347815699037e-03, -1.813316850604431e-03},  {1.701030965093440e-03, 6.691837974867320e-03, 7.943725397923830e-03, -1.122858855098237e-03, -9.628013680571272e-03, 4.594811359003349e-03, 0.000000000000000e+00},  {1.018053316121843e-02, 3.672714278151965e-03, -7.244819660765018e-03, 6.313200012650162e-03, 1.334604311444547e-02, -5.966525604229608e-03, 0.000000000000000e+00},  {2.030114530147141e-02, 3.167421943120657e-02, 3.210578302156224e-02, 3.211642813597061e-05, -1.648658490670257e-02, 7.114137400743913e-03, 0.000000000000000e+00},  {7.474081667641753e-02, 6.560648213564918e-02, 4.423996873193932e-03, 5.227150808764812e-03, 1.908410209701699e-02, -7.341227177209416e-03, 0.000000000000000e+00},  {1.617413214138330e-01, 1.297662008103542e-01, 6.119779010949628e-02, 8.151287424738628e-03, -1.762203378903008e-02, 5.842759697116832e-03, 0.000000000000000e+00},  {3.490773256665088e-01, 2.353413066330261e-01, 3.834704662069997e-02, -3.909250760213392e-03, 1.159176469655407e-02, -2.880926710167779e-03, 0.000000000000000e+00},  {6.275672661464078e-01, 3.322700728291633e-01, 6.736061541770644e-02, 1.364854092432511e-02, -2.812868854284825e-03, 2.073364117223625e-04, 0.000000000000000e+00},  {1.038240962875040e+00, 4.977221330790244e-01, 9.350238918219646e-02, 4.470429624409434e-03, -1.776186795673012e-03, 2.886090026599204e-04, 0.000000000000000e+00},  {1.632448336967657e+00, 6.924764981472531e-01, 9.914264730798589e-02, 2.517724683165884e-04, -3.331417823734106e-04, 1.023411154547020e-03, -1.218948631878079e-03},  {2.423790575631509e+00, 8.879879070201477e-01, 8.984899608599421e-02, -1.522565575326843e-02, -1.350031548780949e-02, 5.139018538174570e-03, 0.000000000000000e+00},  {3.378040526034747e+00, 9.937027626719662e-01, 1.456032128107770e-02, -1.783673232276069e-02, 1.219477720306336e-02, -4.404802065622483e-03, 6.552822083366037e-04},  {4.376912135010808e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_H;
  // ----------------------------------
  package data_C "spline interpolation data for cooling"
    extends Modelica.Icons.Package;
    constant Integer pieces =  14;
    constant Integer[14] order =  {1, 6, 5, 5, 5, 5, 6, 5, 5, 5, 5, 6, 5, 1};
    constant Real[15] breaks = {-4.500000000000000e+01, 5.500000000000000e+01, 5.600000000000000e+01, 5.700000000000000e+01, 5.800000000000000e+01, 5.900000000000000e+01, 6.000000000000000e+01, 6.100000000000000e+01, 6.200000000000000e+01, 6.300000000000000e+01, 6.400000000000000e+01, 6.500000000000000e+01, 6.600000000000000e+01, 6.700000000000000e+01, 1.670000000000000e+02};
    constant Real[14,7] coefs = { {0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00},  {0.000000000000000e+00, -1.177560469260618e-14, -4.298706592021073e-16, -3.197540459561284e-18, 0.000000000000000e+00, 4.153422858121765e-03, -1.935443172327124e-03},  {2.217979685782432e-03, 9.154455256635082e-03, 1.250258099631040e-02, 2.825365134675162e-03, -8.264533294298042e-03, 2.949654312274706e-03, 0.000000000000000e+00},  {2.138550209137974e-02, 2.432585103746601e-02, 8.880197572948109e-04, -7.362249197699484e-04, 6.483738267075483e-03, -1.783574267810382e-03, 0.000000000000000e+00},  {5.056331196563572e-02, 4.091029752199490e-02, 1.974603192233402e-02, 7.362985470428172e-03, -2.434133071976423e-03, 5.017917285712764e-05, 0.000000000000000e+00},  {1.161986729812735e-01, 9.300568135432920e-02, 2.773198163033134e-02, -1.871755088906242e-03, -2.183237207690785e-03, 1.439415945365542e-03, 0.000000000000000e+00},  {2.343207596147025e-01, 1.413185102443376e-01, 2.341145257112331e-02, 3.789455533986028e-03, 5.013842519136925e-03, 1.533114954065979e-03, -1.571300240390782e-03},  {4.078158351969617e-01, 2.178029253930753e-01, 5.662452022270101e-02, 7.749970343377875e-03, -1.089008631639491e-02, 3.393434200189130e-03, 0.000000000000000e+00},  {6.824965990399101e-01, 3.277087026039784e-01, 4.846825535635654e-02, -1.876032920310459e-03, 6.077084684550742e-03, -1.471397040116785e-03, 0.000000000000000e+00},  {1.061403211724369e+00, 4.359684680933794e-01, 6.458869430156176e-02, 7.718335416724655e-03, -1.279900516033185e-03, 2.689450970334558e-03, 0.000000000000000e+00},  {1.571088259990336e+00, 5.966285157342176e-01, 1.069588071588822e-01, 2.949324305593749e-02, 1.216735433563960e-02, -1.266428662143995e-02, 0.000000000000000e+00},  {2.303671893653573e+00, 8.843738434551566e-01, 1.417997961261329e-01, -4.848020581590357e-02, -5.115407877156014e-02, 5.247975665689938e-02, -1.388380909385446e-02},  {3.268807196210443e+00, 9.970124318948814e-01, 5.975136210239666e-03, -5.975136210239620e-03, 2.987568105119810e-03, -5.975136210239622e-04, 0.000000000000000e+00},  {4.268209682589420e+00, 1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00}};
  end data_C;
  // ----------------------------------
  redeclare function extends density_solid "Returns solid density"
  algorithm
    rho := 8.800000000000000e+02;
  end density_solid;
  // ----------------------------------
  redeclare function extends density_liquid "Returns liquid density"
  algorithm
    rho := 7.800000000000000e+02;
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
  This package contains solid and liquid properties for the PCM:  <strong>Rubitherm RT65</strong>.<br><br>
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
end RT65;
