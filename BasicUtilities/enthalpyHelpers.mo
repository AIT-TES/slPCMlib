within slPCMlib.BasicUtilities;
package enthalpyHelpers
  "Functions for evaluation of heat capacities and enthalpies"

  replaceable package PCM =
    slPCMlib.Media_generic.generic_7thOrderSmoothStep;
// ------------------------------------------------------------
  function spHeatCap_solid "Returns solid specific heat capacity"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Real  cp;
  protected
    Real  TT;
  algorithm
    TT := T-PCM.propData.Tref;
    cp := PCM.propData.cpS_linCoef[1]
        + PCM.propData.cpS_linCoef[2]*TT;
  end spHeatCap_solid;
// ------------------------------------------------------------
  function spHeatCap_liquid "Returns liquid specific heat capacity"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Real  cp;
  protected
    Real  TT;
  algorithm
    TT  := T-PCM.propData.Tref;
    cp  := PCM.propData.cpL_linCoef[1]
         + PCM.propData.cpL_linCoef[2]*TT;
  end spHeatCap_liquid;
// ------------------------------------------------------------
  function spHeatCap_baselineMelting
    "Returns baseline heat capacity for melting"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.SpecificHeatCapacity cp;
  protected
    Real  xi_M;
  algorithm
    (xi_M,) := PCM.phaseFrac_complMelting(T);
    cp  := xi_M*spHeatCap_liquid(T)
         + (1.0 - xi_M)*spHeatCap_solid(T);
  end spHeatCap_baselineMelting;
// ------------------------------------------------------------
  function enthalpy_solid "Returns solid enthalpy"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.SpecificEnthalpy h;
  protected
    Modelica.Units.SI.Temperature TT;
  algorithm
    TT := T-PCM.propData.Tref;
     h := PCM.propData.cpS_linCoef[1]*TT
        + PCM.propData.cpS_linCoef[2]*TT*TT/2
        + PCM.propData.href;
  end enthalpy_solid;
// ------------------------------------------------------------
  function enthalpy_liquid "Returns liquid enthalpy"
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.SpecificEnthalpy h;
  protected
    Modelica.Units.SI.Temperature TT;
  algorithm
    // do not use Tref, use Tmax !
    TT := T-PCM.propData.rangeTmelting[2];
    // do not use href !
     h := PCM.propData.cpL_linCoef[1]*TT
        + PCM.propData.cpL_linCoef[2]*TT*TT/2;
  end enthalpy_liquid;
// ------------------------------------------------------------

  function Tcp2pscWrapper
    "converts a function cp =f(T) to an instance of partialScalarFunction"
      extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
    input Real u;
    output Real y;
    input slPCMlib.BasicUtilities.enthalpyHelpers.cpTFunction cptFunction ;
    //  replaceable function cptFunction = slPCMlib.BasicUtilities.enthalpyHelpers.temperatureFunction  "cp(T) function";
  protected
    Real  T;
    Real cp;
  algorithm
    T := u;
    cp := cptFunction(T);
    y :=cp;
  end Tcp2pscWrapper;

  function temperatureFunction
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Real y;
  end temperatureFunction;
  
  partial function cpTFunction
    extends Modelica.Icons.Function;
    input Modelica.Units.SI.Temperature T;
    output Modelica.Units.SI.SpecificHeatCapacity cp;
  end cpTFunction;
annotation (Documentation(
 revisions="<html>
   <ul>
   <li>2022-06-01; initial version; by Tilman Barz </li>
   </ul>
   </html>"));
end enthalpyHelpers;
