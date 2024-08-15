within slPCMlib.Interfaces;
model phTransMod_CurveScaleHysteresis_Differentiated
  "Curve scale hysteresis model, static (differentiated equations)"

  extends basicPhTransModel;

protected
  Real scalerM, scalerS;
  final constant Real eps = Modelica.Constants.small;
  Real xiH_at_T, dxiH_at_T;
  Real xiC_at_T, dxiC_at_T;
  discrete Boolean heatingOn(start=true);

  constant Modelica.Units.SI.Temperature lowLimPhTrange=min(PCM.propData.rangeTmelting[
      1], PCM.propData.rangeTsolidification[1]);
  constant Modelica.Units.SI.Temperature uppLimPhTrange=max(PCM.propData.rangeTmelting[
      2], PCM.propData.rangeTsolidification[2]);

initial algorithm
  assert((PCM.propData.modelForMelting == true)
     and (PCM.propData.modelForSolidification == true),
         "There should be a model for melting AND solidification 
           to simulate hysteresis!", AssertionLevel.error);

// --- the model is robust against following inconsistent parameters
//  assert(PCM.propData.rangeTmelting[1] >= PCM.propData.rangeTsolidification[1],
//         "PCM.propData.rangeTmelting[1] < PCM.propData.rangeTsolidification[1].
//          Phase transition function for complete melting should give always smaller values
//          compared with the function for solidification!",
//         AssertionLevel.error);
//   assert(PCM.propData.rangeTmelting[2] >= PCM.propData.rangeTsolidification[2],
//          "PCM.propData.rangeTmelting[2] < PCM.propData.rangeTsolidification[2].
//           Phase transition function for complete melting should give always smaller values
//           compared with the function for solidification!",
//          AssertionLevel.error);

  if  (indVar.der_T >= 0) then
    (xi,)  := PCM.phaseFrac_complMelting(indVar.T);
    heatingOn := true;
  else
    (xi,)  := PCM.phaseFrac_complSolidification(indVar.T);
    heatingOn := false;
  end if;

algorithm
  when (indVar.T >= lowLimPhTrange)
    and (indVar.T <= uppLimPhTrange)
    and (indVar.der_T > 0) then
      heatingOn :=true;
  end when;
  when (indVar.T >= lowLimPhTrange)
    and (indVar.T <= uppLimPhTrange)
    and (indVar.der_T < 0) then
      heatingOn :=false;
  end when;

equation
  (xiH_at_T, dxiH_at_T) = PCM.phaseFrac_complMelting(indVar.T);
  (xiC_at_T, dxiC_at_T) = PCM.phaseFrac_complSolidification(indVar.T);

  when (indVar.T <= lowLimPhTrange) then
    reinit(xi, 0.0);
  elsewhen (indVar.T >= uppLimPhTrange) then
    reinit(xi, 1.0);
  end when;

  // alternative hysteresis
  // scalerM = ( xiC_at_T - Xi)
  //         / ( xiC_at_T - xiH_at_T+eps);
  // scalerS = ( Xi - xiH_at_T)
  //         / ( xiC_at_T - xiH_at_T + eps);

  if noEvent(xiH_at_T >= xi) then
    scalerM = 1.0;
  elseif noEvent(xi > 1.0) then
    scalerM = 0.0;
  else
    scalerM  = (1.0 - xi)/max((1.0 - xiH_at_T), eps);
  end if;
  if noEvent(xiC_at_T <= xi) then
    scalerS = 1.;
  elseif noEvent(xi < 0.0) then
    scalerS = 0.0;
  else
    scalerS  = xi/max(xiC_at_T, eps);
  end if;

  if noEvent(heatingOn == true) then
    der(xi) = scalerM*dxiH_at_T*indVar.der_T;
    dxi     = scalerM*dxiH_at_T;
  else
    der(xi) = scalerS*dxiC_at_T*indVar.der_T;
    dxi     = scalerS*dxiC_at_T;
  end if;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false),
              graphics={Text(lineColor={108,88,49},
              extent={{-90.0,-90.0},{90.0,90.0}},
              textString="↔")}),
              Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
          <p>
          This specific phase transition model computes the 
          liquid mass phase fraction
          <var> xi </var> as a function of temperature <var> T </var>. <br>
          It is a static (so-called <strong>curve scale</strong>) 
          hysteresis model, e.g. described in  
          </p>
          <blockquote>          
          <p>
           Barz, T., Emhofer, J., Marx, K., Zsembinszki, G., & Cabeza, L. F. 
           (2019). Phenomenological modelling of phase transitions with 
           hysteresis in solid/liquid PCM. Journal of Building Performance 
           Simulation, 12(6), 770-788. 
          <a href>doi.org/10.1080/19401493.2019.1657953</a>. 
          </p>          
          </blockquote>          
          <p> 
          <img src=\"modelica://slPCMlib/Resources/Images/curveScaleHysteresisModel.png\">
          <br>
          </p>
          </html>",
    revisions="<html>
          <ul>
          <li>2022-06-01; initial version; by Tilman Barz </li>
          </ul>
          </html>"));
end phTransMod_CurveScaleHysteresis_Differentiated;
