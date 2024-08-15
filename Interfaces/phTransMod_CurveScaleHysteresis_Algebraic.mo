within slPCMlib.Interfaces;
model phTransMod_CurveScaleHysteresis_Algebraic
  "Curve scale hysteresis model, static (algebraic equations) - obsolete -"

  extends Modelica.Icons.ObsoleteModel;
  extends basicPhTransModel;

protected
  discrete Real T0(start=0.0);
  discrete Real xi0(start=1.0);
  discrete Real xi_at_T0(start=1.0);
  discrete Real scaler(start=1.0);
  discrete Boolean heatingOn(start=true);
  final constant Real eps = Modelica.Constants.small;
  Real xi_at_T(start=1.0), dxi_at_T(start=1.0);

initial algorithm
  assert((PCM.propData.modelForMelting == true)
     and (PCM.propData.modelForSolidification == true),
         "There should be a model for melting AND solidification 
           to simulate hysteresis!", AssertionLevel.error);
  assert(PCM.propData.rangeTmelting[1] >= PCM.propData.rangeTsolidification[1],
         "PCM.propData.rangeTmelting[1] < PCM.propData.rangeTsolidification[1]. 
          Phase transition function for complete melting should give always smaller values 
          compared with the function for solidification!",
         AssertionLevel.error);
   assert(PCM.propData.rangeTmelting[2] >= PCM.propData.rangeTsolidification[2],
          "PCM.propData.rangeTmelting[2] < PCM.propData.rangeTsolidification[2]. 
           Phase transition function for complete melting should give always smaller values 
           compared with the function for solidification!",
          AssertionLevel.error);

  if  (indVar.der_T >= 0) then
    (xi,)  := PCM.phaseFrac_complMelting(indVar.T);
    T0  :=pre(indVar.T);
    xi0 :=pre(xi);
    xi_at_T0  := xi0;
    scaler := 1.0;
    heatingOn := true;
  else
    (xi,)  := PCM.phaseFrac_complSolidification(indVar.T);
    T0  :=pre(indVar.T);
    xi0 :=pre(xi);
    xi_at_T0 := xi0;
    scaler := 1.0;
    heatingOn := false;
  end if;

algorithm
  when (indVar.der_T > 0) then
    T0   :=pre(indVar.T);
    xi0  :=pre(xi);
    (xi_at_T0,)  :=PCM.phaseFrac_complMelting(T0);
    if (xi_at_T0 >= xi0) then
      scaler :=1.0;
    else
      scaler :=(1.0 - xi0)/max((1.0 - xi_at_T0),eps);
    end if;
    heatingOn :=true;
  end when;
  when  (indVar.der_T < 0) then
    T0   :=pre(indVar.T);
    xi0  :=pre(xi);
    (xi_at_T0,)  :=PCM.phaseFrac_complSolidification(T0);
    if (xi_at_T0 <= xi0) then
      scaler :=1.0;
    else
      scaler :=xi0/max(xi_at_T0, eps);
    end if;
    heatingOn :=false;
  end when;

  if noEvent(heatingOn) then
    (xi_at_T, dxi_at_T) :=PCM.phaseFrac_complMelting(indVar.T);
    xi     :=1.0 - scaler*(1.0 - xi_at_T);
    dxi    :=scaler*dxi_at_T;
  else
    (xi_at_T, dxi_at_T) :=PCM.phaseFrac_complSolidification(indVar.T);
    xi     :=scaler*xi_at_T;
    dxi    :=scaler*dxi_at_T;
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
end phTransMod_CurveScaleHysteresis_Algebraic;
