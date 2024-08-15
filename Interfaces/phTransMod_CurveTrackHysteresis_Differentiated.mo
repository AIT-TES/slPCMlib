within slPCMlib.Interfaces;
model phTransMod_CurveTrackHysteresis_Differentiated
  "Curve track hysteresis model, static (differentiated equations)"

  extends basicPhTransModel;

protected
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
    heatingOn := true;
  else
    (xi,)  := PCM.phaseFrac_complSolidification(indVar.T);
    heatingOn := false;
  end if;

algorithm
   when  (indVar.T < uppLimPhTrange) then
      heatingOn :=false;
   end when;
   when (indVar.T > lowLimPhTrange) then
      heatingOn :=true;
   end when;

//equation

   if heatingOn then
     //(xi, dxi) :=PCM.phaseFrac_complMelting(indVar.T);
     (,   dxi)  :=PCM.phaseFrac_complMelting(indVar.T);
     //   print("pos", der(indVar.T), time);
   else
     //(xi, dxi) :=PCM.phaseFrac_complSolidification(indVar.T);
     (,  dxi) :=PCM.phaseFrac_complSolidification(indVar.T);
     //   print("  neg", der(indVar.T), time);
   end if;

   der(xi) :=dxi*indVar.der_T;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false),
              graphics={Text(lineColor={108,88,49},
              extent={{-90.0,-90.0},{90.0,90.0}},
        textString="⇆")}),
              Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
          <p>
          This specific phase transition model computes the 
          liquid mass phase fraction
          <var> xi </var> as a function of temperature <var> T </var>. <br>
          It is a static (so-called <strong>curve track</strong>) 
          hysteresis model, see e.g. 
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
          <img src=\"modelica://slPCMlib/Resources/Images/curveTrackHysteresisModel.png\">
          <br>
          </p>
          </html>",
    revisions="<html>
          <ul>
          <li>2022-06-01; initial version; by Tilman Barz </li>
          </ul>
          </html>"));
end phTransMod_CurveTrackHysteresis_Differentiated;
