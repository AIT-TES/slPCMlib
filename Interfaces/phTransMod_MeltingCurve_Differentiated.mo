within slPCMlib.Interfaces;
model phTransMod_MeltingCurve_Differentiated
  "Melting curve model, static (default model, no hysteresis, differentiated equations)"

  extends basicPhTransModel;


initial algorithm

  (xi,)  := PCM.phaseFrac_complMelting(indVar.T);

equation

  (,  dxi)  = PCM.phaseFrac_complMelting(indVar.T);

  when (indVar.T <= PCM.propData.rangeTsolidification[1]) then
    reinit(xi, 0.0);
  elsewhen (indVar.T >= PCM.propData.rangeTmelting[2]) then
    reinit(xi, 1.0);
  end when;

  der(xi) = dxi*indVar.der_T;

annotation (Icon(coordinateSystem(preserveAspectRatio=false),
    graphics={Text(lineColor={108,88,49},
    extent={{-90.0,-90.0},{90.0,90.0}},
    textString="→")}),
    Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
  Documentation(info="<html>
        <p>
        This specific phase transition model computes the 
        liquid mass phase fraction
        <var> xi </var> as a function of temperature <var> T </var>, 
        see the following example.
        </p>          
        </blockquote>          
        <p> 
        <img src=\"modelica://slPCMlib/Resources/Images/meltingCurveModel.png\">
        <br>
        </p>
        </html>",
  revisions="<html>
<ul>
<li>2022-06-01; initial version; by Tilman Barz </li>
</ul>
</html>"));
end phTransMod_MeltingCurve_Differentiated;
