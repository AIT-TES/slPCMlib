within slPCMlib.Interfaces;
model phTransMod_MeltingCurve_Algebraic
  "Melting curve model, static (no hysteresis, algebraic equations) - obsolete -"

  extends Modelica.Icons.ObsoleteModel;
  extends basicPhTransModel;

equation
  (xi, dxi)  = PCM.phaseFrac_complMelting(indVar.T);

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
end phTransMod_MeltingCurve_Algebraic;
