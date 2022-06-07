within slPCMlib.BasicUtilities;
function smoothStep7thOrder "7th-order smoothstep"
  extends Modelica.Icons.Function;

  input  Real x;
  input  Real startStepX, endStepX;
  output Real ff, dff;
protected
    Real xx, dx, tmp2;
algorithm
    // Scale, and clamp x to 0..1 range, helpers
    xx   :=(x - startStepX)/(endStepX - startStepX);
    dx   :=1.0/(endStepX - startStepX);
    tmp2 :=endStepX - startStepX;
    if     xx < 0.0 then
        // iff :=0.0;
          ff :=0.0;
         dff :=0.0;
    elseif xx > 1.0 then
//          iff := (x - endStepX)
//                 - 20.0*tmp2 /8.0
//                 + 70.0*tmp2 /7.0
//                 - 84.0*tmp2 /6.0
//                 + 35.0*tmp2 /5.0;
          ff :=1.0;
         dff :=0.0;
    else
//            iff := - 20.0*tmp2*xx^(8.0) /8.0
//                 + 70.0*tmp2*xx^(7.0) /7.0
//                 - 84.0*tmp2*xx^(6.0) /6.0
//                 + 35.0*tmp2*xx^(5.0) /5.0;
          ff := -  20.0*xx^7.0 +  70.0*xx^6.0 -  84.0.*xx^5.0 +  35.0*xx^4.0;
         dff :=(- 140.0*xx^6.0 + 420.0*xx^5.0 - 420.0.*xx^4.0 + 140.0*xx^3.0)*dx;
         end if;

  annotation (Icon(graphics,
      coordinateSystem(preserveAspectRatio=false)), Diagram(graphics,
      coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
    <p>
    Function for the evaluation of 7th-order smoothstep described in  
    <a href>https://en.wikipedia.org/wiki/Smoothstep</a>. 
    </p>
    <p>
    The function is 3 times continuously differentiable (C^3 smooth). <br>
    (Note that the effective (or apparent) heat capacity model uses the
    first derivative of the smoothstep function.) 
    </p>
    </html>",
  revisions="<html>
      <ul>
      <li>2022-06-01; initial version; by Tilman Barz </li>
      </ul>
      </html>"));
end smoothStep7thOrder;
