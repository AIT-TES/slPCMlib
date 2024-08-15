within slPCMlib.Interfaces;
model phTransMod_CurveSwitchHysteresis_Algebraic
  "Curve switch hysteresis model, static (algebraic equations) - obsolete -"

  extends Modelica.Icons.ObsoleteModel;
  extends basicPhTransModel;

protected
  discrete Integer modelInd(start=1);
  discrete Real xi0(start=0.5);
  Real xiC_at_T, xiH_at_T;

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

  if  (der(indVar.T) >= 0) then
    modelInd := 1;
    (xi,)  := PCM.phaseFrac_complMelting(indVar.T);
    xi0 := pre(xi);
  else
    modelInd := -1;
    (xi,)  := PCM.phaseFrac_complSolidification(indVar.T);
    xi0 := pre(xi);
  end if;

algorithm
  when        (indVar.T < PCM.propData.rangeTsolidification[1]) then
    modelInd := 1; // activate melting curve
    // print("  ->Low", modelInd, time);
  end when;
  when    (indVar.T > PCM.propData.rangeTmelting[2]) then
    modelInd := -1; // activate solidification curve
    // print("  ->Up", modelInd, time);
  end when;

  when (der(indVar.T) > 0) then
    // if during phase transition
    if ( (PCM.propData.rangeTsolidification[1] < indVar.T) and (indVar.T < PCM.propData.rangeTmelting[2])) then
      modelInd  := 0; // switch: transition from melting to solidification curve
      (xi0,)  := PCM.phaseFrac_complSolidification(indVar.T); // unclocked discrete time variables
      //print("  ->pos sign . change in the phTrRg", modelInd, time);
    end if;
  end when;
  when (der(indVar.T) <= 0) then
    // if during phase transition
    if ( (PCM.propData.rangeTsolidification[1] < indVar.T) and (indVar.T < PCM.propData.rangeTmelting[2])) then
      modelInd  := 0; // switch: transition from cooling to heating curve
      (xi0,)  := PCM.phaseFrac_complMelting(indVar.T); // unclocked discrete time variables
      //print("  ->neg sign . change in the phTrRg", modelInd, time);
    end if;
  end when;

  // --- finish the curve switch,
  //     i.e. apporach the melting/solidification curve during phase switch ---
  (xiC_at_T,)   := PCM.phaseFrac_complSolidification(indVar.T);
  (xiH_at_T,)   := PCM.phaseFrac_complMelting(indVar.T);
  when (xi0 < xiH_at_T) then
    // if during switch transition
    if ( modelInd==0) then
      modelInd := 1;  // activate melting curve
      //print("  ->finish curve switch - now heating curve", modelInd, time);
    end if;
  end when;
  when (xi0 > xiC_at_T) then
      // if during switch transition
    if ( modelInd==0) then
      modelInd := -1;  // activate cooling curve
      //print("  ->finish curve switch - now cooling curve", modelInd, time);
    end if;
  end when;

  if noEvent(modelInd==1) then
    (xi, dxi) := PCM.phaseFrac_complMelting(indVar.T);
    //   print("pos", der(indVar.T), time);
  elseif noEvent(modelInd==-1) then
    (xi, dxi) := PCM.phaseFrac_complSolidification(indVar.T);
    //    print("  neg", der(indVar.T), time);
  elseif noEvent(modelInd==0) then
    xi  := xi0;
    dxi := 0.0;
  end if;

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
          It is a static (so-called <strong>curve switch</strong>) 
          hysteresis model, see e.g. 
          </p>
          <blockquote>          
          <p>
           Barz, T., Emhofer, J., Marx, K., Zsembinszki, G., & Cabeza, L. F. 
           (2019). Phenomenological modelling of phase transitions with 
           hysteresis in solid/liquid PCM. Journal of Building Performance 
           Simulation, 12(6), 770-788. 
          <a href>doi.org/10.1080/19401493.2019.1657953 </a>. 
          </p>          
          </blockquote>          
          <p> 
          <img src=\"modelica://slPCMlib/Resources/Images/curveSwitchHysteresisModel.png\">
          <br>
          </p>
          </html>",
    revisions="<html>
          <ul>
          <li>2022-06-01; initial version; by Tilman Barz </li>
          </ul>
          </html>"));
end phTransMod_CurveSwitchHysteresis_Algebraic;
