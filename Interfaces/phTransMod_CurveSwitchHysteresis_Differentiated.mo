within slPCMlib.Interfaces;
model phTransMod_CurveSwitchHysteresis_Differentiated
  "Curve switch hysteresis model, static (differentiated equations)"

  extends basicPhTransModel;

protected
  discrete Integer modelInd(start=1);
  Modelica.Units.SI.MassFraction xiC_at_T(start=0.5);
  Modelica.Units.SI.MassFraction xiH_at_T(start=0.5);
  Real dxiC_at_T(start=0.), dxiH_at_T(start=0.);

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
  else
    modelInd := -1;
    (xi,)  := PCM.phaseFrac_complSolidification(indVar.T);
  end if;

equation
  (xiC_at_T, dxiC_at_T)   = PCM.phaseFrac_complSolidification(indVar.T);
  (xiH_at_T, dxiH_at_T)   = PCM.phaseFrac_complMelting(indVar.T);

algorithm

  // --- leave the temperature range ---
  when    (indVar.T > PCM.propData.rangeTsolidification[1]) then
    modelInd :=1;  // activate melting curve
    // print("  ->Low", modelInd, time);
  end when;
  when    (indVar.T < PCM.propData.rangeTmelting[2]) then
    modelInd :=-1;  // activate solidification curve
    // print("  ->Up", modelInd, time);
  end when;

  // --- switch between heating and cooling ---
  //     during phase transition
  when (der(indVar.T) > 0) then
    if (( (PCM.propData.rangeTsolidification[1] < indVar.T)
      and (indVar.T < PCM.propData.rangeTmelting[2]))) then
      modelInd  :=0;  // switch: transition from melting to solidification curve
      // print("  ->pos sign . change in the phTrRg", modelInd, time);
    end if;
  end when;
  when (der(indVar.T) <= 0) then
    if ( (PCM.propData.rangeTsolidification[1] < indVar.T)
      and (indVar.T < PCM.propData.rangeTmelting[2])) then
      modelInd  :=0;  // switch: transition from solidification to melting curve
      // print("  ->neg sign . change in the phTrRg", modelInd, time);
    end if;
  end when;

  // --- finish the curve switch,
  //     i.e. apporach the melting/solidification curve during phase switch ---
  when (xi < xiH_at_T) then
    // if during switch transition
    if (modelInd==0) then
      modelInd :=1;   // activate melting curve
      //print("  ->finish curve switch - now heating curve", modelInd, time);
    end if;
  end when;
  when (xi > xiC_at_T) then
    // if during switch transition
    if (modelInd==0) then
      modelInd :=-1;   // activate cooling curve
      //print("  ->finish curve switch - now cooling curve", modelInd, time);
    end if;
  end when;

equation

   when (modelInd==1) then
    // when (modelInd>0) then
      reinit(xi, xiH_at_T);
    elsewhen (modelInd==-1) then
    // elsewhen (modelInd<0) then
     reinit(xi, xiC_at_T);
    end when;

  if (modelInd==1) then
    dxi = dxiH_at_T;
    der(xi) = dxiH_at_T*indVar.der_T;
    // print("pos", der(indVar.T), time);
  elseif (modelInd==-1) then
    dxi = dxiC_at_T;
    der(xi) = dxiC_at_T*indVar.der_T;
    // print("  neg", der(indVar.T), time);
  else
    der(xi) = 0.0;
    dxi = 0.0;
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
end phTransMod_CurveSwitchHysteresis_Differentiated;
