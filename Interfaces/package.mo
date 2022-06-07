within slPCMlib;
package Interfaces "Interfaces for PCM properties and models"
    extends Modelica.Icons.InterfacesPackage;

  replaceable model phTransModMeltingCurve
  "Melting curve model, static (default model, no hysteresis)"
  extends basicPhTransModel;
  equation
    (xi_m, dxi_m)  = PCM.phaseFrac_complMelting(indVar.T);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false),
      graphics={Text(lineColor={108,88,49},
      extent={{-90.0,-90.0},{90.0,90.0}},
      textString="→")}),
      Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
	<p>
	The model predicts the liquid mass phase fraction as function of 
	temperature, see the following example.
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
  end phTransModMeltingCurve;

  replaceable model phTransModCurveTrackHysteresis
    "Curve track hysteresis model, static"

    extends basicPhTransModel;

    discrete Boolean heatingOn(start=true);

  initial algorithm
    if  (indVar.der_T >= 0) then
      heatingOn := true;
    else
      heatingOn := false;
    end if;

  algorithm
     when (heatingOn==true)  and (indVar.T > PCM.propData.rangeTmelting[2]) then
        heatingOn :=false;
     end when;
     when (heatingOn==false) and (indVar.T < PCM.propData.rangeTsolidification[1]) then
        heatingOn :=true;
     end when;

  // --- conditional if-equations
     if noEvent(heatingOn) then
       (xi_m, dxi_m) :=PCM.phaseFrac_complMelting(indVar.T);
       //   print("pos", der(indVar.T), time);
     else
       (xi_m, dxi_m) :=PCM.phaseFrac_complSolidification(indVar.T);
       //   print("  neg", der(indVar.T), time);
     end if;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Text(lineColor={108,88,49},
                extent={{-90.0,-90.0},{90.0,90.0}},
          textString="⇆")}),
                Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
          <p>
          The model predicts the liquid mass phase fraction. 
          It is a static (so-called curve track) hysteresis model, see e.g. 
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
  end phTransModCurveTrackHysteresis;

  replaceable model phTransModCurveScaleHysteresisAlgebraic
    "Curve scale hysteresis model, static (algebraic equations)"
  extends Modelica.Icons.ObsoleteModel;
    extends basicPhTransModel;
  //  extends Modelica.Icons.ObsoleteModel;

  //protected
    discrete Real T0(start=0.0);
    discrete Real Xi0(start=1.0); //(fixed=true, start=0.5);
    discrete Real Xi_at_T0(start=1.0); //(fixed=true, start=0.5);
    discrete Real scaler(start=1.0); //(fixed=true, start=0.5);
    discrete Boolean posRate(start=true); //(,fixed=true, start=true);
    final constant Real eps = Modelica.Constants.small;
    Real Xi_at_T(start=1.0), dXi_at_T(start=1.0);

  //  Real der_T=indVar.der_T;
  initial algorithm

    if  (indVar.der_T >= 0) then
      (xi_m,)  := PCM.phaseFrac_complMelting(indVar.T);
      T0  :=pre(indVar.T);
      Xi0 :=pre(xi_m);
      Xi_at_T0  := Xi0;
      scaler := 1.0; //(1.0 - Xi0)/max((1.0 - Xi_at_T0),eps);
      posRate := true;
    else
      (xi_m,)  := PCM.phaseFrac_complSolidification(indVar.T);
      T0  :=pre(indVar.T);
      Xi0 :=pre(xi_m);
      Xi_at_T0 := Xi0;
      scaler := 1.0; //Xi0/max(Xi_at_T0, eps);
      posRate := false;
    end if;

  algorithm
    when (indVar.der_T > 0) then
      T0   :=pre(indVar.T);
      Xi0  :=pre(xi_m); // unclocked discrete time variables
                        // use pre() if it is algebraic (and not calculated before?),
                        // not previous() which is for clocked variables
      (Xi_at_T0,)  :=PCM.phaseFrac_complMelting(T0);
      if (Xi_at_T0 >= Xi0) then
        scaler :=1.0;
  //     elseif (Xi0 > 1.0) then
  //      scaler := 0.0;
      else
        scaler :=(1.0 - Xi0)/max((1.0 - Xi_at_T0),eps);
      end if;
      posRate :=true;
    end when;
    when  (indVar.der_T < 0) then
      T0   :=pre(indVar.T);
      Xi0  :=pre(xi_m);
      (Xi_at_T0,)  :=PCM.phaseFrac_complSolidification(T0);
      if (Xi_at_T0 <= Xi0) then // + 1e-6
        scaler :=1.0;
      // elseif noEvent(Xi0 < 0.0) then
      //   scalerS = 0.0;
      else
        scaler :=Xi0/max(Xi_at_T0, eps);
      end if;
      posRate :=false;
    end when;

    if noEvent(posRate) then //noEvent(posRate)
      (Xi_at_T, dXi_at_T) :=PCM.phaseFrac_complMelting(indVar.T);
      xi_m     :=1.0 - scaler*(1.0 - Xi_at_T);
      dxi_m    :=scaler*dXi_at_T;
        //   print("pos", indVar.der_T, time);
        //   shifth := 0.0; //(iXiM - ((T-Tref) - scaler*(T-Tref) + scaler*iXiM))*(cp_liquid-cp_solid);
    else
      (Xi_at_T, dXi_at_T) :=PCM.phaseFrac_complSolidification(indVar.T);
      xi_m     :=scaler*Xi_at_T;
      dxi_m    :=scaler*dXi_at_T;
      //  PCMlib.print("  neg XiAtT0 = ", Xi_at_T0, T0);
      //    print("  neg", indVar.der_T, time);
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Text(lineColor={108,88,49},
                extent={{-90.0,-90.0},{90.0,90.0}},
                textString="↔")}),
                Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
          <p>
          The model predicts the liquid mass phase fraction. 
          It is a static (so-called curve scale) hysteresis model as described in  
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
  end phTransModCurveScaleHysteresisAlgebraic;


  replaceable model phTransModCurveScaleHysteresisDifferentiated
  "Curve scale hysteresis model, static (differentiated equations)"

    extends basicPhTransModel;

  //protected
    Real scalerM, scalerS;
   // Real scalerM1, scalerS1;
    final constant Real eps = Modelica.Constants.small;
    Real XiM_at_T, dXiM_at_T;
    Real XiS_at_T, dXiS_at_T;
  //  discrete Boolean isInRange;

  initial algorithm
     if  (indVar.der_T >= 0) then
       // unclocked discrete time variables
       // use pre() if it is algebraic (and not calculated before?),
       // not previous() which is for clocked variables

       (xi_m,)  := PCM.phaseFrac_complMelting(indVar.T);
     else
       (xi_m,)  := PCM.phaseFrac_complSolidification(indVar.T);
     end if;

  equation
    (XiM_at_T, dXiM_at_T) = PCM.phaseFrac_complMelting(indVar.T);
    (XiS_at_T, dXiS_at_T) = PCM.phaseFrac_complSolidification(indVar.T);

  //  if (indVar.T >= PCM.propData.rangeTsolidification[1])
  //      and (indVar.T <= PCM.propData.rangeTmelting[2]) then
  //    isInRange = true;
  //  else
  //    isInRange = false;
  //  end if;

     when (indVar.T <= PCM.propData.rangeTsolidification[1]) then
       reinit(xi_m, 0.0);
     elsewhen (indVar.T >= PCM.propData.rangeTmelting[2]) then
       reinit(xi_m, 1.0);
     end when;

    // alternative hysteresis
    // scalerM = ( XiS_at_T - Xi)
    //         / ( XiS_at_T - XiM_at_T+eps);
    // scalerS = ( Xi - XiM_at_T)
    //         / ( XiS_at_T - XiM_at_T + eps);

  // -> das wieder rein
  //       scalerM  = (1.0 - xi_m)/max((1.0 - XiM_at_T), eps);
  //       scalerS  = xi_m/max(XiS_at_T, eps);

    // -> ab hier zurück
    if noEvent(XiM_at_T >= xi_m) then
      scalerM = 1.0;
    elseif noEvent(xi_m > 1.0) then
      scalerM = 0.0;
    else
      scalerM  = (1.0 - xi_m)/max((1.0 - XiM_at_T), eps);
    end if;
    if noEvent(XiS_at_T <= xi_m) then
      scalerS = 1.;
    elseif noEvent(xi_m < 0.0) then
      scalerS = 0.0;
    else
      scalerS  = xi_m/max(XiS_at_T, eps);
    end if;
    // <- ab hier zurück

    // --- direction changes ---
    if (indVar.der_T >= 0) then
      der(xi_m) = scalerM*dXiM_at_T*indVar.der_T; // - max(xi_m-1.0,0.0); // + max(XiM_at_T-xi_m,0.0);
      dxi_m     = scalerM*dXiM_at_T;
      //   print("pos", indVar.der_T, time);
    else
      der(xi_m) = scalerS*dXiS_at_T*indVar.der_T;
      dxi_m     = scalerS*dXiS_at_T;
      //  PCMlib.print("  neg XiAtT0 = ", Xi_at_T0, T0);
      //    print("  neg", der(T), time);
    end if;
  //      end when;
  //      else
  //
  //        der(xi_m) :=0.0;
  //        dxi_m :=0.0;
  //
  //      end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Text(lineColor={108,88,49},
                extent={{-90.0,-90.0},{90.0,90.0}},
                textString="↔")}),
                Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
          <p>
          The model predicts the liquid mass phase fraction 
          as described in  
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
  end phTransModCurveScaleHysteresisDifferentiated;

  replaceable model phTransModCurveSwitchHysteresis
                    "Curve switch hysteresis model, static"
  extends Modelica.Icons.ObsoleteModel;
    extends basicPhTransModel;
    discrete Integer modelInd(start=1);
    discrete Real Xi0(start=0.5);
protected
    Real Xi_s, Xi_l;

  initial algorithm
    if  (der(indVar.T) >= 0) then
      modelInd := 1;
      (xi_m,)  := PCM.phaseFrac_complMelting(indVar.T);
      Xi0 := pre(xi_m);

    else
      modelInd := -1;
      (xi_m,)  := PCM.phaseFrac_complSolidification(indVar.T);
      Xi0 := pre(xi_m);
    end if;


    // --- on leaving the phase transition temperature range ---
    //    when      (modelInd <> 1) and  (indVar.T < PCM.propData.rangeTsolidification[1]) then

  algorithm
    when        (indVar.T < PCM.propData.rangeTsolidification[1]) then
    //  when     (switch2heating) then
      modelInd := 1; // activate melting curve
    //   print("  ->Low", modelInd, time);
      //    elsewhen  (modelInd <>-1) and  (indVar.T > PCM.propData.rangeTmelting[2]) then
    end when;
    when    (indVar.T > PCM.propData.rangeTmelting[2]) then
    //  else when (switch2cooling) then
      modelInd := -1; // activate solidification curve
    //   print("  ->Up", modelInd, time);
    end when;

    // --- switch between melting/solidification curve ---
    //  when (switch2posRate) then
    when (der(indVar.T) > 0) then
      // if during phase transition
      if ( (PCM.propData.rangeTsolidification[1] < indVar.T) and (indVar.T < PCM.propData.rangeTmelting[2])) then
        modelInd  := 0; // switch: transition from melting to solidification curve
        (Xi0,)  := PCM.phaseFrac_complSolidification(indVar.T); // unclocked discrete time variables
        //print("  ->pos sign . change in the phTrRg", modelInd, time);
      end if;
    end when;
    //  elsewhen (switch2negRate) then
    when (der(indVar.T) <= 0) then
      // if during phase transition
        if ( (PCM.propData.rangeTsolidification[1] < indVar.T) and (indVar.T < PCM.propData.rangeTmelting[2])) then
        modelInd  := 0; // switch: transition from cooling to heating curve
        (Xi0,)  := PCM.phaseFrac_complMelting(indVar.T); // unclocked discrete time variables
        //print("  ->neg sign . change in the phTrRg", modelInd, time);
      end if;
    end when;

    // --- finish the curve switch,
    //     i.e. apporach the melting/solidification curve during phase switch ---

    //  when (finishSwitchPos) then
    (Xi_s,)   := PCM.phaseFrac_complSolidification(indVar.T);
    (Xi_l,)   := PCM.phaseFrac_complMelting(indVar.T);
    when (Xi0 < Xi_l) then
      // if during switch transition
      if ( modelInd==0) then
        modelInd := 1;  // activate melting curve
        //print("  ->finish curve switch - now heating curve", modelInd, time);
      end if;
    end when;
    //  elsewhen (finishSwitchNeg) then
    when (Xi0 > Xi_s) then
        // if during switch transition
      if ( modelInd==0) then
        modelInd := -1;  // activate cooling curve
        //print("  ->finish curve switch - now cooling curve", modelInd, time);
      end if;
    end when;

    // --- now update Xi ---
    // --- conditional if-equations
    //if noEvent(
    if noEvent(modelInd==1) then
      (xi_m, dxi_m) := PCM.phaseFrac_complMelting(indVar.T);
      //   print("pos", der(indVar.T), time);
    elseif noEvent(modelInd==-1) then
      (xi_m, dxi_m) := PCM.phaseFrac_complSolidification(indVar.T);
      //    print("  neg", der(indVar.T), time);
    elseif noEvent(modelInd==0) then
      xi_m  := Xi0;
      dxi_m := 0.0;
    end if;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Text(lineColor={108,88,49},
                extent={{-90.0,-90.0},{90.0,90.0}},
                textString="⇆")}),
                Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
          <p>
          The model predicts the liquid mass phase fraction. 
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
  end phTransModCurveSwitchHysteresis;


  replaceable model phTransModCurveSwitchHysteresisDifferentiated
    "Curve switch hysteresis model, static"

    extends basicPhTransModel;
    discrete Integer modelInd(start=1);
   // discrete Real Xi0(start=0.5);
    Modelica.SIunits.MassFraction xi_SC(start=0.5), xi_MC(start=0.5);
    Real dxi_SC(start=0.), dxi_MC(start=0.);

  initial algorithm
    if  (der(indVar.T) >= 0) then
      modelInd := 1;
      (xi_m,)  := PCM.phaseFrac_complMelting(indVar.T);
   //   Xi0 := pre(xi_m);

    else
      modelInd := -1;
      (xi_m,)  := PCM.phaseFrac_complSolidification(indVar.T);
    //  Xi0 := pre(xi_m);
    end if;
  // --- on leaving the phase transition temperature range ---
  //    when      (modelInd <> 1) and  (indVar.T < PCM.propData.rangeTsolidification[1]) then

  equation
  //algorithm
    (xi_SC, dxi_SC)   = PCM.phaseFrac_complSolidification(indVar.T);
    (xi_MC, dxi_MC)   = PCM.phaseFrac_complMelting(indVar.T);

  algorithm

    // leave the temperature range
    when        (indVar.T < PCM.propData.rangeTsolidification[1]) then
      modelInd :=1;  // activate melting curve
    //   print("  ->Low", modelInd, time);
      //    elsewhen  (modelInd <>-1) and  (indVar.T > PCM.propData.rangeTmelting[2]) then
    end when;
    when    (indVar.T > PCM.propData.rangeTmelting[2]) then
      modelInd :=-1;  // activate solidification curve
    //   print("  ->Up", modelInd, time);
    end when;

    // --- switch between heating and cooling ---
    //     during phase transition
    when (der(indVar.T) > 0) then
      if noEvent(( (PCM.propData.rangeTsolidification[1] < indVar.T)
        and (indVar.T < PCM.propData.rangeTmelting[2]))) then
        modelInd  :=0;  // switch: transition from melting to solidification curve
        //print("  ->pos sign . change in the phTrRg", modelInd, time);
      end if;
    end when;
    when (der(indVar.T) <= 0) then
      if noEvent( (PCM.propData.rangeTsolidification[1] < indVar.T)
        and (indVar.T < PCM.propData.rangeTmelting[2])) then
        modelInd  :=0;  // switch: transition from solidification to melting curve
        //print("  ->neg sign . change in the phTrRg", modelInd, time);
      end if;
    end when;

    // --- finish the curve switch,
    //     i.e. apporach the melting/solidification curve during phase switch ---
     when (xi_m < xi_MC) then
      // if during switch transition
      if noEvent(modelInd==0) then
        modelInd :=1;   // activate melting curve
        //print("  ->finish curve switch - now heating curve", modelInd, time);
      end if;
    end when;
    when (xi_m > xi_SC) then
        // if during switch transition
      if noEvent(modelInd==0) then
        modelInd :=-1;   // activate cooling curve
        //print("  ->finish curve switch - now cooling curve", modelInd, time);
      end if;
    end when;

  equation

     when noEvent(modelInd==1) then
       reinit(xi_m, xi_MC);
     elsewhen noEvent(modelInd==-1) then
       reinit(xi_m, xi_SC);
     end when;

    // --- now update Xi ---
    if noEvent(modelInd==1) then
      dxi_m = dxi_MC;
      der(xi_m) = dxi_MC*indVar.der_T; // - max(xi_m-1.0,0.0); // + max(XiM_at_T-xi_m,0.0);
      //   print("pos", der(indVar.T), time);
    elseif noEvent(modelInd==-1) then
      dxi_m = dxi_SC;
      der(xi_m) = dxi_SC*indVar.der_T;
      //    print("  neg", der(indVar.T), time);
    else //if noEvent(modelInd==0) then
      der(xi_m) = 0.0;
      dxi_m = 0.0;
    end if;


    annotation (Icon(coordinateSystem(preserveAspectRatio=false),
                graphics={Text(lineColor={108,88,49},
                extent={{-90.0,-90.0},{90.0,90.0}},
                textString="⇆")}),
                Diagram(graphics, coordinateSystem(preserveAspectRatio=false)),
      Documentation(info="<html>
          <p>
          The model predicts the liquid mass phase fraction. 
          It is a static (so-called curve switch) 
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
  end phTransModCurveSwitchHysteresisDifferentiated;

annotation (Documentation(info="<html>
      <p>
      This package contains the interfaces for the PCM media and phase 
      transition models, and the implementation of 
      different phase transition models.
      </p>
      <p>   
      </html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end Interfaces;
