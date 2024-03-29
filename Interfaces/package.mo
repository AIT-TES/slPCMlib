﻿within slPCMlib;
package Interfaces "Interfaces for PCM properties and models"
    extends Modelica.Icons.InterfacesPackage;

  model phTransModMeltingCurve
    "Melting curve model, static (default model, no hysteresis)"

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
  end phTransModMeltingCurve;

  model phTransModCurveTrackHysteresis
    "Curve track hysteresis model, static"

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
           to simulate hysteresis!",   AssertionLevel.error);
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
      heatingOn := true;
    else
      heatingOn := false;
    end if;

  algorithm
     when  (indVar.T < uppLimPhTrange) then
        heatingOn :=false;
     end when;
     when (indVar.T > lowLimPhTrange) then
        heatingOn :=true;
     end when;

     if heatingOn then
       (xi, dxi) :=PCM.phaseFrac_complMelting(indVar.T);
       //   print("pos", der(indVar.T), time);
     else
       (xi, dxi) :=PCM.phaseFrac_complSolidification(indVar.T);
       //   print("  neg", der(indVar.T), time);
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
  end phTransModCurveTrackHysteresis;





  model phTransModCurveSwitchHysteresisAlgebraic
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
           to simulate hysteresis!",   AssertionLevel.error);
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
  end phTransModCurveSwitchHysteresisAlgebraic;



  model phTransModCurveSwitchHysteresisDifferentiated
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
          to simulate hysteresis!",   AssertionLevel.error);
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
  end phTransModCurveSwitchHysteresisDifferentiated;

  model phTransModCurveScaleHysteresisAlgebraic
    "Curve scale hysteresis model, static (algebraic equations)"

    // extends Modelica.Icons.ObsoleteModel;
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
           to simulate hysteresis!",   AssertionLevel.error);
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
  end phTransModCurveScaleHysteresisAlgebraic;

  model phTransModCurveScaleHysteresisDifferentiated
    "Curve scale hysteresis model, static (differentiated equations)"

    extends basicPhTransModel;

protected
    Real scalerM, scalerS;
    final constant Real eps = Modelica.Constants.small;
    Real xiH_at_T, dxiH_at_T;
    Real xiC_at_T, dxiC_at_T;
    discrete Boolean heatingOn(start=true);

  initial algorithm
    assert((PCM.propData.modelForMelting == true)
       and (PCM.propData.modelForSolidification == true),
           "There should be a model for melting AND solidification 
           to simulate hysteresis!",   AssertionLevel.error);
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
    when (indVar.T >= PCM.propData.rangeTsolidification[1])
      and (indVar.T <= PCM.propData.rangeTmelting[2])
      and (indVar.der_T > 0) then
        heatingOn :=true;
    end when;
    when (indVar.T >= PCM.propData.rangeTsolidification[1])
      and (indVar.T <= PCM.propData.rangeTmelting[2])
      and (indVar.der_T < 0) then
        heatingOn :=false;
    end when;

  equation
    (xiH_at_T, dxiH_at_T) = PCM.phaseFrac_complMelting(indVar.T);
    (xiC_at_T, dxiC_at_T) = PCM.phaseFrac_complSolidification(indVar.T);

    when (indVar.T <= PCM.propData.rangeTsolidification[1]) then
      reinit(xi, 0.0);
    elsewhen (indVar.T >= PCM.propData.rangeTmelting[2]) then
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
  end phTransModCurveScaleHysteresisDifferentiated;

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
