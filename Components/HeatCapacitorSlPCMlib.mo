within slPCMlib.Components;
model HeatCapacitorSlPCMlib "Lumped thermal PCM element storing heat"

  parameter Modelica.Units.SI.Mass m(displayUnit="kg") = 0.100
    "Mass of PCM element";
  Modelica.Units.SI.Temperature T(start=293.15, displayUnit="degC")
    "Temperature of element";
  Modelica.Units.SI.TemperatureSlope der_T(start=0)
    "Time derivative of temperature (= der(T))";
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port
    annotation (Placement(transformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=90)));

  replaceable package PCM =
    slPCMlib.Media_generic.generic_7thOrderSmoothStep
    constrainedby slPCMlib.Interfaces.partialPCM
    annotation (Dialog(group="PCM and phase transition model"),
    choicesAllMatching=true);

  replaceable slPCMlib.Interfaces.phTransMod_MeltingCurve_Differentiated phTrModel(
      redeclare package PCM = PCM) constrainedby
    slPCMlib.Interfaces.basicPhTransModel(redeclare package PCM = PCM)
    annotation (Dialog(group="PCM and phase transition model"),
      choicesAllMatching=true);

equation
  // assign temperatures to HeatCapacitor
  T = port.T;
  der_T = der(T);

  // input temperature signal to the model
  phTrModel.indVar.T = T;
  phTrModel.indVar.der_T = der_T;

  phTrModel.cp*m*der(port.T) = port.Q_flow;

    annotation (
      Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
              100,100}}), graphics={
          Text(
            extent={{-150,110},{150,70}},
            lineColor={0,0,255},
          textString="%name"),
          Polygon(
            points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,
                13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},
                {-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,
                -89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{
                70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,
                33},{44,41},{36,57},{26,65},{0,67}},
            lineColor={160,160,164},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{
                -76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,
                -83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,
                -77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,
                -73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},
                {-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,
                27},{-48,35},{-44,45},{-40,57},{-58,35}},
            lineColor={0,0,0},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-69,7},{71,-24}},
            lineColor={0,0,0},
            textString="m * cp(T)
(m = %m)")}),
      Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
              {100,100}}), graphics={
          Polygon(
            points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,
                13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},
                {-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,
                -89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{
                70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,
                33},{44,41},{36,57},{26,65},{0,67}},
            lineColor={160,160,164},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{
                -76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,
                -83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,
                -77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,
                -73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},
                {-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,
                27},{-48,35},{-44,45},{-40,57},{-58,35}},
            lineColor={0,0,0},
            fillColor={160,160,164},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-6,-1},{6,-12}},
            lineColor={255,0,0},
            fillColor={191,0,0},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{11,13},{50,-25}},
            lineColor={0,0,0},
            textString="T"),
          Line(points={{0,-12},{0,-96}}, color={255,0,0})}),
    Documentation(info="<html>
        <p>
        This is a generic model for the heat capacity of a PCM. 
        </p><p>    
        It is a modification of 
        <a href=\"modelica://Modelica.Thermal.HeatTransfer.Components.HeatCapacitor\">
        Modelica.Thermal.HeatTransfer.Components.HeatCapacitor</a>        
        and can be used with materials and 
        phase transition models contained in <strong>slPCMlib</strong> . 
        </p>
        <p>
        The model and assumptions (lumped-element model with uniform temperature, etc.) 
        are the same as for the <a href>HeatCapacitor</a>. 
        <br>
        The difference is, that the specific heat capacity is <strong>NOT</strong> 
        constant. 
        The so-called <strong>apparent (or effective)</strong> specific heat capacity 
        <strong>cp</strong> is used, which is a function of temperature. 
        </p>
        <p>
        Select a &lt;PCM&gt; from  
        <a href=\"modelica://slPCMlib.Media\">
        slPCMlib.Media</a>, 
        and a phase transition model 
        &lt;phTrModel&gt; from <a href>slPCMlib.Interfaces</a>. 
        Changes in cp are induced by temperature, use: 
        &lt; phTrModel.indVar.T = T; &gt;. 
        </html>",
  revisions="<html>
       <ul>
       <li>2022-06-01; initial version; by Tilman Barz </li>
       </ul>
       </html>"));
end HeatCapacitorSlPCMlib;
