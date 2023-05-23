within slPCMlib.Components;
model linearT "Generates temperature trajectory of linear segments"

  parameter Modelica.Units.SI.Time tValues[:]={0,2000,4000,6000};
  parameter Modelica.Units.SI.Temperature TValues[:](each displayUnit="degC") = {290,
    312,307,310};
  Modelica.Units.SI.Temperature Tinterpol(displayUnit="degC")
    "Port temperature";

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port
  annotation (Placement(transformation(
        origin={0,-100},
        extent={{-10,-10},{10,10}},
        rotation=90), iconTransformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={110,0})));

equation

  (Tinterpol,) =  Modelica.Math.Vectors.interpolate(tValues,TValues,time);
  port.T =  Tinterpol;

annotation (
  Icon(coordinateSystem(
      preserveAspectRatio=true,
      extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={1,1,1},
          fillColor={242,246,255},
          fillPattern=FillPattern.Solid),
      Line(points={{-74,68},{-74,-94}}, color={192,192,192},
        thickness=1),
      Polygon(
        points={{-74,90},{-82,68},{-66,68},{-74,90}},
        lineColor={192,192,192},
        fillColor={192,192,192},
        fillPattern=FillPattern.Solid),
      Line(points={{-86,-82},{72,-82}},
                                    color={192,192,192},
        thickness=1),
      Polygon(
        points={{94,-82},{72,-74},{72,-90},{94,-82}},
        lineColor={192,192,192},
        fillColor={192,192,192},
        fillPattern=FillPattern.Solid),
      Text(
        extent={{-50,14},{50,-14}},
        lineColor={200,200,200},
        origin={-88,16},
        rotation=90,
        fontSize=40,
        textString="temperature"),
      Text(
        extent={{-42,10},{42,-10}},
        lineColor={200,200,200},
        origin={52,-92},
        rotation=0,
        textString="time",
        fontSize=40),
      Text(
        extent={{-102,100},{192,62}},
        lineColor={0,0,255},
        textString="%name"),
    Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200}),
    Line(points={{-64,-20},{-22,46},{0,-44},{38,8},{76,-20}},    color={28,108,
              200})}),
  Diagram(coordinateSystem(
      preserveAspectRatio=true,
      extent={{-100,-100},{100,100}}), graphics={
      Line(points={{-80,-90},{-80,84}}, color={95,95,95}),
      Polygon(
        points={{-80,97},{-84,81},{-76,81},{-80,97}},
        lineColor={95,95,95},
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid),
      Line(points={{-99,-40},{85,-40}}, color={95,95,95}),
      Polygon(
        points={{97,-40},{81,-36},{81,-45},{97,-40}},
        lineColor={95,95,95},
        fillColor={95,95,95},
        fillPattern=FillPattern.Solid),
      Text(
        extent={{75,-47},{100,-60}},
        lineColor={0,0,0},
        textString="time"),
      Text(
        extent={{-80,99},{-40,82}},
        lineColor={0,0,0},
        textString="y"),
    Line(points={{-76,-10},{-34,56},{-12,-34},{26,18},{64,-10}}, color={28,108,
          200})}),
  Documentation(info="<html>
    <p>
    The component provides a heat port with variable temperature T. 
    </p>
    </html>",
  revisions="<html>
       <ul>
       <li>2022-06-01; initial version; by Tilman Barz </li>
       </ul>
       </html>"));
end linearT;
