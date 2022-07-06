within slPCMlib.Components;
model SineT "Generates sinusoidal temperature variations with linear drift"

  import Modelica.Constants.pi;
//  parameter Real simulTime = 6000 "Simulation time";
  parameter Modelica.Units.SI.Temperature ampl(displayUnit="K") = 2.5
    "Amplitude of sine wave";
  parameter Modelica.Units.SI.Time timePeriod=1000 "Time to complete one cycle";
  parameter Modelica.Units.SI.TemperatureSlope rateT=1/1000
    "Linear temperature drift";
  parameter Modelica.Units.SI.Temperature startT(displayUnit="K") = 24 + 273.15
    "Start temperature";

    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port
    annotation (Placement(transformation(
          origin={0,-100},
          extent={{-10,-10},{10,10}},
          rotation=90), iconTransformation(
            extent={{-10,-10},{10,10}},
            rotation=90,
            origin={110,0})));

equation

     port.T = startT + rateT*time + ampl*Modelica.Math.sin(2*pi/timePeriod*time);
              // (t*frequ - pi/2)  - pi/2

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
          Line(
            points={{-67,0},{-57.6,36.2},{-52.1,55.1},{-47.3,68.4},{-43.1,76.6},
              {-38.9,81.1},{-34.64,81.8},{-30.42,78.6},{-26.201,71.7},{-21.98,61.4},
              {-17.16,46.1},{-11.1,23.2},{1.5,-28.8},{7,-48.2},{11.8,-62.2},{16,
              -71.1},{20.2,-76.4},{24.5,-78},{28.7,-75.6},{32.9,-69.5},{37.1,-59.9},
              {41.9,-45.2},{48,-22.8},{54,2}},
            color={28,108,200},
            smooth=Smooth.Bezier)}),
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
          Line(
            points={{-67,-10},{-57.6,26.2},{-52.1,45.1},{-47.3,58.4},{-43.1,
                66.6},{-38.9,71.1},{-34.64,71.8},{-30.42,68.6},{-26.201,61.7},{
                -21.98,51.4},{-17.16,36.1},{-11.1,13.2},{1.5,-38.8},{7,-58.2},{
                11.8,-72.2},{16,-81.1},{20.2,-86.4},{24.5,-88},{28.7,-85.6},{
                32.9,-79.5},{37.1,-69.9},{41.9,-55.2},{48,-32.8},{54,-8}},
            color={28,108,200},
            smooth=Smooth.Bezier),
          Text(
            extent={{75,-47},{100,-60}},
            lineColor={0,0,0},
            textString="time"),
          Text(
            extent={{-86,99},{-46,82}},
            lineColor={0,0,0},
          textString="T")}),
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
end SineT;
