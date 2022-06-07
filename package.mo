within ;
package slPCMlib "Solid/liquid PCM thermophysical property 
  and phase transition models library"

  annotation (uses(Modelica(version="3.2.3"), ModelicaServices(version="3.2.3")),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
      graphics={              Rectangle(
              extent={{-100,80},{100,-80}},
              lineColor={253,212,198},
              radius=25,
              fillColor={235,238,219},
              fillPattern=FillPattern.Solid,
              lineThickness=0.5,
              pattern=LinePattern.None),
        Ellipse(
          extent={{-100,80},{100,-80}},
          lineColor={108,88,49},
          lineThickness=1),
        Polygon(
          points={{68,72},{84,44},{54,54},{68,72}},
          lineColor={108,88,49},
          lineThickness=1,
          fillColor={108,88,49},
          fillPattern=FillPattern.Solid),
      Line(
        points={{-64,-30},{-46,-30},{-30,-28},{-14,-18},{-6,4},{-2,32},{2,40},{
            6,42},{12,36},{16,22},{18,-16},{22,-28},{32,-28},{56,-28}},
        color={108,88,49},
        thickness=1,
        smooth=Smooth.Bezier),
        Polygon(
          points={{-68,-72},{-84,-44},{-54,-54},{-68,-72}},
          lineColor={108,88,49},
          lineThickness=1,
          fillColor={108,88,49},
          fillPattern=FillPattern.Solid)}));
end slPCMlib;
