within ;
package slPCMlib "Solid/liquid PCM thermophysical property 
  and phase transition models library"

  annotation (uses(Modelica(version="4.0.0"), ModelicaServices(version="4.0.0"),
    Buildings(version="9.0.0")),
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
      graphics={              Rectangle(
        extent={{-100,80},{100,-80}},
        lineColor={200,200,200},
        radius=25,
        fillColor={248,248,248},
        fillPattern=FillPattern.HorizontalCylinder,
        lineThickness=0.5),
      Line(
        points={{-56,-36},{-44,-36},{-28,-34},{-12,-24},{-4,-2},{0,26},{4,34},{
            8,36},{14,30},{18,16},{20,-22},{24,-34},{34,-34},{58,-34}},
        color={238,46,47},
        thickness=1,
        smooth=Smooth.Bezier),
      Line(
        points={{-80,66},{-80,-66},{78,-66}},
        color={28,108,200},
        thickness=1,
        arrow={Arrow.Filled,Arrow.Filled}),
      Line(
        points={{-80,36},{-68,36}},
        color={28,108,200},
        thickness=1),
      Line(
        points={{-80,-6},{-68,-6}},
        color={28,108,200},
        thickness=1),
      Line(
        points={{-80,-48},{-68,-48}},
        color={28,108,200},
        thickness=1)}),
  version="1",
  conversion(noneFromVersion=""));
end slPCMlib;
