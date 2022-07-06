within slPCMlib;
package Media_Rubitherm_SP "Contains different PCM with thermophysical properties for solid and liquid 
   phase and phase transition functions for complete melting/solidification 
   processes"
    extends Modelica.Icons.MaterialProperty;

  annotation (Documentation(info="<html>
      <p>
      This package contains thermophysical properties of various solid/liquid 
      PCM. It holds information on solid and liquid properties 
      (as function of temperature) and provides basic function to compute the 
      liquid mass (and based on this also volume) phase fractions for 
      complete phase transition processes 
      (complete melting or solidification).  
      </p>
      <p>
      It contains 
      <ul>
      <li>generic PCM for which the values can be adapted to match certain 
      behavior, </li>
      <li>specific PCM, which are based on experimental data, e.g. data taken 
      from manufacturer's data sheets. </li>
      </ul>
      </p>
      <p>      </p>
      </html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end Media_Rubitherm_SP;
