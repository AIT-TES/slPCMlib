within slPCMlib;
package Media_Rubitherm_GR "Contains PCM properties for solid and liquid 
   phases, and phase transition functions for complete melting/solidification processes"
    extends Modelica.Icons.MaterialProperty;

  annotation (Documentation(info="<html> 
      <p> 
      This package contains thermophysical properties of various solid/liquid  
      PCM. It holds information on solid and liquid properties  
      and provides phase transition functions to model the  
      liquid mass (and liquid volume) phase fractions during  
      complete melting and/or solidification.   
      </p> 
      <p> 
      The Media_* packages contain either  
      <ul> 
        <li>generic PCM for which the values can be adapted to match certain  
            behavior, or </li> 
        <li>specific PCM which are based on experimental data, e.g. data taken  
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
end Media_Rubitherm_GR;
