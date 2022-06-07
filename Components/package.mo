within slPCMlib;
package Components "Some component models for testing slPCMlib"

  annotation (Documentation(info="<html>
      <p>
      This package contains exemplarily thermal components which use PCM and 
      phase transition models. 
      </p>  
      <p>
        Select a &lt;PCM&gt; from <a href>slPCMlib.Media</a>, and a phase transition model 
        &lt;phTrModel&gt; from <a href>slPCMlib.Interfaces</a>. 
        Changes in cp are induced by temperature, use: 
        &lt; phTrModel.indVar.T = T; &gt;.      
      </html>",
    revisions="<html>
        <ul>
        <li>2022-06-01; initial version; by Tilman Barz </li>
        </ul>
        </html>"));
end Components;
