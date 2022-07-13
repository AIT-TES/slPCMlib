within slPCMlib;
package UsersGuide "User's Guide"
extends Modelica.Icons.Information;

annotation (
Documentation(info="<html>
<p>
The library <strong><u>slPCMlib</u> - 
solid/liquid PCM Modelica library -</strong> contains property models for 
<strong>solid/liqid phase change materials (PCM)</strong> 
showing a 
<strong>non-isothermal phase transition behavior</strong>. 
</p>
<p>
The library contains <strong>generic PCM</strong> and 
<strong>specific commercial PCM</strong>
(media). 
Different phenomenological <strong>phase transition models</strong> 
are implemented 
to account for temperature shifts in latent transition changes, 
e.g. due to multi-step transitions and thermal hysteresis. 
The models predict <strong>effective properties</strong> 
which are valid over the
PCM functional temperature range where latent heat is absorbed and released. 
Based on the properties and adopting the apparent heat capacity method, 
heat transfer problems can be solved 
numerically.  
<br>
</p>
<blockquote> 
  <p>
  <img src=\"modelica://slPCMlib/Resources/Images/slPCMlib.png\"> 
  Selecting material data (either generic or specific PCM), 
  and a phase transition model, heat transfer problems in PCM can be solved.
  </p>
  </blockquote>
  <p>
Assumptions for modeling effective PCM properties:
<ul>
  <li>There are only two phases (two-phase model): a
solid and a liquid phase.</li> 
  <li>Phase transitions are induced by temperature and are
independent of pressure.</li> 
  <li>Phase transitions extend over a temperature range
(non-isothermal phase transitions) and are continuous.</li> 
  <li>Within the phase transition temperature range the
  solid and liquid phases coexist as a homogenous mixture (macroscopic view). 
  The material is then in a semi-solid or semi-liquid state which produces a
mushy zone in the PCM domain.</li> 
  <li>Properties of the mushy state are local
effective (also apparent) mixture properties, which
are defined by a weighting of contributions from
solid and liquid phases. The weighting is based on
the phase change progress, i.e. the mass (or volume) phase fraction.  </li> 
        </ul>
        </p>
      <p>
      Temperature is the input to the model using the 
      &lt; inductionAtNode &gt; connector. 
      For a given temperature input <var> T </var>  
      the liquid mass phase fraction <var> xi </var> is 
      computed, 
      and the following variables are derived:
      </p>  
      <ul>
      <li>liquid volume phase fractions <var> phi </var> </li> 
      <li>effective density <var> rho </var> </li>
      <li>effective thermal conductivity <var> lambda </var> </li>  
      <li>effective specific heat capacity <var> cp </var></li>   
      <li>and specific enthalpy <var> h </var></li>
      <li>baseline heat capacity <var> c_BL </var>, which 
      describes the mixture heat capacity (without the effect of the 
      phase transition enthalpy)</li>
      <li>solid <var> h_S </var> and liquid <var> h_L </var> enthalpies, 
      where the difference <var> h_L - h_S </var> 
      is the temperature-dependent phase transition enthalpy. 
      </ul>
      </p> 
      <p>
      The <a href>slPCMlib.Media</a> package contains 
      <ul>
      <li>generic PCM, for which the user can make copies and adapt 
      the parameters to match certain behavior, </li>
      <li>specific PCM, which are based on experimental data, e.g. data taken 
      from manufacturer's data sheets. </li>
      </ul>
      </p>
      <p>
      The <a href>slPCMlib.Interfaces</a> package contains 
      a PCM media template <a href>slPCMlib.Interfaces.partialPCM</a>
      with partial PCM models, parameter, functions for the solid and liquid 
      properties (single phases) and the evolution of the liquid mass 
      fraction <var> xi </var> 
      for complete solid/liquid or liquid/solid transitions.
      </p>
      <p>
      The <a href>slPCMlib.Interfaces</a> package also contains 
      a basic phase transition model 
      <a href>slPCMlib.Interfaces.basicPhTransModel</a>, which computes the 
      variables given above.  
      The following models 
      extend the basic model: 
      <ul>
      <li>melting curve model 
      <a href>slPCMlib.Interfaces.phTransModMeltingCurve</a> 
      - rate-independent, no hysteresis </li>
      <li>curve track hysteresis model 
      <a href>slPCMlib.Interfaces.phTransModCurveTrackHysteresis</a> 
      - rate-independent, hysteresis</li>
      <li>curve switch hysteresis model 
      <a href>slPCMlib.Interfaces.phTransModCurveSwitchHysteresis</a> 
      - rate-independent, hysteresis</li>      
      <li>curve scale hysteresis model
      <a href>slPCMlib.Interfaces.phTransModCurveScaleHysteresisXXX</a>  with 
      <a href>XXX = {Algebraic, Differentiated}</a>
      - rate-independent, hysteresis</li>      
      </ul>      
      </p>
      <p>
      These phase transition models predict a different behavior for 
      complete and incomplete (interrupted) transitions, depending on the 
      temperature history. The behavior is predicted for <var> xi </var> 
      and affects all properties (variables above). 
      </p>
      <p>
      <strong>Important note on the application:</strong> 
      The library <u>slPCMlib</u> computes effective (apparent) PCM properties 
      for a given temperature input <var> T </var>. 
      The properties can be used to model the heat transfer in a PCM domain.  
      The heat transfer differential (balance) equations 
      can be formulated adopting the 
      apparent heat capacity method. 
      <br>
      The library is not designed/tested to be used with equations 
      formulated adopting the enthalpy model, which would need inverse 
      relation <var>T=f(h)</var>.
      </p>
      <p>      
      The package <a href>slPCMlib.Components</a> contains  
      thermal storage models using PCM.  
      Examples for their use are contained in 
      <a href>slPCMlib.Examples</a>.
      </p>
  </html>",
revisions="<html>
   <ul>
   <li>2022-06-01; initial version; by Tilman Barz </li>
   </ul>
   </html>"));
end UsersGuide;
