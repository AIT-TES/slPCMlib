within slPCMlib;
package UsersGuide "User's Guide"
extends Modelica.Icons.Information;

annotation (
Documentation(info="<html>
<p>
This library <u>slPCMlib</u> contains property models for 
<strong>solid/liqid phase change materials (PCM)</strong>. 
</p>
<blockquote> 
  <p>
  <img src=\"modelica://slPCMlib/Resources/Images/slPCMlib.png\"> 
  </p>
  </blockquote>
<p>
Assumptions:
  <ul>
  <li>Phase transitions are induced by temperature.</li> 
  <li>Phase transitions are pressure independent.</li>  
  <li>Only two phases exist (two-phase model), a solid and a liquid phase.</li> 
  <li>These phases co-exist as homogenous mixture 
  (macroscopic view) over an extended phase transition temperature 
  range (non-isothermal phase transitions).</li> 
  <li>Effective (also apparent) mixture properties can be computed by 
  a simple linear weighting of contributions from solid and liquid phase. 
  The weighting factor is the mass (or volume) liquid phase fraction.</li>   
        </ul>
        </p>
      <p>
      Temperature is the input to the model using the 
      &lt; inductionAtNode &gt; connector. 
      For a given temperature input signal, liquid phase fractions are 
      computed, 
      and the following variables are derived, which are functions of 
      temperature <var> T </var> and liquid mass phase fraction <var> xi_m </var>. 
      </p>  
      <ul>
          <li>liquid mass <var> xi_m </var> and volume<var> xi_v </var> 
              phase fractions </li>         
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
      fraction <var> xi_m </var> 
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
      temperature history. The behavior is predicted for <var> xi_m </var> 
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
      <a href>slPCMlib.Examples</a>
      </p>
  </html>",
revisions="<html>
   <ul>
   <li>2022-06-01; initial version; by Tilman Barz </li>
   </ul>
   </html>"));
end UsersGuide;
