function [Wm,Lm,Wmt,Lmt,Wk,Lk,Wkt,Lkt] = RadialVertexSolutions(t,xi,Cp,Ep,Kp,mup,Q0)

  %M vertex, storage viscosity
  Wm = 1.1901*(mup^2*Q0^3*t(end)/Ep^2)^(1/9)*(1+xi).^0.487.*(1-xi).^(2/3);
  Lm = 0.6944*(Q0^3*Ep*t.^4/mup).^(1/9);

  %Mt vertex, leak-off viscosity
  Wmt = 1.0574*(mup^4*Q0^6*t(end)/Ep^4/Cp^2)^(1/16)*(1+xi).^0.397.*(1-xi).^(5/8);
  Lmt = 0.4502*(Q0^2*t./Cp^2).^(1/4);

  %K vertex, storage toughness
  Wk = 0.6537*(Kp^4*Q0*t(end)/Ep^4)^(1/5)*(1-xi.^2).^(1/2);
  Lk = 0.8546*(Ep*Q0*t/Kp).^(2/5);

  %Kt vertex, leak-off toughness
  Wkt = 0.4744*(Kp^8*Q0^2*t(end)/Ep^8/Cp^2)^(1/8)*(1-xi.^2).^(1/2);
  Lkt = 0.4502*(Q0^2*t/Cp^2).^(1/4);
  
end