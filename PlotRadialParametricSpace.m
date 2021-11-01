function PlotRadialParametricSpace(t,Cp,Ep,Kp,mup,Q0)

  %dimensionless parameters
  tmk = (mup^5*Ep^(13)*Q0^3/Kp^(18))^(1/2);
  tau = t(end)/tmk;
  phi = mup^3*Ep^(11)*Cp^4*Q0/Kp^(14);

  %size of the parameteric space
  taumin = -10;
  taumax = 20;
  phimin = -30;
  phimax = 20;

  %tau mk
  taumkmin = log10(4.54*1e-2);
  taumkmax = log10(2.59*1e6);

  %tau*phi^(5/6) kkt
  taukktmin = log10(5.96*1e-8);
  taukktmax = log10(4.81*1e2);

  %tau*phi^(-1/2) mtkt
  tauktmtmin = log10(4.18);
  tauktmtmax = log10(2.01*1e11);

  %tau*phi^(9/14) mmt
  taummtmin = log10(7.41*1e-6);
  taummtmax = log10(7.20*1e2);


  figure;
  hold on;

  %k
  phik = (taukktmin-taumkmax)/(5/6);
  plot([taumkmax taumkmax],[phimin,phik],'r-','linewidth',2);
  plot([taumkmax taukktmin-5/6*phimin],[phik,phimin],'r-','linewidth',2);
  text((3*taumkmax+taukktmin-5/6*phimin)/4,(2*phimin+phik)/3,'K','fontsize',24);

  %m
  phim = (taummtmin-taumkmin)/(9/14);
  plot([taumkmin taumkmin],[phimin,phim],'b-','linewidth',2);
  plot([taumin taumkmin],[(taummtmin-taumin)/(9/14),phim],'b-','linewidth',2);
  text((2*taumin+taumkmin)/3,(phimin+phim)/2,'M','fontsize',24);

  %kt 
  taukt = (5/6*tauktmtmax+1/2*taukktmax)/(5/6+1/2);%tau-1/2*phi = tauktmtmax
  phikt = (-tauktmtmax+taukktmax)/(5/6+1/2);%tau+5/6*phi = taukktmax
  plot([taukt taumax],[phikt,(tauktmtmax-taumax)/(-1/2)],'m-','linewidth',2);
  plot([taukt taumax],[phikt,(taukktmax-taumax)/(5/6)],'m-','linewidth',2);
  text((taumax+taukt)/2,((taukktmax-taumax)/(5/6)+(tauktmtmax-taumax)/(-1/2))/2,'~K','fontsize',24);

  %mt
  taumt = (1/2*taummtmax+9/14*tauktmtmin)/(1/2+9/14);%tau+9/14*phi = taummtmax
  phimt = (taummtmax-tauktmtmin)/(9/14+1/2);%tau - 1/2*phi = tauktmtmin
  plot([taummtmax-9/14*phimax,taumt],[phimax,phimt],'g-','linewidth',2);
  plot([taumt,tauktmtmin+1/2*phimax],[phimt,phimax],'g-','linewidth',2);
  text((taummtmax-9/14*phimax+tauktmtmin+1/2*phimax)/2,(phimt+2*phimax)/3,'~M','fontsize',24);

  %location of the selected parameters inside the parametric space
  logtau = log10(tau);
  logphi = log10(phi);
  logtau(logtau<taumin) = taumin;
  logtau(logtau>taumax) = taumax;
  logphi(logphi<phimin) = phimin;
  logphi(logphi>phimax) = phimax;   
  plot(logtau,logphi,'ko','markersize',8,'markerfacecolor','k');

  xlim([taumin taumax]);
  ylim([phimin phimax]);
  xlabel('\tau','fontsize',16);
  ylabel('\phi','fontsize',16);
   
end