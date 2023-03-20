function dx = secretion_bicarb(~,x,par,PrCl,PrKa,PrKb,apicalarea,basalarea)

Nal = x(1);
Kl = x(2);
Cll = x(3);
w = x(4);
Na = x(5);
K = x(6);
Cl = x(7);
HCO3 = x(8);
H = x(9);
Va = x(10);
Vb = x(11);
HCOl = x(12);
Hl = x(13);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NaKbasalfactor = 0.7;                                                    % Fraction of NaK ATPase in the basal membrane
JNaKb = NaKbasalfactor*basalarea * par.aNaK * ( par.r * par.Ke^2 * Na^3 ...
                  / ( par.Ke^2 + par.alpha1 * Na^3 ) );
JNaKa = (1-NaKbasalfactor)*apicalarea * par.aNaK * ( par.r * Kl^2 * Na^3 ...
                  / ( Kl^2 + par.alpha1 * Na^3 ) ); 

%%%%%%%%%%%%%%%%%
VCl = par.RTF * log( Cll / Cl );                                          % Nernst Potentials
VHCO = par.RTF * log( HCOl / HCO3 );
VKb = par.RTF * log( par.Ke / K );                                        % basal K 
VKa = par.RTF * log( Kl / K ) ;                                           % apical K 
VtNa = par.RTF * log( Nal / par.Nae );                                    % tight junction Na
VtK = par.RTF * log( Kl / par.Ke );                                       % tight junction K

Vt = Va - Vb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ca2+ activated channels.

JCl = par.GCl * PrCl * ( Va + VCl ) / par.F;                              % fS.micro-metres^2.mV.mol.C^-1
JKb = par.GK * PrKb * ( Vb - VKb ) / par.F;                               % basal KCa flux
JKa = 0.5*par.GK * PrKa * ( Va - VKa ) / par.F;                           % apical KCa flux. 0.5 is the change from Shan for new model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tight Junction Na+ and K+ currents

JtNa = par.GtNa * par.St * ( Vt - VtNa ) / par.F;                         % fS.micro-metres^2.mV.mol.C^-1
JtK  = par.GtK * par.St * ( Vt - VtK ) / par.F;                           % fS.micro-metres^2.mV.mol.C^-1 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water fluxes (apical, basal and tight junction)

Qa =  par.La * ( 2 * ( Nal + Kl - Na - K - H ) - par.CO20 + par.Ul );     % micro-metres^3.s^-1
Qb =  par.Lb * ( 2 * ( Na + K + H ) + par.CO20 - par.Ie);
Qt =  par.Lt * ( 2 * ( Nal + Kl ) + par.Ul - par.Ie);                     % micro-metres^3.s^-1
Qtot=(Qa+Qt);                                                             % micro-metres^3.s^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ K+ 2Cl- co-transporter (Nkcc1)

JNkcc1 = par.aNkcc1 * basalarea * ( par.a1 - par.a2 * Na * K * Cl^2 ) ...
                                             / ( par.a3 + par.a4 * Na * K * Cl^2 );     
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)

JAe4 = basalarea * par.G4 * ( ( par.Cle / ( par.Cle + par.KCl ) ) * ( Na / ( Na + par.KNa ) ) ...
             * ( HCO3 / ( HCO3 + par.KB ) )^2 );       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ / H+ Anion exchanger (Nhe1)

JNhe1 = basalarea * par.G1 * ( ( par.Nae / ( par.Nae + par.KNa ) ) * ( H / ( par.KH + H ) )...
                          - ( Na / ( Na + par.KNa ) ) * ( par.He / ( par.KH + par.He ) ) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bicarbonate Buffer (Reaction Term) and apical current
% JBB is a reaction inside the cell, note how it depends on the
% cellular volume

JBB = w * par.GB * ( par.kp * par.CO20 - par.kn * HCO3 * H );               % cellular buffering
JBBA = par.wl * par.GB * ( par.kp *0.8* par.CO20 - par.kn * HCOl * Hl );    % apical buffering
JHCO = 0.3 * par.GCl * PrCl * ( Va + VHCO ) / par.F;                        % apical HCO current

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations

% Nal = x(1);
% Kl = x(2);
% Cll = x(3);
% w = x(4);
% Na = x(5);
% K = x(6);
% Cl = x(7);
% HCO3 = x(8);
% H = x(9);
% Va = x(10);
% Vb = x(11);
% HCOl = x(12);
% Hl = x(13);

Jw = Qb - Qa;

dx(1) = ( JtNa - Qtot * Nal + 3*JNaKa )/par.wl;
dx(2) = ( JtK - Qtot * Kl + JKa - 2*JNaKa)/par.wl;
dx(3) = ( - JCl - Qtot * Cll )/par.wl;
dx(4) = Jw;                              
dx(5) = ( JNkcc1 - 3 * (JNaKa+JNaKb) + JNhe1 - JAe4 - dx(4) * Na ) / w;
dx(6) = ( JNkcc1 + 2 * (JNaKa+JNaKb) - JKb - JKa - dx(4) * K ) / w;
dx(7) = ( 2 * JNkcc1 + JAe4 + JCl - dx(4) * Cl ) / w;
dx(8) = ( JBB - 2 * JAe4 + JHCO - dx(4) * HCO3 ) / w;
dx(9) = ( JBB - JNhe1 - dx(4) * H ) / w;
dx(10) = -JCl - JHCO - JNaKa - JKa - JtK - JtNa;                             
dx(11) = -JKb - JNaKb       + JtK + JtNa;
dx(12) = ( JBBA - JHCO - Qtot * HCOl ) / par.wl; 
dx(13) = ( JBBA - Qtot * Hl) / par.wl;

dx = dx';

end