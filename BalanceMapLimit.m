function [y_lim, yd_lim,par]=BalanceMapLimit(par)
%Å@This function calculate the boundary of stable region on the balance map

if nargin < 1
    par.g=9.8;  % gravity
    par.LI=0.8; % leg length
    par.LS=0.3; % swing COM length
    par.MI=60-12;   % total mass of stance leg and body
    par.MS=12;      % mass of swing leg
    par.OmegaI = sqrt(par.g/par.LI); % natural frequency of the inverted pendulum
    par.OmegaS = sqrt(par.g/par.LS); % natural frequency of the simple pendulum
    par.Omega0 = par.OmegaI/par.OmegaS; % natural frequency ratio
    coefK = par.MI/(par.MI+par.MS); % parameter K in ICRA paper
    par.A=[par.OmegaI^2/coefK, (1-coefK)/coefK*par.OmegaI^2;
        -par.OmegaS^2/coefK, -par.OmegaS^2/coefK];
    [V,S]=eig(par.A);
    dummy = inv(V);
    par.T(1,1:2)=dummy(1,1:2)/(dummy(1,1)+dummy(1,2));
    par.T(2,1:2)=dummy(2,1:2)/(dummy(2,1)+dummy(2,2));
    par.Tinv=inv(par.T);    % diagonalization matrix
    par.ChiOmegaI = sqrt(S(1,1)); 
    par.ChiOmegaS = sqrt(-S(2,2));
%Å@Analytical solution of the eigen value
%     par.ChiOmegaI = sqrt(1/(2*coefK)*(DeltaOmega+sqrt(DeltaOmega^2+4*coefK*par.OmegaI^2*par.OmegaS^2)));
%     par.ChiOmegaS = sqrt(1/(2*coefK)*(-DeltaOmega+sqrt(DeltaOmega^2+4*coefK*par.OmegaI^2*par.OmegaS^2)));
    par.ChiOmega0 = par.ChiOmegaI/par.ChiOmegaS;
    %X1_dot0 = 0.5; % minimum veloicty of inverted pendulum
    par.ES = 0.6; %Orbital Energy of swing pendulum
end

Omega0 = par.ChiOmega0;

E0_Positive=[5:-0.01:0.01];
E0_Negative=[-0.001:-0.01:-0.75];

% Calculation of stability limit of the phase difference by using numerical calculation with
% non-linear dynamics model of a compass model
Delay_Positive=StabilityBoundaryPositive_Model(E0_Positive,Omega0,'nlin',par); % positive enregy ratio
Delay_Negative=StabilityBoundaryNegative_Model(E0_Negative,Omega0,'nlin',par); % negative energy ratio
y_lim_p=length(E0_Positive);
yd_lim_p=length(E0_Positive);
for cnt1=1:length(E0_Positive)
    y_lim_p(cnt1)=sqrt(E0_Positive(cnt1))/Omega0*sinh(Delay_Positive(cnt1));
    yd_lim_p(cnt1)=sqrt(E0_Positive(cnt1))*cosh(Delay_Positive(cnt1));
end
y_lim_n=length(E0_Negative);
yd_lim_n=length(E0_Negative);
for cnt1=1:length(E0_Negative)
    y_lim_n(cnt1)=sqrt(-E0_Negative(cnt1))/Omega0*cosh(Delay_Negative(cnt1));
    yd_lim_n(cnt1)=sqrt(-E0_Negative(cnt1))*sinh(Delay_Negative(cnt1));
end
y_lim=[y_lim_p,y_lim_n];
yd_lim=[yd_lim_p,yd_lim_n];

save('balancemap_limit.mat','y_lim', 'yd_lim','par');
