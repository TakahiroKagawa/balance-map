function []=TwoLinkModel_comp(E0,psi,par)
% This function compares the simulation results between the linear and
% nonlinear compass gait models
close all;
if nargin < 1
    E0=0.25;
end
if nargin < 2
    psi=0;
end
if nargin < 3
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
%@Analytical solution of the eigen value
%     par.ChiOmegaI = sqrt(1/(2*coefK)*(DeltaOmega+sqrt(DeltaOmega^2+4*coefK*par.OmegaI^2*par.OmegaS^2)));
%     par.ChiOmegaS = sqrt(1/(2*coefK)*(-DeltaOmega+sqrt(DeltaOmega^2+4*coefK*par.OmegaI^2*par.OmegaS^2)));
    par.ChiOmega0 = par.ChiOmegaI/par.ChiOmegaS;
    par.ES = 0.6; %Orbital Energy of swing pendulum
end

global gPar;
gPar = par;
par.EI = par.ES*E0; %Orbital Energy of inverted pendulum
par.E0 = par.EI/par.ES;    % Energy ratio
par.psi = psi; % phase of swing angle related to stance angle
T0=0;
chi0=DefInitChi(T0,par.E0,par.psi,par.ChiOmega0);

duration=0.8;
sampling=0.01;
%iteration=duration/sampling;

%opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
x0=chi2x(chi0,par);
y0=x2y(x0,par);

% for simulation with nonlinear 2link dynamics, use '@f_nlin'. for linear
% model, use '@f_lin'

opts = odeset('RelTol',1e-12,'AbsTol',1e-16,'Refine',30,'Events',@IntersectionPoint_backward);
h=0:-sampling:-duration;
[t_n,y_n] =  ode113(@(t,y)f_lin(t,y,par),h,y0,opts); % Integrate for one stride

opts = odeset('RelTol',1e-12,'AbsTol',1e-16,'Refine',30,'Events',@IntersectionPoint_forward);
h=0:sampling:duration;
[t_p,y_p] = ode113(@(t,y)f_lin(t,y,par),h,y0,opts); % Integrate for one stride
t = [flipud(t_n);t_p(2:end)];
y = [flipud(y_n);y_p(2:end,:)];
%Normal trajectory
opts = odeset('RelTol',1e-12,'AbsTol',1e-16,'Refine',30,'Events',@IntersectionPoint_forward);
h=0:sampling:2*duration;
y_init=y_n(end-1,:);

[t_lin,y_lin] = ode113(@(t,y)f_lin(t,y,par),h,y_init,opts); % Integrate for one stride
x_lin=y2x(y_lin,par);
chi_lin = x2chi(x_lin,par);
[t_nlin,y_nlin] = ode113(@(t,y)f_nlin(t,y,par),h,y_init,opts); % Integrate for one stride
xn_nlin=y2x_nlin(y_nlin,par);
chi_nlin = x2chi(y2x(y_nlin,par),par);

x_nlin=y2x(y_nlin,par);
chi_nlin = x2chi(x_nlin,par);


% u=[-100,0;0,0];ParturbationOnset=1*sampling; ParturbationDuration=20*sampling; % backward falling
  u=[0,0;-12,0];ParturbationOnset=10*sampling; ParturbationDuration=10*sampling;% 0.01 interval




%Disturbed trajectory
duration = 0.6;
h=0:sampling:ParturbationOnset;
[t1,y1] = ode113(@(t,y)f_lin(t,y,par),h,y_init,opts); % Integrate for swing phase

y0=y1(end,:);
h=ParturbationOnset:sampling:ParturbationOnset+ParturbationDuration;
[t2,y2] = ode113(@(t,y)f_lin_dist(t,y,par,u),h,y0,opts); % Integrate for swing phase
y0=y2(end,:);
h=t2(end):sampling:duration;
[t3,y3] = ode113(@(t,y)f_lin(t,y,par),h,y0,opts);
yp_lin=[y1;y2(2:end,:);y3(2:end,:)];
xp_lin=y2x(yp_lin,par);
chip_lin = x2chi(xp_lin,par);
tp_lin=[t1;t2(2:end);t3(2:end)];


%Disturbed trajectory
duration = 0.6;
h=0:sampling:ParturbationOnset;
[t1,y1] = ode113(@(t,y)f_nlin(t,y,par),h,y_init,opts); % Integrate for swing phase
if length(h) ==2
    t1 = h';
    y1 =[y1(1,:);y1(end,:)];
end
y0=y1(end,:);
h=ParturbationOnset:sampling:ParturbationOnset+ParturbationDuration;
[t2,y2] = ode113(@(t,y)f_nlin_dist(t,y,par,u),h,y0,opts); % Integrate for swing phase
y0=y2(end,:);
h=t2(end):sampling:duration;
[t3,y3] = ode113(@(t,y)f_nlin(t,y,par),h,y0,opts);
yp_nlin=[y1;y2(2:end,:);y3(2:end,:)];
xp_nlin=y2x_nlin(yp_nlin,par);
chip_nlin = x2chi(y2x(yp_nlin,par),par);
tp_nlin=[t1;t2(2:end);t3(2:end)];

% normal movement
figure(1);
plot(t_lin,x_lin(:,1),'-r');hold on;
plot(t_lin,x_lin(:,2),'--r');
plot(t_nlin,x_nlin(:,1),'-b');
plot(t_nlin,x_nlin(:,2),'--b');
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Position (m)');
grid on
figure(31);
plot(t_lin,chi_lin(:,1),'-r');hold on;
plot(t_lin,chi_lin(:,2),'--r');
plot(t_nlin,chi_nlin(:,1),'-b');
plot(t_nlin,chi_nlin(:,2),'--b');
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Position (m)');
grid on


BalanceMap(2);
DataSize=length(t_lin);
chi_lin = x2chi(x_lin,par);
EquivInitialChi_lin = VirtualInitialState(chi_lin,par);
plot(par.ChiOmega0*EquivInitialChi_lin(1,1),EquivInitialChi_lin(1,2),'o','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChi_lin(DataSize,1),EquivInitialChi_lin(DataSize,2),'s','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChi_lin(:,1),EquivInitialChi_lin(:,2),'-r');

DataSize=length(t_nlin);
chi_nlin = x2chi(x_nlin,par);
EquivInitialChi_nlin = VirtualInitialState(chi_nlin,par);
plot(par.ChiOmega0*EquivInitialChi_nlin(1,1),EquivInitialChi_nlin(1,2),'ob');
plot(par.ChiOmega0*EquivInitialChi_nlin(DataSize,1),EquivInitialChi_nlin(DataSize,2),'sb');
plot(par.ChiOmega0*EquivInitialChi_nlin(:,1),EquivInitialChi_nlin(:,2),'-b');

StickPicture(y_nlin,par,5,11)
% perturbed movement

figure(3);
plot(tp_lin,xp_lin(:,1),'-r');hold on;
plot(tp_lin,xp_lin(:,2),'--r');
plot(tp_nlin,xp_nlin(:,1),'-b');
plot(tp_nlin,xp_nlin(:,2),'--b');
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Position (m)');
grid on

figure(33);
plot(tp_lin,chip_lin(:,1),'-r');hold on;
plot(tp_lin,chip_lin(:,2),'--r');
plot(tp_nlin,chip_nlin(:,1),'-b');
plot(tp_nlin,chip_nlin(:,2),'--b');
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Position (m)');
grid on

BalanceMap(4);
DataSize=length(tp_lin);
chip_lin = x2chi(xp_lin,par);
EquivInitialChip_lin = VirtualInitialState(chip_lin,par);
plot(par.ChiOmega0*EquivInitialChip_lin(1,1),EquivInitialChip_lin(1,2),'o','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChip_lin(DataSize,1),EquivInitialChip_lin(DataSize,2),'s','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChip_lin(:,1),EquivInitialChip_lin(:,2),'-r');

DataSize=length(tp_nlin);
chip_nlin = x2chi(xp_nlin,par);
EquivInitialChip_nlin = VirtualInitialState(chip_nlin,par);
plot(par.ChiOmega0*EquivInitialChip_nlin(1,1),EquivInitialChip_nlin(1,2),'ob');
plot(par.ChiOmega0*EquivInitialChip_nlin(DataSize,1),EquivInitialChip_nlin(DataSize,2),'sb');
plot(par.ChiOmega0*EquivInitialChip_nlin(:,1),EquivInitialChip_nlin(:,2),'-b');

figure(5);
[PosMargin_lin, NegMargin_lin]=MarginFunction(chip_lin,par);
plot(tp_lin,NegMargin_lin,'-r');hold on
[PosMargin_nlin, NegMargin_nlin]=MarginFunction(chip_nlin,par);
plot(tp_nlin,NegMargin_nlin,'-b');
plot(tp_lin,zeros(length(tp_lin),1),'-k','Linewidth',2);
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Margin for backward balance loss (m/s)');
grid on;

figure(6);
plot(tp_lin,PosMargin_lin,'-r');hold on
plot(tp_nlin,PosMargin_nlin,'-b');
plot(tp_lin,zeros(length(tp_lin),1),'-k','Linewidth',2);
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Margin for forward balance loss (m/s)');
grid on;
StickPicture(yp_nlin,par,5,12);




function [PosMargin, NegMargin]=MarginFunction(chi,par)
dummy = VirtualInitialState(chi,par);
EquivInitialChi_p(:,1)=par.ChiOmega0*dummy(:,1);
EquivInitialChi_p(:,2)=dummy(:,2);
load('balancemap_limit.mat');
[ChiLength,~]=size(chi);
[BorderLength]=length(y_lim);
PosMargin=100000*ones(ChiLength,1);
NegMargin=100000*ones(ChiLength,1);
for cnt1 = 1:BorderLength
    NormLim(cnt1)=norm([par.ChiOmega0*y_lim(cnt1),yd_lim(cnt1)]);
end

for cnt1 = 1:ChiLength
    normChi=norm(EquivInitialChi_p(cnt1,:));
    % positive border
    for cnt2 = 1:BorderLength
        Distance(cnt2)=norm(EquivInitialChi_p(cnt1,:)-[par.ChiOmega0*y_lim(cnt2),yd_lim(cnt2)]);
        
        if Distance(cnt2)<PosMargin(cnt1);
            if normChi < NormLim(cnt2)
                PosMargin(cnt1)=Distance(cnt2);
            else
                PosMargin(cnt1)=-Distance(cnt2);                
            end
        end
    end
    % negative border
    EI=0.5*EquivInitialChi_p(cnt1,2)^2-0.5*EquivInitialChi_p(cnt1,1)^2;
    if EI<0
        NegMargin(cnt1)=-abs(EquivInitialChi_p(cnt1,1)+EquivInitialChi_p(cnt1,2))/sqrt(2);
    else
        NegMargin(cnt1)=abs(EquivInitialChi_p(cnt1,1)+EquivInitialChi_p(cnt1,2))/sqrt(2);
    end
end
function T0state =  VirtualInitialState(chi,par)
[DataNum,~]=size(chi);
T0state=zeros(DataNum,2);
for cnt1 =1 : DataNum
    T=-atan2(chi(cnt1,2),chi(cnt1,4));
    InvMatB(1,1)=cosh(par.ChiOmega0*T);
    InvMatB(1,2)=1/par.ChiOmega0*sinh(par.ChiOmega0*T);
    InvMatB(2,1)=par.ChiOmega0*sinh(par.ChiOmega0*T);
    InvMatB(2,2)=cosh(par.ChiOmega0*T);
    T0state(cnt1,:)=(InvMatB* chi(cnt1,[1,3])')';
end


function [Et,Ek,Ep]=MechEnergy(y,par)
% Mechanical energy of double pendulum

Ep=(par.MI+par.MS)*par.g*par.LI*cos(y(:,1))-...
    par.MS*par.g*par.LS*cos(y(:,2));

Ek=1/2*(par.MI+par.MS)*par.LI^2*y(:,3).^2+...
   1/2*par.MS*(par.LS^2*y(:,4).^2 -2*par.LI*par.LS*y(:,3).*y(:,4).*cos(y(:,2)-y(:,1)));
Et=Ep+Ek;


function x = y2x(y,par)
% state variables to position
x(:,1)=-par.LI*y(:,1);
x(:,2)= par.LI*y(:,2);
x(:,3)=-par.LI*y(:,3);
x(:,4)= par.LI*y(:,4);

function x = y2x_nlin(y,par)
% state variables to position
x(:,1)=-par.LI*sin(y(:,1));
x(:,2)= par.LI*sin(y(:,2));
x(:,3)=-par.LI*y(:,3).*cos(y(:,1));
x(:,4)= par.LI*y(:,4).*cos(y(:,2));

function y = x2y(x,par)
% state variables to position
y(:,1)=-1/par.LI*x(:,1);
y(:,2)= 1/par.LI*x(:,2);
y(:,3)=-1/par.LI*x(:,3);
y(:,4)= 1/par.LI*x(:,4);

function x = prime2x(x_prime,par)
% 
% Sigma1 = par.ChiOmegaI^2;
% Sigma2 = -par.ChiOmegaS^2;
% % Alpha = 1/sqrt(par.OmegaI^4+(par.OmegaS^2+coefK*Sigma2)^2);
% % Beta = 1/sqrt(par.OmegaI^4+(par.OmegaS^2+coefK*Sigma1)^2);
% Alpha = 1/(par.OmegaI^2+(par.OmegaS^2+coefK*Sigma2)*par.LS/par.LI);
% Beta = 1/(par.OmegaI^2+(par.OmegaS^2+coefK*Sigma1)*par.LS/par.LI);
% P=[Alpha*par.OmegaI^2, Alpha*(par.OmegaS^2+coefK*Sigma2)*par.LS/par.LI;
% Beta*par.OmegaI^2, Beta*(par.OmegaS^2+coefK*Sigma1)*par.LS/par.LI];
% Pinv=inv(P);
% Normaized position to state variables

x(:,1)=par.Tinv(1,1)*x_prime(:,1)+par.Tinv(1,2)*x_prime(:,2);
x(:,2)=par.Tinv(2,1)*x_prime(:,1)+par.Tinv(2,2)*x_prime(:,2);
x(:,3)=par.Tinv(1,1)*x_prime(:,3)+par.Tinv(1,2)*x_prime(:,4);
x(:,4)=par.Tinv(2,1)*x_prime(:,3)+par.Tinv(2,2)*x_prime(:,4);

function x_prime = x2prime(x,par)
% coefK = par.MI/(par.MI+par.MS);
% Sigma1 = par.ChiOmegaI^2;
% Sigma2 = -par.ChiOmegaS^2;
% % Alpha = 1/sqrt(par.OmegaI^4+(par.OmegaS^2+coefK*Sigma2)^2);
% % Beta = 1/sqrt(par.OmegaI^4+(par.OmegaS^2+coefK*Sigma1)^2);
% Alpha = 1/(par.OmegaI^2+(par.OmegaS^2+coefK*Sigma2)*par.LS/par.LI);
% Beta = 1/(par.OmegaI^2+(par.OmegaS^2+coefK*Sigma1)*par.LS/par.LI);
% P=[Alpha*par.OmegaI^2, Alpha*(par.OmegaS^2+coefK*Sigma2)*par.LS/par.LI;
% Beta*par.OmegaI^2, Beta*(par.OmegaS^2+coefK*Sigma1)*par.LS/par.LI];

x_prime(:,1)=par.T(1,1)*x(:,1)+par.T(1,2)*x(:,2);
x_prime(:,2)=par.T(2,1)*x(:,1)+par.T(2,2)*x(:,2);
x_prime(:,3)=par.T(1,1)*x(:,3)+par.T(1,2)*x(:,4);
x_prime(:,4)=par.T(2,1)*x(:,3)+par.T(2,2)*x(:,4);

function x = chi2x(chi,par,ES)
if nargin < 3
    ES = par.ES;
end

% Normaized position to state variables
x_prime(:,1)=sqrt(2*ES)/par.ChiOmegaS*chi(:,1);
x_prime(:,2)=sqrt(2*ES)/par.ChiOmegaS*chi(:,2);
x_prime(:,3)=sqrt(2*ES)*chi(:,3);
x_prime(:,4)=sqrt(2*ES)*chi(:,4);
x = prime2x(x_prime,par);

function chi = x2chi(x,par)
x_prime = x2prime(x,par);
[DataSize,~]=size(x_prime);
ES = zeros(DataSize,1);
chi = zeros(DataSize,4);
for cnt1 = 1:DataSize
    ES(cnt1) = 1/2*(x_prime(cnt1,4)^2)+1/2*(par.ChiOmegaS*x_prime(cnt1,2))^2;
    
    % State variables to Normalized position
    chi(cnt1,1)=par.ChiOmegaS/sqrt(2*ES(cnt1))*x_prime(cnt1,1);
    chi(cnt1,2)=par.ChiOmegaS/sqrt(2*ES(cnt1))*x_prime(cnt1,2);
    chi(cnt1,3)=1/sqrt(2*ES(cnt1))*x_prime(cnt1,3);
    chi(cnt1,4)=1/sqrt(2*ES(cnt1))*x_prime(cnt1,4);
end

function ydot=f_nlin(t,y,par)    %#ok<INUSL>
% ODE definition
% y1: theta1
% y2: theta2
% y3: theta1dot
% y4: theta2dot

% First order differential equations for Simplest Walking Model
M=zeros(2,2);
V=zeros(2,1);
G=zeros(2,1);
M(1,1)=(par.MI+par.MS)*par.LI^2;
M(1,2)=-par.MS*par.LI*par.LS*cos(y(2)-y(1));
M(2,1)=M(1,2);
M(2,2)=par.MS*par.LS^2;
V(1)= par.MS*par.LI*par.LS*y(4)^2*sin(y(2)-y(1));
V(2)=-par.MS*par.LI*par.LS*y(3)^2*sin(y(2)-y(1));
G(1)=-(par.MI+par.MS)*par.g*par.LI*sin(y(1));
G(2)=par.MS*par.g*par.LS*sin(y(2));

yddot=M\(-V-G);
ydot=zeros(4,1);
ydot(1) = y(3);
ydot(2) = y(4);
ydot(3:4)=yddot;

function ydot=f_nlin_dist(t,y,par,u)    %#ok<INUSL>
% ODE definition
% y1: theta1
% y2: theta2
% y3: theta1dot
% y4: theta2dot

% First order differential equations for Simplest Walking Model
M=zeros(2,2);
V=zeros(2,1);
G=zeros(2,1);
M(1,1)=(par.MI+par.MS)*par.LI^2;
M(1,2)=-par.MS*par.LI*par.LS*cos(y(2)-y(1));
M(2,1)=M(1,2);
M(2,2)=par.MS*par.LS^2;
V(1)= par.MS*par.LI*par.LS*y(4)^2*sin(y(2)-y(1));
V(2)=-par.MS*par.LI*par.LS*y(3)^2*sin(y(2)-y(1));
G(1)=-(par.MI+par.MS)*par.g*par.LI*sin(y(1));
G(2)=par.MS*par.g*par.LS*sin(y(2));

Jh = [-par.LI*cos(y(1)),0;
    -par.LI*sin(y(1)),0];
Jt = [-par.LI*cos(y(1)),par.LI*cos(y(2));
    -par.LI*sin(y(1)),par.LI*sin(y(2))];
f=Jh'*u(1,:)'+Jt'*u(2,:)';

yddot=M\(-G-V+f);

ydot=zeros(4,1);
ydot(1) = y(3);
ydot(2) = y(4);
ydot(3:4)=yddot;

function ydot=f_lin(t,y,par)    %#ok<INUSL>
% ODE definition
% y1: theta1
% y2: theta2
% y3: theta1dot
% y4: theta2dot

% First order differential equations for Simplest Walking Model
M=zeros(2,2);
G=zeros(2,1);
M(1,1)=(par.MI+par.MS)*par.LI^2;
% M(1,2)=0;
M(1,2)=-par.MS*par.LI*par.LS;
M(2,1)=-par.MS*par.LI*par.LS;
M(2,2)=par.MS*par.LS^2;
G(1)=-(par.MI+par.MS)*par.g*par.LI*y(1);
G(2)=par.MS*par.g*par.LS*y(2);


yddot=M\(-G);
ydot=zeros(4,1);
ydot(1) = y(3);
ydot(2) = y(4);
ydot(3:4)=yddot;


function ydot=f_lin_dist(t,y,par,u)    %#ok<INUSL>
% ODE definition
% y1: theta1
% y2: theta2
% y3: theta1dot
% y4: theta2dot

% First order differential equations for Simplest Walking Model
M=zeros(2,2);
G=zeros(2,1);
M(1,1)=(par.MI+par.MS)*par.LI^2;
% M(1,2)=0;
M(1,2)=-par.MS*par.LI*par.LS;
M(2,1)=-par.MS*par.LI*par.LS;
M(2,2)=par.MS*par.LS^2;
G(1)=-(par.MI+par.MS)*par.g*par.LI*y(1);
G(2)=par.MS*par.g*par.LS*y(2);
Jh = [-par.LI*cos(y(1)),0;
    -par.LI*sin(y(1)),0];
Jt = [-par.LI*cos(y(1)),par.LI*cos(y(2));
    -par.LI*sin(y(1)),par.LI*sin(y(2))];
f=Jh'*u(1,:)'+Jt'*u(2,:)';

yddot=M\(-G+f);
ydot=zeros(4,1);
ydot(1) = y(3);
ydot(2) = y(4);
ydot(3:4)=yddot;

function StickPicture(y,par,interval,fignum)
figure(fignum);clf;
axis([-0.5 0.5 0 1]);
axis square;
grid on
hold on
Hip=[-par.LI*sin(y(:,1)),par.LI*cos(y(:,1))];
Foot=Hip+[par.LI*sin(y(:,2)),-par.LI*cos(y(:,2))];
for cnt1=1:interval:length(Hip)
    x=[0,Hip(cnt1,1),Foot(cnt1,1)];
    y=[0,Hip(cnt1,2),Foot(cnt1,2)];
    plot(x,y,'-k');
end
set(gca,'FontSize',16);
set(gca,'Ytick',[0 0.5 1]);

function GraphicPendula(y,par,fignum)

figure(fignum);clf;
axis([-1 1 -1 1]);
axis square;
grid on

Hip=[-par.LI*sin(y(:,1)),par.LI*cos(y(:,1))];
Foot=Hip+[par.LI*sin(y(:,2)),-par.LI*cos(y(:,2))];
StickPic=line('XData',[0,Hip(1,1),Foot(1,1)],'YData',[0,Hip(1,2),Foot(1,2)],'EraseMode','xor');
for cnt1=1:1:length(Hip)
    x=[0,Hip(cnt1,1),Foot(cnt1,1)];
    y=[0,Hip(cnt1,2),Foot(cnt1,2)];
    set(StickPic,'XData',x,'YData',y);
    drawnow;
    pause(0.03);
end
set(gca,'FontSize',16);
set(gca,'Ytick',[0 0.5 1]);



function [val,ist,dir]=IntersectionPoint_backward(t,y)   %#ok<INUSL>
% Check for heelstrike collision using zero-crossing detection
global gPar;
par  = gPar;
x=y2x(y',par);
chi=x2chi(x,par);
% val = chi(2)-chi(1);  % Geometric collision condition, when = 0
% val = y(1)+y(2);  % Geometric collision condition, when = 0
val = x(2)-x(1);  % Geometric collision condition, when = 0
ist = 1;            % Stop integrating if collision found
dir = 0;            % Condition only true when passing from - to +

function [val,ist,dir]=IntersectionPoint_forward(t,y)   %#ok<INUSL>
% Check for heelstrike collision using zero-crossing detection
global gPar;
par  = gPar;
x=y2x(y',par);
chi=x2chi(x,par);
val = x(2)-x(1);  % Geometric collision condition, when = 0
% val = y(1)+y(2);  % Geometric collision condition, when = 0
ist = 1;            % Stop integrating if collision found
dir = -1;            % Condition only true when passing from - to +