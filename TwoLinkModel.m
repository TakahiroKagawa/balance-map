function [chi,t,par]=TwoLinkModel(E0,psi,mode,par)
if nargin < 1
    E0=0.5;
end
if nargin < 2
    psi=0;
end
if nargin < 3
    mode = 'lin';
end
if nargin < 4
    par.g=9.8;
    par.LI=0.8;
    par.LS=0.3;
    par.MI=60-12;
    par.MS=12;
    par.OmegaI = sqrt(par.g/par.LI); % inverted pendulum
    par.OmegaS = sqrt(par.g/par.LS); % simple pendulum
    par.Omega0 = par.OmegaI/par.OmegaS;
    coefK = par.MI/(par.MI+par.MS);
    par.A=[par.OmegaI^2/coefK, (1-coefK)/coefK*par.OmegaI^2;
        -par.OmegaS^2/coefK, -par.OmegaS^2/coefK];
    [V,S]=eig(par.A);
    dummy = inv(V);
    par.T(1,1:2)=dummy(1,1:2)/(dummy(1,1)+dummy(1,2));
    par.T(2,1:2)=dummy(2,1:2)/(dummy(2,1)+dummy(2,2));
    par.Tinv=inv(par.T);
    par.ChiOmegaI = sqrt(S(1,1));
    par.ChiOmegaS = sqrt(-S(2,2));

%     par.ChiOmegaI = sqrt(1/(2*coefK)*(DeltaOmega+sqrt(DeltaOmega^2+4*coefK*par.OmegaI^2*par.OmegaS^2)));
%     par.ChiOmegaS = sqrt(1/(2*coefK)*(-DeltaOmega+sqrt(DeltaOmega^2+4*coefK*par.OmegaI^2*par.OmegaS^2)));
    par.ChiOmega0 = par.ChiOmegaI/par.ChiOmegaS;
    %X1_dot0 = 0.5; % minimum veloicty of inverted pendulum
    ES0 = 0.6;
    par.ES = ES0; %Orbital Energy of swing pendulum
end
par.EI = par.ES*E0; %Orbital Energy of inverted pendulum
par.E0 = par.EI/par.ES;    % Energy ratio
par.psi = psi; % phase of swing angle related to stance angle

T0=0;
chi0=DefInitChi(T0,par.E0,par.psi,par.ChiOmega0);

duration=1.5;
sampling=0.01;
%iteration=duration/sampling;
opts = odeset('RelTol',1e-12,'AbsTol',1e-16,'Refine',30,'Events',@IntersectionPoint);
%opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
x0=chi2x(chi0,par);
y0=x2y(x0,par);


% for simulation with nonlinear 2link dynamics, use '@f_nlin'. for linear
% model, use '@f_lin'
if strcmp(mode,'lin')
    h=0:sampling:duration;
    [t,y] = ode113(@f_lin,h,y0,opts,par); % Integrate for one stride    
end

if strcmp(mode,'nlin')
    h=0:sampling:duration;
    [t,y] = ode113(@f_nlin,h,y0,opts,par); % Integrate for one stride    
end
x=y2x(y,par);
chi=x2chi(x,par);
% figure(101);
% plot(t,[-y(:,1),y(:,2)]);
% figure(102);
% plot(t,[-y(:,3),y(:,4)]);
% figure(1);
% BalanceMapPendula(t,y,par,20)
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

function GraphicPendula(y,par,fignum)
close all;
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
function BalanceMapPendula(t,y,par,fignum)
figure(fignum);clf;
[DataSize,~]=size(y);
x = y2x(y,par);
chi = x2chi(x,par);
EquivInitialChi = VirtualInitialState(chi,par);

fig = openfig('balancemap.fig','new');hold on;ax1=gca;
figure(fignum);set(gcf,'position', [50   500   960   420]);
s1=subplot(1,2,1);
figobj1=get(ax1,'children');
copyobj(figobj1,s1);axis([-1 1 0 1.5]);
close(fig)

balancemapPoint=line('XData',par.ChiOmega0*EquivInitialChi(1,1),'YData',EquivInitialChi(1,2),'EraseMode','xor','Marker','o','MarkerEdgeColor','r','MarkerSize',10);
subplot(1,2,2);
plot([-10;10], [0,0]);
axis([-par.LI par.LI -par.LI par.LI]);
axis equal;
grid on;
Hip=[-par.LI*sin(y(:,1)),par.LI*cos(y(:,1))];
Foot=Hip+[par.LI*sin(y(:,2)),-par.LI*cos(y(:,2))];
StickPic=line('XData',[0,Hip(1,1),Foot(1,1)],'YData',[0,Hip(1,2),Foot(1,2)],'EraseMode','xor');
for cnt1=2:1:min(100,DataSize)
    subplot(1,2,1);
    set(balancemapPoint,'XData',par.ChiOmega0*EquivInitialChi(cnt1,1),'YData',EquivInitialChi(cnt1,2));
    drawnow;
    subplot(1,2,2);    
    x=[0,Hip(cnt1,1),Foot(cnt1,1)];
    y=[0,Hip(cnt1,2),Foot(cnt1,2)];
    set(StickPic,'XData',x,'YData',y);
    drawnow;
    pause(0.01);
end

function DrawGraph(y,ydef,par,t,tdef,fignum)
figure(fignum);clf;
[DataSize,~]=size(y);
x = y2x(y,par);
chi = x2chi(x,par);
EquivInitialChi = VirtualInitialState(chi,par);
xdef = y2x(ydef,par);
chidef = x2chi(xdef,par);
EquivInitialChidef = VirtualInitialState(chidef,par);
subplot(1,2,1);
% tdef=0:sampling:sampling*(length(ydef)-1);
plot(tdef,[-par.LI*ydef(:,1),par.LI*ydef(:,2)],'--');hold on;
% tp=0:sampling:sampling*(length(y)-1);
plot(t,[-par.LI*y(:,1),par.LI*y(:,2)]);
grid on;
xlabel('Time (s)'); ylabel('Position (m)');
set(gca,'YLim',[-0.4 0.4]);set(gca,'XLim',[0 0.65]);


fig = openfig('balancemap.fig','new');hold on;ax1=gca;
figure(fignum);set(gcf,'position', [800   500   960   420]);
s1=subplot(1,2,2);
figobj1=get(ax1,'children');
copyobj(figobj1,s1);axis([-1 1 0 1.5]);
close(fig)
hold on;
plot(par.ChiOmega0*EquivInitialChi(1,1),EquivInitialChi(1,2),'o','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChi(DataSize,1),EquivInitialChi(DataSize,2),'s','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChi(:,1),EquivInitialChi(:,2),'-r');
plot(par.ChiOmega0*EquivInitialChidef(1,1),EquivInitialChidef(1,2),'ob');
plot(par.ChiOmega0*EquivInitialChidef(end,1),EquivInitialChidef(end,2),'sb');
plot(par.ChiOmega0*EquivInitialChidef(:,1),EquivInitialChidef(:,2),'-b');



function [val,ist,dir]=IntersectionPoint(t,y,par)   %#ok<INUSL>
% Check for heelstrike collision using zero-crossing detection
%x=y2x(y',par);
x=y2x(y',par);
chi=x2chi(x,par);
% val = chi(2)-chi(1);  % Geometric collision condition, when = 0

val = x(2)-x(1);  % Geometric collision condition, when = 0
ist = 1;            % Stop integrating if collision found
dir = -1;            % Condition only true when passing from - to +ng from - to +