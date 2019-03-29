function []=TwoLinkModel_sim(E0,psi,mode,par)
% This function simulate walking with compass gait model
% nonlinear compass gait models
close all;
if nargin < 1
    E0=0.25;
end
if nargin < 2
    psi=0;
end
if nargin < 3
    mode = 'lin';
end
if nargin < 4
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
chi0=DefInitChi(T0,par.E0,par.psi,par.ChiOmega0); % initial condition

duration=1.0;
sampling=1/300;



x0=chi2x(chi0,par);
y0=x2y(x0,par);

% for simulation with nonlinear 2link dynamics, use '@f_nlin'. for linear
% model, use '@f_lin'
if strcmp(mode,'lin')
    opts = odeset('RelTol',1e-12,'AbsTol',1e-16,'Refine',30,'Events',@IntersectionPoint_backward);
    % Backward integration to calculate the state at toe-off
    h=0:-sampling:-duration;
    [t_n,y_n] =  ode113(@(t,y)f_lin(t,y,par),h,y0,opts); % Integrate for one stride
    
    %Normal trajectory
    opts = odeset('RelTol',1e-12,'AbsTol',1e-16,'Refine',30,'Events',@IntersectionPoint_forward);
    h=0:sampling:2*duration;
    y0=y_n(end-1,:);
    [tn,yn] = ode113(@(t,y)f_lin(t,y,par),h,y0,opts); % Integrate for one stride
    xn=y2x(yn,par);
    chin=x2chi(xn,par);

end
if strcmp(mode,'nlin')
    opts = odeset('RelTol',1e-8,'AbsTol',1e-12,'Refine',30,'Events',@IntersectionPoint_backward);
    % Backward integration to calculate the state at toe-off
    h=0:-sampling:-duration;
    [t_n,y_n] = ode113(@(t,y)f_nlin(t,y,par),h,y0,opts); % Integrate for one stride

    %Normal trajectory
    opts = odeset('RelTol',1e-8,'AbsTol',1e-12,'Refine',30,'Events',@IntersectionPoint_forward);
    h=0:sampling:2*duration;
    y0=y_n(end-1,:);
    [tn,yn] = ode113(@(t,y)f_nlin(t,y,par),h,y0,opts); % Integrate for one stride
    xn=y2x_nlin(yn,par);
    chin=x2chi(y2x(yn,par),par);
end

% disturbance
u=[-100,0;0,0];ParturbationOnset=sampling*round(0.033333/sampling); ParturbationDuration=sampling*round(0.2/sampling); % backward falling
%u=[0,0;-12,0];ParturbationOnset=sampling*round(0.1/sampling); ParturbationDuration=sampling*round(0.1/sampling);% 0.01 interval


if strcmp(mode,'lin')
    %Disturbed trajectory
    h=0:sampling:ParturbationOnset;
    [t1,y1] = ode113(@(t,y)f_lin(t,y,par),h,y0,opts); % Integrate for swing phase
   if length(h)==2
       t1=[t1(1);t1(end)];
       y1=[y1(1,:);y1(end,:)];
   end
    y0=y1(end,:);
    h=ParturbationOnset:sampling:ParturbationOnset+ParturbationDuration;
    [t2,y2] = ode113(@(t,y)f_lin_dist(t,y,par,u),h,y0,opts); % Integrate for swing phase
    y0=y2(end,:);
    h=t2(end):sampling:duration;
    [t3,y3] = ode113(@(t,y)f_lin(t,y,par),h,y0,opts);
    yp=[y1;y2(2:end,:);y3(2:end,:)];
    tp=[t1;t2(2:end);t3(2:end)];
    xp=y2x(yp,par);
    chip=x2chi(xp,par);

end

if strcmp(mode,'nlin')
    %Disturbed trajectory
    h=0:sampling:ParturbationOnset;
    [t1,y1] = ode113(@(t,y)f_nlin(t,y,par),h,y0,opts); % Integrate for swing phase
    
    y0=y1(end,:);
    h=ParturbationOnset:sampling:ParturbationOnset+ParturbationDuration;
    [t2,y2] = ode113(@(t,y)f_nlin_dist(t,y,par,u),h,y0,opts); % Integrate for swing phase
    y0=y2(end,:);
    h=t2(end):sampling:duration;
    [t3,y3] = ode113(@(t,y)f_nlin(t,y,par),h,y0,opts);
    yp=[y1;y2(2:end,:);y3(2:end,:)];
    tp=[t1;t2(2:end);t3(2:end)];
    xp=y2x_nlin(yp,par);
    chip = x2chi(y2x(yp,par),par);

end

t=[0:sampling:sampling*(length(yp)-1)];
BalanceMapPendula(tp,yp,par,20)
% skip=1;
% Movie_BalanceMapPendula(tp,yp,skip,par,20);
DrawGraph(yp(1:end,:),yn(1:end,:),par,tp(1:end),tn(1:end),10)
[PosMargin, NegMargin]=MarginFunction(chip,par);
figure(30);
plot(t,[PosMargin, NegMargin]);grid on;
figure(30);

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

function ydot=f_nlin(t,y,par)    
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

function ydot=f_nlin_dist(t,y,par,u)    
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

function ydot=f_lin(t,y,par)   
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


function ydot=f_lin_dist(t,y,par,u)   
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


function BalanceMapPendula(t,y,par,fignum)
% This function shows animations of trajectory on the balance map and stick
% picture of the compass gait model

figure(fignum);clf;
[DataSize,~]=size(y);
x = y2x(y,par);
chi = x2chi(x,par);
EquivInitialChi = VirtualInitialState(chi,par);
BalanceMap(999);hold on;ax1=gca;
figure(fignum);set(gcf,'position', [50   500   960   420]);
s1=subplot(1,2,1);
figobj1=get(ax1,'children');
copyobj(figobj1,s1);axis([-0.6 0.6 0 1.2]);axis square;

xlabel('Position');
ylabel('Velocity');
set(gca,'Fontsize', 16);
close(999)
h_BalancemapPoint=animatedline(par.ChiOmega0*EquivInitialChi(1,1),EquivInitialChi(1,2),'Marker','o','MarkerEdgeColor','r','MarkerSize',10);

subplot(1,2,2);
plot([-10;10], [0,0]);
axis([-0.8 0.8 -0.4 1.2]);
set(gca,'Ytick',[0 0.5 1]);
axis square;
grid on;
xlabel('Horizotal axis (m)');ylabel('Vertical axis (m)');set(gca,'Fontsize', 16);

Hip=[-par.LI*sin(y(:,1)),par.LI*cos(y(:,1))];
Foot=Hip+[par.LI*sin(y(:,2)),-par.LI*cos(y(:,2))];
x=[0,Hip(1,1),Foot(1,1)];
y=[0,Hip(1,2),Foot(1,2)];
h_stick = animatedline(x,y,'Color','k','LineWidth',2);
for cnt1=2:1:DataSize
    subplot(1,2,1);
    clearpoints(h_BalancemapPoint);
    addpoints(h_BalancemapPoint,par.ChiOmega0*EquivInitialChi(cnt1,1),EquivInitialChi(cnt1,2));
    %set(balancemapPoint,'XData',par.ChiOmega0*EquivInitialChi(cnt1,1),'YData',EquivInitialChi(cnt1,2));
    drawnow;
    subplot(1,2,2);    
    x=[0,Hip(cnt1,1),Foot(cnt1,1)];
    y=[0,Hip(cnt1,2),Foot(cnt1,2)];
    clearpoints(h_stick);
    addpoints(h_stick,x,y);
    drawnow;
    pause(0.01);
end

function DrawGraph(y,ydef,par,t,tdef,fignum)
% This function shows the trajectories of Xst and Xsw and balance map (normal and perturbed conditions)

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
plot(tdef,[-ydef(:,1),ydef(:,2)],'--');hold on;
% tp=0:sampling:sampling*(length(y)-1);
plot(t,[-y(:,1),y(:,2)]);
grid on;
xlabel('Time (s)'); ylabel('Position (m)');
set(gca,'YLim',[-0.4 0.4]);set(gca,'XLim',[0 0.65]);


BalanceMap(999);hold on;ax1=gca;;hold on;ax1=gca;
figure(fignum);set(gcf,'position', [800   500   960   420]);
s1=subplot(1,2,2);
figobj1=get(ax1,'children');
copyobj(figobj1,s1);axis([-1 1 0 1.2]);
close(999)
hold on;
plot(par.ChiOmega0*EquivInitialChi(1,1),EquivInitialChi(1,2),'o','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChi(DataSize,1),EquivInitialChi(DataSize,2),'s','MarkerEdgeColor','r','MarkerSize',10);
plot(par.ChiOmega0*EquivInitialChi(:,1),EquivInitialChi(:,2),'-r');
plot(par.ChiOmega0*EquivInitialChidef(1,1),EquivInitialChidef(1,2),'ob');
plot(par.ChiOmega0*EquivInitialChidef(end,1),EquivInitialChidef(end,2),'sb');
plot(par.ChiOmega0*EquivInitialChidef(:,1),EquivInitialChidef(:,2),'-b');




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