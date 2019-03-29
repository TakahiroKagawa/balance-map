function Te=TwoLinkFallDetection(E0,psi,X2_dot0)
if nargin < 1
    E0negative=[-0.25:0.05:-0.05,-0.025];
    E0positive=[0.025,0.05:0.05:1.25];
    E0=[E0negative,E0positive];
%      E0=0;
end
if nargin < 2
     psi=[-0.5:0.25:1.25];
%     psi=0.1:0.1:1;
end

if nargin <3
    X2_dot0=2;
end

Te=size(length(E0),length(psi));
for cnt1=1:length(E0)
    for cnt2=1:length(psi)
        %Computer simulation of swing movement with two-link model
        [chi,t,par]=TwoLinkModel(E0(cnt1),psi(cnt2),X2_dot0);
         Te(cnt1,cnt2)=t(end);
    end
end


BalanceMap(1,par);
figure(1)
for cnt1=1:length(E0)
    for cnt2=1:length(psi)
        if E0(cnt1)>0
            y_p=sqrt(E0(cnt1))/par.Omega0*sinh(psi(cnt2));
            yd_p=sqrt(E0(cnt1))*cosh(psi(cnt2));
        else
            if E0(cnt1)<0
                y_p=sqrt(-E0(cnt1))/par.Omega0*cosh(psi(cnt2));
                yd_p=sqrt(-E0(cnt1))*sinh(psi(cnt2));
            else
                % for par.E0=0, par.psi indicates chi_dot_0
                y_p=psi(cnt2)/par.Omega0;
                yd_p=psi(cnt2);
            end
        end
        if Te(cnt1,cnt2)<1
            %Stable touchdown
            plot(par.Omega0*y_p,yd_p,'ob','Markersize',10);
        else
            %Unstable touchdown
            plot(par.Omega0*y_p,yd_p,'xr','Markersize',10);
        end
    end
end

