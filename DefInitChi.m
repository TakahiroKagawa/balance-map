function  chi0 = DefInitChi(T0,E0,psi,Omega0)
% define initial normalized position
if E0>0
    chi0(1)=sqrt(E0)/Omega0*sinh(Omega0*T0+psi);
    chi0(2)=sin(T0);
    chi0(3)=sqrt(E0)*cosh(Omega0*T0+psi);
    chi0(4)=cos(T0);
else
    if E0<0
        chi0(1)=sqrt(-E0)/Omega0*cosh(Omega0*T0+psi);
        chi0(2)=sin(T0);
        chi0(3)=sqrt(-E0)*sinh(Omega0*T0+psi);
        chi0(4)=cos(T0);
    else
        % for par.E0=0, par.psi indicates chi_0
        chi0(1)=psi/Omega0*exp(Omega0*T0);
        chi0(2)=sin(T0);
        chi0(3)=psi*exp(Omega0*T0);
        chi0(4)=cos(T0);
    end
end