function [delta]=StabilityBoundaryPositive_Model(E0_array,omega0,mode,par)
delta=zeros(1,length(E0_array));

for cnt1=1:length(E0_array)
    %BinarySearch 
    deltaMax = 3;
    deltaMin = 0;
    
    ipsi = 1;
    while ipsi>0.001
        %BinarySearch 
        [chi,t,par]=TwoLinkModel(E0_array(cnt1),(deltaMax+deltaMin)/2,mode,par);
        if par.ChiOmega0 ~= omega0
            error('something wrong with parameter of omega.');
        end
        if t(end)<1
            deltaMin=(deltaMax+deltaMin)/2;
        else
            deltaMax=(deltaMax+deltaMin)/2;
        end
        ipsi=deltaMax-deltaMin;
    end
    delta(cnt1)=deltaMin;
end


