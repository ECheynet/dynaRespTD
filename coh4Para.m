    function [coh] = coh4Para(Cy,Cz,dy,dz,U,f)
        
        % check vertical coeff
        if numel(Cy)>4;    error('C must be a vector with  max 4 columns or rows');end
        if numel(Cy)==1,
            Cy = [Cy(:)',0,1,0];
        elseif  numel(Cy)==2,
            Cy = [Cy(:)',1,0];
        elseif  numel(Cy)==3,
            Cy = [Cy(:)',0];
        end
        
        % check vertical coeff
        if numel(Cz)>4;    error('C must be a vector with  max 4 columns or rows');end
        if numel(Cz)==1,
            Cz = [Cz(:)',0,1,0];
        elseif  numel(Cz)==2,
            Cz = [Cz(:)',1,0];
        elseif  numel(Cz)==3,
            Cz = [Cz(:)',0];
        end
        
        % lateral separation
        a1 = (Cy(1).*f).^2;
        a2 = Cy(2).^2;
        AA = (sqrt(a1+a2).*dy./U).^(Cy(3));
        
        % vertical separation
        b1 = (Cz(1).*f).^2;
        b2 = Cz(2).^2;
        BB = (sqrt(b1+b2).*dz./U).^(Cz(3));
        
        % combination
        coh = exp(-sqrt(AA.^2+BB.^2));
        coh = coh.*cos(Cy(4).*dy.*f./U).*cos(Cz(4).*dz.*f./U);
        
        
    end