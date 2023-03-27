function [Do] = dynaRespTD_noCoupling(Bridge,u,w,t,DOF,varargin)
% function [Do] = dynaRespTD_noCoupling(Bridge,u,w,t,DOF,varargin) Computes
% the bridge displacement response from the bridge and turbulent wind
% characteristics. The input force are linearised and the modal responses
%  are uncoupled.
% 
% INPUTS
% Bridge: structure variables that contains the information about the bridge
% u : matrix [Nyy x N]: Horizontal wind velocity histories normal to the deck
% w : matrix [Nyy x N]: Vertical wind velocity histories normal to the deck
% t:  vector 1xN]: time vector
% varargin
%  rho: air density (default value is 1.25)
%  k: term introducing an aerodynamic damping to the torsional deck motion.
%     default value is 0 (0.25 for flat plate aerodynamic)
% 
% OUTPOUT
%  Do : matrix [Nyy x N] of the bridge displacement histories along Nyy
%  nodes, for the selected degrees of freedom, at each time step (N steps in
%  total)
%
% Author: E. Cheynet - UiS/UiB - 23.11.2019
%

p = inputParser();
p.CaseSensitive = false;
p.addOptional('rho',1.25);
p.addOptional('k',0);
p.parse(varargin{:});
% shorthen the variables name
rho  = p.Results.rho ;
k  = p.Results.k ;

% shortened variable
D = Bridge.D;
B = Bridge.B;
Cd = Bridge.Cd;
Cl = Bridge.Cl;
Cm = Bridge.Cm;
dCd = Bridge.dCd;
dCl = Bridge.dCl;
dCm =Bridge.dCm;
phi = Bridge.phi;
wn = Bridge.wn;
zetaStruct = Bridge.zetaStruct;
dt = median(diff(t));

N = numel(t);

% PREALLOCATION
% Main loop
meanU = nanmean(u);
[Nyy]= size(phi,3);

Do = zeros(Nyy,N);
Vo = zeros(Nyy,N);
dummyVo = zeros(Nyy,1);
dummyDo = zeros(Nyy,1);

for idt=1:N
    [Fmodal,M,K,C] = getLoad_1DOF(DOF,Bridge,meanU,u(idt,:)',w(idt,:)',dummyDo,dummyVo);
    if idt ==1 % initial acceleration
        DoM = zeros(size(M));
        VoM = zeros(size(M));
        AoM = M\(Fmodal-C.*VoM-K.*DoM);
    end
    [DoM,VoM,AoM,Do(:,idt),Vo(:,idt)] = Newmark_1D(DOF,dt,DoM,VoM,AoM,Fmodal,M,K,C);
    
    dummyVo = Vo(:,idt);
    dummyDo = Do(:,idt);
    
end


    function [Fmodal,M,K,C] = getLoad_1DOF(DOF,Bridge,meanU,u1,w1,Do,Vo)
        
        if round(abs(mean(u1(:))))== abs(mean(meanU(:)))
            u1 = u1(:)'-meanU(:)';
        end
        X = Bridge.y.*Bridge.L;
        CST = 1/2*B*meanU*rho;
        switch DOF
            % -------------------------------------------------------------------------
            case 'lateral'
                myZeta = zetaStruct(1,:);
                myPhi = squeeze(phi(1,:,:));
                myWn = wn(1,:);
                Fmean = CST.*meanU.*(D/B)*Cd;
                Mtot = Bridge.m+2*Bridge.mc; % total mass of deck along y axis in kg/m
                Bq1 = CST(:)*[2.*(D/B)*Cd, (D/B*dCd-Cl)];% 1 x 2
                Cae1 = CST(:)*[-2*D/B*Cd];% 1 x 1
                Kae1 = zeros(size(Cae1));% 1 x 1
                % -------------------------------------------------------------------------
            case 'vertical'
                myZeta = zetaStruct(2,:);
                myPhi = squeeze(phi(2,:,:));
                myWn = wn(2,:);
                Fmean = CST.*meanU.*[Cl];
                Mtot = Bridge.m+2*Bridge.mc; % total mass of deck along z axis in kg/m
                Bq1 = CST(:)*[2*Cl, (dCl+D/B*Cd)];% 1 x 2
                Cae1 = -CST(:)*[dCl + D/B*Cd];% 1 x 1
                Kae1 = zeros(size(Cae1));% 1 x 1
                % -------------------------------------------------------------------------
            case 'torsional'
                myZeta = zetaStruct(3,:);
                myPhi = squeeze(phi(3,:,:));
                myWn = wn(3,:);
                Fmean = CST.*meanU.*[B*Cm];
                Mtot = Bridge.m_theta; % kg.m^2/m
                Bq1 = CST(:).*[2*B*Cm B*dCm];% 1 x 2
                Cae1 = -CST(:).*k*dCm*B*B;% 1 x 1
                Kae1 = 1/2*rho*B*meanU.^2.*dCm*B;% 1 x 1
        end
        
        % MODAL MASS AND STIFNESS CALCULATION
        M_modal  = trapz(X,Mtot.*myPhi.^2,2)';
        K_modal = myWn.^2.*M_modal; % frequency must be in rad/Hz !!!
        C_modal = 2.*myWn.*M_modal.*myZeta;% frequency must be in rad/Hz !!!
        % no coupling:
        M = diag(M_modal);
        K = diag(K_modal);
        C = diag(C_modal);
        
        
        % wind load calculation
        F1 = sum(Bq1'.*[u1(:)';w1(:)']);        
        FCae1 = (Cae1(:).*Vo(:))';
        FKae1 = (Kae1(:).*Do(:))';
        
        Ftot = Fmean+FCae1+FKae1+F1;
        
        % modal load
        Fmodal  = diag(trapz(X,myPhi.*Ftot,2));
        
    end
    function [x1,dx1,ddx1,Do,Vo] = Newmark_1D(DOF,dt,x0,dx0,ddx0,F,M,K,C,varargin)
        % options: default values
        inp = inputParser();
        inp.CaseSensitive = false;
        inp.addOptional('alpha',1/12);
        inp.addOptional('beta',1/2);
        inp.parse(varargin{:});
        % shorthen the variables name
        alphaCoeff = inp.Results.alpha ;
        beta = inp.Results.beta;
        
        switch DOF
            case 'lateral'
                myPhi = squeeze(phi(1,:,:));
            case 'vertical'
                myPhi = squeeze(phi(2,:,:));
            case 'torsional'
                myPhi = squeeze(phi(3,:,:));
        end
        
        aDT2 = (alphaCoeff.*dt.^2);
        aDT = (alphaCoeff.*dt);
        
        A = (1./aDT2.*M+beta/aDT*C+K);
        B1 = F+M.*(1./aDT2*x0+1./aDT*dx0+(1/(2*alphaCoeff)-1)*ddx0);
        B2 = C.*(beta/aDT*x0+(beta/alphaCoeff-1).*dx0);
        B3 = C.*(beta/alphaCoeff-2)*dt/2*ddx0;
        
        x1 = A\(B1+B2+B3);
        ddx1= 1/aDT2.*(x1-x0)-1/aDT.*dx0-(1/(2*alphaCoeff)-1).*ddx0;
        dx1= dx0+(1-beta).*dt*ddx0+beta.*dt*ddx1;
        
        x2=diag(x1);
        dx2=diag(dx1);
        Do = myPhi'*x2;
        Vo = myPhi'*dx2;
    end

end
