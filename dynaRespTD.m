function [Do] = dynaRespTD(Bridge,u,w,t,varargin)
% function [Do] = dyna_NL_TD(Bridge,u,w,t,varargin) computes the bridge
% wind-induced lateral, vertical and torsional displacement response in the
% time domain. 
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
%  Do : matrix [3 x Nyy x N] of the bridge displacement histories along Nyy
%  nodes, for the three degrees of freedom, at each time step (N steps in
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

if abs(mean(u(1,:),'omitnan'))<1e-4,    warning('u has a zero mean value');end

% Maximal incidence angle accepted (radians)
alphaMaxD = 20*pi/180;
alphaMinD = -20*pi/180;
alphaMax = 20*pi/180;
alphaMin = -20*pi/180;


dt = median(diff(t));
N = numel(t);

% Definition of Nyy, and Nmodes:
[~,Nmodes,Nyy]= size(Bridge.phi);

if size(u,2)==Nyy && size(u,1)==N
    u=u';
elseif size(u,2)~=Nyy && size(u,1)~=Nyy
    error(' u is not defined in Nyy nodes')
elseif  size(u,2)~=N && size(u,1)~=N
    error(' u is not defined in N time steps')
end

if size(w,2)==Nyy && size(w,1)==N
    w=w';
elseif size(w,2)~=Nyy && size(w,1)~=Nyy
    error(' w is not defined in Nyy nodes')
elseif  size(w,2)~=N && size(w,1)~=N
    error(' w is not defined in N time steps')
end


if max(Bridge.y)== 1
    y = Bridge.y.*Bridge.L;
else
    warning('Bridge.y is not normalized');
end

%% MODAL MASS AND STIFNESS CALCULATION
Mtot = diag([Bridge.m+2*Bridge.mc,Bridge.m+2*Bridge.mc,Bridge.m_theta]);
phi0 = reshape(phi,[],Nyy);
phi0_N = bsxfun(@times,phi0,1./max(abs(phi0),[],2));
Nm = size(phi0,1);

Mtot = repmat(Mtot,round(Nm/3),round(Nm/3));
M = zeros(Nm);

for pp=1:Nm
    for qq=1:Nm/3
        if (pp+3*qq)<=Nm
            Mtot(pp,pp+3*qq) = 0;
            Mtot(pp+3*qq,pp) = 0;
        end
    end
    for qq=1:Nm
        M(pp,qq)  = trapz(y,phi0_N(pp,:).*phi0_N(qq,:).*Mtot(pp,qq));
    end
end

K = diag(wn(:)).^2*M;
C = 2.*diag(wn(:))*M*diag(zetaStruct(:));

%% PREALLOCATION
% Check for second derivative of aerodynamic coefficients
if ~isfield(Bridge,'ddCd'),    ddCd=0;else    ddCd = Bridge.ddCd;end
if ~isfield(Bridge,'ddCl'),    ddCl=0;else    ddCl = Bridge.ddCl;end
if ~isfield(Bridge,'ddCm'),    ddCm=0;else    ddCm = Bridge.ddCm;end

if ~isfield(Bridge,'d3Cd'),    d3Cd=0;else    d3Cd = Bridge.d3Cd;end
if ~isfield(Bridge,'d3Cl'),    d3Cl=0;else    d3Cl = Bridge.d3Cl;end
if ~isfield(Bridge,'d3Cm'),    d3Cm=0;else    d3Cm = Bridge.d3Cm;end

Cdrag=[Cd,dCd,ddCd,d3Cd];
Clift=[Cl,dCl,ddCl,d3Cl];
Cmoment=[Cm,dCm,ddCm,d3Cm];

%% INITIALISATION
Do = zeros(3,Nyy,N);
rt=zeros(Nyy,1);
dry=zeros(Nyy,1);
drz=zeros(Nyy,1);
drt=zeros(Nyy,1);

for idt=1:N
    % get  modal forces
    [Fmodal]= getFmodal(u(:,idt),w(:,idt),dry,drz,drt,D,B,rt);
    
    % Numerical solver
    if idt ==1
        % initial acceleration
        DoM = zeros(size(M));
        VoM = zeros(size(M));
        AoM = M\(Fmodal-C.*VoM-K.*DoM);
    end
    
    [DoM,VoM,AoM,Do(:,:,idt),Vo] = Newmark(dt,DoM,VoM,AoM,Fmodal,M,K,C);
    
    rt = Do(3,:,idt)';
    dry = Vo(1,:)';
    drz = Vo(2,:)';
    drt = Vo(3,:)';
end

%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%
    function [Fmodal]= getFmodal(u,w,dry,drz,drt,D,B,rt)
        Wtot = (w-drz-k*B.*drt);
        Vrel = sqrt((u-dry).^2+Wtot.^2); % is [Nyy x 1]
        beta = atan2(Wtot,(u-dry)); % is [Nyy x 1]
        alpha = rt+beta; % is [Nyy x N]
        COEFF = 1/2*rho*B.*Vrel.^2;
        
        Fdrag = D/B.*getCoeffAero(Cdrag,alpha,alphaMinD,alphaMaxD);
        Flift = getCoeffAero(Clift,alpha,alphaMin,alphaMax);
        Fmoment = B.*getCoeffAero(Cmoment,alpha,alphaMin,alphaMax);
        
        Ftot(:,1)=COEFF.*(Fdrag.*cos(beta) - Flift.*sin(beta));
        Ftot(:,2)=COEFF.*(Fdrag.*sin(beta) + Flift.*cos(beta));
        Ftot(:,3)=COEFF.*(Fmoment);
        
        F = repmat(Ftot,1,Nmodes)';
        Fmodal = trapz(y,F.*phi0_N,2);
        Fmodal = diag(Fmodal(:));
    end
    function [x1,dx1,ddx1,Do,Vo] = Newmark(dt,x0,dx0,ddx0,F,M,K,C,varargin)
        % options: default values
        inp = inputParser();
        inp.CaseSensitive = false;
        inp.addOptional('alpha',1/12);
        inp.addOptional('beta',1/2);
        inp.parse(varargin{:});
        % shorthen the variables name
        alphaCoeff = inp.Results.alpha ;
        beta = inp.Results.beta;
        
        aDT2 = (alphaCoeff.*dt.^2);
        aDT = (alphaCoeff.*dt);
        
        A = (1./aDT2.*M+beta/aDT*C+K);
        B1 = F+M.*(1./aDT2*x0+1./aDT*dx0+(1/(2*alphaCoeff)-1)*ddx0);
        B2 = C.*(beta/aDT*x0+(beta/alphaCoeff-1).*dx0);
        B3 = C.*((beta/alphaCoeff-2)*dt/2*ddx0);
        
        x1 = A\(B1+B2+B3);
        ddx1= 1/aDT2.*(x1-x0)-1/aDT.*dx0-(1/(2*alphaCoeff)-1).*ddx0;
        dx1= dx0+(1-beta).*dt*ddx0+beta.*dt*ddx1;
        
        x2=diag(x1);
        dx2=diag(dx1);
        Vo = zeros(3,Nyy);
        Do = zeros(3,Nyy);
        for oo = 1:3
            Do(oo,:) = phi0(oo:3:end,:)'*x2(oo:3:end);
            Vo(oo,:) = phi0(oo:3:end,:)'*dx2(oo:3:end);
        end
    end
    function Caero = getCoeffAero(C,alpha,alphaMin,alphaMax)
        alpha(alpha>alphaMax)=alphaMax;
        alpha(alpha<alphaMin)=alphaMin;
        Caero = C(1) + C(2)*alpha + 1/2*C(3)*alpha.^2 + 1/6*C(4)*alpha.^3; % Polynomial approximation
    end


end