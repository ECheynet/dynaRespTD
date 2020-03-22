
function [Sr] = dynaRespFD(Bridge,Wind,indR,varargin)
% [Sr,f] = dynaRespFD(Bridge,Wind,r) computes the displacement
% response of a structure subjected to a turbulent wind load in the
% frequency domain. The response is therefore a power spectral density
% (PSD) Sr that expresses as a function of the frequency vector f.
%
%% Syntax
%
%  [Sr,f] = dynaRespFD(Bridge,Wind,r)
%  [Sr,f] = dynaRespFD(Bridge,Wind,r,'cohType','empirical')
%  [Sr,f] = dynaRespFD(Bridge,Wind,r,'cohType','ESDU')
%  [Sr,f] = dynaRespFD(Bridge,Wind,r,'Nfreq',1024)
%
% Input
%  * Bridge:structure variable containing structure properties
%  * Wind : structure variable containing wind turbulence info
%  * indR: integer : indice indicating the position of the calculated response
%  Optional inputs:
%       - 'cohType': 'empirical' (default) or 'ESDU' models
%       - 'Nfreq': number of frequency points (default:Nfreq = 2048). Nfreq
%       is not used is the measured time series are used as input.
%       - term introducing an aerodynamic damping to the torsional deck motion.
% %     default value is 0 (0.25 for flat plate aerodynamic)
%
% Output
%
% * Sr : PSD of the bridge response at a single position y(indR)
% * f :  frequency vector associated with Sr
%
%% Example
%
% %  For this example,we use the suspension bridge model for the Lysefjord
% bridge in Norway, and we plot the PSD of the vertical response
% displacement
%
% load bridgeData
% [Srz,f] = dynaRespFD(Bridge,Wind,10); % compute PSD of response
% figure
% loglog(f,Srz);
% xlabel(' frequency (Hz)');
% ylabel('S_{rz} (m^2/Hz)');
%
%% Author Info
% Author: Etienne Cheynet, Univeristy of Stavanger, Norway.
% last modified: 10/07/2019
%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('cohType','empirical');
p.addOptional('pos',[]);
p.addOptional('k',0);
p.addOptional('rho',1.25);  % air density
p.addOptional('AA','none'); % Aero Admittance (AA) lateral
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
cohType = p.Results.cohType ;
AA = p.Results.AA;
k = p.Results.k;
rho = p.Results.rho;
%% Preprocessing
y = Bridge.y*Bridge.L; % Nyy = numel(y);
f = Wind.f;
Su = Wind.Su;
Sw = Wind.Sw;
% Cross-spectrum between u and w included or not
if ~isfield('Wind','Co_uw')
    Co_uw = zeros(size(Su));
else
    Co_uw = Wind.Co_uw;
end
Nfreq = numel(f);
N=size(Su,2);
if N~=Nfreq
    error('size(Su,2) ~=numel(f) is not acceptable. You should have an  numel(f)== size(Su,2)');
end
%% Get basic wind and bridge properties
% get matrix of mass, buffeting forces, aerodynamic damping and stifness
[M,Bq,Cae,Kae] = getBasicProperties(Wind,Bridge);
% Get aerodynamic cross-admittance functions
K = f(:)*Bridge.B./Wind.U;% reduced frequency
[Xiu,Xiw] = getAeroAdmittance(AA,K);
if size(Xiu,1)>= size(Xiu,2)
    Xiu=Xiu';
    Xiw=Xiw';
end
%% MODAL MASS AND STIFNESS CALCULATION
[N1,N2]=size(Bridge.phi); % must be [Nmodes,Nyy]
if N1>N2,    Bridge.phi=Bridge.phi';end
Cae_modal = Cae*trapz(y,Bridge.phi.^2,2)';
Kae_modal = Kae*trapz(y,Bridge.phi.^2,2)';
Mtot  = trapz(y,M.*Bridge.phi.^2,2)';
K_modal = Bridge.wn.^2.*Mtot;
C_modal = 2.*Bridge.wn.*Mtot.*Bridge.zetaStruct;
Ktot = K_modal-Kae_modal;
Ctot  = C_modal-Cae_modal;
if any(Ktot<=0)
    warning('torsional divergence reached !')
    Sr = nan(1,Nfreq);
    return
end
%% BRIDGE RESPONSE CALCULATION
% matrix distance between each nodes:
dy = abs(bsxfun(@plus,y,-y'));
% matrix wind velocity for cohernece ( if different altitudes)
U_coh = bsxfun(@plus,Wind.U,Wind.U')*0.5;
Sr = zeros(1,Nfreq);
for ii=1:Nfreq % at each frequency step:
    % get spectral matrix of wind  fluctuations with coherence
    [Svv] = getSvv(U_coh,dy,f(ii),Wind,Su(:,ii),Sw(:,ii),Co_uw(:,ii),cohType);
    % Get PSD of wind load (physical base)
    [Sqq] = getSq(Svv,Bq,Xiu(:,ii),Xiw(:,ii));
    % switch from physical to modal base
    [SQ] = JointAcceptance(Sqq,Bridge.phi,y);
    % mechanical admittance function
    omega_i = 1i*2*pi*f(ii);
    Hmeca = 1./(Ktot + Ctot.*omega_i+Mtot.*(omega_i).^2);
    % bridge response at each frequency step in the modal base
    SR = abs(Hmeca.^2.*SQ);
    % come back to physical base
    Sr(ii) = SR*Bridge.phi(:,indR).^2;
end
%%
% *************************************************************************
% ********************       NESTED FUNCTIONS      ************************
% *************************************************************************
    function [M,Bq,Cae,Kae] = getBasicProperties(Wind,Bridge)
        %  No modal coupling
        CST = 1/2*Bridge.B*Wind.U*rho;
        switch Bridge.DOF
            case 'lateral'
                M = Bridge.m+2*Bridge.mc; % total mass of deck along y axis in kg/m
                Bq = CST(:)*[2.*(Bridge.D/Bridge.B)*Bridge.Cd, (Bridge.D/Bridge.B*Bridge.dCd-Bridge.Cl)];% 1 y 2
                Cae = -CST(:)*(2*Bridge.D/Bridge.B*Bridge.Cd);% 1 y 1
                Kae = zeros(size(Cae));% 1 y 1
            case 'vertical'
                M = Bridge.m+2*Bridge.mc; % total mass of deck along z axis in kg/m
                Bq = CST(:)*[2*Bridge.Cl, (Bridge.dCl +Bridge.D/Bridge.B*Bridge.Cd)];% 1 y 2
                Cae = -CST(:)*(Bridge.dCl+Bridge.D/Bridge.B*Bridge.Cd);% 1 y 1
                Kae = zeros(size(Cae));% 1 y 1
            case 'torsional'
                M = Bridge.m_theta; % kg.m^2/m
                Bq = CST(:)*[2*Bridge.B*Bridge.Cm, Bridge.B*Bridge.dCm];% 1 y 2
                Cae = -CST(:)*(k*Bridge.dCm*Bridge.B.^2);% 1 y 1
                Kae = CST(:).*Wind.U*(Bridge.dCm*Bridge.B);% 1 y 1
            otherwise
                error(['You have written: Bridge.DOF = ''',Bridge.DOF,'''. You have to write ''lateral'', ''vertical or ''torsional'' ']);
        end
    end
    function [Xiu,Xiw] = getAeroAdmittance(AA,K)
        if strcmpi(AA,'Liepmann')
            Xiu= sqrt(1./(1+2*pi.*K));
            Xiw= sqrt(1./(1+2*pi.*K));
        elseif strcmpi(AA,'Holmes')
            Xiu= sqrt(1./(1+4*K));
            Xiw= sqrt(1./(1+4*K));
        elseif strcmpi(AA,'none')
            Xiu = ones(size(K));
            Xiw = ones(size(K));
        else
            error('Aerodynamic admittance specified is unknown');
        end
    end
    function [Svv] = getSvv(U_coh,dy,f,Wind,Su,Sw,Co_uw,cohType)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get coherence for u and w
        if strcmpi(cohType,'empirical')
            [cohU] = coh1Para(U_coh,dy,f,Wind.Cuy);% cocoherence for u-component
            [cohW] = coh1Para(U_coh,dy,f,Wind.Cwy);% cocoherence for w-component
            if max(Co_uw)==0,
                cohUW = zeros(size(cohU));
            else
                cohUW = 1/2.*(Co_uw(:)*Co_uw(:)')/(Su(:)*Sw(:)').*(cohU+cohW);
            end
        elseif strcmpi(cohType,'ESDU')
            [cohU,~,cohW] = cohESDU86010(dy,f,Wind.U);
            cohUW = 1/2.*(Co_uw(:)*Co_uw(:)')/(Su(:)*Sw(:)').*(cohU+cohW);
        else
            error('cohType is unknown. it must be ''ESDU'', or''empirical'' ')
        end
        % beuild the spectral matrix
        % Svv is a 3D matrix [2 y Nyy y Nyy]
        Svv(1,1,:,:) = sqrt(Su(:)*Su(:)').*cohU;
        Svv(2,2,:,:) = sqrt(Sw(:)*Sw(:)').*cohW;
        Svv(1,2,:,:) =  sqrt(Co_uw(:)*Co_uw(:)').*cohUW;
        Svv(2,1,:,:) = Svv(1,2,:,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    function [Sqq] = getSq(Svv,Bq,Xiu,Xiw)
        A = squeeze((Bq(1,1)*Xiu).^2.*Svv(1,1,:,:));
        B = squeeze((Bq(1,2)*Xiw).^2.*Svv(2,2,:,:));
        C = squeeze((Bq(1,1)*Xiu.*Bq(1,2)*Xiw).*Svv(1,2,:,:));
        D = squeeze((Bq(1,2)*Xiw.*Bq(1,1)*Xiu).*Svv(2,1,:,:));
        Sqq = squeeze(A+B+C+D); % is [ Nyy y Nyy ]
    end
    function [SQ] = JointAcceptance(Sqq,phi,Y)
        % [SQ] = JointAcceptance(Sqq,phi,Y) calculates the PSD of the wind load
        % in the modal base, based on the PSD  of the wind load in the physical
        % base.
        %% Inpout:
        % phi: mode shapes defined in eigenBridge.m
        % Sqq : PSD of the partially correlated wind load . It size is [Nyy y Nyy]
        % Y : vector [1 y Nyy] : discretisation of bridge girder into Nyy points .
        % !!!! Y is not normalized, and this is excpected.
        %% Output
        %  SQ : PSD of the wind load in the modal base, is a [Nmodes y Nmodes]
        %  matrix
        %%
        [Nmodes,~]=size(phi);
        SQ = zeros(1,Nmodes);
        for pp =1:Nmodes
            phiphi = squeeze(phi(pp,:))'*squeeze(phi(pp,:));
            dummy = phiphi.*Sqq;
            SQ(pp) = trapz(Y,trapz(Y,dummy,2),1);
        end
        %  clear dummy* phiphi
        SQ = abs(SQ);
    end
    function [coh] = coh1Para(U,dy,freq,C)
        % force horizontal vector.
        if numel(C)>1
            error('too many coefficient are specified for the coherence function');
        end
        coh = exp(-C.*dy.*freq/U);
    end
    function [cohU,cohV,cohW] = cohESDU86010(D,f,U)
        % [coh] = cohESDU86010(z,h,Lu,D,f,U,component,option) Calculate the
        % simplified ESDU coherence model.
        % INPUTS
        % D: scalar:  spatial separation
        % f: scalar:  frequency
        % U: scalar:  mean wind speed
        % reduced frequency
        eta1 = 2*pi.*f.*D./U;
        % root-coherence
        cohU = exp(-1.15.*eta1.^(1.5));
        cohV = exp(-0.65.*eta1.^(1.3));
        cohW = exp(-0.65.*eta1.^(1.3));
    end
% *************************************************************************
% ********************             END             ************************
% *************************************************************************
end
