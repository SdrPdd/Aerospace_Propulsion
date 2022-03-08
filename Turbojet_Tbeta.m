clear all
close all
clc

%% TURBOJET VARYING OF T_4 AND beta_c

%% DATI INIZIALI

V0 = 0;             %Fixed point condition [m/s]
pa = 101325;        %[Pa]
Ta = 288;           %[K]
Qf = 43e6;          %[J/kg]
Au = 0.126;         %[m^2]
cp = 1005;          %[J/kg K]
gamma = 1.4;
R = cp*(gamma-1)/gamma;


beta_c = 1:50; 
T04_v = 1200:200:1600;

% Efficiency
eps_d = 0.95;
eta_c = 0.88;
eta_t = 0.94;

%% INITIALIZATIONS
% Working with Array it is better to initialize 
% the loop outputs as a vector of zeros

S = zeros(length(T04_v),length(beta_c));
TSFC = zeros(length(T04_v),length(beta_c));
Ssp = zeros(length(T04_v),length(beta_c));
eta_th = zeros(length(T04_v),length(beta_c));
eta_p = zeros(length(T04_v),length(beta_c));
eta = zeros(length(T04_v),length(beta_c));
p7 = zeros(size(beta_c));
T7 = zeros(size(beta_c));
Vu = zeros(size(beta_c));
strozz_r = zeros(size(beta_c));

%% REAL CASE

for i = 1:length(T04_v)
    T04 = T04_v(i);
    
    %DIFFUSER
    p02 = pa.*eps_d;
    T02 = Ta;           %Hp adiabatic diffuser
    
    %COMPRESSOR
    p03 = p02.*beta_c;
    T03 = T02.*(1+(beta_c.^((gamma-1)./gamma)-1)./eta_c);
    
    %COMBUSTION CHAMBER
    % Hp we have no losses
    p04 = p03;
    f = cp.*(T04-T03)./Qf;
        
    %TURBINE
    T05 = T04-(T03-T02)./(1+f);
    T05_i = T04-(T04-T05)./eta_t;
    p05 = p04.*(T05_i./T04).^(gamma./(gamma-1));
    
    %NOZZLE
    % Hp adapted nozzle
    % Hp pressure and temperature at the nozzle inlet are in stagnant conditions
    p07 = p05;
    T07 = T05;
    p_cr_p0 = -1; %((gamma+1)/2^(gamma/(1-gamma)));
    
    for ib = 1:length(beta_c)
        p07_b = p07(ib);
        T07_b = T07(ib);
        
        if pa/p07_b <=p_cr_p0
            p7(ib) = p07_b.*p_cr_p0;
            T7(ib) = 2./(gamma+1).*T07_b;
            Vu(ib) = (gamma.*R.*T7(ib)).^0.5;
            strozz_r(ib) = 1;
        else 
            p7(ib) = pa;
            T7(ib) = T07_b.*(p7(ib)./p07_b).^((gamma-1)./gamma);
            if p7(ib)>p07(ib)
                Vu(ib) = 0;
            else
                Vu(ib) = (2.*cp.*T07_b.*(1-(pa./p07_b)^((gamma-1)/gamma)))^0.5;
            end
            strozz_r(ib) = 0;
        end
    end
    
    rho7 = p7./R./T7;
    
    
    mu = rho7.*Vu.*Au;
    ma = mu./(1+f);
    mf = ma.*f;
    
    % EFFICIENCY
    S(i,:) = ma.*((1+f.*Vu-V0)+p7-pa).*Au;
    TSFC(i,:) = mf./S(i,:);
    Ssp(i,:) = S(i,:)./ma;
    eta_th(i,:) = (Vu.^2-V0.^2)./(2.*f.*Qf);
    eta_p(i,:) = 2.*V0./Vu./(1+V0./Vu);
    eta(i,:) = (Vu-V0).*V0./f.*Qf;
    
    %PLOT
    figure(1)
    pl = plot(beta_c,S(i,:),'DisplayName',['T04 = ',num2str(T04),'K']);
    hold all
    grid on
    grid minor
    xlabel('beta_c','Interpreter','Tex');
    ylabel('S[N]','Interpreter','Tex');
    legend show
    
    figure(2)
    pl = plot(beta_c,TSFC(i,:),'DisplayName',['T04 = ',num2str(T04),'K']);
    hold all
    grid on
    grid minor
    xlabel('beta_c','Interpreter','Tex');
    ylabel('TSFC','Interpreter','Tex');
    legend show
    
end
    
    
