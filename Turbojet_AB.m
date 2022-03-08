clear all
close all
clc

%% TURBOJET WITH AFTER BURNER

%% INITIAL DATA

V0 = 805/3.6;       %[m/s]
pa = 0.458e5;       %[Pa]
Ta = 248;           %[K]
Qf = 43e6;          %[J/kg]
Au = 0.0935;        %[m^2]
cp = 1005;          %[J/kg K]
gamma = 1.4;
R = cp*(gamma-1)/gamma;
beta_c = 4; 
T04 = 1100;         %[K]

% Efficiency
eps_d = 0.95;
eta_c = 0.88;
eta_t = 0.94;
eta_pb = 0.95;
eta_b = 0.95;
eta_pab = 0.92;
eta_ab = 0.95;
T07_ab = 1800;      %[K]


%% IDEAL CASE

% ASYNOTIC CONDITION
M0 = V0/sqrt(gamma*R*Ta);

% DIFFUSER
p1 = pa;
T1 = Ta;
p02_id = p1.*(1+(gamma-1)./2.*M0.^2).^(gamma./(gamma-1));
T02_id = T1.*(1+(gamma-1)./2.*M0.^2);

%COMPRESSOR
p03_id = beta_c*p02_id;
T03_id = beta_c^((gamma-1)/gamma)*T02_id;

% COMBUSTION CHAMBER
p04_id = p03_id;
f_id = cp*(T04-T03_id)/Qf;

%TURBINE
T05_id = T04-(T03_id-T02_id)/(1+f_id);
p05_id = p04_id*(T05_id/T04)^(gamma/(gamma-1));

%NOZZLE
p07_id = p05_id;
T07_id = T05_id;
p_cr_p0 = ((gamma+1)/2)^(gamma/(1-gamma));

if pa/p07_id <= p_cr_p0
    p7_id = p07_id*p_cr_p0;
    T7_id = 2/(gamma+1)*T07_id;
    Vu_id = (gamma*R*T7_id)^0.5;
    strozz_id = 1;
else
    p7_id = pa;
    T7_id = T07_id*(p7_id/p07_id)^((gamma-1)/gamma);
    Vu_id = (2*cp*T07_id*(1-pa/p07_id));
    strozz_id = 0;
end

rho7_id = p7_id/R/T7_id;


mu_id = rho7_id*Vu_id*Au;
ma_id = mu_id/(1+f_id);
mf_id = ma_id*f_id;

% EFFICIENCY
S_id = ma_id*((1+f_id)*Vu_id-V0)+(p7_id-pa)*Au;
TSFC_id = mf_id/S_id;
Ssp_id = S_id/ma_id;
eta_th_id = (Vu_id^2-V0^2)/(2*f_id*Qf);
eta_p_id = 2*V0/Vu_id/(1+V0/Vu_id);
eta_id = (Vu_id-V0)*V0/f_id*Qf;


%% REAL CASE
% DIFFUSER
p02 = p02_id*eps_d;
T02 = T02_id;

%COMPRESSOR
p03 = beta_c*p02;
T03 = (1+(beta_c^((gamma-1)/gamma)-1)/eta_c)*T02;

% COMBUSTION CHAMBER
p04 = p03*eta_pb;
f = cp*(T04-T03)/eta_b/Qf;

%TURBINE
T05 = T04-(T03-T02)/(1+f);
T05_i = T04-(T04-T05)/eta_t;
p05 = p04*(T05_i/T04)^(gamma/(gamma-1));

%NOZZLE
p07 = p05;
T07 = T05;

if pa/p07 <= p_cr_p0
    p7 = p07*p_cr_p0;
    T7 = 2/(gamma+1)*T07;
    Vu = (gamma*R*T7)^0.5;
    strozz_id = 1;
else
    p7 = pa;
    T7 = T07*(p7/p07)^((gamma-1)/gamma);
    Vu = (2*cp*T07*(1-pa/p07))^0.5;
    strozz_id = 0;
end

rho7 = p7/R/T7;


mu = rho7*Vu*Au;
ma = mu/(1+f);
mf = ma*f;

% EFFICIENCY
S= ma*((1+f)*Vu-V0)+(p7-pa)*Au;
TSFC = mf/S;
Ssp = S/ma;
eta_th = (Vu^2-V0^2)/(2*f*Qf);
eta_p = 2*V0/Vu/(1+V0/Vu);
eta = (Vu-V0)*V0/f*Qf;

%% AFTER BURNER
p07_ab = eta_pab*p05;
f_ab = (1+f)*cp*(T07_ab-T05)/(eta_ab*Qf-cp*T07_ab);

if pa/p07_ab <= p_cr_p0
    p7_ab = p07_ab*p_cr_p0;
    T7_ab = 2/(gamma+1)*T07_ab;
    Vu_ab = (gamma*R*T7_ab)^0.5;
    strozz_r = 1;
else
    p7_ab = pa;
    T7_ab = T07_ab*(p7_ab/p07_ab)^((gamma-1)/gamma);
    Vu_ab = (2*cp*T07_ab*(1-pa/p07_ab))^0.5;
    strozz_r = 0;
end

rho7_ab = p7_ab/R/T7_ab;

% NOZZLE WITH AFTER-BURNER
Au_ab = rho7*Vu*Au/rho7_ab/Vu_ab;


mu_ab = mu;
ma_ab = mu/(1+f+f_ab);
mf1_ab = ma_ab*f;
mf2_ab = ma_ab*f_ab;
mf_tot = mf1_ab+mf2_ab;

% EFFICIENCY
S_ab= ma_ab*((1+f+f_ab)*Vu_ab-V0)+(p7_ab-pa)*Au_ab;
TSFC_ab = mf_tot/S_ab;
Ssp_ab = S_ab/ma_ab;
eta_th_ab = (Vu_ab^2-V0^2)/(2*(f+f_ab)*Qf);
eta_p_ab = 2*V0/Vu_ab/(1+V0/Vu_ab);
eta_ab = (Vu_ab-V0)*V0/(f+f_ab)*Qf;

%% COMPARISON
Rapp_Au = Au_ab/Au;
Rapp_S = S_ab/S;
Rapp_TSFC = TSFC_ab/TSFC;
Rapp_Ssp = Ssp_ab/Ssp;
Rapp_eta_th = eta_th_ab/eta_th;
Rapp_eta_p = eta_p_ab/eta_p;
Rapp_eta = eta_ab/eta;











