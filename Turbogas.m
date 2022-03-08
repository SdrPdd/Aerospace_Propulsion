clear all
close all
clc

%% TURBOGAS


%DATA
p1 = 101325;        %[Pa]
T1 = 288;           %[K]
beta_c = 1:200;
T3_v = 1200:200:1600;%[K]
gamma = 1.4;
R = 287;             %[J/kg k]
eta_c = 0.8;
eta_t = 0.9;
cp = gamma*R/(gamma -1);

% Initializations
eta_th_id = zeros(length(T3_v),length(beta_c));
eta_th_r = zeros(length(T3_v),length(beta_c));
Lu_id = zeros(length(T3_v),length(beta_c));
Lu_r = zeros(length(T3_v),length(beta_c));

for i = 1:length(T3_v)
    T3 = T3_v(i);
% IDEAL
p2 = beta_c.*p1;
T2_id = T1.*beta_c.^((gamma-1)./gamma);
eta_th_id(i,:) = 1-1./(beta_c.^((gamma-1)./gamma));
tau_c = T2_id./T1;
tau = T3./T1;

Lu_id(i,:) = cp.*T1.*(tau-tau_c).*(1-1./tau_c);

T4_id = T3./tau_c;

% REAL
T2 = T1 + (T2_id - T1)./eta_c;
eta_th_r(i,:) = eta_th_id(i,:).*(tau.*eta_c.*eta_t-tau_c)./ ...
    ((tau.*eta_c-tau_c)+(1+eta_c));
Lu_r(i,:) = cp.*T1.*(tau.*eta_t.*(1-1./beta_c.^((gamma-1)./gamma))- ...
    (beta_c.^((gamma-1./gamma)-1)./eta_c));

figure(1)
pl = plot(beta_c, eta_th_id(i,:),'DisplayName', ...
    ['Ideal, T3 = ',num2str(T3),'K']);
hold all
col = pl.Color;
plot(beta_c,eta_th_r(i,:),'Color',col,'LineStyle','--','DisplayName', ...
    ['Real, T3 = ',num2str(T3),'K']);
ylim([0 1])
grid on
grid minor
xlabel('\beta_c','Interpreter','Tex');
ylabel('\eta_t_h)','Interpreter','Tex');
legend show

figure(2)
pl = plot(beta_c, Lu_id(i,:)./cp./T1,'DisplayName', ...
    ['Ideale, T3 = ',num2str(T3),'K']);
hold all
col = pl.Color;
plot(beta_c,Lu_r(i,:)./cp./T1,'Color',col,'LineStyle','--', ...
    'DisplayName',['Reale, T3 = ',num2str(T3),'K']);
grid on
grid minor
xlabel('\beta_c','Interpreter','Tex');
ylabel('$\frac(L_u)()c_p T_1$','Interpreter','Latex');
ylim([0 2])

end
legend show


