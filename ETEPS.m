clear, clc, close
% -----------------------------------------------------------
% Electrolyzer Performance Analysis
% Authors:
% Jalal Heidari – 10946612
% Maryam Pasandidehnia – 10951180
% MohammadSina GhanbariPakdehi – 10952079
% Seyedmehdi Khalkhali – 10959631
% -----------------------------------------------------------
%% Data
% Reactions and Electrodes
del_E_std = 1.23;                   % [V]
P_bar = 5;                          % [bar] - Operating pressure
P = P_bar/1.0125;                   % [atm] 
P_mmHg = P*760;                     % [mmHg] 
T_rel = 80;                         % [Celsius] - Operating temperature
T_abs = T_rel+273.15;               % [K] 
F = 96500;                          % [C/eq] - Faraday constant
ne = 2;                             % [-] Exchanged Electron
nu = [-1, -0.5, 1];                 % Stochiometric Coefficient [H2 O2 H2O]
ba = 0.083;                         % [V/dec] - Anodic Tafel slope
bc = 0.055;                         % [V/dec] - Cathodic Tafel slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE THESE VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%

i0a = 2.778;                        % [A/m2] - Anodic exchange current density - Acta 3030
% i0a = 0.134535;                   % [A/m2] - Anodic exchange current density - Practical
i0c = 4.769;                        % [A/m2] - Cathodic exchange current density - Acta 4030
% i0c = 0.740568;                   % [A/m2] - Cathodic exchange current density - Practical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE THESE VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%

i_lim = 400;                        % [A/m2] - Limit current density
% Cell Size and Membrane
W = 0.3;                            % [m] - Electrodes width
H = 0.3;                            % [m] - Electrodes height
d = 5e-4;                           % [m] - Membrane thickness
S = H*W;                            % [m2] - Electrode Area
wg = 3.5e-3;                        % [m] - Cell thickness
VC = wg*S;                          % [m3] - Cell volume
sigma = 15.1;                       % [S/m] - Membrane conductivity
% Electrolyte Properties and Calculations
CKOH = 30;                          % [%] - KOH content
m = (CKOH*(183.1221-0.56845*T_abs+984.5679*exp(CKOH/115.96277)))/5610.5; % [M] - KOH molarity
sigma_KOH = -2.04*m-0.0028*(m^2)+0.005332*m*T_abs+207.2*(m/T_abs)+0.001043*(m^3) ... 
    -0.0000003*(m^2)*(T_abs^2);     % [S/m] - Electrolyte conductivity
R_bubble_free = (1/sigma_KOH)*((wg-d)/S);       % [ohm] - Bubble-free electrolyte resistance
P0_H2O = exp(37.047 - (6275.7 / T_abs) - (3.4159 * log(T_abs)));    % [atm] - Equilibrium water pressure
P_H2O = exp(0.016214 - 0.13802 * m + 0.1933 * sqrt(m) + 1.0239 * log(P0_H2O)); % [atm] - Partial pressure of water vapor

%% Calaculations
i_vector = 5:5:400;            % [A/m2] - Current density vector
delE_Nernst = del_E_std-((8.314*T_abs)/(2*F)*log(((P - P_H2O)^1.5) / (P_H2O / P0_H2O)));    % [V] - Nernst Potential
% Preallocation for calculation variables
theta = ones(1,length(i));      % [%] - Bubble Coverage
eta_a = ones(1,length(i));      % [V] - OverPotential Anode
eta_c = ones(1,length(i));      % [V] - OverPotential Cathode
R_bubble_layer = ones(1,length(i)); % [ohm] - Bubble layer electrolyte resistance
sigma_el = ones(1,length(i));   % [S/m] - Electrode conductivity
delV = ones(1,length(i));       % [V] Cell Potential
Power = ones(1,length(i));      % [Watt] Cell Power
V_H2 = ones(1,length(i));       % [Nl/h] Produced H2 
V_O2 = ones(1,length(i));       % [Nl/h] Produced O2 
V_H2O = ones(1,length(i));      % [Nl/h] Required H2O

% Main Loop
for i=1:length(i_vector)
    theta(i) = (-97.25 + 182 * (T_abs / 298) - 84 * (T_abs / 298)^2) ...
               * (i_vector(i) / 300000)^0.3 * (P / (P - P_H2O));
    eta_a(i) = ba*log(i_vector(i)/i0a)-ba*log(1-theta(i));
    eta_c(i) = bc*log(i_vector(i)/i0c)-bc*log(1-theta(i));
    R_bubble_layer(i) = R_bubble_free * ((1 / ((1 - 0.6666 * theta(i)) ^ 1.5)) - 1);
    sigma_el(i) = (wg-d)/((R_bubble_free+R_bubble_layer(i))*S);
    delV(i) = delE_Nernst+eta_a(i)+eta_c(i)+((wg-d)/sigma_el(i))*i_vector(i)+(d/sigma)*i_vector(i);
    Power(i) = delV(i)*S*i_vector(i);
    V_H2(i) = (((i_vector(i)*S)*0.00201588)/(2*F*0.08375))*3600*1000;
    V_O2(i) = (((i_vector(i)*S)*0.0319988)/(4*F*1.354))*3600*1000;
    V_H2O(i) = (((i_vector(i)*S)*0.01801528)/(2*F*997))*3600*1000;
end

%% Graphics
figure(1)
subplot(1,2,1)
plot(i_vector,delV)
xlabel('Current Density [A/m^{2}]'); ylabel('\DeltaV [V]'); title('Cell Voltage versus Current Density') ; grid on

subplot(1,2,2)
plot(i_vector,eta_a)
hold on 
plot(i_vector,eta_c)
xlabel('Current Density [A/m^{2}]'); ylabel('OverPotential [V]'); title('Overpotential versus Current Density') ; grid on
legend('Anode','Cathode')

figure(2)
plot(i_vector,Power)
xlabel('Current Density [A/m^{2}]'); ylabel('Power [W]'); title('Cell Power versus Current Density') ;grid on

figure(3)
plot(i_vector,V_H2)
hold on
plot(i_vector,V_O2)
xlabel('Current Density [A/m^{2}]'); ylabel('Volume Rate [Nl/h]'); grid on
legend('H_{2}','O_{2}')

%% Post Processing
%               Scenario 1 Scenario 2 Scenario 3
Total_gas_flow = [1.54733E-05 1.71033E-05 2.3136E-05];          % [kmol/s]
H2_CO_rate = [0.745543068 0.77100698 1.643711836];              % [-]
H2_rate = [4.05189E-06 3.85476E-06 5.85926E-06];                % [kmol/s]
CO_rate = [5.43482E-06 4.99964E-06 3.56465E-06];                % [kmol/s]
H2_CO_rate_req = 2.1;                                           % [-]
H2_rate_req = CO_rate* H2_CO_rate_req;                          % [kmol/s]
H2_rate_req_el = (H2_rate_req -H2_rate)*3600;                   % [kmol/h]
H2_single_cell_molar = V_H2(end)*((1*T_abs)/(P*298))*(P/(0.0821*T_abs))/1000;                  % [kmol/h]

N_cells = ceil(H2_rate_req_el/H2_single_cell_molar)';            % [-]
Power_cells = (N_cells*Power(end)/1000);                         % [kW]               
Flux_of_req_H2O = (V_H2O(end)*N_cells);                          % [L/h]

%% Output
Scenario = [1 2 3]';
output = table(Scenario, N_cells, Power_cells, Flux_of_req_H2O);
disp(output)