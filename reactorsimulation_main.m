% This script is the main reactor simulation file, which plots the
% concentration and component molar flow rates. To run, simply enter
% 'reactorsimulation_main' in the command window. The runtime is
% approximately 10 seconds. 

% (PAtm,D,L,QG,QL)
PAtm = 4.7415;
D = 2.79;
L = 5.58;
QG = 4.30;
QL = 0.228;
% PAtm = 10, D = 3, L = 15, QL = 0.24, QG = 5
% Target Hydrogenation Rate = 182.1 mol/s
%% Catalyst Parameters
dp = 0.0019; % Catalyst Diameter (m)
CatalystBulkDensity = 650;

%% Reactor Dimensions
A = (pi/4)*(D^2); % Cross Sectional Area
V = A*L;

%% Composition of Input Streams
% Liquid Stream: 9 wt% EAQ, 21 wt% THEAQ, 35 wt% p-Xylene, 35 wt% 2-Octanol
% Input Mass Fraction (Row 1 = EAQ, 2 = THEAQ, 3 = EAQH2, 4 = THEAQH2, 5 = p-Xylene, 6 = 2-Octanol)
% (Column 1 = Mass Fraction, 2 = MW)
LiquidMassComp = [0.09 236.3; 0.189 240.3;  0.00 238.3; 0.021 242.3; 0.35 106.2; 0.35 130.231];
% Calculation of Mole Fraction
[NumComp,~] = size(LiquidMassComp);
TotalMole = 0;
TotalFraction = 0;
for i = 1:NumComp
    TotalFraction = TotalFraction + LiquidMassComp(i,1);
end
if TotalFraction < 0.999 || TotalFraction > 1.001
    disp('WARNING: Total Mass Fraction Does Not Equal Unity')
end
for i = 1:NumComp
    TotalMole = TotalMole + LiquidMassComp(i,1)/LiquidMassComp(i,2);
end
LiquidMoleComp = zeros(1,NumComp);
for i = 1:NumComp
    LiquidMoleComp(i) = (LiquidMassComp(i,1)/LiquidMassComp(i,2))/TotalMole;
end
MWavg = zeros(1,NumComp);
for i = 1:NumComp
    MWavg(i) = LiquidMoleComp(i).*LiquidMassComp(i,2);
end
MWavg = sum(MWavg);
clear TotalMole
clear TotalFraction

%% Liquid Flow Rate (Based on Stoichiometry and Target Conversion)
%QL = 0.231;

%% Temperature and Pressure
TCel = 40; % Temperature in Celsius
PAtmOriginal = PAtm;

%% Molar Flow Rate

TotalMoleFlow = QL*LiquidDensity(TCel)/(0.001*MWavg);
%EAQ 
EAQMoleFlow = LiquidMoleComp(1)*TotalMoleFlow;
EAQMoleFlowInitial = EAQMoleFlow;
%THEAQ 
THEAQMoleFlow = LiquidMoleComp(2)*TotalMoleFlow;
TotalQuinone = EAQMoleFlow + THEAQMoleFlow;
%EAQH2
EAQH2MoleFlow = LiquidMoleComp(3)*TotalMoleFlow;
%THEAQH2
THEAQH2MoleFlow = LiquidMoleComp(4)*TotalMoleFlow;
%H2
H2MoleFlow = (101325*QG)/(8.314*(TCel+273.15));
AQStarting = THEAQMoleFlow+EAQMoleFlow;
%% Reactor Simulation
nDifferentialElements = 1000;
dL = L/nDifferentialElements; % Creating Differential Element of L
dV = dL*A;

DegradationPdtsFlow = 0;
EAQ = zeros(1,nDifferentialElements);
THEAQ = zeros(1,nDifferentialElements);
EAQH2 = zeros(1,nDifferentialElements);
THEAQH2 = zeros(1,nDifferentialElements);
DegradationPdts = zeros(1,nDifferentialElements);
zcoord = zeros(1,nDifferentialElements);
TemperatureProfile = zeros(1,nDifferentialElements);
PressurePlot = zeros(1,nDifferentialElements);
CriticalFailure = false;
ErrorFlag = -1;
troubleshooter = zeros(1,nDifferentialElements);
CumulativeHeatRemoved = 0;
Iter = 0;

for i = 1 :nDifferentialElements
    PressurePlot(i) = PAtm;
    LiquidConc = LiquidConcentration(TCel);
    EAQConc = (EAQMoleFlow/TotalMoleFlow)*LiquidConc;
    THEAQConc = (THEAQMoleFlow/TotalMoleFlow)*LiquidConc;
    troubleshooter(i) = H2MoleFlow;
    [hT(i), kDeg,dpdz, k, kLaL(i), kSaS(i)] = tricklekinetics(D,QL,QG,dp,TCel,PAtm,EAQConc,THEAQConc);
    %% No of Moles in Each Reaction
    LiquidMoleComp(1) = EAQMoleFlow/TotalMoleFlow;
    LiquidMoleComp(2) = THEAQMoleFlow/TotalMoleFlow;  
    if LiquidMoleComp(1) > 0
        mainReaction1 = dV*PAtm*((k(1)*LiquidMoleComp(1))/(LiquidMoleComp(1)+LiquidMoleComp(2)));
    else
        mainReaction1 = 0;
    end
    if LiquidMoleComp(2) > 0
        mainReaction2 = dV*PAtm*((k(2)*LiquidMoleComp(2))/(LiquidMoleComp(1)+LiquidMoleComp(2)));
    else
        mainReaction2 = 0;
    end
    degReaction1 = kDeg(1)*dV;
    degReaction2 = kDeg(2)*dV;
    degReaction3 = kDeg(3)*dV;
    %% Quinone Mass Balance
    EAQMoleFlowTemp = EAQMoleFlow - mainReaction1 - degReaction1 - degReaction3;
    if EAQMoleFlowTemp > 0
        EAQMoleFlow = EAQMoleFlowTemp;
        EAQH2MoleFlow = EAQH2MoleFlow + mainReaction1;
    elseif EAQMoleFlowTemp <= 0
        mainReaction1 = EAQMoleFlow;
        degReaction1 = 0;
        degReaction3 = 0;
        EAQH2MoleFlow = EAQH2MoleFlow + EAQMoleFlow;
        EAQMoleFlow = 0;
    end
    THEAQMoleFlowTemp = THEAQMoleFlow - mainReaction2 + degReaction1 - degReaction2;
    if THEAQMoleFlowTemp > 0
        THEAQMoleFlow = THEAQMoleFlowTemp;
        THEAQH2MoleFlow = THEAQH2MoleFlow + mainReaction2; 
    elseif THEAQMoleFlowTemp <= 0
        THEAQH2MoleFlow = THEAQH2MoleFlow + THEAQMoleFlow; 
        mainReaction2 = THEAQMoleFlow;
        THEAQMoleFlow = 0;
        degReaction2 = 0;
    end
    clear EAQMoleFlowTemp
    clear THEAQMoleFlowTemp
    %% Execution of Equilibrium Reaction
    MolesTransferred = EAQH2MoleFlow;
    THEAQMoleFlowTemp = THEAQMoleFlow - MolesTransferred;
    if THEAQMoleFlowTemp > 0
        THEAQMoleFlow = THEAQMoleFlowTemp;
        EAQH2MoleFlow = EAQH2MoleFlow - MolesTransferred;
        EAQMoleFlow = EAQMoleFlow + MolesTransferred;
        THEAQH2MoleFlow = THEAQH2MoleFlow + MolesTransferred;
    elseif THEAQMoleFlowTemp <= 0
        MolesTransferred = THEAQMoleFlow;
        THEAQMoleFlow = THEAQMoleFlow - MolesTransferred;
        EAQH2MoleFlow = EAQH2MoleFlow - MolesTransferred;
        EAQMoleFlow = EAQMoleFlow + MolesTransferred;
        THEAQH2MoleFlow = THEAQH2MoleFlow + MolesTransferred;
    end
    clear THEAQMoleFlowTemp
    %% Mass Balance for Remaining Components
    DegradationPdtsFlow = DegradationPdtsFlow + degReaction2 + degReaction3;
    % Note: For hydrogen gas, 1 mole of H2 is used up per mole of main reactions, 2 moles for degradation reactions
    H2MoleFlowTemp = H2MoleFlow - mainReaction1 - mainReaction2 - 2*degReaction1 - 2*degReaction2 - 2*degReaction3;
    if H2MoleFlowTemp > 0
        H2MoleFlow= H2MoleFlowTemp;
    elseif H2MoleFlowTemp <= 0
        break
    end
    QG = H2MoleFlow*8.314*(TCel+273.15)/101325; %Gas Flow Rate at 1 atm
    %% Energy Balance (Adiabatic Operation)
    MolarReactionHeat = 104*1000; %J/mol of reaction for main reactions 1 and 2
    HeatEvolved = (mainReaction1 + mainReaction2)*MolarReactionHeat;
    HeatRemoved = pi*((D/2)^2)*dL*hT(i)*(TCel-30);
    CumulativeHeatRemoved = CumulativeHeatRemoved + HeatRemoved;
    dT = (HeatEvolved-HeatRemoved)/(H2MoleFlow*GasHeatCapacity(TCel)+TotalMoleFlow*LiquidHeatCapacity(TCel));
    TemperatureProfile(i) = TCel;
    TCel = TCel+dT;
    %% Momentum Balance (Pressure Drop)
    PAtm = PAtm - (dpdz/101325)*dL;
    if PAtm < 1
        CriticalFailure = true;
        break
    end
    %% For Plotting of Molar Flow Rates
    EAQFlow(i) = EAQMoleFlow;
    THEAQFlow(i) = THEAQMoleFlow;
    EAQH2Flow (i) = EAQH2MoleFlow;
    THEAQH2Flow (i) = THEAQH2MoleFlow;
    H2MoleFlowProfile(i) = H2MoleFlow;
    DegradationFlow(i) = DegradationPdtsFlow;
    %% For Plotting of Concentration Curves
    EAQ(i) = LiquidConc*EAQMoleFlow/TotalMoleFlow;
    THEAQ(i) = LiquidConc*THEAQMoleFlow/TotalMoleFlow;
    EAQH2(i) = LiquidConc*EAQH2MoleFlow/TotalMoleFlow;
    THEAQH2(i) = LiquidConc*THEAQH2MoleFlow/TotalMoleFlow;
    DegradationPdts(i) = LiquidConc*DegradationPdtsFlow/TotalMoleFlow;
    zcoord(i) = i*dL;
    Iter = Iter + 1;
end
if Iter < nDifferentialElements
    EAQ(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*EAQ(Iter);
    THEAQ(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*THEAQ(Iter);
    EAQH2(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*EAQH2(Iter);
    THEAQH2(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*THEAQH2(Iter);
    DegradationPdts(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*DegradationPdts(Iter);
    EAQFlow(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*EAQFlow(Iter);
    THEAQFlow(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*THEAQFlow(Iter);
    EAQH2Flow(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*EAQH2Flow(Iter);
    THEAQH2Flow(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*THEAQH2Flow(Iter);
    DegradationFlow(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*DegradationFlow(Iter);
    H2MoleFlowProfile(Iter+1:nDifferentialElements) = ones(1,numel(Iter+1:nDifferentialElements)).*H2MoleFlowProfile(Iter);
    zcoord(Iter+1:nDifferentialElements) = dL.*(Iter+1:1:nDifferentialElements);
end
%% Creation of Output
if CriticalFailure ~= true
    MoleFlowPlot = figure;
    ConcPlot = figure;
    % Plotting of Concentration Curves
    figure(ConcPlot)
    plot(zcoord,EAQ,'Linewidth',1.5)
    hold on
    plot(zcoord,THEAQ,'Linewidth',1.5)
    plot(zcoord,EAQH2,'Linewidth',1.5)
    plot(zcoord,THEAQH2,'Linewidth',1.5)
    xlim([0 L])
    set(gca,'FontSize',12)
    xlabel('Reactor Length (m)','FontSize',14,'FontWeight','bold')
    ylabel('Concentration (mol/m3)','FontSize',14,'FontWeight','bold')

    % Plotting of Molar Flow Rates
    figure(MoleFlowPlot)
    plot(zcoord,EAQFlow,'Linewidth',1.5)
    hold on
    plot(zcoord,THEAQFlow,'Linewidth',1.5)
    plot(zcoord,EAQH2Flow,'Linewidth',1.5)
    plot(zcoord,THEAQH2Flow,'Linewidth',1.5)
    plot(zcoord,H2MoleFlowProfile,'Linewidth',1.5)
    xlim([0 L])
    set(gca,'FontSize',12)
    xlabel('Reactor Length (m)','FontSize',14,'FontWeight','bold')
    ylabel('Molar Flow Rate (mol/s)','FontSize',14,'FontWeight','bold')
    
    % Calculation of Process Cooling Water
    WaterMassRequired = CumulativeHeatRemoved/(4184*5); %Water Flow Required for 5C rise in water temp
    WaterRequired = WaterMassRequired/1000; %To calculate Volumetric Flow Rate
    
    % To print output
    HydrogenationRate = EAQH2MoleFlow+THEAQH2MoleFlow;
    HydrogenationRateStr = strcat('Hydrogenation Rate: ', string(HydrogenationRate), ' mol/s');
    AQEnding = THEAQMoleFlow+EAQMoleFlow;
    Conversion = ((AQStarting-AQEnding)/AQStarting)*100;
    ConversionStr = strcat('Conversion: ', string(Conversion), '%');
    EAQH2Final = strcat('Outlet Flow Rate of EAQH2: ', string(EAQH2MoleFlow), ' mol/s');
    THEAQH2Final = strcat('Outlet Flow Rate of THEAQH2: ', string(THEAQH2MoleFlow), ' mol/s');
    EAQFinal = strcat('Outlet Flow Rate of EAQ: ', string(EAQMoleFlow), ' mol/s');
    THEAQFinal = strcat('Outlet Flow Rate of THEAQ: ', string(THEAQMoleFlow), ' mol/s');
    OutletPressure = strcat('Outlet Pressure: ', string(PAtm), ' atm');
    TotalPressureDrop = strcat('Total Pressure Drop: ', string(PAtmOriginal - PAtm), ' atm');
    TotalDegradationProducts = strcat('Total Degradation Products: ', string(DegradationPdtsFlow), ' mol/s');
    EAQTopUpRequired = EAQMoleFlowInitial - EAQMoleFlow;
    EAQTopUpDisp = strcat('EAQ Top-Up Required: ',string(EAQTopUpRequired), ' mol/s');
    WaterRequiredStr = strcat('Process Water Flow: ', string(WaterRequired), ' m3/s');
    disp(HydrogenationRateStr); disp(ConversionStr); disp(EAQH2Final); disp(THEAQH2Final)
    disp(EAQFinal); disp(THEAQFinal); disp(TotalDegradationProducts); disp(OutletPressure); disp(TotalPressureDrop)
    disp(EAQTopUpDisp); disp(WaterRequiredStr);
elseif CriticalFailure == true
    disp('Critical Failure: Outlet Pressure Less than 1. Reaction Cannot Proceed.')
end

%% Reaction Kinetics and Mass Transfer Correlations
function [hT, kDeg,dpdz, k, kLaL, kSaS] = tricklekinetics(D,QL,QG,dp,TCel,PAtm,EAQConc,THEAQConc)
A = (pi/4)*(D^2); % Cross Sectional Area
UL = QL/A; % Superficial Liquid Velocity
QGActual = QG/PAtm;
UG = QGActual/A; % Superficial Gas Velocity

%% Unit Conversion for Inputs
P = 101325 * PAtm;
TK = 273.15 + TCel;

%% Physical Property Constants
rhoL = LiquidDensity(TCel); % Density of Liquid Phase (kg/m3)
rhoG = P*2/(8.314*TK*1000); % Density of Gas Phase (kg/m3)
L = UL*rhoL; % Superficial Mass Flux of Liquid Phase
G = UG*rhoG; % Superficial Mass Flux of Gas Phase
E = 0.4; % Bed Porosity (unitless)
uL = LiquidViscosity(TCel); % Dynamic Viscosity of Liquid Phase
uH2 = GasViscosity(TCel); % Dynamic Viscosity of Gas Phase
g = 9.81; % Gravitational Constant
sigma = SurfaceTension(TCel); % Surface Tension of Liquid Phase
ac = 6*(1-E)/dp; % Specific Surface Area of Catalyst (m2/m3)
DH2 = 1.05*(10^-6)*exp(-1520/TK); % Diffusivity of Hydrogen Gas in Xylene
%ac = ((1/(1-E))/((4/3)*pi*((dp/2)^3)))*(4*pi*((dp/2)^2)); %Catalyst area m2/m3
HH2 = HenryConstant(TCel); %(Pa/Concentration)

%% Dimensionless Numbers
ReSL = (dp*UL*rhoL)/(uL); % Reynolds Number for Liquid Phase
ReSG = (dp*UG*rhoG)/(uH2); % Reynolds Number for Gas Phase
WeSL = ((UL)^2)/(rhoL*sigma); % Weber Number for Liquid Phase
%% Determination of Flow Regime
phiC = ac/(D^2);
dpdT = dp/D;
spraypulse = (phiC^(-0.1))*(ReSL^0.52)*(ReSG^-0.47)-0.34;
pulsetrickle = (phiC^(-0.2))*(ReSL^0.27)*(ReSG^0.2)*(dpdT^(-0.5))-18;
spraytrickle = (phiC^(-0.1))*(ReSL^-0.25)*(ReSG^0.67)*(dpdT^(-0.5))-52;
spraywavy = (ReSL^0.45)*(ReSG^0.13)-10;
wavytrickle = (phiC^(-0.1))*(ReSL^0.2)*(ReSG^0.3)*(dpdT^(-0.5))-520;
pulsebubble = (phiC^(-0.6))*(ReSL^1.13)*(ReSG^-0.2)*(dpdT^(-0.8))-790;
flowtransitions = [spraypulse pulsetrickle spraytrickle spraywavy wavytrickle pulsebubble];

%% Liquid Phase Mass Transfer
XG = (UG/UL)*((rhoL/rhoG)^0.5);
dk = dp*(16*(E^3))/((9*pi*((1-E)^2))^(1/3));
Sc = uL/(rhoL*DH2);
% Low Interaction Regime
kLaLLow = (DH2*(2*(10^(-4)))/(dk^2))*(((XG^0.25)*(ReSL^0.2)*(WeSL^0.2)*(Sc^0.5)*((ac*dk/(1-E))^0.25))^3.4);
% High Interaction Regime
kLaLHigh = (DH2*0.45/(dk^2))*(((XG^0.5)*(ReSL^0.8)*(WeSL^0.2)*(Sc^0.5)*((ac*dk/(1-E))^0.25))^1.3);
% Transition Flow
kLaLTrans = (DH2*0.091/(dk^2))*(((XG^0.25)*(ReSL^0.2)*(WeSL^0.2)*(Sc^0.3)*((ac*dk/(1-E))^0.25))^3.8);
kLaL = kLaLHigh; % Correlation for high interaction regime is selected, since pulse flow is occurring
%% Pressure Drop Correlations
XL = (G/L)*((rhoL/rhoG)^0.5);
dpLdz = (ReSL/(1-E))*(150+1.75*ReSL/(1-E))*(uL^2)*(((1-E)/E)^3)/(rhoL*(dp^3)); %Ergun Equation for Liquid Phase
dpGdz = (ReSG/(1-E))*(150+1.75*ReSG/(1-E))*(uH2^2)*(((1-E)/E)^3)/(rhoG*(dp^3)); %Ergun Equation for Gas Phase
chi = (dpLdz/dpGdz)^0.5;
dpdz = (dpLdz+dpGdz)*(10^(0.416/(((log10(chi))^2)+0.666)));

%% Liquid to Solid Mass Transfer
% Wetting Factor, Mills and Dudukovic (1981)
FrSL = (UL^2)/(9.81*dp);
nCE = 1.0 - exp(-1.35*(ReSL^0.333)*(FrSL^0.235)*(WeSL^-0.17)*((ac*dp/((1-E)^2))^-0.0425));
Bd = (E*0.66*(chi^0.81))/(1+0.66*(chi^0.81)); % Dynamic Liquid Hold-Up
ReIL = (dp*UL*rhoL)/(uL*Bd); % Reynolds Number for Liquid Phase based on interstitial velocity
ReIG = ((dp*UG*rhoG)/(uH2*(1-E)))*(E/Bd); % Reynolds Number for Gas Phase based on interstitial velocity
GaL = (g*(dp^3)*(rhoL^2)*(E^3))/((uL^2)*((1-E)^3)); % Galilei Number for Liquid Phase based on interstitial velocity
Sc = uL/(rhoL*DH2);
ReL = ReSL*(E/Bd)/(1-E);
ShS = 0.847*(ReL^0.495)*(Sc^(1/3))/nCE;
kS = ShS*(DH2/dp);
kSaS = kS*ac;
%% Intrinsic Kinetics
%EAQ
kAPB1 = (0.003*650*1000*65.58/60)*exp(-17041/(8.314*TK));
%THEAQ
kAPB2 = (0.003*650*1000*2.62/60)*exp(-10629/(8.314*TK));

%% Effective Kinetics
% EAQ
keff1 = ((HH2*(1/kLaL + 1/kSaS)) + 1/kAPB1)^-1;
kA1 = 101325*keff1;
k(1) = kA1;

%THEAQ
keff2 = ((HH2*(1/kLaL + 1/kSaS)) + 1/kAPB2)^-1;
kA2 = 101325*keff2;
k(2) = kA2;

%% Kinetics of Degradation Reactions
% 50C
k1_50 = 2.70 * (10^-2) * (1000/(101325*3600)); % mol/(kg.Pa.s)
k2_50 = 1.27 * (10^-4) * (1000/(101325*3600)); % mol/(kg.Pa.s)
k3_50 = 3.57 * (10^-4) * (1000/(101325*3600)); % mol/(kg.Pa.s)
b1_50 = 80.7/1000; % m3/mol
b2_50 = 26.5/1000; % m3/mol
alpha_50 = 0.98;

% 70C
k1_70 = 3.94 * (10^-2) * (1000/(101325*3600)); % mol/(kg.Pa.s)
k2_70 = 3.62 * (10^-4) * (1000/(101325*3600)); % mol/(kg.s)
k3_70 = 1.05 * (10^-3) * (1000/(101325*3600)); % mol/(kg.Pa.s)
b1_70 = 52.3/1000; % m3/mol
b2_70 = 26.2/1000; % m3/mol
alpha_70 = 1.0;

% Performing Linear Interpolation
pk1 = polyfit([50 70],[k1_50 k1_70],1);
pk2 = polyfit([50 70],[k2_50 k2_70],1);
pk3 = polyfit([50 70],[k3_50 k3_70],1);
pb1 = polyfit([50 70],[b1_50 b1_70],1);
pb2 = polyfit([50 70],[b2_50 b2_70],1);
palpha = polyfit([50 70],[alpha_50 alpha_70],1);

k1 = pk1(1)*TCel + pk1(2);
k2 = pk2(1)*TCel + pk2(2);
k3 = pk3(1)*TCel + pk3(2);
b1 = pb1(1)*TCel + pb1(2);
b2 = pb2(1)*TCel + pb2(2);
alpha = palpha(1)*TCel + palpha(2);
m = 650; % Catalyst Holdup

%% Interpolated Intrinsic Kinetics of Degradation Reactions
kDegIn1 = m*k1*b1*EAQConc*P/((1+b1*EAQConc+b2*THEAQConc)^2); %in mol/(m3.s)
kDegIn2 = m*k2*b2*THEAQConc*P/((1+b1*EAQConc+b2*THEAQConc)^2); %in mol/(m3.s)
kDegIn3 = m*k3*EAQConc^alpha; %in mol/(m3.s)

%% Effective Kinetics of Degradation Reactions
kDeg1 = (HH2*(1/kLaL + 1/kSaS + 1/kDegIn1))^-1;
kDeg2 = (HH2*(1/kLaL + 1/kSaS + 1/kDegIn2))^-1;
kDeg3 = (HH2*(1/kLaL + 1/kSaS + 1/kDegIn3))^-1;
kDeg = [kDeg1 kDeg2 kDeg3];

%% Heat Transfer Correlations: Nusselt Number
kCond = LiquidConductivity(TCel);
aspect = dp/D;
Nu = 2.51*(1-exp(-4.71/aspect))*(ReSL^0.68);
hT = Nu*kCond/dp;
end

%% Functions for Estimation of Physical Properties
% Henrys Constant
function H = HenryConstant(TCel)
% This function gives henry's constant as Concentration/Pressure (mol/m3Pa)
% T is in Celsius, valid from 20 to 80 C. 
% Liquid is 35% p-Xylene, 35% 2-Octanol and 30% EAQ
HA = [51991.17518 50230.27863 48464.29551 46850.46935 45228.43414 43758.05165 42282.51953 40938.1205 39602.26278 38374.56578 37244.89491];
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a curve, a second order polynomial is fitted. 
p = polyfit(TA,HA,2);
H = p(1)*TCel^2 + p(2)*TCel + p(3);
if TCel<min(TA) || TCel>max(TA)
    disp('WARNING: Temperature has exceeded interpolation range for a nonlinear model')
end
end

% Liquid Density
function rhoL = LiquidDensity(TCel)
rho = [931.1 922.5 913.8 905.1 896.2 887.3 878.3 869.2 860.1 850.9 841.7];
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,rho,1);
rhoL = p(1)*TCel + p(2);
end

% Liquid Viscosity
function uL = LiquidViscosity(TCel)
u = [0.6523 0.5275 0.4349 0.3650 0.3111 0.2689 0.2353 0.2083 0.1862 0.1680 0.1528].*0.001;
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a curve, a second order polynomial is fitted. 
p = polyfit(TA,u,2);
uL = p(1)*(TCel^2) + p(2)*TCel + p(3);
if TCel<min(TA) || TCel>max(TA)
    disp('WARNING: Temperature has exceeded interpolation range for a nonlinear model')
end
end

% Gas Viscosity
function uG = GasViscosity(TCel)
u = [0.009046 0.009244 0.009442 0.00964 0.00984 0.01004 0.01025 0.01046 0.01068 0.01091 0.01115].*0.001;
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,u,1);
uG = p(1)*TCel + p(2);
end

% Surface Tension
function sigma = SurfaceTension(TCel)
sigmaA = [0.04787 0.04627 0.04469 0.04311 0.04155 0.04001 0.03848 0.03698 0.03550 0.03404 0.03262];
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,sigmaA,1);
sigma = p(1)*TCel + p(2);
end

% Liquid Concentration (Assume 30% Quinones, 35% p-Xylene, 35% Methanol)
function CL = LiquidConcentration(TCel)
CLA = [6943 6888 6831 6774 6717 6658 6599 6539 6477 6415 6352];
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,CLA,1);
CL = p(1)*TCel + p(2);
end

% Liquid Heat Capacity (Assume 30% Quinones, 35% p-Xylene, 35% Methanol)
function HeatCapL = LiquidHeatCapacity(TCel)
HeatCap = [0.22731 0.23183 0.23651 0.24117 0.24576 0.25054 0.25514 0.25980 0.26450 0.26920 0.27380].*1000; %Units in J/molK
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,HeatCap,1);
HeatCapL = p(1)*TCel + p(2);
end

function HeatCapG = GasHeatCapacity(TCel)
HeatCap = [0.02937 0.02939 0.02942 0.02945 0.02949 0.02954 0.02961 0.02969 0.02981 0.02996 0.03015].*1000; %Units in J/molK
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a curve, a second order polynomial is fitted. 
p = polyfit(TA,HeatCap,2);
HeatCapG = p(1)*(TCel^2) + p(2)*TCel + p(3);
if TCel<min(TA) || TCel>max(TA)
    disp('WARNING: Temperature has exceeded interpolation range for a nonlinear model')
end
end

% Liquid Conductivity
function kCond = LiquidConductivity(TCel)
Conductivity = [0.1452 0.1430 0.1408 0.1385 0.1362 0.1340 0.1317 0.1294 0.1271 0.1249 0.1226]; %Units in W/(m.K)
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,Conductivity,1);
kCond = p(1)*TCel + p(2);
end