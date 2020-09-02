 % Approximate runtime: 30 minutes
% This function performs the optimisation process to minimise process cost.
% To run, simply key in 'cost_optimisation_main' in the command window. 

%% Expected Results: 
% 5 years = [4.74146229918933,2.78935179056159,5.57870358112163,4.21763300131390,0.228230508256410]
% Total Cost for 5 years excl feed material cost = [7126813.29352452]
% Total Annualised Cost = [16500867.6080469]

%% Main Code Starts Below
XGuess = [5 2.5 6.2 5 0.23];
A = []; B = []; Aeq = []; Beq = []; LB = [4 2 4 3 0.2]; UB = [];
options = optimoptions('fmincon','Algorithm','sqp');
output = fmincon(@(X)costobj(X),XGuess,A,B,Aeq,Beq,LB,UB,@nonlcon,options);
TotCost = costobj(output); %Optimised Cost, Excluding Hydrogen Feed Cost

%Cost of Feed
QG = output(4);
PricePerKgUSD = 4.50; %(NREL, 2009)
MoleFlowH2 = QG*101325/(8.314*(313.15));
MassFlowH2 = MoleFlowH2/2000;
AnnualH2Cost = 1.4*PricePerKgUSD*MassFlowH2*3600*8000;
ProcessWaterCost = 184431.542; % Calculated Separately
TotalAnnualisedCost = TotCost/5 + AnnualH2Cost + ProcessWaterCost;

function TotCost = costobj(X)
PAtm = X(1); D = X(2); L = X(3); QG = X(4); QL = X(5); 
TCel = 40;

%% Calculation of Operating Costs
OperationYears = 5; %Years of Assumed Operation
% Energy Required for Compression
QGinlet = QG*1.01325/2.6;
QGoutlet = QG/PAtm;
if (QGinlet-QGoutlet)>0
    CompressorPowerOut = CompressorPowerRequired(QGinlet,PAtm); %Interpolated from HYSYS
else
    CompressorPowerOut = 0;
end
MWhCost = 94.3; %5-year Average Cost of electricity (Wholesale, USEP) per MWh, Singapore
kWhCost = MWhCost/1000; %Cost of electricity per kWh, Singapore
Efficiency = 0.85; %Assume Pump Efficiency
CompressorPowerIn = CompressorPowerOut./Efficiency;

% Energy Required for Pumping Feed
ExchangerPdrop = (10/14.7);
PumpPowerOut = QL*101325*(PAtm+ExchangerPdrop-(1.2/1.01325)); %Units in Watts, Isothermal Assumption
PumpPowerIn = PumpPowerOut./Efficiency;

% Energy Required for Pumping Product
ProductPumpOut = QL*101325*(4-2.5); %Units in Watts, Isothermal Assumption
ProductPumpIn = ProductPumpOut./Efficiency;
% Hourly Cost of Compression and Pumping
HourlyCost = kWhCost*(CompressorPowerIn + PumpPowerIn + ProductPumpIn)/1000;

%% Calculation of Capital Costs
% Cost of Reaction Vessel
S = 17100; %MAWP (stress), units in PSI
Ppsig = PAtm.*14.6959 - 14.6959; %Operating pressure in psig
Pd = exp(0.60608+0.91615.*(log(Ppsig))) + 0.0015655.*((log(Ppsig)).^2); %Design Pressure, psig
Di = D./0.3048; %Internal Reactor Diameter in feet
Lft = L./0.3048; %Reactor Length in feet
Lft = Lft*1.2; %Add 20% to account for column internals
E = 1; %fractional weld efficiency
tp = Pd.*Di.*12./(2.*S.*E-1.2.*Pd); %Pre-adjusted shell wall thickness
tw = 0.22.*(Di.*12+18).*((Lft.*12).^2)./(S.*(Di.^2)); %additional thickness for wind, as vertical vessel
tv = 0.5.*tw+tp; %Thickness after adjusting for wind protection
tc = 0.125; %corrosion allowance
ts = tc + tv; %Final shell thickness
rho = 503.8; %in lb/ft3
W = pi.*(Di+ts./12).*(Lft+0.8.*Di).*(ts./12).*rho;
MatCostUSD = 1.64.*W; %US$1.64 per kg of 316 Steel

% Cost of Internals (3 x Sieve Trays for Liquid Redistribution)
CEPCI_2019 = 610.5;
CEPCI_2010 = 532.9;
CostPerTray2010 = 130 + 440*(D^1.8);
TotalInternalsCostUSD = 3*(CEPCI_2019/CEPCI_2010)*CostPerTray2010;

% Cost of Compressor (Centrifugal)
CompressorkW = CompressorPowerOut/1000;
CompressorCost2010 = 580000 + 20000*(CompressorkW^0.6); % From Towler
CompressorCostUSD = (CEPCI_2019/CEPCI_2010)*CompressorCost2010;

% Cost of Pump (Centrifugal) - 2 required
FlowRateLitres = QL*1000;
FlowPerPump = FlowRateLitres/2;
CostPerPump2010 = 8000 + 240*(FlowPerPump^0.9); % From Towler
TotalPumpCostUSD = 2*(CEPCI_2019/CEPCI_2010)*CostPerPump2010; % Cost of 2 pumps

% Cost of Heat Exchanger
hT = 1000; %Heat transfer coefficient, (W/m2K) from Smith, 2005
LiquidMoleFlow = QL*LiquidConcentration(TCel); %in mol/s
TempRiseRequired = 10; %Temperature rise required, in K
HeatCapL = LiquidHeatCapacity(TCel); %Liquid Heat Capacity, in J/mol.K
HeatTransferRequired = LiquidMoleFlow*HeatCapL*TempRiseRequired; %Units in Watts
LMTD = 25; %Log mean temperature difference
HeatExchangerArea = HeatTransferRequired/(LMTD*hT);
ExchangerCost2010 = 28000 + 54*(HeatExchangerArea^1.2); % From Towler
ExchangerCostUSD = (CEPCI_2019/CEPCI_2010)*ExchangerCost2010;

% Cost of Catalyst
CostPerKgUSD = 90; %Cost per kg 
CatalystCost = (650.*(pi.*(D.^2)./4).*L).*CostPerKgUSD; %Assume US$58/kg
CapCost = 1.4.*(MatCostUSD+CatalystCost+CompressorCostUSD+TotalPumpCostUSD+TotalInternalsCostUSD+ExchangerCostUSD);
OpCost = HourlyCost.*8000*OperationYears;
TotCost = CapCost+OpCost;
end

function [c, ceq] = nonlcon(X)
PAtm = X(1);
D = X(2);
L = X(3);
QG = X(4); 
QL = X(5); 
[EAQH2MoleFlow,Pdrop,HydrogenationRate,ErrorFlag] = helper_reactorsimulation(PAtm,D,L,QG,QL);
QGmin = (182.1*8.314*(50+273.15))/(101325*PAtm);
c = [-(HydrogenationRate-182.1); -(PAtm-Pdrop-1); -(L/D-2); QGmin-QG; ErrorFlag;abs(EAQH2MoleFlow)-0.001];
ceq = [];
end

function CompressorPowerOut = CompressorPowerRequired(QGinlet,PAtm)
% This function calculates the compressor power required per m3/s of flow,
% for an inlet pressure of 2.6 bar
% CompressorPowerOut in Watts
% mole flow basis = 371.4 kgmol/hr per 1m3/s of flow at 2.6 Bar, 30C
Power = [0 48.96 145.1 225.3 356.5 510.0 594.8 704.8 857.4]*1000;
Pressure = [2.5660 3 4 5 7 10 12 15 20];
PowerPerM3 = spline(Pressure,Power,PAtm);
CompressorPowerOut = PowerPerM3*QGinlet;
end

% Liquid Heat Capacity (Assume 30% Quinones, 35% p-Xylene, 35% Methanol)
function HeatCapL = LiquidHeatCapacity(TCel)
HeatCap = [0.22731 0.23183 0.23651 0.24117 0.24576 0.25054 0.25514 0.25980 0.26450 0.26920 0.27380].*1000; %Units in J/molK
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,HeatCap,1);
HeatCapL = p(1)*TCel + p(2);
end

% Liquid Concentration (Assume 30% Quinones, 35% p-Xylene, 35% Methanol)
function CL = LiquidConcentration(TCel)
CLA = [6943 6888 6831 6774 6717 6658 6599 6539 6477 6415 6352];
TA = [20 30 40 50 60 70 80 90 100 110 120];
% Since this is a straight line, a linear fit will suffice. 
p = polyfit(TA,CLA,1);
CL = p(1)*TCel + p(2);
end