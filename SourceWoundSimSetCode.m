% Runs a set of model simulations to test the threshold between simple wound response and tumor formation.

% Susan Massey
% July 2016

clear all 
% close all
% clc

% DECLARE GLOBAL PARAMETERS:
global EC50 Qr ninitx Dp K km DrBar RhoRBar p0 r0 S_max decay

% SET TIME AND SPACE VARIABLES: 
w = 2; % pdepe symmetry parameter; spherical = 2 (cylinder=1, slab = 0)
x = linspace(0,1,375); % cm; 1 cm spherically symmetric space domain vector with 375 values 
%1/375 is in cm, 10/375 gives mm, 100/375 = tenths of millimeters or about 1/4 tenth of a mm
t0=0;  % days; start time 0 days
tf=250; % days; end time 50 days; vs. lifespan = 3 years, so 365 days/year*3years = 1095
tsteps=tf*2+1; %stepping by every 1/2 day
t = linspace(t0,tf,tsteps); % days; time vector with 101 values from t0 to tf
ninitx=.03; %cm; radius of initial sphere of influence of injection 
nt=t;
nx=x;

% DEFINE UNVARIED PARAMETERS:
Dp = 0.0005; % cm^2/day; 
K = (2.3)*10^8; % cells/cm^3; carrying capacity. Computed assuming cell radius = 13.365 microns (ie, volume of 10^4 microns^3)
km = 30; % ng/mL; half-max receptor activation;
EC50=10^0.5; %ng/mL; half max PDGF for downstream effects; Pringle 1989.
DrBar = (2/3)*8.7e-5; %cm2/day; max possible rate of diffusion of r cells
RhoRBar = log(2)/(18/24); % 1/day; max possible rate of proliferation of r cells
Qr=10^(-5.15); % consumption rate or PDGF by r cells

% DEFINE INITIAL CONDITIONS:  
O2a= 250/((4/3)*pi*ninitx^3); %cells/mL = cells/cc = cells/cm^3; background level of r cells
r0 = O2a; % cells/cm^3 = cells/mL
p0 = 0.8145; % ng/mL

% DEFINE STRUCTURE TO STORE DATA
WoundData=struct;

% DEFINE OPTIONS FOR pdepe SOLVER
options=odeset('RelTol',1e-10,'AbsTol',1e-14);

% START TIMER
tic

% BEGIN LOOPING TO DEFINE VARIED PARAMETERS:

% NUMBER OF TIMES TO ITERATE LOOPS
I=10;  % S_MAX LOOP
J=10;  % DECAY LOOP

% INITIATE VECTORS TO STORE VARIED PARAMETER VALUES
Smaxes=zeros(1,I);
decays=zeros(1,J);

% S_MAX LOOP
for i=1:I
    S_max=10*i;
    Smaxes(i)=S_max;
    
% DECAY LOOP
    for j=1:J
        decay=0.01*j;
        decays(j)=decay;
        
% SET INDEX FOR SIMULATION SAVE:
        idx=(i-1)*J+j;

% SOLVE PDES
sol = pdepe(w,@SourceWoundPDEs,@WoundICs,@WoundBCs,nx,nt,options);

% SAVE OUTPUT IN STRUCTURE
WoundData(idx).r = sol(:,:,1); %r; 
WoundData(idx).p = sol(:,:,2); %p; 
% WoundData(idx).t = t; 

% REPORT TIME OF SOLVE
fprintf(['Simulation ',num2str(idx),' of ',num2str(I*J),' complete. \n']);
toc

% END THE LOOPS
    end;
end;

% SAVE WORKSPACE VARIABLES AND STRUCTURE

save SourceWoundSimSet

