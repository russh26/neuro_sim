%% function DEMO
%  RUNS EXAMPLE OF MILEUSNIC IA AFFERENT MODEL IMPLEMENTATION
%
%  Reference: Mileusnic, M. P., Brown, I. E., Lan, N., & Loeb, G. E. 
%  (2006). Mathematical models of proprioceptors. I. Control and 
%  transduction in the muscle spindle. Journal of neurophysiology, 96(4), 
%  1772-1788. doi:10.1152/jn.00868.2005

function DEMO

% creates model coefficients mat file
CreateInxCoefficients

% sampling frequency
nRate = 1000;
% simulation time
tSim = 3.3;

% temporal profile of muscle length relative to optimal muscle length (L0)
L = [0.95*ones(1,1400),0.95:0.0005:1.08,1.08*ones(1,1639)];

% dynamic and static gamma input
Gamma_dyn = 0;
Gamma_st = 70;

nIaFireRate = MileusnicModel(L,tSim,nRate,Gamma_dyn,Gamma_st);

%% Plotting
time = linspace(0,tSim,numel(nIaFireRate));
% time = 0:0.0010001:tSim;
figure; 
subplot(2,1,1);
plot(time,nIaFireRate,'k','LineWidth',1.5);
ylim([0 450]);
title('Ia Firing Rate');
ylabel('Pulses per second (ppm)');

subplot(2,1,2);
plot(time,L,'k','LineWidth',1.5);
title('Muscle Length');
ylabel('Length (L0)');
xlabel('Time (s)');

disp(num2str(max(nIaFireRate)));
