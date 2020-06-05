%% function nIaFireRate = MileusnicModel(L,tSim,nRate,Gamma_dyn,Gamma_st)
%  CALCULATES IA PRIMARY AFFERENT FIRING FREQUENCY USING THE MILEUSNIC
%  MODEL
%
%  INPUTS -----------------------------------------------------------------
%  L            <numeric> muscle length relative to optimal length [1 x n]
%  tSim         <numeric> simulation time in seconds [n]
%  nRate        <numeric> sampling frequency in hertz [n]
%  Gamma_dyn    <numeric> dynamic gamma frequency [n]
%  Gamma_st     <numeric> static gamma frequency [n]
%
%  OUTPUTS ----------------------------------------------------------------
%  nIaFireRate  <numeric> Ia primary afferent firing frequency in hertz 
%                   [1 x n]
%
%  DEPENDENCIES -----------------------------------------------------------
%  TorsionODE   <function> embedded function
%  Coefficients <mat file> mat file of model parameters (see:
%                   CreateInxCoefficients.m)
%
%  REFERENCES -------------------------------------------------------------
%  Mileusnic, M. P., Brown, I. E., Lan, N., & Loeb, G. E. 
%  (2006). Mathematical models of proprioceptors. I. Control and 
%  transduction in the muscle spindle. Journal of neurophysiology, 96(4), 
%  1772-1788. doi:10.1152/jn.00868.2005

function nIaFireRate = MileusnicModel(L,tSim,nRate,Gamma_dyn,Gamma_st)

% fraction of smaller input between either bag1 or bag2 + chain
S = 0.156;

% Bag1
[nIa1_Bag1,~] = TorsionODE(1,L,Gamma_dyn,Gamma_st,nRate,tSim);

% Bag2
[nIa1_Bag2,~] = TorsionODE(2,L,Gamma_dyn,Gamma_st,nRate,tSim);

% Chain
[nIa1_Ch,~] = TorsionODE(3,L,Gamma_dyn,Gamma_st,nRate,tSim);

% Add Bag2 and Chain
nIaBC = nIa1_Bag2 + nIa1_Ch;

% Larger + S*Smaller
for i = 1:numel(nIa1_Bag1)
    if nIa1_Bag1 > nIaBC
        nIaFireRate(i) = nIa1_Bag1(i)+S*nIaBC(i);
    else
        nIaFireRate(i) = nIaBC(i)+S*nIa1_Bag1(i);
    end
end

% prevents negative firing frequencies
nIaFireRate(nIaFireRate<0) = 0;
end

%% CALCULATES AFFERENT POTENTIAL FOR A SINGLE INTRAFUSAL FIBER OR CHAIN)
function [nIa1,nIa2] = TorsionODE(inx,L,Gamma_dyn,Gamma_st,nRate,tSim)
% load model coefficients
Params = load('Coefficients.mat');

%% Relate fusimotor firing frequency to activation (assuming static firing frequency)

% % inx denotes which structure is being solved (1=bag1, 2=bag2, 3=chain)

% relationship of fusimotor frequency to activation
fdyn = Gamma_dyn/(Gamma_dyn+Params.freq(1));
fst = Gamma_st/(Gamma_st+Params.freq(3));

%initial conditions for position and velocity
x(1) =0;
v(1) = 0;

Beta = (Params.Beta0(inx) + Params.Beta1(inx)*fdyn + Params.Beta2(inx)*fst);
Gamma = Params.g1(inx)*fdyn + Params.g2(inx)*fst;

T = linspace(0,tSim,tSim*nRate);

% Calulate dL and ddL and resample to maintain matrix size.
dL = diff(L).*nRate;
dL = resample(dL,numel(L),numel(dL));

ddL = diff(L,2).*nRate^2;
ddL = resample(ddL,numel(L),numel(ddL));

% figure; 
% subplot(3,1,1);
% plot(L);
% title('Length');
% subplot(3,1,2);
% plot(dL);
% title('Change in Length');
% subplot(3,1,3);
% plot(ddL);
% ('Second derivative of length');

for i = 1:numel(T)
    % Coeff of asymmetry in F-V curve during shortening or lengthening
    if dL < 0
        C = Params.cs(inx);
    else
        C = Params.cl(inx);
    end
    
    % x = Torsion
    % v = dx/dt
    % z = dv/dt
    
    odeFunV = @(t,x,v) v;
    
    % equation for tension
    odeFuncZ = @(t,x,v)(Params.ksr(inx)/Params.m(inx)*(C.*Beta.*sign(dL(i)-v/Params.ksr(inx))*...
        abs(dL(i)-v/Params.ksr(inx))^Params.a(inx)*...
        (L(i)-Params.l0sr(inx)-x./Params.ksr(inx)-Params.R(inx))+...
        Params.kpr(inx)*(L(i)-Params.l0sr(inx)-x/Params.ksr(inx)...
        -Params.l0pr(inx))+Params.m(inx)*ddL(i)+Gamma-x));
    
    %% 4th Order Runge-Kutta for 2nd Order Equation
    
    h = tSim/nRate;
    
    k_1 = odeFunV(T(i),x(i),v(i));
    L_1 = odeFuncZ(T(i),x(i),v(i));
    
    k_2 = odeFunV(T(i)+0.5*h,x(i)+0.5*h*k_1,v(i)+0.5*h*L_1);
    L_2 = odeFuncZ(T(i)+0.5*h,x(i)+0.5*h*k_1,v(i)+0.5*h*L_1);
    
    k_3 = odeFunV((T(i)+0.5*h),(x(i)+0.5*h*k_2),(v(i)+0.5*h*L_2));
    L_3 = odeFuncZ((T(i)+0.5*h),(x(i)+0.5*h*k_2),(v(i)+0.5*h*L_2));
    
    k_4 = odeFunV((T(i)+h),(x(i)+k_3*h),(v(i)+L_3*h)); % Corrected        
    L_4 = odeFuncZ((T(i)+h),(x(i)+k_3*h),(v(i)+L_3*h));

    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
    v(i+1) = v(i) + (1/6)*(L_1+2*L_2+2*L_3+L_4)*h;  % main equation
    
    %% Solve afferent potential for each structure (bag1, bag2, or chain)
    if inx ==1
        nIa1(i) = Params.G(inx).*(x(i)./Params.ksr(inx)-(Params.lsr(inx)-Params.l0sr(inx)));
        nIa2(i) = zeros(numel(nIa1,1));
    elseif inx ==2 || inx ==3
        nIa1(i) = Params.G(inx).*(x(i)./Params.ksr(inx)-(Params.lsr(inx)-Params.l0sr(inx)));
        nIa2(i) = Params.G(inx).*(Params.X(inx).*(Params.Lsec(inx)./Params.l0sr(inx).*(Params.lsr(inx)-Params.l0sr(inx))...
            +(1-Params.X(inx)).*(Params.Lsec(inx)/Params.l0pr(inx))*(L(i)-x(i)/Params.ksr(inx)-Params.l0sr(inx)-Params.lpr(inx))));
    end
    
end
% figure; plot(x);
% title('x');
% % ylim([0 100]);
% figure; plot(v);
% title('V');

% nIa1 = Params.G.*(x/Params.ksr-(Params.lsr-Params.l0sr));
end
