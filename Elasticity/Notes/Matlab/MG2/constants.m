% constants.m %
% http://mt.seas.upenn.edu/Stuff/e/main.pdf %
% Wang, Jin, Khachaturyan, J. Appl. Phys. 92 (2002) 1351 %

clear;
global H0 volume0 mesh voxel0 S X0 K KS KK EpsAvg X
global hole
global Mu Nu Lambda DLambda DMu
global Sigma 
global Eps Eps0
global muo nuo lambdao alphao betao Epso
clear;

% Vitreloy1: Zr: 41.2 Be: 22.5 Ti: 13.8 Cu: 12.5 Ni: 10 in GPa %
E = 92.1;
nu = 0.387;

% Vitreloy105: Zr: 52.5 Ti: 5 Cu: 17.9 Ni: 14.6 Al:10 in GPa %
% Applied Physics Letters 81 (2002) 4739 %
E = 88.6;
nu = 0.371;

% for 2D %
nu = nu*2;
% 2D isotropic medium - not 3D %
mu = E/(2+2*nu);
lambda = E*nu/(1-nu^2);
alpha = (1+nu)/2;
beta  = nu;

mesh = [256 128];
% mesh = [128 64]; %
% mesh = [64 32]; %

VOXEL_HEIGHT_IN_NM = 2;
% supercell: in nm %
H0 = [(1+sqrt(5))/2 0 
      0 1] * mesh(2) * VOXEL_HEIGHT_IN_NM;
G0 = 2*pi*inv(H0)';

volume0 = abs(det(H0));
voxel0 = volume0;
for i = 1 : length(mesh)
  s{i} = 0 : 1/mesh(i) : 1-1/mesh(i)/2;
  voxel0 = voxel0 / mesh(i);
end;

[S{1},S{2}] = meshgrid(s{1},s{2});
X0{1} = S{1} * H0(1,1) + S{2} * H0(2,1);
X0{2} = S{1} * H0(1,2) + S{2} * H0(2,2);
K{1}  = (S{1}-round(S{1}))*mesh(1)*G0(1,1)+(S{2}-round(S{2}))*mesh(2)*G0(2,1);
K{2}  = (S{1}-round(S{1}))*mesh(1)*G0(1,2)+(S{2}-round(S{2}))*mesh(2)*G0(2,2);
K{3}  = sqrt(K{1}.^2 + K{2}.^2);
% maxK = max(max(K{3})) * 0.5; %
maxK = max(max(K{3}));
KS = (K{3} <= maxK);
K{3}(1,1) = 1;
K{1}  = K{1} ./ K{3};
K{1}(1,1) = 0;
K{2}  = K{2} ./ K{3};
K{2}(1,1) = 0;
% Voight: xx yy xy %
KK{1} = K{1} .* K{1};
KK{2} = K{2} .* K{2};
KK{3} = K{1} .* K{2};

hole = [
%     0.5 0.5   0.5 0.1  0.1 %
%     0.5 0.5   0.5 0.05  0.001 %
%     0.5 0.5   0.5 0.04  0.001 %
%     0.41 0.41 0.5 0.1  0.1 %
       ];

Nu = ones(mesh(2),mesh(1)) * nu;
Mu = ones(mesh(2),mesh(1)) * mu;
for n = 1 : size(hole,1),
  DS = S{1} - hole(n,1);
  DS1 = DS - round(DS);
  DS = S{2} - hole(n,2);
  DS2 = DS - round(DS);
  DX1 = DS1 * H0(1,1) + DS2 * H0(2,1);
  DX2 = DS1 * H0(1,2) + DS2 * H0(2,2);
  DR2 = DX1.^2 + DX2.^2;
  factor1to0 = 0.5 + tanh((DR2/hole(n,4)^2-1)/hole(n,5))/2;
  Mu = Mu .* factor1to0;
end;

figure(1); clf;
pcolor(X0{1}, X0{2}, Mu); shading flat; axis equal;
% 2D isotropic medium - not 3D %
Lambda = 2*Nu./(1-Nu) .* Mu;

muo = mu;
nuo = nu;
lambdao = lambda;
alphao = alpha;
betao = beta;
DLambda = lambdao - Lambda;
DMu     = muo     - Mu;

Eps0{1} = X0{1}*0;
Eps0{2} = X0{1}*0;
Eps0{3} = X0{1}*0;

EPSo = Eps0;
Eps = Eps0;
Epso = Eps0;


units;
T = 300;
kT = T * BOLZ * J_IN_EV;
trialfreq = voxel0 * 1e13;
Q0amp = 5;               % eV
OmegaF = voxel0 / 10;    % Magnitude of transformation volume

% strainrate = 10^2; %
% strainrate = 1; %
% strainrate = 1e-2; %
strainrate = 1e-4;
% strainrate = 1e-6; %
% strainrate = 1e-8; %
% 1/s %
