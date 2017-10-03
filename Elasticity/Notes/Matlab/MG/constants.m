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

mu = E/(2+2*nu);
lambda = 2*nu/(1-2*nu) * mu;
alpha = 1/2/(1-nu);
beta  = nu/(1-nu);

% mesh = [256 128]; %
% mesh = [128 64]; %
mesh = [64 32];
% supercell: in nm %
H0 = [(1+sqrt(5))/2 0 0
      0 1 0
      0 0 1] * mesh(2) / 2;
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
K{3}  = K{2} * 0;
K{4}  = sqrt(K{1}.^2 + K{2}.^2);
maxK = max(max(K{4})) * 0.1;
% maxK = max(max(K{4})) * 1; %
KS = K{4} < maxK;
K{4}(1,1) = 1;
K{1}  = K{1} ./ K{4};
K{1}(1,1) = 0;
K{2}  = K{2} ./ K{4};
K{2}(1,1) = 0;
% Voight: xx yy zz yz xz xy %
KK{1} = K{1} .* K{1};
KK{2} = K{2} .* K{2};
KK{3} = K{3} .* K{3};
KK{4} = K{2} .* K{3};
KK{5} = K{1} .* K{3};
KK{6} = K{1} .* K{2};

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
Lambda = 2*Nu./(1-2*Nu) .* Mu;

muo = mu;
nuo = nu;
lambdao = lambda;
alphao = alpha;
betao = beta;
DLambda = lambdao - Lambda;
DMu     = muo     - Mu;

NonZeroComponents = [1 2 6];

ResidualStrainSTD = 0;
Eps0{1} = randn(mesh(2),mesh(1)) * ResidualStrainSTD; 
Eps0{2} = -Eps0{1};
Eps0{3} = X0{1}*0; 
Eps0{4} = X0{1}*0; 
Eps0{5} = X0{1}*0; 
Eps0{6} = randn(mesh(2),mesh(1)) * ResidualStrainSTD; 
Eps0{1} = Eps0{1} - mean(mean(Eps0{1}));
Eps0{2} = Eps0{2} - mean(mean(Eps0{2}));
Eps0{3} = Eps0{3} - mean(mean(Eps0{3}));
Eps0{4} = Eps0{4} - mean(mean(Eps0{4}));
Eps0{5} = Eps0{5} - mean(mean(Eps0{5}));
Eps0{6} = Eps0{6} - mean(mean(Eps0{6}));

EPSo = Eps0;
Eps = Eps0;
Epso = Eps0;

