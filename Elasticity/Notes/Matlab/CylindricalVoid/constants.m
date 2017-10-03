% constants.m %

clear;
global H0 volume0 mesh voxel0 S X0 K KS KK EpsAvg X
global hole
global Mu Nu Lambda DLambda DMu
global Sigma 
global Eps Eps0
global muo nuo lambdao alphao betao Epso
clear;

mu = 1000;
nu = 0.3;
lambda = 2*nu/(1-2*nu) * mu;
alpha = 1/2/(1-nu);
beta  = nu/(1-nu);

H0 = [1 0 0
      0 1 0
      0 0 1];
G0 = 2*pi*inv(H0)';
% mesh = [1024 1024]; %
% mesh = [512 512]; %
mesh = [256 256];
% mesh = [128 128]; %

volume0 = abs(det(H0));
voxel0 = volume0;
for i = 1 : length(mesh)
  s{i} = 0 : 1/mesh(i) : 1-1/mesh(i)/2;
  voxel0 = voxel0 / mesh(i);
end;

[S{1}, S{2}] = meshgrid(s{1},s{2});
X0{1} = S{1} * H0(1,1) + S{2} * H0(2,1);
X0{2} = S{1} * H0(1,2) + S{2} * H0(2,2);
K{1}  = (S{1}-round(S{1}))*mesh(1)*G0(1,1)+(S{2}-round(S{2}))*mesh(2)*G0(2,1);
K{2}  = (S{1}-round(S{1}))*mesh(1)*G0(1,2)+(S{2}-round(S{2}))*mesh(2)*G0(2,2);
K{3}  = K{2} * 0;
K{4}  = sqrt(K{1}.^2 + K{2}.^2);
maxK = max(max(K{4})) * 0.4;
% maxK = max(max(K{4})); %
KS = (K{4} <= maxK);
K{4}(1,1) = 1;
K{1}  = K{1} ./ K{4};
K{1}(1,1) = 0;
K{2}  = K{2} ./ K{4};
K{2}(1,1) = 0;
% Voight: xx yy  zz yz xz xy %
KK{1} = K{1} .* K{1};
KK{2} = K{2} .* K{2};
KK{3} = K{3} .* K{3};
KK{4} = K{2} .* K{3};
KK{5} = K{1} .* K{3};
KK{6} = K{1} .* K{2};

EpsAvg = 1/2/mu/(1+nu) * [-nu 1 -nu 0 0 0];
H = H0 * (eye(3) + T2VoightToFull(EpsAvg));
X{1} = S{1} * H(1,1) + S{2} * H(2,1);
X{2} = S{1} * H(1,2) + S{2} * H(2,2);

hole = [
%     0.5 0.5   0.5 0.1  0.1 %
%     0.5 0.5   0.5 0.05  0.001 %
    0.5 0.5   0.5 0.02  0.01
%         0.41 0.41 0.5 0.1  0.1 %
       ];

for i = 1 : 6,
  Eps0{i} = X0{1}*0; 
  Epso{i} = Eps0{i};
  EPSo{i} = X0{1}*0;
  Eps{i}  = X0{1}*0 + EpsAvg(i); 
  VirtualStrain{i} = Eps{i} - Epso{i};
  TrueStrain{i} = Eps{i} - Eps0{i};
end;
Mu = ones(mesh(2),mesh(1)) * mu;
Nu = ones(mesh(2),mesh(1)) * nu;
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
pcolor(X{1}, X{2}, Mu); shading flat; axis equal;
Lambda = 2*Nu./(1-2*Nu) .* Mu;

muo = mu;
nuo = nu;
lambdao = lambda;
alphao = alpha;
betao = beta;
DLambda = lambdao - Lambda;
DMu     = muo     - Mu;

% PlaneStrainComponents = [1 2 3 6]; %
PlaneStrainComponents = [1 2 6];

Tau0 = IsotropicMul(Lambda,   Mu, Eps0);
Tau1 = IsotropicMul(DLambda, DMu, Eps);
for i = 1 : 6,
  Tau{i} = Tau0{i} + Tau1{i};
end;
Epso = IsotropicInv(Tau, lambdao, muo);
