% Analysis.m %

constants; 
load EpsAvg.out;
load Eps01.out; Eps0{1}=Eps01;
Eps0{2}=-Eps0{1};
load Eps03.out; Eps0{3}=Eps03;
HomogeneousStressSolver;

H = H0 * (eye(2) + T2VoightToFull(EpsAvg));
X{1} = S{1} * H(1,1) + S{2} * H(2,1);
X{2} = S{1} * H(1,2) + S{2} * H(2,2);
for i = 1 : 2,
  U{i} = ( F{i} - alphao * G .* K{i} ) / muo / complex(0,1) ./ K{3};
  u{i} = real(ifft2(EPS{i})) / voxel0;
  X{i} = X{i} + u{i};
end

% figure(1); clf; %
% pcolor(X{1}, X{2}, sqrt(Eps0{1}.^2+Eps0{2}.^2+2*Eps0{3}.^2)); %
% shading flat; axis equal;  axis off; drawnow; %

% load SS.out; %
% figure(2); clf; plot(SS(:,1), SS(:,2), '-o');  %
% xlabel('\epsilon_{11}'); ylabel('\sigma_{11} [GPa]'); drawnow; %

% figure(3); clf; %
% pcolor(X{1}, X{2}, sqrt(Sigma{1}.^2+Sigma{2}.^2+2*Sigma{3}.^2));  %
% shading flat; axis equal;  axis off; drawnow; %

Mises = reshape( sqrt( Eps0{1}.^2 + Eps0{2}.^2 + 2*Eps0{3}.^2 ), ...
                 mesh(2)*mesh(1), 1);
MeshMises = 100;
MisesMin = 0;
MisesMax = max(Mises);
MisesDel = (MisesMax - MisesMin) / MeshMises;
edges = MisesMin + ( 0 : MeshMises )' * MisesDel;
Edges = edges; Edges(1) = Edges(1)-10*eps; 
Edges(MeshMises+1) = Edges(MeshMises+1)+10*eps;

den = histc(Mises,Edges)/length(Mises)/MisesDel;
denmax = max(den(2:MeshMises+1));
plot(0,0,'ro', edges(2:MeshMises+1), den(2:MeshMises+1), 'd');
axis([MisesMin MisesMax 0 denmax]);
xlabel('net von Mises transformation strain of a voxel');
ylabel('Probability density');
title(sprintf('%d \\times %d nm: T = %d K, strain rate = %.0e/s, strain = %.2g', round(H0(2,2)), round(H0(1,1)),T,strainrate,EpsAvg(1)))
printepspdf(sprintf('den%d',round(H0(2,2))));
