% sanity.m %

constants;

for m = 1 : 6,
  Epso{m} = X0{1}*0; 
end
% Epso{2}(mesh(2)/4,mesh(1)*3/4) = 1; %
Epso{2}(1,1) = 1;
for m = 1 : 6,
  EPSo{m} = fft2(Epso{m}) * voxel0;
end;
TUT = X0{1}*0; 
for m = 1 : mesh(2),
  for n = 1 : mesh(1),
    phase = K{4}(m,n) * ( K{1}(m,n) * X0{1} + K{2}(m,n) * X0{2} );
    TUT(m,n) = TUT(m,n) + sum(sum(exp(-j*phase).*Epso{2}))*voxel0;
  end;
end;
EPSo{2}
TUT 

Tut = X0{1}*0; 
for m = 1 : mesh(2),
  for n = 1 : mesh(1),
    phase = K{4}(m,n) * ( K{1}(m,n) * X0{1} + K{2}(m,n) * X0{2} );
    Tut = Tut + TUT(m,n)*exp(j*phase)/volume0;
  end;
end;
Epso{2}
Tut

SIGMAo = IsotropicMul(lambdao, muo, EPSo);

F{1} = SIGMAo{1} .* K{1} +  SIGMAo{6} .* K{2};
F{2} = SIGMAo{6} .* K{1} +  SIGMAo{2} .* K{2};
F{3} = SIGMAo{5} .* K{1} +  SIGMAo{4} .* K{2};
G = K{1}.*F{1} + K{2}.*F{2} + K{3}.*F{3};

EPS{1} = (F{1}.*K{1} - alphao*G.*KK{1})/2/muo;
EPS{2} = (F{2}.*K{2} - alphao*G.*KK{2})/2/muo;
EPS{3} = (F{3}.*K{3} - alphao*G.*KK{3})/2/muo;
EPS{4} = (F{2}.*K{3} + F{3}.*K{2} - 2*alphao*G.*KK{4})/4/muo;
EPS{5} = (F{1}.*K{3} + F{3}.*K{1} - 2*alphao*G.*KK{5})/4/muo;
EPS{6} = (F{1}.*K{2} + F{2}.*K{1} - 2*alphao*G.*KK{6})/4/muo;

EPS_analytic{1} = X0{1}*0; 
EPS_analytic{2} = K{2}.*K{2}*voxel0;
EPS_analytic{3} = X0{1}*0; 
EPS_analytic{4} = K{2}.*K{3}/2*voxel0;
EPS_analytic{5} = X0{1}*0; 
EPS_analytic{6} = K{1}.*K{2}/2*voxel0;

EPS{1} - EPS_analytic{1}
EPS{2} - EPS_analytic{2}
EPS{3} - EPS_analytic{3}
EPS{4} - EPS_analytic{4}
EPS{5} - EPS_analytic{5}
EPS{6} - EPS_analytic{6}

for i = 1 : 6,
  Eps{i} = real(ifft2(EPS{i})) / voxel0 + EpsAvg(i);
  VirtualStrain{i} = Eps{i} - Epso{i};
end;
Sigma = IsotropicMul(lambdao, muo, VirtualStrain);

% check stress self equilibrium %
XX=Sigma{1}; YY=Sigma{2}; XY=Sigma{6};
% R2=(X0{1}+1).^2+(X0{2}+1).^2; YY=(X0{2}+1)./R2; XY=(X0{1}+1)./R2; %
divergence = ...
    (circshift(YY,[-1 0])-circshift(YY,[1 0]))*mesh(2)+...
    (circshift(XY,[0 -1])-circshift(XY,[0 1]))*mesh(1);
fprintf (1, '%e\n', norm(divergence,'fro') / norm(YY,'fro') )
