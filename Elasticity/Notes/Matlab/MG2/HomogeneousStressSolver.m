for i = 1:3,
  EPSo{i} = fft2(Eps0{i}) * voxel0;
  EPSo{i} = EPSo{i} .* KS;
end;
SIGMAo = IsotropicMul(lambdao, muo, EPSo);
F{1} = SIGMAo{1} .* K{1} +  SIGMAo{3} .* K{2};
F{2} = SIGMAo{3} .* K{1} +  SIGMAo{2} .* K{2};
G = K{1}.*F{1} + K{2}.*F{2};
EPS{1} = (F{1}.*K{1} - alphao*G.*KK{1})/muo;
EPS{2} = (F{2}.*K{2} - alphao*G.*KK{2})/muo;
EPS{3} = (F{1}.*K{2} + F{2}.*K{1} - 2*alphao*G.*KK{3})/2/muo;
for i = 1:3,
  Eps{i} = real(ifft2(EPS{i})) / voxel0 + EpsAvg(i);
  TrueStrain{i} = Eps{i} - Eps0{i};
end;
Sigma = IsotropicMul(Lambda,  Mu,  TrueStrain);
Stress = [ mean(mean(Sigma{1})) mean(mean(Sigma{2})) mean(mean(Sigma{3})) ];
