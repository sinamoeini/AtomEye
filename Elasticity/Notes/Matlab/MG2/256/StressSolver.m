iter = 0;  
while (1)
  EpsoLast = Epso;
  for i = 1:3,
    EPSo{i} = fft2(Epso{i}) * voxel0;
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
  end;
  Tau0 = IsotropicMul(Lambda,   Mu, Eps0);
  Tau1 = IsotropicMul(DLambda, DMu, Eps);
  for i = 1 : 3,
    Tau{i} = Tau0{i} + Tau1{i};
  end;
  Epso = IsotropicInv(Tau, lambdao, muo);
  maxchange = MaxChange(EpsoLast, Epso);
  iter = iter + 1;
  if ( maxchange < 1e-2 * norm(EpsAvg) ) 
    break;
  end
end

for i = 1 : 3,
  TrueStrain{i} = Eps{i} - Eps0{i};
end;
Sigma = IsotropicMul(Lambda,  Mu,  TrueStrain);

X{1} = S{1} * H(1,1) + S{2} * H(2,1);
X{2} = S{1} * H(1,2) + S{2} * H(2,2);
for i = 1 : 2,
  U{i} = ( F{i} - alphao * G .* K{i} ) / muo / complex(0,1) ./ K{3};
  u{i} = real(ifft2(EPS{i})) / voxel0;
  X{i} = X{i} + u{i};
end

Stress = [ mean(mean(Sigma{1})) mean(mean(Sigma{2})) mean(mean(Sigma{3})) ];
