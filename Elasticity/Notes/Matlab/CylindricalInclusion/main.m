% main.m %
% http://mt.seas.upenn.edu/Stuff/e/main.pdf %
% Wang, Jin, Khachaturyan, J. Appl. Phys. 92 (2002) 1351 %

constants;

step = 0;  printfreq = 1;
while (1)
  EpsoLast = Epso;
  for i = 1 : 6,
    EPSo{i} = fft2(Epso{i}) * voxel0;
    EPSo{i} = EPSo{i} .* KS;
  end;
  SIGMAo = IsotropicMul(lambdao, muo, EPSo);
  F{1} = SIGMAo{1} .* K{1} +  SIGMAo{6} .* K{2};
  F{2} = SIGMAo{6} .* K{1} +  SIGMAo{2} .* K{2};
  F{3} = SIGMAo{5} .* K{1} +  SIGMAo{4} .* K{2};
  G = K{1}.*F{1} + K{2}.*F{2} + K{3}.*F{3};
  EPS{1} = (F{1}.*K{1} - alphao*G.*KK{1})/muo;
  EPS{2} = (F{2}.*K{2} - alphao*G.*KK{2})/muo;
  EPS{3} = (F{3}.*K{3} - alphao*G.*KK{3})/muo;
  EPS{4} = (F{2}.*K{3} + F{3}.*K{2} - 2*alphao*G.*KK{4})/2/muo;
  EPS{5} = (F{1}.*K{3} + F{3}.*K{1} - 2*alphao*G.*KK{5})/2/muo;
  EPS{6} = (F{1}.*K{2} + F{2}.*K{1} - 2*alphao*G.*KK{6})/2/muo;
  for i = 1 : 6,
    Eps{i} = real(ifft2(EPS{i})) / voxel0 + EpsAvg(i);
  end;
  
  Tau0 = IsotropicMul(Lambda,   Mu, Eps0);
  Tau1 = IsotropicMul(DLambda, DMu, Eps);
  for i = 1 : 6,
    Tau{i} = Tau0{i} + Tau1{i};
  end;
  Epso = IsotropicInv(Tau, lambdao, muo);
  maxchange = MaxChange(EpsoLast, Epso);
  if mod(step, printfreq) == 0
    fprintf (1, 'step=%d, max strain change = %e\n', step, maxchange);
  end;
  step = step + 1;
  if ( maxchange < 1e-5 ) 
    break;
  end
end

% break; %

for i = 1 : 6,
  TrueStrain{i}    = Eps{i} - Eps0{i};
end;
Sigma = IsotropicMul(Lambda,  Mu,  TrueStrain);

figure(1); clf;
pcolor(X{1}, X{2}, Sigma{2}); shading flat; axis equal;

figure(2); clf;
for m = [ 
%     round(hole(1,2)*mesh(2))+1  %
%     round((hole(1,2)+0.675*hole(1,4)/H(2,2))*mesh(2))+1  %
%     round((hole(1,2)+ 1.35*hole(1,4)/H(2,2))*mesh(2))+1  %
    round((hole(1,2)+ 5*hole(1,4)/H(2,2))*mesh(2))+1 
        ]
  x = X0{1}(m,:) - hole(1,1)*H(1,1);
  y = X0{2}(m,:) - hole(1,2)*H(2,2);
  sigmaxy = mu * strain * hole(1,4)^2 * (x.^4 - 6*x.^2.*y.^2 + y.^4 ) ./ ...
            (x.^2+y.^2).^3 / (1-nu) .* (x.^2 + y.^2 > hole(1,4)^2);
  plot(X0{1}(m,:), sigmaxy, 'k', X0{1}(m,:), Sigma{6}(m,:), 'ro');
  xlabel('x'); ylabel(sprintf('\\sigma_{xy}(x,y=%g)',y(1)));
  title(sprintf('H_{11}=%g, H_{22}=%g, R=%g, \\epsilon_{xy}^0=%g, mesh_1=%d, mesh_2=%d', H(1,1), H(2,2), hole(1,4), strain, mesh(1), mesh(2)) );
  hold on;
end

legend('Analytical', 'Fourier');
print(gcf, '-depsc', 'CylindricalInclusion.eps');
