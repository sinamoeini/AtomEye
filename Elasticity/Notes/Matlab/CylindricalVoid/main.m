% main.m %
% http://mt.seas.upenn.edu/Stuff/e/main.pdf %
% Wang, Jin, Khachaturyan, J. Appl. Phys. 92 (2002) 1351 %

constants;

step = 0;  printfreq = 1; more off;
while (1)
  EpsoLast = Epso;
  for i = PlaneStrainComponents,
    EPSo{i} = fft2(Epso{i}) * voxel0;
    EPSo{i} = EPSo{i} .* KS;
%     Epso{i} = real(ifft2(EPSo{i})) / voxel0; %
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
  for i = PlaneStrainComponents,
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
  fprintf (1, '%e\n', norm(divergence,'fro') / norm(YY,'fro') );

  Tau0 = IsotropicMul(Lambda,   Mu, Eps0);
  Tau1 = IsotropicMul(DLambda, DMu, Eps);
  for i = 1 : 6,
    Tau{i} = Tau0{i} + Tau1{i};
  end;

  mix = 1;
  EpsoMix = IsotropicInv(Tau, lambdao, muo);
  for i = PlaneStrainComponents,
    Epso{i} = mix*EpsoMix{i} + (1-mix)*Epso{i};
  end
  
  avgchange = AvgChange(EpsoLast, EpsoMix);
  if mod(step, printfreq) == 0
    fprintf (1, 'step=%d, avg strain change = %e\n', step, avgchange);
  end;
  step = step + 1;
  if ( avgchange < 1e-4 * norm(EpsAvg) ) 
    break;
  end
end

for i = PlaneStrainComponents,
  TrueStrain{i} = Eps{i} - Eps0{i};
end;
Sigma = IsotropicMul(Lambda,  Mu,  TrueStrain);

figure(1); clf;
pcolor(X{1}, X{2}, Sigma{2}); shading flat; axis equal;
pcolor(X{1}, X{2}, Epso{3}); shading flat; axis equal;
pcolor(X{1}, X{2}, Eps{3}); shading flat; axis equal;

for m = [ 
    round(hole(1,2)*mesh(2))+1 
%     round((hole(1,2)+0.675*hole(1,4)/H(2,2))*mesh(2))+1  %
%     round((hole(1,2)+ 1.35*hole(1,4)/H(2,2))*mesh(2))+1  %
%     round((hole(1,2)+ 3*hole(1,4)/H(2,2))*mesh(2))+1  %
        ]
  figure(2); clf;
  x = X0{1}(m,:) - hole(1,1)*H(1,1);
  y = X0{2}(m,:) - hole(1,2)*H(2,2);
  sigmayy = (1 + ( 3 * hole(1,4)^4 * (x.^4 - 6*x.^2 .* y.^2 + y.^4) + ...
                   hole(1,4)^2 * (x.^6+13*x.^4.*y.^2+7*x.^2.*y.^4-5*y.^6) )/...
             2 .* (x.^2+y.^2).^-4) .* (x.^2+y.^2 > hole(1,4)^2);
  plot(X0{1}(m,:), sigmayy, 'k', X0{1}(m,:), Sigma{2}(m,:), 'ro');
  xlabel('x'); 
  ylabel(sprintf('\\sigma_{yy}(x,y=%1.0gR)/\\sigma_{yy}^\\infty', ...
                 y(1)/hole(1,4)));
  title(sprintf('H^0_{11}=%g, H^0_{22}=%g, R=%g, mesh_1=%d, mesh_2=%d', ...
                H0(1,1), H0(2,2), hole(1,4), mesh(1), mesh(2)) );
  legend('Analytical', 'Fourier');
  print(gcf, '-depsc', sprintf('CylindricalVoid%1.0g.eps',y(1)/hole(1,4)));
end
