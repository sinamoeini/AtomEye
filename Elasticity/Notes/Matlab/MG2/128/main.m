% http://mt.seas.upenn.edu/Stuff/e/main.pdf %
% Wang, Jin, Khachaturyan, J. Appl. Phys. 92 (2002) 1351 %

constants; 

step = 0; time = 0; printperiod = 5; npic = 0; more off;
strainincrement = 1e-4; 
Modes = ceil(voxel0 * 4);
% How much elevation, tilt, and roughness of energy landscape in strain space %
RandomRatio = 0.2;
Q0 = Q0amp * (1-RandomRatio + RandomRatio*rand(mesh(2),mesh(1),Modes));

% pure shear transformations %
DEps0{1} = randn(mesh(2),mesh(1),Modes)/2*OmegaF/voxel0;
DEps0{2} = -DEps0{1};
DEps0{3} = randn(mesh(2),mesh(1),Modes)/2*OmegaF/voxel0;
% diffusional relaxation time: Acta Materialia 57 (2009) 881­892 %
tTemp = 1 / (trialfreq * exp(-0.37 / kT));
% time of last transformation %
tlast = Eps0{1} * 0 - 100*tTemp;
unix('rm -rf Jpg/ Eps/; mkdir Jpg/ Eps/');
strainhistory = [0]; stresshistory = [0];

% supercell strain %
EpsAvg = [0 0 0];
% residual stress generator: satisfies self equilibration %
ResidualStrainSTD = OmegaF/voxel0 * 0.1;
Eps0{1} = randn(mesh(2),mesh(1)) * ResidualStrainSTD / 2; 
Eps0{2} = -Eps0{1};
Eps0{3} = randn(mesh(2),mesh(1)) * ResidualStrainSTD / 2; 
Eps0{1} = Eps0{1} - mean(mean(Eps0{1}));
Eps0{2} = Eps0{2} - mean(mean(Eps0{2}));
Eps0{3} = Eps0{3} - mean(mean(Eps0{3}));
HomogeneousStressSolver;
ResidualSigma = Sigma;
% Damage trackers %
PermSoft = 10;
TempSoft = 30;
PermSoftening = X0{1}*0;
PermSofteningCap = -log(0.8);
TempSoftening = X0{1}*0;

EpsAvg = 0.01*[1 -nu 0];
% incremental strain in reference to the present %
Eps0{1} = X0{1}*0;
Eps0{2} = X0{1}*0;
Eps0{3} = X0{1}*0;
% incremental stress distribution Sigma and supercell Stress %
HomogeneousStressSolver;

while (1)
  % softening - mimicking free volume creation and annihilation %
  Softening = PermSoftening + TempSoftening .* exp((tlast-time)/tTemp);
  SIGMA{1} = Sigma{1} + ResidualSigma{1};
  SIGMA{2} = Sigma{2} + ResidualSigma{2};
  SIGMA{3} = Sigma{3} + ResidualSigma{3};
  for m = 1 : Modes,
    % bias due to stress distribution %
    DWork = ( SIGMA{1} .* DEps0{1}(:,:,m) + ...
              SIGMA{2} .* DEps0{2}(:,:,m) + ...
          2 * SIGMA{3} .* DEps0{3}(:,:,m) ) * voxel0 * 1e-18 * J_IN_EV;
    Q(:,:,m) = Q0(:,:,m) .* exp(-Softening) - DWork/2;
  end
  
  if ( min(min(min(Q))) < 0 )
    % athermal plasticity %
    [Qmin,m] = min(Q,[],3);
    [Row,Col] = find(Qmin<0);
%     break %
    for p = 1 : length(Row)
      row = Row(p); col = Col(p); mode = m(row,col);
      VoxelShearTransform;
    end
    HomogeneousStressSolver;
    EpsAvg = EpsAvg - isotropicinv([0 Stress(2) Stress(3)], lambda, mu);
  else
    q = reshape(Q, mesh(2)*mesh(1)*Modes, 1);
    rate = trialfreq * exp( - q / kT );
    timeincrement = 1 / sum(rate);
    if ( timeincrement < strainincrement/strainrate )
      % thermal plasticity %
      dice = rand(1);
      roll = cumsum(rate) / sum(rate);
      idx = min(find( (roll > dice)));
      mode = ceil( idx / mesh(2) / mesh(1) );
      col = ceil( idx / mesh(2) - (mode - 1)*mesh(1) );
      row = idx - (mode-1)*mesh(2)*mesh(1) - (col-1) * mesh(2);
      VoxelShearTransform;
      HomogeneousStressSolver;
      EpsAvg = EpsAvg - isotropicinv([0 Stress(2) Stress(3)], lambda, mu);
      time = time + timeincrement;
      EpsAvg = EpsAvg + strainrate*timeincrement*[1 -nu 0];
      Stress(1) = Stress(1) + E * strainrate*timeincrement;
      Sigma{1}  = Sigma{1}  + E * strainrate*timeincrement;
%       break %
    else
      % pure elasticity %
      time = time + strainincrement/strainrate;
      EpsAvg = EpsAvg + strainincrement*[1 -nu 0];
      Stress(1) = Stress(1) + E * strainincrement;
      Sigma{1}  = Sigma{1}  + E * strainincrement;
    end
  end
  
  if ( EpsAvg(1)-strainhistory(length(strainhistory)) >= ...
       printperiod * strainincrement * 0.99 )
    strainhistory = [strainhistory EpsAvg(1)];
    stresshistory = [stresshistory Stress(1)];
    H = H0 * (eye(2) + T2VoightToFull(EpsAvg));
    X{1} = S{1} * H(1,1) + S{2} * H(2,1);
    X{2} = S{1} * H(1,2) + S{2} * H(2,2);
    for i = 1 : 2,
      U{i} = ( F{i} - alphao * G .* K{i} ) / muo / complex(0,1) ./ K{3};
      u{i} = real(ifft2(EPS{i})) / voxel0;
      X{i} = X{i} + u{i};
    end
    fprintf (1, 'step=%d, time=%e s, strain = %g\n', step, time, EpsAvg(1));
    fprintf (1,'stress = [%g %g %g] GPa\n', Stress(1),Stress(2),Stress(3)); 
    npic = npic + 1;
    figure(1); clf;
    pcolor(X{1}, X{2}, sqrt(Eps0{1}.^2+Eps0{2}.^2+2*Eps0{3}.^2)); 
    text(H(1,1)*0.18, -H(2,2)*0.1, ...
         sprintf ('t = %.2e sec: %g \\times %g nm\n', time,H(2,2),H(1,1)));
    shading flat; axis equal;  axis off; drawnow;
%     title('\epsilon_{\rm transform}'); %
    print(gcf, '-djpeg', sprintf('Jpg/Strain%03d.jpg',npic))
    print(gcf, '-depsc', sprintf('Eps/Strain%03d.eps',npic))
    figure(2); clf; plot(strainhistory, stresshistory, '-o'); 
    xlabel('\epsilon_{11}'); ylabel('\sigma_{11} [GPa]'); drawnow;
    legend(sprintf('T = %d K, strain rate = %.0e/s',T,strainrate),4);
    print(gcf, '-djpeg', sprintf('Jpg/SS%03d.jpg',npic))
    print(gcf, '-depsc', sprintf('Eps/SS%03d.eps',npic))
    figure(3); clf;
    pcolor(X{1}, X{2}, sqrt(Sigma{1}.^2+Sigma{2}.^2+2*Sigma{3}.^2)); 
    shading flat; axis equal;  axis off; drawnow;
%     title('\sigma [GPa]'); %
    print(gcf, '-djpeg', sprintf('Jpg/Stress%03d.jpg',npic))
    print(gcf, '-depsc', sprintf('Eps/Stress%03d.eps',npic))
    unix(sprintf('convert -trim -resize 600 Jpg/Strain%03d.jpg -trim -resize 600 Jpg/Stress%03d.jpg -trim -resize 600 Jpg/SS%03d.jpg -append -resize 600x400! Jpg/%03d.jpg',npic,npic,npic,npic));
    unix(sprintf('rm Jpg/Strain%03d.jpg Jpg/Stress%03d.jpg Jpg/SS%03d.jpg',npic,npic,npic));
  end;
  step = step + 1;  

  if (abs(EpsAvg(1)) > 0.1) 
    break 
  end
  
end

save -ascii EpsAvg.out EpsAvg;
A=Eps0{1}; save -ascii Eps01.out A;
A=Eps0{3}; save -ascii Eps03.out A;
A=[strainhistory' stresshistory']; save -ascii SS.out A;
