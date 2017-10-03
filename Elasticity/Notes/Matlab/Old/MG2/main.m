constants;

units;
kT = 300 * BOLZ * J_IN_EV;
% kT =  77 * BOLZ * J_IN_EV; %
% kT = 600 * BOLZ * J_IN_EV; %
trialfreq = voxel0 * 1e13;
Q0amp = 2;               % eV
OmegaF = voxel0 / 10;    % Magnitude of transformation volume 
GeneralSoftening = 30;
SpecificSoftening = GeneralSoftening * 0;

% 1/s %
% strainrate = 1;   %
strainrate = 1e-4;  
% strainrate = 1e-6;  %
% strainrate = 1e-9; %

step = 0; time = 0; more off;
strainincrement = 1e-4; 
Modes = 12;
RandomRatio = 0.3;
% residual stress effect %
Q0 = Q0amp * (1-RandomRatio + RandomRatio*rand(mesh(2),mesh(1),Modes));
% pure shear transformations %
DEps0{1} = randn(mesh(2),mesh(1),Modes)/2*OmegaF/voxel0;
DEps0{2} = -DEps0{1};
DEps0{3} = randn(mesh(2),mesh(1),Modes)/2*OmegaF/voxel0;
% diffusional relaxation time: Acta Materialia 57 (2009) 881­892 %
trela = 1 / (trialfreq * exp(-0.37 / kT));
% time of last transformation %
tlast = Eps0{1} * 0 - 100*trela;
unix('rm -rf Jpg/; mkdir Jpg/');
strainhistory = [0]; stresshistory = [0];

% supercell strain %
EpsAvg = 0*[1 -nu 0];
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
% incremental strain in reference to the present %
Eps0{1} = X0{1}*0;
Eps0{2} = X0{1}*0;
Eps0{3} = X0{1}*0;
% incremental stress distribution Sigma and supercell Stress %
HomogeneousStressSolver;

while (1)
  for m = 1 : Modes,
    % bias due to stress distribution %
    SIGMA{1} = Sigma{1} + ResidualSigma{1};
    SIGMA{2} = Sigma{2} + ResidualSigma{2};
    SIGMA{3} = Sigma{3} + ResidualSigma{3};
    DWork = ( SIGMA{1} .* DEps0{1}(:,:,m) + ...
              SIGMA{2} .* DEps0{2}(:,:,m) + ...
          2 * SIGMA{3} .* DEps0{3}(:,:,m) ) * voxel0 * 1e-18 * J_IN_EV;
    % softening - mimicking free volume creation %
    Softening = GeneralSoftening * ...
        ( Eps0{1}.^2        + Eps0{2}.^2        + 2*Eps0{3}.^2) +...
                      SpecificSoftening * ...
        ( Eps0{1}.*DEps0{1}(:,:,m) + Eps0{2}.*DEps0{2}(:,:,m) + ...
        2*Eps0{3}.*DEps0{3}(:,:,m) ).^2 ./ ...
        (DEps0{1}(:,:,m).^2 + DEps0{2}(:,:,m).^2 + 2*DEps0{3}(:,:,m).^2);
    Q(:,:,m) = ( Q0(:,:,m) + Q0amp*(1-exp((tlast-time)/trela)) - DWork/2 ) ...
        .* exp(-Softening);
  end
  
  if ( min(min(min(Q))) < 0 )
    % athermal plasticity %
    [Qmin,m] = min(Q,[],3);
    [row,col] = find(Qmin<0);
%     break %
    for p = 1 : length(row)
      for i = 1:3,
        Eps0{i}(row(p),col(p)) = Eps0{i}(row(p),col(p)) + ...
                                DEps0{i}(row(p),col(p),m(row(p),col(p)));
      end
      tlast(row(p),col(p)) = time;
      Q0(row(p),col(p),:) = Q0amp*(1-RandomRatio + RandomRatio*rand(Modes,1));
      DEps0{1}(row(p),col(p),:) = randn(Modes,1)/2*OmegaF/voxel0;
      DEps0{3}(row(p),col(p),:) = randn(Modes,1)/2*OmegaF/voxel0;
      DEps0{2}(row(p),col(p),:) = -DEps0{1}(row(p),col(p),:);
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
      for i = 1:3,
        Eps0{i}(row,col) = Eps0{i}(row,col) + DEps0{i}(row,col,mode);
      end
      tlast(row,col) = time;
      Q0(row,col,:) = Q0amp*(1-RandomRatio + RandomRatio*rand(Modes,1));
      DEps0{i}(row,col,:) = randn(Modes,1)/2*OmegaF/voxel0;
      DEps0{3}(row,col,:) = randn(Modes,1)/2*OmegaF/voxel0;
      DEps0{2}(row,col,:) = -DEps0{1}(row,col,:);
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
  end;
  
  if ( EpsAvg(1)-strainhistory(length(strainhistory)) >= strainincrement*0.99 )
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
    figure(1); clf;
    pcolor(X{1}, X{2}, sqrt(Eps0{1}.^2+Eps0{2}.^2+2*Eps0{3}.^2)); 
    shading flat; axis equal;  axis off; 
%     title('\epsilon_{\rm transform}'); %
    fprintf (1, 'step=%d, time=%e s, strain = %g\n', step, time, EpsAvg(1));
    fprintf (1,'stress = [%g %g %g] GPa\n', Stress(1),Stress(2),Stress(3)); 
    print(gcf, '-djpeg', sprintf('Jpg/Strain%06d.jpg',step))
    figure(2); clf; plot(strainhistory, stresshistory, '-o'); drawnow;
    xlabel('\epsilon_{11}'); ylabel('\sigma_{11} [GPa]'); 
    print(gcf, '-djpeg', sprintf('Jpg/SS%06d.jpg',step))
    figure(3); clf;
    pcolor(X{1}, X{2}, sqrt(Sigma{1}.^2+Sigma{2}.^2+2*Sigma{3}.^2)); 
    shading flat; axis equal;  axis off; 
%     title('\sigma [GPa]'); %
    print(gcf, '-djpeg', sprintf('Jpg/Stress%06d.jpg',step))
    unix(sprintf('convert Jpg/Strain%06d.jpg Jpg/Stress%06d.jpg Jpg/SS%06d.jpg -append Jpg/all%06d.jpg',step,step,step,step));
    unix(sprintf('rm Jpg/Strain%06d.jpg Jpg/Stress%06d.jpg Jpg/SS%06d.jpg',step,step,step));
  end;
  step = step + 1;  

  if (abs(EpsAvg(1)) > 0.3)
    break;
  end;
  
end
