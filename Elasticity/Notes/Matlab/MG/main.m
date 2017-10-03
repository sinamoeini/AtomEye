constants;

units;
kT = BOLZ * 300;
trialfreq = voxel0 * 1e12;
Q0 = 9 * EV_IN_J;
OmegaF = 1;    % Magnitude of transformation volume 
GeneralSoftening = 4;
SpecificSoftening = GeneralSoftening * 0.3;

strain = 0.015;
% 1/s %
strainrate = 1;  
% strainrate = 1e-4;   %
% strainrate = 1e-6;  %
% strainrate = 1e-9; %

step = 0;  time = 0; printfreq = 1;  more off;
strainincrement = 1e-3; 
Modes = 12;
DEps0{1} = randn(mesh(2),mesh(1),Modes) / 2 * OmegaF / voxel0;
DEps0{6} = randn(mesh(2),mesh(1),Modes) / 2 * OmegaF / voxel0;
DEps0{2} = -DEps0{1};
unix('rm -rf Jpg/; mkdir Jpg/');

strainhistory = [0];
stresshistory = [0];
while (1)
  strainhistory = [strainhistory strain];
  EpsAvg = strain * [1 -nu -nu 0 0 0];
  H = H0 * (eye(3) + T2VoightToFull(EpsAvg));

  StressSolver;
  
  stresshistory = [stresshistory mean(mean(Sigma{1}))];
  mean(mean(Sigma{2}))
  mean(mean(Sigma{3}))
  break
  
  for m = 1 : Modes,
    DWork(:,:,m) = ( Sigma{1} .* DEps0{1}(:,:,m) + ...
                     Sigma{2} .* DEps0{2}(:,:,m) + ...
                 2 * Sigma{6} .* DEps0{6}(:,:,m) ) * voxel0 * 1e-18;
    Softening(:,:,m) = GeneralSoftening * ...
        (Eps0{1}.^2        + Eps0{2}.^2        + 2*Eps0{6}.^2) +...
        SpecificSoftening * ...
        ( Eps0{1}.*DEps0{1}(:,:,m) + Eps0{2}.*DEps0{2}(:,:,m) + ...
        2*Eps0{6}.*DEps0{6}(:,:,m) ).^2 ./ ...
        (DEps0{1}(:,:,m).^2 + DEps0{2}(:,:,m).^2 + 2*DEps0{6}(:,:,m).^2);
  end
  DQ = DWork / 2;
  Q = (Q0 - DQ) .* exp(-Softening);
  
  if ( min(min(min(Q))) < 0 )
    [Qmin,idx] = min(Q,[],3);
    idx = idx .* (Qmin<0);
    for m = 1 : Modes,
      for i = NonZeroComponents,
        Eps0{i} = Eps0{i} + DEps0{i}(:,:,m) .* (idx==m);
        DEps0{i}(:,:,m) = DEps0{i}(:,:,m) .* (idx~=m) + ...
            randn(mesh(2),mesh(1))/2*OmegaF/voxel0 .* (idx==m);
      end
      DEps0{2}(:,:,m) = -DEps0{1}(:,:,m);
    end;
  else
    q = reshape(Q,mesh(2)*mesh(1)*Modes,1);
    rate = trialfreq * exp( - q / kT );
    timeincrement = 1 / sum(rate);
    if ( timeincrement * strainrate < strainincrement )
      time = time + timeincrement;
      dice = rand(1);
      roll = cumsum(rate) / sum(rate);
      idx = min(find( (roll > dice)));
      mode = ceil( idx / mesh(2) / mesh(1) );
      col = ceil( idx / mesh(2) - (mode - 1)*mesh(1) );
      row = idx - (mode-1)*mesh(2)*mesh(1) - (col-1) * mesh(2);
      for i = NonZeroComponents,
        Eps0{i}(row,col) = Eps0{i}(row,col) + DEps0{i}(row,col,mode);
%         DEps0{i}(row,col,:) = [ -DEps0{i}(row,col,mode) %
%                     randn(Modes-1,1)/2*OmegaF/voxel0 ]; %
      end
%       DEps0{2}(row,col,:) = -DEps0{1}(row,col,:); %
    else
      time = time + strainincrement/strainrate;
    end
  end;
  
  if mod(step, printfreq) == 0
    figure(1); clf;
    pcolor(X{1}, X{2}, Eps0{1}.^2+Eps0{2}.^2+2*Eps0{6}.^2); 
    shading flat; axis equal;  axis off; 
    fprintf (1, 'step=%d, time=%e s, strain = %e\n', step, time, strain);
    print(gcf, '-djpeg', sprintf('Jpg/%06d.jpg',step))
    figure(2); clf; plot(strainhistory, stresshistory); drawnow;
    xlabel('\epsilon_{11}'); ylabel('\sigma_{11} [GPa]'); 
  end;
  step = step + 1;  
  strain = strain + strainrate * time;

  if (abs(strain) > 0.3)
    break;
  end;
  
end


% mesh(2)=4; mesh(1)=3; Modes=3; Q=randn(mesh(2),mesh(1),Modes)+1 %
% idx=find(reshape(Q,mesh(2)*mesh(1)*Modes,1)==Q(3,2,3));       %
% mode = ceil( idx / mesh(2) / mesh(1) ) %
% col = ceil( idx / mesh(2) - (mode - 1)*mesh(1) ) %
% row = idx - (mode-1)*mesh(2)*mesh(1) - (col-1) * mesh(2) %
% Q(3,2,3) - Q(row,col,mode) %
