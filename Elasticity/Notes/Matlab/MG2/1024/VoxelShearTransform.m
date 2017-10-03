Eps0{1}(row,col) = Eps0{1}(row,col) + DEps0{1}(row,col,mode);
Eps0{2}(row,col) = Eps0{2}(row,col) + DEps0{2}(row,col,mode);
Eps0{3}(row,col) = Eps0{3}(row,col) + DEps0{3}(row,col,mode);

DeltaSoftening = DEps0{1}(row,col,mode)^2 + ...
    DEps0{2}(row,col,mode)^2 + 2*DEps0{3}(row,col,mode)^2;

PermSoftening(row,col) = PermSoftening(row,col) + DeltaSoftening * PermSoft;
if (PermSoftening(row,col) > PermSofteningCap) 
    PermSoftening(row,col) = PermSofteningCap;
end

TempSoftening(row,col) = DeltaSoftening * TempSoft;
tlast(row,col) = time;

Q0(row,col,:) = Q0amp * (1-RandomRatio + RandomRatio*rand(Modes,1));
DEps0{1}(row,col,:) = 0.5 * randn(Modes,1)/2*OmegaF/voxel0;
DEps0{3}(row,col,:) = 0.5 * randn(Modes,1)/2*OmegaF/voxel0;
DEps0{2}(row,col,:) = -DEps0{1}(row,col,:);
