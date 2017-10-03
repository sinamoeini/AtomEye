% main.m %

N = 10000;

R = SphericallyIsotropicRotationMatrix;

eta5 = [];
for i = 1 : N
  eta = SphericallyIsotropicRandomStrain;
  Eta = R' * eta * R;
  eta5(i,:) = [eta(1,3) Eta(1,3)];
end

figure(1); hist(eta5(:,1),100);
figure(2); hist(eta5(:,2),100);
