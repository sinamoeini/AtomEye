% SphericallyIsotropicRandomStrain.m %

function eta = SphericallyIsotropicRandomStrain()
  h = randn(1,3);
  k = h - sum(h)/3;
  Rhat = SphericallyIsotropicRotationMatrix;
  eta = Rhat' * diag(k) * Rhat;
  