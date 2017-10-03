function Stress = isotropicmul(lambda, mu, Strain)
  tra = lambda .* ( Strain(1) + Strain(2) );
  Stress(1) = tra + 2*mu .* Strain(1);
  Stress(2) = tra + 2*mu .* Strain(2);
  Stress(3) = 2*mu .* Strain(3);
