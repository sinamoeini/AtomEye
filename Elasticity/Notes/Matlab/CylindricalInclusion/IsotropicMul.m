function Stress = IsotropicMul(lambda, mu, Strain)
  tra = lambda .* ( Strain{1} + Strain{2} + Strain{3} );
  Stress{1} = tra + 2*mu .* Strain{1};
  Stress{2} = tra + 2*mu .* Strain{2};
  Stress{3} = tra + 2*mu .* Strain{3};
  Stress{4} = 2*mu .* Strain{4};
  Stress{5} = 2*mu .* Strain{5};
  Stress{6} = 2*mu .* Strain{6};
