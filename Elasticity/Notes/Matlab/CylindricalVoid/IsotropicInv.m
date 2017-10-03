function Strain = IsotropicInv(Stress, lambda, mu)
  tra = -lambda/2/mu/(3*lambda+2*mu) * (Stress{1} + Stress{2} + Stress{3});
  Strain{1} = tra + Stress{1}/2/mu;
  Strain{2} = tra + Stress{2}/2/mu;
  Strain{3} = tra + Stress{3}/2/mu;
  Strain{4} = Stress{4}/2/mu;
  Strain{5} = Stress{5}/2/mu;
  Strain{6} = Stress{6}/2/mu;
  