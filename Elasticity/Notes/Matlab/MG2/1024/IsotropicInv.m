function Strain = IsotropicInv(Stress, lambda, mu)
  tra = -lambda/2/mu/(2*lambda+2*mu) * (Stress{1} + Stress{2});
  Strain{1} = tra + Stress{1}/2/mu;
  Strain{2} = tra + Stress{2}/2/mu;
  Strain{3} = Stress{3}/2/mu;
  