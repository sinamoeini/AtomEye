
eps = 0.00001;
eta = 0.5;

S_search = eta*log(eta)+(1-eta)*log(1-eta);

S_total = log(eps) - 1 - (1-eta)*log(1-eta);

S_total / S_search