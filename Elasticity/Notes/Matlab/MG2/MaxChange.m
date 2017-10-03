function maxchange = MaxChange(EpsoLast, Epso)
  TotalChange = (Epso{1}-EpsoLast{1}).^2 + (Epso{2}-EpsoLast{2}).^2 + ...
            2 * (Epso{3}-EpsoLast{3}).^2;
  maxchange = sqrt(max(max(TotalChange)));
