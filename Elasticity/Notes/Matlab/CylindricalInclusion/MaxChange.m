function maxchange = MaxChange(EpsoLast, Epso)
  TotalChange = Epso{1} * 0;
  for i = 1 : 6
    Change{i} = Epso{i} - EpsoLast{i};
    TotalChange = TotalChange + Change{i}.^2;
  end;
  maxchange = sqrt(max(max(TotalChange)));
