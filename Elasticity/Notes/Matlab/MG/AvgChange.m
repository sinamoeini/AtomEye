function avgchange = AvgChange(EpsoLast, Epso)
  TotalChange = Epso{1} * 0;
  for i = 1 : 6
    Change{i} = Epso{i} - EpsoLast{i};
    TotalChange = TotalChange + Change{i}.^2;
  end;
  avgchange=sqrt(sum(sum(TotalChange))/size(TotalChange,1)/size(TotalChange,2));
