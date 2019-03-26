function B=basisFunction(type,order,a,b)
if type==3
  for i=1:3
    if i==1
      B(i)=1;
    elseif i==2
      B(i)=0.5*(3*b+1);
    elseif i==3
      B(i)=a*(1-b);
    end
  end
elseif type==4
    for i=1:3
      if i==1
      B(i)=1;
    elseif i==2
      B(i)=a;
    elseif i==3
      B(i)=b;
    end
    end
end