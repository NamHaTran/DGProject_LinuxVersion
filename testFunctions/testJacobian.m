%test jacobian
xGauss=[-0.774596669241483377035853079956 0 0.774596669241483377035853079956];
wGauss=[0.555555555555555555555 0.888888888888888888888 0.55555555555555555555555];
GaussPt=zeros(3,3,2);
wGaussPt=zeros(3,3,2);
F=ones(3,3);
J2D=zeros(3,3);
%X=[-2 1 -1];
%Y=[3 3 6];
X=[-1 1 -1];
Y=[-1 -1 1];
for na=1:3
    for nb=1:3
        GaussPt(na,nb,1)=xGauss(na);
        GaussPt(na,nb,2)=xGauss(nb);
        wGaussPt(na,nb,1)=wGauss(na);
        wGaussPt(na,nb,2)=wGauss(nb);
    end
end

  for na=1:3
    for nb=1:3
        aG=GaussPt(na,nb,1);
        bG=GaussPt(na,nb,2);
        J2D(na,nb)=J2DCal(3,X,Y,aG,bG);
        B=basisFunction(3,3,aG,bG);
        F(na,nb)=B(1);
    end
end

inte=calcIntegral(F,J2D,wGaussPt)