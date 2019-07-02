function integral=calcIntegral(F,Jacobi,weight)
integral=0;
for na=1:3
    for nb=1:3
        w1=weight(na,nb,1);
        w2=weight(na,nb,2);
        integral=integral+w1*w2*Jacobi(na,nb)*F(na,nb);
    end
end