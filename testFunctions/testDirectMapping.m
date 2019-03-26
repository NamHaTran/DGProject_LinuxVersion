clear;
%direct mapping
X=[1 2 3 1];
Y=[0.5 -6 -2 0.5];

a=[-1 1 1 -1 -1];
b=[-1 -1 1 1 -1];

xCo=[];
yCo=[];

aCo=[-1]; bCo=-0.9:0.1:0.9;
for i=1:18
    aCo=[aCo 1];
    xCo = [xCo 0.25*(1 - aCo)*(1 - bCo(i))*X(1) + 0.25*(1 + aCo)*(1 - bCo(i))*X(2) + 0.25*(1 + bCo(i))*X(3)];
    yCo = [yCo 0.25*(1 - aCo)*(1 - bCo(i))*Y(1) + 0.25*(1 + aCo)*(1 - bCo(i))*Y(2) + 0.25*(1 + bCo(i))*Y(3)];
end

figure;
hold on;
plot(X,Y,'-');
plot(xCo,yCo,'r*');


figure;
hold on;
plot(a,b,'-');
plot(aCo,bCo,'r*');