clear;
clc;
Xinput=0.0020023143709651178;
Yinput=3.3895735607835068e-18;
C=[0.0019000000000000000 2.8151099999999999e-18];
A=[0.0020125500000000001 -0.00010547300000000000];
B=[0.0020153100000000002 3.4625400000000000e-18];

aG=[-0.7745 -0.7745 -0.7745;
    0 0 0;
    0.7745 0.7745 0.7745];
bG=[-0.7745 0 0.7745;
    -0.7745 0 0.7745;
    -0.7745 0 0.7745];
aGL=[-1 -1 -1;
    0 0 0;
    1 1 1];
bGL=[-1 0 1;
    -1 0 1;
    -1 0 1];

D=[0.0018974000000000000 -9.9438300000000000e-05];
X=[A(1) B(1) C(1) D(1) A(1)];
Y=[A(2) B(2) C(2) D(2) A(2)];
Ax=(A(1)-B(1)-D(1)+C(1))/4;
Ay=(A(2)-B(2)-D(2)+C(2))/4;
Bx=(-A(1)+B(1)-D(1)+C(1))/4;
By=(-A(2)+B(2)-D(2)+C(2))/4;
Cx=(-A(1)-B(1)+D(1)+C(1))/4;
Cy=(-A(2)-B(2)+D(2)+C(2))/4;
Dx=(A(1)+B(1)+D(1)+C(1))/4-Xinput;
Dy=(A(2)+B(2)+D(2)+C(2))/4-Yinput;

%{
Ax=(A(1)-B(1))/4;
Ay=(A(2)-B(2))/4;
Bx=(-A(1)+B(1))/4;
By=(-A(2)+B(2))/4;
Cx=(-A(1)-B(1)+2*C(1))/4;
Cy=(-A(2)-B(2)+2*C(2))/4;
Dx=(A(1)+B(1)+2*C(1))/4-Xinput;
Dy=(A(2)+B(2)+2*C(2))/4-Yinput;
%}
AA=-Ax*By+Ay*Bx;
BB=-Ax*Dy+Bx*Cy-Cx*By+Ay*Dx;
CC=Dx*Cy-Cx*Dy;
p=[AA BB CC];
No=roots(p);
if length(No)==1
    aOut=No(1);
elseif length(No)==2
    for i=1:2
        if (abs(No(i))<=1 && abs(No(i))>=0)
            aOut=No(i);
            break;
        else
            aOut=2;
        end
    end
else
    disp('phuong trinh vo nghiem');
    aOut=2;
end
bOut=-(aOut*By+Dy)/(aOut*Ay+Cy);

%cach 2
AA1=-Ay*Bx+Ax*By;
BB1=-Ay*Dx+By*Cx-Cy*Bx+Ax*Dy;
CC1=Dy*Cx-Cy*Dx;
p1=[AA1 BB1 CC1];
No1=roots(p1);
if length(No1)==1
    aOut1=No1(1);
elseif length(No1)==2
    for i=1:2
        if (abs(No1(i))<=1 && abs(No1(i))>=0)
            aOut1=No1(i);
            break;
        else
            aOut1=2;
        end
    end
else
    disp('phuong trinh vo nghiem');
    aOut1=2;
end
bOut1=-(aOut1*Bx+Dx)/(aOut1*Ax+Cx);

figure;
hold on;
plot(X,Y,'-b');
plot(Xinput,Yinput,'or');
grid on;
axis equal;
figure;
hold on;
plot([-1 1 1 -1 -1],[-1 -1 1 1 -1],'-b');
plot(aOut,bOut,'or');
plot(aOut1,bOut1,'*c');
for nG1=1:3
    for nG2=1:3
        plot(aG(nG1,nG2),bGL(nG1,nG2),'ok');
        plot(aGL(nG1,nG2),bG(nG1,nG2),'ok');
    end
end
grid on;
axis equal;