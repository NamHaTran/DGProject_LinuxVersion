clear;
clc;
inpoed=dlmread('inpoed.txt',' ',13,1);
points=dlmread('Points.txt');
nedge=size(inpoed,2);
figure;
hold on;
for iedge=1:nedge
    point1=inpoed(1,iedge);
    point2=inpoed(2,iedge);
    XCoor=[points(point1+1,2) points(point2+1,2)];
    YCoor=[points(point1+1,3) points(point2+1,3)];
    plot(XCoor,YCoor,'-b');
end
axis equal;