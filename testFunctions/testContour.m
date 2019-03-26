clc;
close;
clear;

%read data
toado=dlmread('cellCentroid.txt');
rhoArr=dlmread('rho20100.txt');
rho=rhoArr(:,4);
rhouArr=dlmread('rhou20100.txt');
rhou=rhouArr(:,4);
rhovArr=dlmread('rhov20100.txt');
rhov=rhovArr(:,4);
rhoEArr=dlmread('rhoE20100.txt');
rhoE=rhoEArr(:,4);
p=(1.4-1).*(rhoE-0.5.*rhou.^2./rho-0.5.*rhov.^2./rho);
T=p./(268.*rho);
vMag=sqrt((rhou./rho).^2+(rhov./rho).^2);
%data=importdata('static_pressure');
x=toado(:,1);
y=toado(:,2);
%p=data.data(:,4);
%enforce y==0 where y<1e-9 (to avoid NaN at centerline)
y(y<1e-9)=0;
%make rectangular grid
Nx = length(rho);
Ny = Nx;
xx = linspace(min(x),max(x),Nx);
yy = linspace(min(y),max(y),Ny);
[xi,yi]=meshgrid(xx,yy);
%interpolate results on rectangular grid
zi=griddata(x,y,rho,xi,yi);
zi2=griddata(x,y,p,xi,yi);
zi3=griddata(x,y,T,xi,yi);
zi4=griddata(x,y,vMag,xi,yi);
%mirror
%xi = [ xi(end:-1:2,: ); xi];
%yi = [-yi(end:-1:2,: ); yi];
%zi = [ zi(end:-1:2,: ); zi];

%plot
%figure()
%colormap jet;
%hcb=colorbar;
%colorTitleHandle = get(hcb,'Title');
%titleString = '[Density]';
%set(colorTitleHandle ,'String',titleString);
%title('Density Distribution');
%xlabel('X coordinate');
%ylabel('Y coordinate');
%contourf(xi,yi,zi,100, 'edgecolor', 'none');
%axis equal;
%axis([0 0.0508 -0.015 0.015]);

figure('Name','Density contour','NumberTitle','off');
[CC,hh]=contourf(xi,yi,zi,20,'LineColor','w');
set(hh,'LineColor','none');
axis equal;
colormap jet;
xlabel('z (m)');
ylabel('y (m)');
h = colorbar;
ylabel(h, 'kg/m3');

inpoed=dlmread('inpoed.txt',' ',13,1);
points=dlmread('Points.txt');
nedge=size(inpoed,2);
%figure;
hold on;
for iedge=1:nedge
    point1=inpoed(1,iedge);
    point2=inpoed(2,iedge);
    XCoor=[points(point1+1,2) points(point2+1,2)];
    YCoor=[points(point1+1,3) points(point2+1,3)];
    plot(XCoor,YCoor,'-b');
end
axis equal;


figure('Name','Pressure contour','NumberTitle','off');
[CC,hh]=contourf(xi,yi,zi2,20,'LineColor','w');
set(hh,'LineColor','none');
axis equal;
colormap jet;
xlabel('z (m)');
ylabel('y (m)');
h = colorbar;
ylabel(h, 'kg/m3');

inpoed=dlmread('inpoed.txt',' ',13,1);
points=dlmread('Points.txt');
nedge=size(inpoed,2);
%figure;
hold on;
for iedge=1:nedge
    point1=inpoed(1,iedge);
    point2=inpoed(2,iedge);
    XCoor=[points(point1+1,2) points(point2+1,2)];
    YCoor=[points(point1+1,3) points(point2+1,3)];
    plot(XCoor,YCoor,'-b');
end
axis equal;


figure('Name','Temperature contour','NumberTitle','off');
[CC,hh]=contourf(xi,yi,zi3,20,'LineColor','w');
set(hh,'LineColor','none');
axis equal;
colormap jet;
xlabel('z (m)');
ylabel('y (m)');
h = colorbar;
ylabel(h, 'kg/m3');

inpoed=dlmread('inpoed.txt',' ',13,1);
points=dlmread('Points.txt');
nedge=size(inpoed,2);
%figure;
hold on;
for iedge=1:nedge
    point1=inpoed(1,iedge);
    point2=inpoed(2,iedge);
    XCoor=[points(point1+1,2) points(point2+1,2)];
    YCoor=[points(point1+1,3) points(point2+1,3)];
    plot(XCoor,YCoor,'-b');
end
axis equal;

figure('Name','Temperature contour','NumberTitle','off');
[CC,hh]=contourf(xi,yi,zi4,20,'LineColor','w');
set(hh,'LineColor','none');
axis equal;
colormap jet;
xlabel('z (m)');
ylabel('y (m)');
h = colorbar;
ylabel(h, 'kg/m3');

inpoed=dlmread('inpoed.txt',' ',13,1);
points=dlmread('Points.txt');
nedge=size(inpoed,2);
%figure;
hold on;
for iedge=1:nedge
    point1=inpoed(1,iedge);
    point2=inpoed(2,iedge);
    XCoor=[points(point1+1,2) points(point2+1,2)];
    YCoor=[points(point1+1,3) points(point2+1,3)];
    plot(XCoor,YCoor,'-b');
end
axis equal;