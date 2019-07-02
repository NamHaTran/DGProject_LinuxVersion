function j2d=J2DCal(type,X,Y,a,b)
if type==3
    xA=X(1);
    xB=X(2);
    xC=X(3);
    yA=Y(1);
    yB=Y(2);
    yC=Y(3);
    dxa = (1 - b)*(xB - xA) / 4.0;
    dxb = a * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
    dya = (1 - b)*(yB - yA) / 4.0;
    dyb = a * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
elseif type==4
    xA=X(1);
    xB=X(2);
    xC=X(3);
    xD=X(4);
    yA=Y(1);
    yB=Y(2);
    yC=Y(3);
    yD=Y(4);
    dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*b;
	dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*a;
	dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*b;
	dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*a;
end
j2d = dxa * dyb - dxb * dya;