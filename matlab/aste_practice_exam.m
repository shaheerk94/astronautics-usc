e = 0.5;
f = acos(e^2);
E = acos((e + cos(f)/(1+e*cos(f))));
u = 396800;
a = 20000;

t = (E - e*sin(e))/(u/a^3);
tf = E * sqrt(a^3/u)

f = 225 * pi/180;
T = 300 * 60; 
e = 0.25; 
a = (T^2* u/(2*pi))^(1/3); 
p = a * (1-e^2); 
di = 50/180 * pi; 

dv = di/(sqrt(p/u) * (1/(1 + e * cos(f))));
dv

vp = sqrt((u/a) * (1+e)/(1-e));
vp
vc = sqrt(u/a)
dv = vp - vc