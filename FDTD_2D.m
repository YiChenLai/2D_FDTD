clc; clear; close all

% Ba 
eps0 = 8.854187817*1e-12;
mu0  = 4*pi*1e-7;
c0 = 1/sqrt(eps0*mu0);
um   = 1e-6;
nm   = 1e-9;

% Parameter Setting
x_boundary = [0, 500*um];
y_boundary = [0, 500*um];
wavelength = 100*nm;
f0 = c0/wavelength;

nt = 1000;

% Enviroment build
x = linspace(x_boundary(1),x_boundary(2),501)*wavelength;
y = linspace(y_boundary(1),y_boundary(2),501)*wavelength;


dx = x(2);dy = y(2);
nx = length(x); ny = length(y);
dt = 1/(c0*sqrt(1/dx^2+1/dy^2));


Hx = zeros(nx,ny);
Hy = zeros(nx,ny);
Ez = zeros(nx,ny);
cex = Ez; cey = Ez; chz = Ez;

const_mu  = -dt/mu0;
const_eps = -dt/eps0;

colormap('jet')



for t = 1:nt
    % CEx & CEy
    cex(:,1:end-1) = diff(Ez,1,2)/dy;
    cex(:,end) = (0-Ez(:,end))/dy;
    Hx = Hx+const_mu.*cex;

    cey(1:end-1,:) = -diff(Ez,1,1)/dx;
    cey(end,:) = -(0-Ez(end,:))/dx;
    Hy = Hy+const_mu.*cey;

    % CHz
    chz(1,1) = (Hy(1,1)-0)/dx-(Hx(1,1)-0)/dy;
    chz(2:end,1) = (diff(Hy(:,1),1,1)/dx) - ((Hx(2:end,1)-0)/dy);
    chz(1,2:end) = (Hy(1,2:end)-0)/dx-diff(Hx(1,:))/dy;
    
    chz(2:end,2:end) = diff(Hy(:,2:end),1,1)/dx-diff(Hx(2:end,:),1,2)/dy;
    
    Ez = Ez-const_eps.*chz;
    

    Ez(20,round(ny/2)) = Ez(20,round(ny/2))...
        +sin(2*pi*f0*dt*t).*exp(-.5*((t-20)/8)^2);
    
    pcolor(x,y,Ez')
    xlabel('X axis')
    ylabel('Y axis')
    shading interp
    axis equal; axis tight
    colorbar
    drawnow;
    
    
end






