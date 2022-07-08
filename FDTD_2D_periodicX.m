clc; clear; close all
%%
% Basic physical parmeter
eps0 = 8.854187817*1e-12;
mu0  = 4*pi*1e-7;
c0 = 1/sqrt(eps0*mu0);
um   = 1e-6;
nm   = 1e-9;
n_air = 1;
%%
% Parameter Setting
x_boundary = [0, 1*um];
y_boundary = [0, 4*um];
wavelength = 500*nm;

substrate = [0,0,diff(x_boundary),1*um];
n_s = 1.7742;   % Al2O3
grating = [250*nm,substrate(2)+substrate(4),500*nm,1*um];
n_g = 2.4431;   % GaN
dx = wavelength/10;
dy = wavelength/20;

f0 = c0/wavelength;
nt = 1000;
%% Enviroment build
x = x_boundary(1):dx:x_boundary(2);
y = y_boundary(1):dy:y_boundary(2);

nx = length(x); ny = length(y);
Hx = zeros(nx,ny);
Hy = zeros(nx,ny);
Ez = zeros(nx,ny);
cex = Ez; cey = Ez; chz = Ez;

dt = 1/(c0*sqrt(1/dx^2+1/dy^2));

%% Environment Material setting
structure_data = [substrate,n_s;grating,n_g];
environment = n_air*ones(nx,ny);

for i = 1:2
    [~,x_1] = min(abs(x-structure_data(i,1)));
    [~,y_1] = min(abs(y-structure_data(i,2)));
    [~,x_2] = min(abs(x-(structure_data(i,1)+structure_data(i,3))));
    [~,y_2] = min(abs(y-(structure_data(i,2)+structure_data(i,4))));
    environment(x_1:x_2,y_1:y_2) = structure_data(i,5);
end


mu = mu0;
eps = environment.^2*eps0;
%%
const_mu  = -dt./mu;
const_eps = -dt./eps;

figure(1)
pcolor(x,y,environment')
shading interp
colormap('jet')
hold on
rectangle('Position',substrate)
rectangle('Position',grating)
axis equal; axis tight
xlim([min(x),max(x)])
ylim([min(y),max(y)])
colorbar
hold off

figure(2)
colormap('jet')
for t = 1:nt
    % CEx & CEy
    cex(:,1:end-1) = diff(Ez,1,2)/dy;
    cex(:,end) = (0-Ez(:,end))/dy;
    Hx = Hx+const_mu.*cex;

    cey(1:end-1,:) = -diff(Ez,1,1)/dx;
    cey(end,:) = -(Ez(1,:)-Ez(end,:))/dx;
    Hy = Hy+const_mu.*cey;

    % CHz
    chz(1,1) = (Hy(1,1)-Hy(end,1))/dx-(Hx(1,1)-0)/dy;
    chz(2:end,1) = (diff(Hy(:,1),1,1)/dx) - ((Hx(2:end,1)-0)/dy);
    chz(1,2:end) = (Hy(1,2:end)-Hy(end,2:end))/dx-diff(Hx(1,:))/dy;
    
    chz(2:end,2:end) = diff(Hy(:,2:end),1,1)/dx-diff(Hx(2:end,:),1,2)/dy;
    
    Ez = Ez-const_eps.*chz;
    

%     Ez(:,2) = Ez(:,2)+sin(2*pi*f0*dt*(t-1));
    Ez(:,1) = Ez(:,2)+exp(-1i*2*pi*f0*dt*(t-1));
    
    mesh(x,y,real(Ez'))
    xlabel('X axis')
    ylabel('Y axis')
    shading flat
%     axis equal; axis tight
    colorbar
    title(['t = ',num2str(t-1)])
    drawnow
    
    
end