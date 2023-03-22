
%% A Matlab code written by and developed by HUSSEIN A. H. Muhammed March 2023.
%% B.Sc.H and M.Sc. (Honuors) in Geosciences and Computational Geophysics, respectively.
%% This code simulates the 1D decoupled Burguer's Equation with a source term (Q^-), 
%% in Conventional coordinate system.

clear;

% 1D decoupled Burguer's Equation
% Ct = -uCx + KCxx + Q^-

%% Prepare the movie file
    vidObj = VideoWriter('DecoupledBE1D.avi');
    open(vidObj);

%% Domain
% Space
Lx=10;
dx=0.5;
nx=fix(Lx/dx);
x=linspace(0, Lx, nx);
%Re = 0.05;     % Re is the Reynold's Number = v/rho*U*L

% Time
T=7;

CFL=0.1;
dt=CFL*dx/1;
nt=fix(T/dt);

% Field Arrays
% Variables
C=zeros(nx,1);

% Parameters
u=zeros(nx,1);
K=zeros(nx, 1);     % Reynold's number
K(:)=1;

% Intial Conditions
t=0;
C(:)=0;
u(:)=sin(pi*x);
K(:)=1;

% Time stepping Loop
for n=1:nt
    % Boundary Conditions
    C(end) = C(end-1);
    C(1) = C(end);

    % Source
    if n==1
        C(1)=1;
    end

    % CDFDTD Solution
    t=t+dt;
    Cn=C;
    for i=2:nx-1
        % Advection term
        A = u(i)*(Cn(i)-Cn(i-1))/dx;
        % Diffusion term
        D = K(i)*(Cn(i+1)-2*Cn(i)+Cn(i-1))/dx^2;
        % Euler's Method
        C(i) = Cn(i) + dt*(-A + D);
    end

% Visualize at selected steps
clf;
plot(x, C, '-O');
title(sprintf('1D-BurgerEqu. elapsed time = %.2f', t));
axis([0 Lx 0 0.1]);
xlabel('displacement (unit)');
ylabel('Temperature (°)');
shg; pause(0.1);

 % Write each frame to the file
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
end


%% Close the file
close(vidObj);

