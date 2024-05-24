close all
clear all
tic

L = 256; % Grid size
dx=1; % Unit grid spacing
x=[-L/2:dx:L/2-dx]'; % Grid

[row col] = size(x);
%% Parameters
alpha = 1; % Fixed value of alpha used in calculating A
A     = 0.3; % Value of A to allow gain to be nonzero
kappa = 0.0905602; % Fixed value of kappa used in calculating A

initial_condition = A + 0.001*exp(1i*kappa*x);

% Long Range periodic Interaction matrix
LRI_MATRIX = AL_PERIODIC(L/2,alpha); 

% Define the problem: d/dt u_n = -i * ( 1 + |u_n|^2 ) * L_alpha u_n) 
RHS = @(t,v) -1i*( LRI_MATRIX*v + abs(v).^2 .* LRI_MATRIX * v );

[tt,uu] = ode23(RHS, [0,30],initial_condition); % Solve problem

% Solution plotting
figure 
%imagesc(tt,x,abs(uu').^2); 
pcolor(tt,x,abs((uu-A)').^2);
colorbar; 
shading interp;
xlabel("TIME");
ylabel("Grid"), 

title('Initial condition: $A + e^{i\kappa n}$','Interpreter','latex','FontSize',16)
subtitle("\kappa="+kappa+", A="+A+", \alpha="+alpha)

%% Compute and plot PSD
%figure
%subplot(2,1,1)
%n=length(x);
%uhat = fft(uu, n);
%PSD =  uhat.*conj(uhat)/n;
%dt = 0.001;
%freq = 1/(dt*n)*(0:n-1);
%L = 1:floor(n/2);
%plot(freq,PSD, 'LineWidth',2.5)
%set(gca,'FontSize',16)
%xlabel('Frequency')
%ylabel('PSD')

%subplot(2,1,2)
%plot(freq(L),PSD(L,:), 'LineWidth',2.5)
%set(gca,'FontSize',16)
%xlabel('Frequency(L)')
%ylabel('PSD(L)')

sprintf("Simulation completed in: %0.5f", toc)

