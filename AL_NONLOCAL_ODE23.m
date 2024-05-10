tic

L = 256; % Grid size
dx=1; % Unit grid spacing
x=[-L/2:dx:L/2-dx]'; % Grid

[row col] = size(x);
%% Parameters
alpha =1; % Fixed value of alpha used in calculating A
A = 0.05; % Value of A to allow gain to be nonzero
kappa = 0.00654995; % Fixed value of kappa used in calculating A

initial_condition = A + exp(1i*kappa*x);

% Long Range periodic Interaction matrix
LRI_MATRIX = AL_PERIODIC(L/2,alpha); 

% Define the problem: d/dt u_n = -i * ( 1 + |u_n|^2 ) * L_alpha u_n) 
RHS = @(t,v) -1i*( LRI_MATRIX*v + abs(v).^2 .* LRI_MATRIX * v );

[tt,uu] = ode23(RHS, [0 30], initial_condition); % Solve problem

% Solution plotting
figure 
subplot(2,1,1)
%imagesc(tt,x,abs(uu').^2); colorbar; shading interp;
pcolor(tt,x,abs(uu').^2);
xlabel("TIME");
ylabel("Grid"), 

title('Initial condition: $A + e^{i\kappa n}$','Interpreter','latex','FontSize',16)
subtitle("\kappa="+kappa+", A="+A+", \alpha="+alpha)

colorbar;
cmocean('balance')
shading interp
fontsize(16,"points")
%saveas(gcf, sprintf("images/APRIL/CW_thetamax_%0.2f_kmax_%0.2f_A_%0.2f.png",theta,kappa,A))

%% Plot the PSD at the final time
subplot(2,1,2)
n=length(tt);
uhat = fft(uu', n);
PSD =  abs(uhat)^2/n; 

kappa = (2*pi/n)*[-n/2:n/2-1]; % Discrete wavenumbers

plot(kappa, PSD)
title('PSD at final time: $|\hat{u}(t)|^2$','Interpreter','latex','FontSize',16)
colorbar;
cmocean('balance')
shading interp
fontsize(16,"points")


sprintf("Simulation completed in: %0.5f", toc)
