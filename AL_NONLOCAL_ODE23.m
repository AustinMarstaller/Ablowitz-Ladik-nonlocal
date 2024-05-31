close all
clear all
tic

L = 256; % Grid size
dx=1; % Unit grid spacing
x=[-L/2:dx:L/2-dx]'; % Grid

[row col] = size(x);
%% Parameters
alpha = 1; % Fixed value of alpha used in calculating A
A     = 0.25; % Value of A to allow gain to be nonzero
kappa = 0.0653585; % Fixed value of kappa used in calculating A
delta = 0.0001;

initial_condition = A + delta*exp(1i*kappa*x);

% Long Range periodic Interaction matrix
LRI_MATRIX = AL_PERIODIC(L/2,alpha); 

% Define the problem: d/dt u_n = -i * ( 1 + |u_n|^2 ) * L_alpha u_n) 
RHS = @(t,u) -1i*( LRI_MATRIX*u + abs(u).^2 .* (LRI_MATRIX * u) );

[tt,u] = ode23t(RHS, [0,1000],initial_condition); % Solve problem

% Solution plotting
figure 

%imagesc(tt,x,abs(u-A)'.^2); 
pcolor(tt,x,abs((u-A)').^2);
cmocean('balance')
colorbar; 
shading interp;
xlabel("TIME");
ylabel("Grid"), 

title('Initial condition: $A + \delta e^{i\kappa n}$','Interpreter','latex','FontSize',16)
subtitle("\kappa="+kappa+", A="+A+", \alpha="+alpha+", \delta="+delta)
%exportgraphics(gcf,'MI_attempt.pdf','ContentType','vector');


%% Compute and plot spatial FT at final time
figure
subplot(2,1,1)
n=length(x);
uhat = fft(u(end,:)-A, n);
PSD =  uhat.*conj(uhat)/n;
PSD = fftshift(PSD);
freq = (2*pi/n)*[-L/2:L/2-1];
%plot(freq,PSD, 'LineWidth',2.5)

plot(freq(126:133),PSD(126:133), 'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Frequency')
ylabel('Spatial FT')
title('Spatial FT at final time: $|\widehat{(u-A)}_{tf}|^2$','Interpreter','latex','FontSize',16)
%exportgraphics(gcf,'SpatialFT_attempt.pdf','ContentType','vector');
%Dormand-Prince integratorexportgraphics(f,'bartransparent.pdf','ContentType','vector',...
%               'BackgroundColor','none')

subplot(2,1,2)
n=length(x);
uhat = fft(u(end,:), n);
PSD =  uhat.*conj(uhat)/n;
PSD = fftshift(PSD);
freq = (2*pi/n)*[-L/2:L/2-1];
plot(freq(126:133),PSD(126:133), 'LineWidth',2.5)
set(gca,'FontSize',16)
xlabel('Frequency')
ylabel('Spatial FT')
title('Spatial FT at final time: $|\hat{u}_{tf}|^2$','Interpreter','latex','FontSize',16)



%% Coefficient matrix
sprintf("Simulation completed in: %0.5f", toc)
function M = AL_PERIODIC(N,a)
tic

I_N = -N:1:N-1;
r_INDEX = 1:1:2*N-1;

%% Initialize the coefficient matrix for the LRI operator

%% Populate the diagonal 
M_DIAG = 2*zeta(1+a)*(1 / (2*N)^(1+a) ) * eye(2*N);

temp = zeros(1,2*N);

n = 1;
for r = 1:numel(r_INDEX)
    % Go through the r indexed sum.
    
    % Compute the indices n+r, n-r.
    congruence_sum = mod(I_N(n) + r_INDEX(r), 2*N);
    congruence_difference = mod(I_N(n) - r_INDEX(r), 2*N);
    
    % Compute the coefficient values
    znr =  (hurwitzZeta(1+a, r/(2*N)));
       
    znr = 1/(2*N)^(1+a) * znr;
    
    % Increment M_DENSE_firstRow on the corresponding columns
    temp(congruence_sum        + 1) = temp(congruence_sum        + 1) + znr;
    temp(congruence_difference + 1) = temp(congruence_difference + 1) + znr;
end
% Do the same order flipping that the original did
M_DENSE = toeplitz([temp(N+1:end),temp(1:N)]);

M = M_DENSE + M_DIAG;
sprintf("Nonlocal operator constructed in: %0.5f",toc)
end
