tic
%% i ∂_t u_n = -(u_(n+1)+ u_(n-1) - 2u_n) - |u_n|^2 (L_alpha u)_n

% Integer lattice
L = 256; % Grid size
dx=1; % Unit grid spacing
x=[-L/2:dx:L/2-dx]'; % Grid

[row col] = size(x);
%% Parameters
alpha =25; % Fixed value of alpha used in calculating A
theta = -3.14159;
lambda = 1.17228; % This value isn't used in the initial condition. Just for clarity
A = 0.3; % Value of A to allow gain to be nonzero
kappa = 1.01934; % Fixed value of kappa used in calculating A

%initial_condition = A + rand([row,1])*0.01;
%initial_condition = A*exp(1i*(-lambda * 0 + x*theta)) + exp(1i*kappa*x);

initial_condition = exp(-(sin(x)).^2);
% Long Range periodic Interaction matrix
LRI_MATRIX = AL_PERIODIC(L/2,alpha); 

% Define the problem: d/dt u_n = -i * ( 1 + |u_n|^2 ) * L_alpha u_n) 
RHS = @(t,v) -1i*( LRI_MATRIX*v + abs(v).^2 .* LRI_MATRIX * v );

[tt,uu] = ode23(RHS, [0 30], initial_condition); % Solve problem

% Monitor the AL-NLS norm
%{
[row_uu col_uu] = size(uu);
mass = zeros(row_uu, 1);

for i = 1:size(tt)
    mass(i) = sum(log(1 + abs(uu(i,:).^2)));
end
%}
% Solution plotting
figure 
%imagesc(tt,x,abs(uu').^2); colorbar; shading interp;
pcolor(tt,x,abs(uu').^2);
xlabel("TIME");
ylabel("Grid"), 

title('Initial condition: $Ae^{i \theta x} + e^{i \kappa x}$','Interpreter','latex','FontSize',16)
subtitle("\theta="+theta+", \kappa="+kappa+", A="+A+", \alpha="+alpha)

colorbar;
cmocean('balance')
 shading interp
fontsize(16,"points")
%saveas(gcf, sprintf("images/APRIL/CW_thetamax_%0.2f_kmax_%0.2f_A_%0.2f.png",theta,kappa,A))

sprintf("Simulation completed in: %0.5f", toc)
