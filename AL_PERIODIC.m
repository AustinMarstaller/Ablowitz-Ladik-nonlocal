% Modified Long range interaction on the lattice with periodic BC on n in {-N, ..., N-1}
% See Appendix for the exact expression of periodic operator
function M = AL_PERIODIC_TEST(N,a)
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
    
    % Don't need this anymore, since we won't be overwriting the same entry
    % in `data`. We'll be incrementing the same index in M_DENSE_firstRow
    % twice.
%     if congruence_sum == congruence_difference
%         znr = 2*znr;
%     end
    
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