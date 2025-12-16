clear; 

d = 2;
k = 2;
% database
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];
I = eye(2);
Obs = Z;
% initialize for integral
n_U = 2000;

total_dim = d^2 * (d)^(2*k);

cvx_begin sdp
    variable C(total_dim, total_dim) complex semidefinite

    % constraints for comb
    C >= 0;
    PartialTrace(C,6,[d d d d d d]) == kron(PartialTrace(C,[4 5 6],[d d d d d d]),eye(d^2)/d^2);
    PartialTrace(C,[2 3 4 5 6],[d d d d d d]) == eye(d)*d^2;

    % constraints from any unitary
    cost = 0;
    for i = 1:n_U
    % QETLAB Function
    U = RandomUnitary(d);
    U_ketbra = U(:)*U(:)';
    term1 = PartialTrace(PartialTranspose(C, [2 3 4 5], [d d d d d d]) * Tensor(I,PermuteSystems(kron(U_ketbra,U_ketbra), [1 3 2 4], [d d d d]),I), [2 3 4 5], [d d d d d d]);  % C \star |U⟩⟩⟨⟨U|^{⊗2}
    term2 = Tensor(I,Obs.');  % I \otimes O^T
    term3 = U * Obs * U';  % UOU^{-1}

    trace_term = norm(PartialTrace(term1.' * term2, 2, [d d]) - term3,"fro");
    cost = cost + trace_term;
     end
    cost = cost / (n_U);

    minimize cost

cvx_end