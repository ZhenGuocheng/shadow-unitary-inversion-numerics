clear; 

d = 2;

% database
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];
I = eye(d);
Obs = Z;

P23 = [1 0 0 0;
       0 0 1 0;
       0 1 0 0;
       0 0 0 1];
%Schur Matrix
SW_trans = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
0 1/sqrt(2) 0 0 0 0 0 0 0 0 0 0 1/sqrt(2) 0 0 0;
0 0 0 0 1/sqrt(2) 0 0 0 0 0 0 0 0 1/sqrt(2) 0 0;
0 0 0 0 0 0 0 1/sqrt(2) 0 0 0 0 0 0 1/sqrt(2) 0;
0 0 0 0 0 0 0 0 0 0 1/sqrt(2) 0 0 0 0 1/sqrt(2);
0 1/sqrt(2) 0 0 0 0 0 0 0 0 0 0 -1/sqrt(2) 0 0 0;
0 0 0 0 1/sqrt(2) 0 0 0 0 0 0 0 0 -1/sqrt(2) 0 0;
0 0 0 0 0 0 0 1/sqrt(2) 0 0 0 0 0 0 -1/sqrt(2) 0;
0 0 0 0 0 0 0 0 0 0 1/sqrt(2) 0 0 0 0 -1/sqrt(2);
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0];


SWAP23 = kron(kron(I,P23),I);

% initialize for integral
n_U = 2000;

cvx_begin sdp
    variable lamb(8) nonnegative                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    
    % construct C
    CI3 = [lamb(1)*eye(3), zeros(3,3), zeros(3,3), zeros(3,3);
           zeros(3,3), lamb(2)*eye(3), zeros(3,3), zeros(3,3);
           zeros(3,3), zeros(3,3), lamb(3)*eye(3), zeros(3,3);
           zeros(3,3), zeros(3,3), zeros(3,3), lamb(4)*eye(3)];
    CI1 = [lamb(5), 0, 0, 0;
           0, lamb(6), 0, 0;
           0, 0, lamb(7), 0;
           0, 0, 0, lamb(8)];

    C = [CI3, zeros(12, 4);
         zeros(4, 12), CI1];

    C1 = SWAP23' * SW_trans * C * SW_trans' * SWAP23;

    % constraints for comb
    PartialTrace(C1,4,[d d d d]) == kron(PartialTrace(C1,[3 4], [d d d d]),eye(d)/(d));
    PartialTrace(C1,[2 3 4],[d d d d]) == eye(d)*(d);
    % constraints from any unitary
    cost = 0; 
     for i = 1:n_U
        U = RandomUnitary(d);
        U_ketbra = U(:)*U(:)';
        term1 = PartialTrace(PartialTranspose(C1, [2 3], [d d d d]) * kron(kron(I,U_ketbra),I), [2 3], [d d d d]);  % C \star |U⟩⟩⟨⟨U|
        term2 = Tensor(eye(d),Obs.');  % I_{P}\otimes O^T
        term3 = U * Obs * U';  % UOU^{-1}
            
        % Goal Function
        trace_term = norm(PartialTrace(term1.' * term2, 2, [d d]) - term3,'fro');
        cost = cost + trace_term;
     end
    cost = cost / (n_U);
    

    minimize cost

cvx_end

cost
