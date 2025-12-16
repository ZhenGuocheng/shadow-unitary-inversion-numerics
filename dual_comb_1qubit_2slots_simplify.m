clear; 

d = 2;

% database
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];
I = eye(2);
Obs = Z;

% initialize for integral
n_U = 2000;

% construct SW trans
B = load("U_sch.mat"); % B is the Schur Matrix for U^{\otimes 3}, U \in U(2)
B = B.data;

label = [1,4,2,2];
t = 2;

Id = eye(2^t);
I2 = eye(2);

U = {};
for a = 1:numel(label)-1
    col = sum(label(1:a)):sum(label(1:a+1));
    col = col(1:numel(col)-1);
    for l = 1:2
        for k = 1:2^t
            for j = col
                U{end+1} = Tensor(B(:,j), Id(:, k), I2(:, l));
            end
        end
    end
end

U = cell2mat(U);
SW_trans = U;

cvx_begin sdp
    variable lambda(24)

    variable B0203(4, 4) complex
    variable B0607(4, 4) complex

    variable A0109(2, 2) complex
    variable A0210(2, 2) complex
    variable A0211(2, 2) complex
    variable A0310(2, 2) complex
    variable A0311(2, 2) complex
    variable A0412(2, 2) complex
    variable A0513(2, 2) complex
    variable A0614(2, 2) complex
    variable A0615(2, 2) complex
    variable A0714(2, 2) complex
    variable A0715(2, 2) complex
    variable A0816(2, 2) complex
    variable A1415(2, 2) complex
    variable A0203(2, 2) complex
    variable A1011(2, 2) complex
    variable A0607(2, 2) complex

    % construct C
    blocks = cell(8,1);
    for i = 1:8
        blocks{i} = lambda(i) * eye(4);
    end
    C_4_diag = blkdiag(blocks{:});
    
    O = zeros(4,4);
    C_4_off = [O, O, O, O, O, O, O, O;
               O, O, B0203, O, O, O, O, O;
               O, B0203', O, O, O, O, O, O;
               O, O, O, O, O, O, O, O;
               O, O, O, O, O, O, O, O;
               O, O, O, O, O, O, B0607, O;
               O, O, O, O, O, B0607', O, O;
               O, O, O, O, O, O, O, O];

    C_4 = C_4_diag + C_4_off;

    blocks = cell(16,1);
    for i = 9:24
        blocks{i} = lambda(i) * eye(2);
    end
    C_2_diag = blkdiag(blocks{:});
    O = zeros(2,2);
    C_2_off = [O, O, O, O, O, O, O, O, A0109, O, O, O, O, O, O, O;
                 O, O, A0203, O, O, O, O, O, O, A0210, A0211, O, O, O, O, O;
                 O, A0203', O, O, O, O, O, O, O, A0310, A0311, O, O, O, O, O;
                 O, O, O, O, O, O, O, O, O, O, O, A0412, O, O, O, O;
                 O, O, O, O, O, O, O, O, O, O, O, O, A0513, O, O, O;
                 O, O, O, O, O, O, A0607, O, O, O, O, O, O, A0614, A0615, O;
                 O, O, O, O, O, A0607', O, O, O, O, O, O, O, A0714, A0715, O;
                 O, O, O, O, O, O, O, O, O, O, O, O, O, O, O, A0816;
                 A0109', O, O, O, O, O, O, O, O, O, O, O, O, O, O, O;
                 O, A0210', A0310', O, O, O, O, O, O, O, A1011, O, O, O, O, O;
                 O, A0211', A0311', O, O, O, O, O, O, A1011', O, O, O, O, O, O;
                 O, O, O, A0412', O, O, O, O, O, O, O, O, O, O, O, O;
                 O, O, O, O, A0513', O, O, O, O, O, O, O, O, O, O, O;
                 O, O, O, O, O, A0614', A0714', O, O, O, O, O, O, O, A1415, O;
                 O, O, O, O, O, A0615', A0715', O, O, O, O, O, O, A1415', O, O;
                 O, O, O, O, O, O, O, A0816', O, O, O, O, O, O, O, O];

    C_2 = C_2_diag + C_2_off;

    C_BK = {C_4, C_2}; 
    C = blkdiag(C_BK{:});

    %SWAP to the original computation bases
    C1 = PermuteSystems(SW_trans * C * SW_trans', [1 5 3 4 2 6], [d d d d d d]);

    % constraints for comb
    C1 >= 0;
    PartialTrace(C1,6,[d d d d d d]) == kron(PartialTrace(C1,[5 6],[d d d d d d]),eye(d)/d);
    PartialTrace(C1,[4 5 6],[d d d d d d]) == kron(PartialTrace(C1,[3 4 5 6],[d d d d d d]),eye(d)/d);
    PartialTrace(C1,[2 3 4 5 6],[d d d d d d]) == eye(d)*d^2;

    % constraints from any unitary
    cost = 0;
    for i = 1:n_U
    % QETLAB Function
    U = RandomUnitary(d);
    U_ketbra = U(:)*U(:)';
    term1 = PartialTrace(PartialTranspose(C1, [2 3 4 5], [d d d d d d]) * kron(kron(I,kron(U_ketbra,U_ketbra)),I), [2 3 4 5], [d d d d d d]);  % C1 \star |U⟩⟩⟨⟨U|^{⊗2}
    term2 = Tensor(I,Obs.');  % I \otimes O^T
    term3 = U * Obs * U';  % UOU^{-1}
        
    % Goal Function
    trace_term = norm(PartialTrace(term1.' * term2, 2, [d d]) - term3,'fro');
    cost = cost + trace_term;
    end
    cost = cost / (n_U);
    
    minimize cost

cvx_end