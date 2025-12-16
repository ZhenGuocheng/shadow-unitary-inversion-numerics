clear; 

%% initialization
d = 2;

% database
X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];
I = eye(d);
Obs = Z;

% initialize for integral
n_U = 2000;

%% construct SW trans
yd = {{5},{3,3,3},{1,1}};
l1 = [1,2,2,3,2,3,3,4];
l2 = [1,2];
t = 3;
b_dims = [5,3,1];

B = load("U_sch_3slot.mat"); % B is the Schur Matrix for U^{\otimes 4}, U \in U(2)
B = B.data;
SW_trans = ConstructSchurTransformation(B, cell2mat([yd{:}]), [2^t, 2]);

%% generate equivalent blocks
block_store = {};
for yd_e=yd
    tup_block = {};
    for r=1:numel(yd_e{1})
        for l=l2
            for k=l1
                tup_block{end+1} = [yd_e{1}{1},k,l];
            end
        end
    end
    block_store{end+1} = tup_block;
end

equiv_block_coor = {};
for tup=block_store
    local_coor = {};
    for i = 1:numel(tup{1})
        for j = 1:numel(tup{1})
            if isequal(tup{1}{i}, tup{1}{j}) && (i < j)
                local_coor{end+1} = [i;j];
            end
        end
    end
    local_coor = cell2mat(local_coor);
    equiv_block_coor{end+1} = local_coor;
end


%% cvx optimization
cvx_begin sdp

    % construct C
    variable lamb1(numel(block_store{1}))
    variable lamb2(numel(block_store{2}))
    variable lamb3(numel(block_store{3}))
    
    variable A1(b_dims(1), b_dims(1), size(equiv_block_coor{1},2)) complex
    variable A2(b_dims(2), b_dims(2), size(equiv_block_coor{2},2)) complex
    variable A3(b_dims(3), b_dims(3), size(equiv_block_coor{3},2)) complex
    
    expression X1(b_dims(1)*numel(block_store{1}), b_dims(1)*numel(block_store{1}))
    expression X2(b_dims(2)*numel(block_store{2}), b_dims(2)*numel(block_store{2}))
    expression X3(b_dims(3)*numel(block_store{3}), b_dims(3)*numel(block_store{3}))

    % Block-diagonal Matrix
    C1 = ConstructBlockVariable(lamb1, A1, X1, equiv_block_coor{1}, b_dims(1), numel(block_store{1}));
    C2 = ConstructBlockVariable(lamb2, A2, X2, equiv_block_coor{2}, b_dims(2), numel(block_store{2}));
    C3 = ConstructBlockVariable(lamb3, A3, X3, equiv_block_coor{3}, b_dims(3), numel(block_store{3}));

    C_BK = {C1, C2, C3}; 
    C = blkdiag(C_BK{:});

    % SWAP to the original computation bases
    C = PermuteSystems(SW_trans * C * SW_trans', [1 5 4 6 2 7 3 8], [d d d d d d d d]);

    % constraints for comb
    C >= 0;
    PartialTrace(C,8, [d d d d d d d d])==kron(PartialTrace(C,[7 8], [d d d d d d d d]),eye(d)/d);
    PartialTrace(C,[6 7 8], [d d d d d d d d])==kron(PartialTrace(C,[5 6 7 8], [d d d d d d d d]),eye(d)/d);
    PartialTrace(C,[4 5 6 7 8], [d d d d d d d d])==kron(PartialTrace(C,[3 4 5 6 7 8], [d d d d d d d d]),eye(d)/d);
    PartialTrace(C,[2 3 4 5 6 7 8], [d d d d d d d d])==eye(d)*d^3;

    % constraints from any unitary
    cost = 0;
    for i = 1:n_U
    % QETLAB Function
    U = RandomUnitary(d);
    U_ketbra = U(:)*U(:)';
    term1 = PartialTrace(PartialTranspose(C, [2 3 4 5 6 7], [d d d d d d d d]) * Tensor(I, U_ketbra, U_ketbra, U_ketbra, I), [2 5 3 6 4 7], [d d d d d d d d]);  % C \star |U⟩⟩⟨⟨U|^{⊗3}
    term2 = Tensor(I,Obs.');  % I \otimes O
    term3 = U * Obs * U';  % UOU^{-1}
        
    % Goal Function
    trace_term = norm(PartialTrace(term1.' * term2, 2, [d d]) - term3,'fro');
    cost = cost + trace_term;
    end
    cost = cost / (n_U);
    
    minimize cost

cvx_end


%% functions
function U = ConstructSchurTransformation(vec, yd_label, dims)
    label = [1, yd_label];
    Id = eye(dims(1));
    I2 = eye(dims(2));
    U = {};
    for a = 1:numel(label)-1
        col = sum(label(1:a)):sum(label(1:a+1));
        col = col(1:numel(col)-1);
        for l = 1:dims(2)
            for k = 1:dims(1)
                for j = col
                    U{end+1} = Tensor(vec(:,j), Id(:, k), I2(:, l));
                end
            end
        end
    end
    U = cell2mat(U);
end


function C = ConstructBlockVariable(lamb, A, X, coors, b_dim, coor_dim)

    % diagonal blocks
    blocks = cell(coor_dim, 1);
    for i = 1:coor_dim
        blocks{i} = lamb(i) * eye(b_dim);
    end
    C_diag = blkdiag(blocks{:});

    % off-diagonal blocks
    C_off = X;
    for j = 1:size(coors,2)
        coor = coors(:, j);
        precoor = coor - ones(2, 1);
        C_off(precoor(1)*b_dim+1:coor(1)*b_dim, precoor(2)*b_dim+1:coor(2)*b_dim) = A(:, :, j);
        C_off(precoor(2)*b_dim+1:coor(2)*b_dim, precoor(1)*b_dim+1:coor(1)*b_dim) = A(:, :, j)';
    end
    C = C_diag + C_off;
end
