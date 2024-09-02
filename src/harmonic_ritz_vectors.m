function dy = harmonic_ritz_vectors(Fold, G, k, s, V, tol)
    opts.tol = tol;
    opts.v0 = ones(size(Fold, 1), 1);
    E = zeros(s, k);
    D = zeros(k, 1);
    
    % Compute the harmonic Ritz vectors
    [E2, D2] = eigs(Fold, G, k, 'LM', opts);
    for p = 1:k
        D(p, 1) = abs(D2(p, p));
    end
    [~, I] = sort(D, 1);
    for q = 1:k
        E(:, q) = E2(:, I(q, 1));
    end
    dy0 = V * E; % yi = Q * gi???
    
    dy = [];
    
    % If dy0 has complex components, we separate each eigenvector
    % into real and complex parts and treat as two distinct vectors
    if isreal(dy0) == 0
        ij = 1;
        jj = 0;
        while size(dy, 2) <= k && ij <= k
            if isreal(dy0(:, ij)) == 0 && norm(real(dy0(:, ij))) > 0
                dy(:, jj + 1) = real(dy0(:, ij));
                jj = size(dy, 2);
                if ij <= k
                    dy(:, jj + 1) = abs(imag(dy0(:, ij)) * sqrt(1));
                    jj = size(dy, 2);
                    if ij < k
                        if dy0(:, ij) == conj(dy0(:, ij + 1))
                            ij = ij + 2;
                        else
                            ij = ij + 1;
                        end
                    end
                end
            else
                dy(:, jj + 1) = dy0(:, ij);
                ij = ij + 1;
                jj = size(dy, 2);
            end
        end
    else
        dy = dy0;
    end
end