function qrupdate_gs(a::AbstractMatrix, q_in::AbstractMatrix, r_in::AbstractMatrix)
    # Incremental thin QR update via double modified Gram-Schmidt: extends an
    # existing factorisation of a's first n_prev columns to all n_cols columns.
    # q_in may have fewer rows than a (the Hessenberg gains a row each Arnoldi
    # step); the new bottom row is zero for previously-factorised columns
    # (guaranteed by the Hessenberg structure of a's callers), so padding it
    # with zeros is exact, not an approximation.
    n_rows, n_cols = size(a)
    n_prev = size(q_in, 2)
    T = eltype(a)

    q_out = zeros(T, n_rows, n_cols)
    r_out = zeros(T, n_cols, n_cols)

    n_q_rows = size(q_in, 1)
    q_out[1:n_q_rows, 1:n_prev] .= q_in
    r_out[1:n_prev, 1:n_prev] .= r_in

    for j = (n_prev+1):n_cols
        w = a[:, j]
        for _ = 1:2   # two passes of MGS (reorthogonalisation)
            for i = 1:(j-1)
                proj = dot(view(q_out, :, i), w)
                r_out[i, j] += proj
                w .-= proj .* view(q_out, :, i)
            end
        end
        r_out[j, j] = norm(w)
        q_out[:, j] = w / r_out[j, j]
    end

    return q_out, r_out
end
