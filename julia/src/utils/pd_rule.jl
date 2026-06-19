function pd_rule(
    m::Int,
    n::Int,
    m_initial::Int,
    m_min::Int,
    m_max::Int,
    m_step::Int,
    res::AbstractVector,
    iter::Int,
    alpha_p::Real,
    alpha_d::Real,
)::Tuple{Int,Int}
    if iter > 3
        mj =
            m + ceil(
                Int,
                alpha_p * (res[iter] / res[iter-1]) +
                alpha_d * ((res[iter] - res[iter-2]) / (2 * res[iter-1])),
            )
    elseif iter > 2
        mj = m + ceil(Int, alpha_p * (res[iter] / res[iter-1]))
    else
        mj = m_initial
    end

    if mj < m_min
        m_initial += m_step
        mj = m_initial
    end

    if mj > m_max
        mj = m_max
    end

    if m_initial + m_step > n
        mj = n
    end

    return mj, m_initial
end
