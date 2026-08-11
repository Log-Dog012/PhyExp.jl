# 一元线性回归：加权最小二乘 与 York 正交回归
#
# 对应 Python 版的 SLR.py，改进点：
# - 回归结果（斜率/截距）**自带不确定度**（Python 版只返回裸数值）
# - 提供 x、y 都带不确定度时的 York 正交回归
# - 无不确定度数据时自动退化为普通最小二乘，并用残差估计参数不确定度

"""
    LinearFitResult

一元线性回归结果。

字段：
- `a`：截距（带不确定度，与 `y` 同单位）
- `b`：斜率（带不确定度，单位 = `y` 的单位 / `x` 的单位）
- `r`：线性相关系数
- `chi2`：加权残差平方和 χ²
- `dof`：自由度（n − 2）

可用 [`slope`](@ref) / [`intercept`](@ref) / [`correlation`](@ref) 提取，
直接用 `roundsig` 格式化输出。
"""
struct LinearFitResult{A,B}
    a::A
    b::B
    r::Float64
    chi2::Float64
    dof::Int
end

slope(r::LinearFitResult) = r.b
intercept(r::LinearFitResult) = r.a
correlation(r::LinearFitResult) = r.r

function Base.show(io::IO, r::LinearFitResult)
    print(io, "一元线性回归结果\n")
    print(io, "  截距 a = ", roundsig(r.a; str=true), "\n")
    print(io, "  斜率 b = ", roundsig(r.b; str=true), "\n")
    print(io, "  相关系数 r = ", round(r.r; digits=4), "\n")
    print(io, "  χ²/ν = ", round(r.chi2 / max(r.dof, 1); digits=3), "  (ν = ", r.dof, ")")
end

# ---------- 数据分解 ----------

# 统一提取 (数值向量, 单位或nothing, 不确定度向量)
_decomp(v::AbstractVector{<:Real}) = (Float64.(v), nothing, zeros(Float64, length(v)))
function _decomp(v::AbstractVector{<:Measurement})
    (Float64.(value.(v)), nothing, Float64.(uncertainty.(v)))
end
function _decomp(v::AbstractVector{<:Quantity{<:Measurement}})
    u = unit(first(v))
    (Float64.(ustrip.(nominal.(v))), u, Float64.(ustrip.(uncertval.(v))))
end
function _decomp(v::AbstractVector{<:Quantity{<:Real}})
    u = unit(first(v))
    (Float64.(ustrip.(v)), u, zeros(Float64, length(v)))
end
function _decomp(q::Quantity{<:AbstractVector})
    u = unit(q)
    m = q.val
    if m isa AbstractArray{<:Measurement}
        (Float64.(value.(m)), u, Float64.(uncertainty.(m)))
    else
        (Float64.(m), u, zeros(Float64, length(m)))
    end
end

# ---------- 加权最小二乘 ----------

"""
    fit_line(x, y)

带不确定度的一元线性回归 `y = a + b·x`（加权最小二乘）。

- `x`、`y`：长度相等的测量数组。元素可以是普通数值、`Measurement`、
  带单位（`Quantity`）或单位+不确定度组合。
- `y` 带不确定度时，权重取 `w_i = 1/σᵧᵢ²`；返回的 `a`、`b` 带标准不确定度。
- `y` 无不确定度时退化为普通最小二乘，参数不确定度由残差估计。

返回 [`LinearFitResult`](@ref)。
"""
function fit_line(x, y)
    length(x) == length(y) || throw(DimensionMismatch("x 与 y 的长度必须相等。"))
    xv, xu, _ = _decomp(x)
    yv, yu, σy = _decomp(y)
    n = length(xv)

    if all(σy .== 0)
        w = ones(n)                       # 无不确定度 → 等权重
        residual_based = true
    else
        w = 1.0 ./ σy.^2
        residual_based = false
    end

    S = sum(w)
    Sx = sum(w .* xv)
    Sy = sum(w .* yv)
    Sxx = sum(w .* xv .^ 2)
    Sxy = sum(w .* xv .* yv)
    D = S * Sxx - Sx^2

    b = (S * Sxy - Sx * Sy) / D
    a = (Sy - b * Sx) / S

    yfit = a .+ b .* xv
    resid = yv .- yfit
    chi2 = residual_based ? sum(resid .^ 2) : sum(w .* resid .^ 2)
    dof = n - 2

    if residual_based
        s2 = dof > 0 ? sum(resid .^ 2) / dof : 0.0
        σb = sqrt(s2 / D)
        σa = sqrt(s2 * Sxx / D)
    else
        σb = sqrt(S / D)
        σa = sqrt(Sxx / D)
    end

    r = cor(xv, yv)

    a_out = a ± σa
    b_out = b ± σb
    if yu !== nothing                       # 装回单位
        a_out = a_out * yu
        b_out = b_out * (xu === nothing ? yu : yu / xu)
    end

    LinearFitResult(a_out, b_out, r, chi2, dof)
end

# ---------- York 正交回归 ----------

"""
    fit_line_york(x, y; maxiter=50, tol=1e-12)

**York 正交回归**：当 `x` 与 `y` **都**带不确定度时的一元线性拟合
（迭代加权，同时考虑两轴不确定度，假定 x、y 独立测量）。

- 若 x、y 均无不确定度，自动退化为 [`fit_line`](@ref)。
- 返回 [`LinearFitResult`](@ref)，其中 `a`、`b` 带标准不确定度。

参考：York et al., *Am. J. Phys.* 72, 367 (2004)。
"""
function fit_line_york(x, y; maxiter=50, tol=1e-12)
    length(x) == length(y) || throw(DimensionMismatch("x 与 y 的长度必须相等。"))
    xv, xu, σx = _decomp(x)
    yv, yu, σy = _decomp(y)
    n = length(xv)

    if all(σx .== 0) && all(σy .== 0)
        return fit_line(x, y)
    end

    # σ=0 的元素给一个相对小的下限，避免权重无穷大 / 下溢；
    # 单轴无不确定度时 York 迭代自动退化为该轴的加权最小二乘
    scale = maximum(abs.(vcat(xv, yv)))
    σfloor = max(1e-12 * scale, 1e-8)
    sx = max.(σx, σfloor)
    sy = max.(σy, σfloor)

    Wx = 1.0 ./ sx .^ 2
    Wy = 1.0 ./ sy .^ 2
    α = sqrt.(Wx .* Wy)

    # 初始斜率：普通最小二乘（取标称值，York 迭代在纯数值上做）
    xm = xv; ym = yv
    b = nominal(slope(fit_line(xm, ym)))

    # York 迭代（注意：for 循环作用域不向外部泄漏，先声明循环外用到的变量）
    W = zeros(n)
    U = zeros(n)
    V = zeros(n)
    x̄w = 0.0
    ȳw = 0.0
    for _ in 1:maxiter
        W = Wx .* Wy ./ (Wx .+ b^2 .* Wy)          # 假定 x、y 独立（相关系数 r=0）
        x̄w = sum(Wx .* xm) / sum(Wx)
        ȳw = sum(Wy .* ym) / sum(Wy)
        U = xm .- x̄w
        V = ym .- ȳw
        β = W .* (U ./ Wx .+ b .* V ./ Wy)
        b_new = sum(W .* β .* V) / sum(W .* β .* U)
        abs(b_new - b) < tol && (b = b_new; break)
        b = b_new
    end

    a = ȳw - b * x̄w
    σb = sqrt(1 / sum(W .* U .^ 2))
    σa = sqrt(1 / sum(W) + x̄w^2 * σb^2)

    # 相关系数与 χ²
    r = cor(xm, ym)
    chi2 = sum(W .* (ym .- (a .+ b .* xm)) .^ 2)

    a_out = a ± σa
    b_out = b ± σb
    if yu !== nothing
        a_out = a_out * yu
        b_out = b_out * (xu === nothing ? yu : yu / xu)
    end
    LinearFitResult(a_out, b_out, r, chi2, n - 2)
end
