# 不确定度计算：A 类、B 类、合成与异常值检验
#
# 对应 Python 版的 AB_uncert.py。

"""
    ua(x)

计算 **A 类不确定度**（测量列标准误差，贝塞尔公式）：

    u_A = s / √n

其中 s 为样本标准差（自由度 n−1），n 为测量次数。

# 支持输入
- 普通数值数组：`ua([1.01, 1.02, 1.00])`
- 整体带单位：`ua([1.01, 1.02, 1.00]u"V")` → 返回带单位的 A 类不确定度
- 元素带单位：`ua([1.01u"V", 1.02u"V", 1.00u"V"])`
- 元素为测量值（`Measurement`）：自动只取标称值

至少需要 2 个测量值，否则报错。
"""
function ua end

function ua(x::AbstractArray{<:Real})
    length(x) < 2 && throw(ArgumentError("A 类不确定度至少需要 2 个测量值，实际 $(length(x)) 个。"))
    std(x) / sqrt(length(x))
end

function ua(x::AbstractArray{<:Measurement})
    ua(value.(x))
end

function ua(q::Quantity{<:AbstractArray})
    ua(ustrip(q)) * unit(q)
end

function ua(x::AbstractArray{<:Quantity})
    u = unit(first(x))
    ua(ustrip.(x)) * u
end

ua(x::Number) = throw(ArgumentError("A 类不确定度需要一组测量值，输入为标量。"))

# 支持的分布类型及其包含因子 K
# 均匀分布：u_B = Δ/√3；正态分布（3σ 原则）：u_B = Δ/3
function _包含因子(dist)
    if dist in (:uniform, "均匀")
        return √3
    elseif dist in (:normal, "正态")
        return 3.0
    end
    throw(ArgumentError("未知分布类型 $(dist)，目前支持 :uniform(均匀) 与 :normal(正态)。"))
end

"""
    ub(Δ; dist=:uniform, K=nothing)

把仪器误差限 Δ 转换为 **B 类不确定度**：

- `dist=:uniform`（默认）：均匀分布，u_B = Δ/√3
- `dist=:normal`：正态分布（3σ 原则），u_B = Δ/3
- 显式给 `K`：u_B = Δ/K（此时忽略 `dist`）

Δ 可以带单位（`Quantity`），结果自动保留单位。
"""
function ub(Δ; dist=:uniform, K=nothing)
    k = K === nothing ? _包含因子(dist) : K
    k <= 0 && throw(ArgumentError("分布因子 K 必须为正数。"))
    Δ / k
end

"""
    combine(u1, u2, u3...)

**不确定度合成**（各分量平方和的平方根）：

    u = √(u₁² + u₂² + u₃² + …)

输入可以是：
- 不确定度数值（可直接带单位）
- 测量值（`Measurement` / `Quantity{Measurement}`）：自动提取其标准不确定度参与合成
"""
function combine(a, b, cs...)
    us = (_combval(a), _combval(b), map(_combval, cs)...)
    sqrt(sum(u^2 for u in us))
end

# 合成时取值：Measurement 取不确定度，普通数值/带单位数值直接用原值
_combval(x::Measurement) = uncertainty(x)
_combval(q::Quantity{<:Measurement}) = Quantity(uncertainty(q.val), unit(q))
_combval(x) = x

# ---------- 异常值检验 ----------

"""
    zscore(x)

标准化得分：z_i = (x_i − x̄) / s。支持带单位数组与 `Measurement` 数组。
"""
zscore(v::AbstractArray{<:Real}) = (v .- mean(v)) ./ std(v)
zscore(v::AbstractArray{<:Measurement}) = zscore(value.(v))
zscore(q::Quantity{<:AbstractArray}) = zscore(ustrip(q))
zscore(v::AbstractArray{<:Quantity}) = zscore(ustrip.(nominal.(v)))

"""
    outliers(x; method=:grubbs, alpha=0.05)

返回布尔掩码，标记离群值（`true` 表示该点为离群值）。

- `method=:grubbs`（默认）：格拉布斯检验，显著性水平 `alpha`（0.05 或 0.01）
- `method=:zscore`：3σ 原则，|z| > 3

配合 `reject` 使用可直接剔除。
"""
function outliers end

outliers(x; method=:grubbs, alpha=0.05) =
    method == :grubbs ? grubbs_mask(x, alpha) : zscore(_vals(x)) .> 3

"""
    reject(x; kwargs...)

剔除离群值后返回剩余数组。参数同 [`outliers`](@ref)。
"""
reject(x; kwargs...) = x[.!outliers(x; kwargs...)]

"""
    grubbs(x; alpha=0.05)

格拉布斯检验：返回 `true` 表示数据列中存在离群值。
"""
grubbs(x; alpha=0.05) = any(outliers(x; method=:grubbs, alpha=alpha))

_vals(v::AbstractArray{<:Real}) = Float64.(v)
_vals(v::AbstractArray{<:Measurement}) = Float64.(value.(v))
_vals(q::Quantity{<:AbstractArray}) = Float64.(ustrip(q))
_vals(v::AbstractArray{<:Quantity}) = Float64.(ustrip.(nominal.(v)))

# 迭代格拉布斯检验：每次剔除最可疑的点，重新检验，直到无离群值
function grubbs_mask(x, alpha=0.05)
    vals = _vals(x)
    n = length(vals)
    n < 3 && throw(ArgumentError("格拉布斯检验至少需要 3 个数据点，实际 $n 个。"))
    keep = trues(n)
    y = collect(vals)
    idx = collect(1:n)
    while length(y) >= 3
        m = mean(y)
        s = std(y)
        s == 0 && break
        z = abs.((y .- m) ./ s)
        i = argmax(z)
        if z[i] > _grubbs_crit(length(y), alpha)
            keep[idx[i]] = false
            deleteat!(y, i)
            deleteat!(idx, i)
        else
            break
        end
    end
    .!keep
end

# 格拉布斯临界值表（单侧检验），n = 3:30
const _GRUBBS_05 = [1.155, 1.481, 1.715, 1.887, 2.020, 2.126, 2.215, 2.290, 2.355,
                    2.412, 2.462, 2.507, 2.549, 2.585, 2.624, 2.657, 2.681, 2.709,
                    2.733, 2.758, 2.781, 2.802, 2.822, 2.841, 2.859, 2.876, 2.893, 2.908]
const _GRUBBS_01 = [1.155, 1.496, 1.764, 1.973, 2.139, 2.274, 2.387, 2.482, 2.564,
                    2.636, 2.699, 2.755, 2.806, 2.852, 2.894, 2.932, 2.968, 3.001,
                    3.031, 3.060, 3.087, 3.112, 3.135, 3.157, 3.178, 3.199, 3.218, 3.236]

function _grubbs_crit(n::Int, alpha::Float64)
    if n <= 30
        tbl = alpha < 0.05 ? _GRUBBS_01 : _GRUBBS_05
        return tbl[n-2]
    end
    # n > 30：用 30 时的临界值近似（并提示）
    @warn "格拉布斯临界值表仅覆盖 n ≤ 30，n=$n 时采用 n=30 的临界值近似。"
    alpha < 0.05 ? _GRUBBS_01[end] : _GRUBBS_05[end]
end
