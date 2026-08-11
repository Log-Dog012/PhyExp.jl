# 有效数字修约与测量值提取工具
#
# 对应 Python 版的 sigfigs.py / error.py。

"""
    nominal(x)

提取测量值的标称值（不含不确定度）。

- 普通数值：原样返回
- `Measurement`：返回 `value(x)`
- 带单位（`Quantity`）：去掉单位返回数值；带单位且带不确定度时返回带单位的标称值
- 数组：逐元素提取
"""
nominal(x::Real) = x
nominal(m::Measurement) = value(m)
nominal(q::Quantity{<:Real}) = ustrip(q)
nominal(q::Quantity{<:Measurement}) = Quantity(value(q.val), unit(q))
nominal(v::AbstractArray) = nominal.(v)

"""
    uncertval(x)

提取测量值的标准不确定度；若输入本身没有不确定度信息，返回 0。

- `Measurement`：返回 `uncertainty(x)`
- 带单位且带不确定度：返回带单位的标准不确定度
- 带单位但无不确定度：返回带单位的 0
"""
uncertval(x::Real) = 0.0
uncertval(m::Measurement) = uncertainty(m)
uncertval(q::Quantity{<:Real}) = 0.0 * unit(q)
uncertval(q::Quantity{<:Measurement}) = Quantity(uncertainty(q.val), unit(q))
uncertval(v::AbstractArray) = uncertval.(v)

# ---------- 有效数字修约 ----------

"""
    roundsig(m; str=false)

按大学物理实验的修约规范处理测量值：

1. 不确定度取 **1 位**有效数字；若其首位数字为 **1 或 2**，取 **2 位**；
2. 测量值修约到与不确定度的末位**对齐**；
3. 数值过小（末位在小数点后 5 位以上）或过大（末位在百位以上）时，
   自动改用科学计数法 `(a ± b)×10ᵏ` 输出。

输入可以是 `Measurement`、带单位（`Quantity`）或普通数值。
`str=true` 时返回格式化字符串（例如 `"9.84 ± 0.05"`），否则返回修约后的测量值对象。
"""
function roundsig end

function roundsig(m::Measurement; str=false)
    v, e = m.val, m.err
    iszero(e) && return str ? string(v) : m

    k = floor(Int, log10(abs(e)))
    firstdigit = floor(Int, abs(e) / 10.0^k)
    nd = firstdigit <= 2 ? 2 : 1            # 不确定度取几位有效数字
    e_r = round(e, sigdigits=nd)
    k_r = floor(Int, log10(abs(e_r)))
    pos = -k_r + (nd - 1)                   # 不确定度末位对应的小数位

    # 数值过小 / 过大 → 科学计数法（k 取测量值数量级，使标称值落在 1~10）
    if pos > 4 || pos < -2
        k = iszero(v) ? k_r : floor(Int, log10(abs(v)))
        vp = v / 10.0^k
        ep = e_r / 10.0^k
        kk = floor(Int, log10(abs(ep)))
        p2 = max(-kk + (nd - 1), 0)         # 不确定度末位对应 vp 的小数位
        vpp = round(vp, digits=p2)
        epp = round(ep, digits=p2)
        s = "($(_fmt(vpp, p2)) ± $(_fmt(epp, p2)))×10$(_superscript(k))"
        return str ? s : (vpp * 10.0^k) ± (epp * 10.0^k)
    end

    v_r = round(v, digits=pos)
    e_r2 = round(e_r, digits=pos)
    if str
        return "$(_fmt(v_r, pos)) ± $(_fmt(e_r2, pos))"
    end
    v_r ± e_r2
end

roundsig(q::Quantity{<:Measurement}; str=false) =
    str ? "($(roundsig(q.val; str=true))) $(unit(q))" : Quantity(roundsig(q.val), unit(q))
roundsig(q::Quantity{<:Real}; str=false, sigdigits=2) =
    str ? string(round(ustrip(q); sigdigits=sigdigits)) : ustrip(q)

"""
    roundsig(x::Real; sigdigits=2)

对没有不确定度的普通数值按指定位有效数字修约（用于 ua/ub 等不确定度数值的显示）。
"""
roundsig(x::Real; str=false, sigdigits=2) =
    str ? string(round(x; sigdigits=sigdigits)) : round(x; sigdigits=sigdigits)

# 按指定小数位格式化：小数位 ≤ 0 时输出整数形式
_fmt(x::Real, pos::Int) = pos > 0 ? string(round(x, digits=pos)) : string(Int(round(x, digits=pos)))

_superscripts = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
_superscript(k::Integer) = (k < 0 ? "⁻" : "") * join(_superscripts[d+1] for d in reverse(digits(abs(k))))

# ---------- 相对误差 ----------

"""
    relerr(测量值, 真值)

相对误差：(测量值 - 真值) / 真值。单位会自动约掉。
"""
relerr(x, x0) = (x - x0) / x0

"""
    relerr(m::Measurement)

测量值的**相对不确定度**：σ / |标称值|。
"""
relerr(m::Measurement) = abs(m.err / m.val)
relerr(q::Quantity{<:Measurement}) = abs(q.val.err / q.val.val)

# ---------- 多次测量综合 ----------

"""
    mean_measurement(x; ub=0.0)

多次测量结果汇总：返回 `平均值 ± (A类 ⊕ B类) 合成不确定度`。

- `x`：测量值数组，可带单位（`Quantity` 或元素带单位的数组）
- `ub`：仪器 B 类不确定度（数值，单位需与 `x` 一致）
"""
function mean_measurement(x; ub=0.0)
    xv, xu = _decompose(x)
    m = mean(xv)
    σ = combine(ua(xv), _as_ub(ub, xu))
    out = m ± σ
    xu === nothing ? out : out * xu
end

_as_ub(ub::Real, ::Nothing) = ub
_as_ub(ub::Real, ::Any) = ub
_as_ub(ub::Quantity, xu) = ustrip(ub)      # 带单位的 B 类不确定度，统一取数值
_as_ub(ub::Measurement, xu) = value(ub)

# 把输入统一为 (数值向量, 单位或 nothing)
_decompose(v::AbstractVector{<:Real}) = (Float64.(v), nothing)
_decompose(v::AbstractVector{<:Measurement}) = (Float64.(value.(v)), nothing)
function _decompose(v::AbstractVector{<:Quantity})
    u = unit(first(v))
    (Float64.(ustrip.(nominal.(v))), u)
end
function _decompose(q::Quantity{<:AbstractVector})
    u = unit(q)
    m = q.val
    if m isa AbstractArray{<:Measurement}
        (Float64.(value.(m)), u)
    else
        (Float64.(m), u)
    end
end
