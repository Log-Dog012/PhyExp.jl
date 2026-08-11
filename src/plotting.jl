# 绘图封装（基于 Plots.jl）
#
# 对应 Python 版的 plt.py，改进点：
# - `plotfit` 一键出"数据点误差棒 + 回归直线"图
# - 自动从测量值提取误差棒、自动标注单位

# 取元素视图：整体带单位的 Quantity 取内部数组，普通数组原样
_elements(q::Quantity) = q.val
_elements(v::AbstractArray) = v

# 数组的单位（无单位返回 nothing）
_unitof(v::AbstractArray) = eltype(v) <: Quantity ? unit(first(v)) : nothing
_unitof(q::Quantity) = unit(q)

"""
    errorbarplot(x, y; xerr=nothing, kwargs...) -> Plots.Plot

带误差棒的散点图。`x`、`y` 元素为测量值（`Measurement` / `Quantity{Measurement}`）
时自动提取误差棒；也可用 `xerr`、`yerr` 显式指定。

返回 `Plots.Plot` 对象，可继续美化或 `savefig`。
"""
function errorbarplot(x, y; xerr=nothing, yerr=nothing, xlabel="自变量", ylabel="因变量", kwargs...)
    xe = _elements(x)
    ye = _elements(y)
    xv = Float64.(ustrip.(nominal.(xe)))
    yv = Float64.(ustrip.(nominal.(ye)))
    xer = xerr === nothing ? Float64.(ustrip.(uncertval.(xe))) : Float64.(xerr)
    yer = yerr === nothing ? Float64.(ustrip.(uncertval.(ye))) : Float64.(yerr)
    xu = _unitof(x); yu = _unitof(y)
    xl = xu === nothing ? xlabel : "$xlabel / $xu"
    yl = yu === nothing ? ylabel : "$ylabel / $yu"
    Plots.scatter(xv, yv; xerror=xer, yerror=yer, xlabel=xl, ylabel=yl, kwargs...)
end

"""
    plotfit(x, y; title="带不确定度的一元线性回归图", xlabel="自变量", ylabel="因变量", kwargs...)

一键绘制**带不确定度的一元线性回归图**：数据点（含误差棒）+ 拟合直线。

内部调用 [`fit_line`](@ref)，单位自动标注在坐标轴上。
返回 `Plots.Plot` 对象，可用 `Plots.savefig(p, "fit.png")` 保存。
"""
function plotfit(x, y; title="带不确定度的一元线性回归图", xlabel="自变量", ylabel="因变量",
                 linecolor=:red, markercolor=:steelblue, marker=:circle, kwargs...)
    fit = fit_line(x, y)
    xe = _elements(x)
    ye = _elements(y)
    xv = Float64.(ustrip.(nominal.(xe)))
    yv = Float64.(ustrip.(nominal.(ye)))
    σy = Float64.(ustrip.(uncertval.(ye)))
    xu = _unitof(x); yu = _unitof(y)
    xl = xu === nothing ? xlabel : "$xlabel / $xu"
    yl = yu === nothing ? ylabel : "$ylabel / $yu"

    lo, hi = extrema(xv)
    xgrid = range(lo, hi; length=100)
    a = ustrip(nominal(intercept(fit)))     # 数值化后画图（单位已标注在坐标轴）
    b = ustrip(nominal(slope(fit)))
    ygrid = a .+ b .* xgrid

    p = Plots.scatter(xv, yv; yerror=σy, label="数据点", marker=marker,
                      markercolor=markercolor, xlabel=xl, ylabel=yl, title=title, kwargs...)
    Plots.plot!(p, xgrid, ygrid; label="回归线", color=linecolor)
    p
end
