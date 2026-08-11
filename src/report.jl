# 实验报告表格输出：Markdown / LaTeX
#
# 对应 Python 版的 notebook.py table()，改进点：
# - 测量值列自动按有效数字规范修约后显示
# - 带单位列自动在表头附加单位（如 "电压 U / V"）
# - 支持 LaTeX 表格输出（实验报告常用）

# 单元格文本：测量值自动修约，带单位只显示数值（单位进表头）
_cell(x::Measurement) = roundsig(x; str=true)
_cell(q::Quantity{<:Measurement}) = roundsig(q.val; str=true)
_cell(q::Quantity{<:Real}) = string(ustrip(q))
_cell(x) = string(x)

# 列单位（字符串形式），无单位返回 ""
function _colunit(v::AbstractVector)
    T = eltype(v)
    T <: Quantity ? string(unit(first(v))) : ""
end
_colunit(q::Quantity{<:AbstractVector}) = string(unit(q))

function _table_rows(cols::Vector, index::Bool, n::Int)
    rows = String[]
    for i in 1:n
        cells = [_cell(col[i]) for col in cols]
        prefix = index ? "$i | " : ""
        push!(rows, "| $prefix" * join(cells, " | ") * " |")
    end
    rows
end

"""
    mdtable(data::AbstractDict; index=true, title=nothing, caption=nothing, label=nothing)

把列式数据生成 **Markdown 表格**字符串（可用 `Markdown.parse` 渲染）。

- `data`：`Dict`（或有序字典），键为列名，值为列数据（可为 `Measurement` /
  带单位数组，自动修约显示、表头附单位）
- `index`：是否添加序号列
- `title`：表题（自动加"表："前缀）；`caption`：表注；`label`：交叉引用标签
"""
function mdtable(data::AbstractDict; index=true, title=nothing, caption=nothing, label=nothing)
    headers = collect(keys(data))
    cols = collect(values(data))
    n = length(first(cols))
    all(c -> length(c) == n, cols) ||
        throw(DimensionMismatch("所有列长度必须相同。"))

    hdrs = [string(h) * (isempty(_colunit(c)) ? "" : " / " * _colunit(c))
            for (h, c) in zip(headers, cols)]

    header_row = "| " * (index ? "Index | " : "") * join(hdrs, " | ") * " |"
    sep_row = "| " * (index ? "--- | " : "") * join(["---" for _ in hdrs], " | ") * " |"
    rows = _table_rows(cols, index, n)

    out = join(vcat(header_row, sep_row, rows), "\n")
    title !== nothing && (out = "表：$title\n\n" * out)
    caption !== nothing && (out *= "\n\n*$caption*")
    label !== nothing && (out *= "\n\n<!-- $label -->")
    out
end

"""
    latextable(data::AbstractDict; index=true, caption=nothing, label=nothing)

把列式数据生成 **LaTeX 表格**（`tabular`）字符串，测量值自动修约并以
`1.234 ± 0.012` 形式输出，可直接放入报告源文件。
"""
function latextable(data::AbstractDict; index=true, caption=nothing, label=nothing)
    headers = collect(keys(data))
    cols = collect(values(data))
    n = length(first(cols))
    all(c -> length(c) == n, cols) ||
        throw(DimensionMismatch("所有列长度必须相同。"))

    ncol = length(headers) + (index ? 1 : 0)
    spec = join(["c" for _ in 1:ncol], " ")
    hdrs = [string(h) * (isempty(_colunit(c)) ? "" : "~(\\mathrm{" * _colunit(c) * "})")
            for (h, c) in zip(headers, cols)]
    if index
        hdrs = vcat("序号", hdrs)
    end

    lines = String[]
    push!(lines, "\\begin{table}[htbp]")
    caption !== nothing && push!(lines, "  \\caption{$caption}")
    label !== nothing && push!(lines, "  \\label{$label}")
    push!(lines, "  \\centering")
    push!(lines, "  \\begin{tabular}{$spec}")
    push!(lines, "    \\hline")
    push!(lines, "    " * join(hdrs, " & ") * " \\\\")
    push!(lines, "    \\hline")
    for i in 1:n
        cells = [_latex_cell(col[i]) for col in cols]
        index && pushfirst!(cells, "$i")
        push!(lines, "    " * join(cells, " & ") * " \\\\")
    end
    push!(lines, "    \\hline")
    push!(lines, "  \\end{tabular}")
    push!(lines, "\\end{table}")
    join(lines, "\n")
end

_latex_cell(m::Measurement) = replace(roundsig(m; str=true), "±" => "\\pm")
_latex_cell(q::Quantity{<:Measurement}) = _latex_cell(q.val)
_latex_cell(q::Quantity{<:Real}) = string(ustrip(q))
_latex_cell(x) = string(x)
