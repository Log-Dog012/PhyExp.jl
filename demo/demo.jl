# PhyExp.jl 端到端演示：仿"密里根油滴实验"数据处理流程
#
# 运行：julia --project=.. demo.jl   （在 PhyExp.jl/demo 目录下）
# 或：  julia --project=PhyExp.jl demo/demo.jl

using PhyExp
using Unitful
using Measurements
using Statistics
using Plots

Plots.default(fontfamily = "simhei")   # GR 后端可用的中文字体（黑体）；其他平台可换 "SimSun" 等

println("═"^64)
println("  PhyExp.jl 演示：油滴实验数据处理（带不确定度 + 单位）")
println("═"^64)

# ---------- 1. 构造带不确定度的物理量 ----------
println("\n◆ 1. 带不确定度的物理量（原生 ± 语法 + 单位）")
t = (2.12 ± 0.05)u"s"          # 下落时间
d = (6.00 ± 0.02)u"mm"         # 下落距离
v = d / t                       # 速度，误差自动传播、单位自动运算
println("  下落时间 t = ", roundsig(t; str=true))
println("  下落距离 d = ", roundsig(d; str=true))
println("  速度   v = d/t = ", roundsig(v; str=true))

# ---------- 2. 多次测量 → A类/B类/合成 ----------
println("\n◆ 2. 多次测量的不确定度分析")
# 同一油滴 5 次下落时间（秒）
t_list = [2.12, 2.15, 2.09, 2.13, 2.11]
uA = ua(t_list)                          # A 类不确定度
uB = ub(0.02)                            # 仪器误差限 0.02 s → B 类（均匀分布，÷√3）
u_tot = combine(uA, uB)                  # 合成
println("  A 类不确定度 uA = ", roundsig(uA; str=true))
println("  B 类不确定度 uB = ", roundsig(uB; str=true))
println("  合成不确定度 uc = ", roundsig(u_tot; str=true))

# 一键汇总（Python 版的 多次测量结果）
t̄ = mean_measurement(t_list; ub=0.02)
println("  测量结果 t = ", roundsig(t̄; str=true), " s")

# ---------- 3. 异常值检验 ----------
println("\n◆ 3. 异常值检验（格拉布斯）")
data = [10.1, 10.2, 10.0, 10.3, 10.1, 15.7, 10.2]
mask = outliers(data)
println("  离群点标记：", mask)
println("  剔除后数据：", reject(data))

# ---------- 4. 一元线性回归（带不确定度） ----------
println("\n◆ 4. 一元线性回归（加权最小二乘，结果带不确定度）")
# 弗兰克-赫兹实验：栅极电压 U 与板极电流 I 的线性段
U = [45.0, 50.0, 55.0, 60.0, 65.0, 70.0]u"V"
I = ([2.35, 2.78, 3.21, 3.64, 4.08, 4.50] .± 0.05)u"μA"

fit = fit_line(U, I)
println(fit)                              # 自动打印拟合报告
println()
println("  斜率 b = ", roundsig(slope(fit); str=true))
println("  截距 a = ", roundsig(intercept(fit); str=true))
println("  相关系数 r = ", round(correlation(fit); digits=4))

# 预测/插值：U = 52 V 时的电流（误差传播自动完成）
I_pred = intercept(fit) + slope(fit) * 52u"V"
println("  U = 52 V 时预测电流 = ", roundsig(I_pred; str=true))

# ---------- 5. York 正交回归（x、y 都有不确定度） ----------
println("\n◆ 5. York 正交回归（x、y 均有不确定度）")
xf = collect(1.0:8.0)
yf = 1.2 .+ 0.8 .* xf
fit_y = fit_line_york(xf .± 0.05, yf .± 0.08)
println("  斜率 b = ", roundsig(slope(fit_y); str=true))
println("  截距 a = ", roundsig(intercept(fit_y); str=true))

# ---------- 6. 报告表格（Markdown / LaTeX） ----------
println("\n◆ 6. 实验报告表格")
# 列名用裸名 + 带单位数组，库自动在表头附单位、自动修约测量值
md = mdtable(Dict("栅极电压 U" => U,
                  "板极电流 I" => I,
                  "ΔI" => uncertval.(I));
             title="弗兰克-赫兹实验数据", caption="线性段数据", label="tab:fh")
println(md)

lt = latextable(Dict("栅极电压 U" => U, "板极电流 I" => I);
                caption="弗兰克-赫兹实验数据")
println("\n（LaTeX 表格已生成，可写入报告源文件）")

# ---------- 7. 绘图 ----------
println("\n◆ 7. 回归拟合图（保存为 fit_demo.png）")
p = plotfit(U, I; title="弗兰克-赫兹实验：I-U 关系", xlabel="栅极电压 U", ylabel="板极电流 I")
Plots.savefig(p, joinpath(@__DIR__, "fit_demo.png"))
println("  已保存 demo/fit_demo.png")

println("\n" * "═"^64)
println("  演示完成。更多用法见 README.md 与源码 docstring。")
println("═"^64)
