# PhyExp.jl

**大学物理实验数据处理库（Julia 版）** — 简化数据处理，自动生成规范报告要素。

`PhyExp.jl` 是同名 Python 库 [`phyexp`](https://gitee.com/Log-Dog012/phyexp) 的 Julia 重构版，目标是更简洁、更一致、更少踩坑的 API。

---

## 与 Python 版对比：更好用在哪儿？

| 方面 | Python 版（pint + uncertainties） | **Julia 版（Unitful + Measurements）** |
|---|---|---|
| 构造带不确定度的量 | `ufloat(5.2, 0.1)` + `Q_(..., "m")` | **`(5.2 ± 0.1)u"m"`**（原生语法） |
| 多库拼接（单位+不确定度） | 200 行代码处理两个 ureg 冲突 | **多重分派，0 行胶水代码** |
| 算术运算 | 多个 isinstance 分支 | **`(1.0±0.1)m / (2.0±0.2)s` 自动传播** |
| 线性回归结果 | 只返回裸数 `a, b` | **`a ± σa`, `b ± σb` 带不确定度** |
| x、y 都带不确定度的回归 | ❌ 不支持 | ✅ **York 正交回归** |
| 异常值检验 | README 提到未实现 | ✅ **格拉布斯 + 3σ** |
| 报告表格 | Markdown only | **Markdown + LaTeX** |
| API 一致性 | 中文函数名混搭 | **英文短名 + 中文别名（兼容旧代码）** |

---

## 安装

```julia
using Pkg
Pkg.develop(path="path/to/PhyExp.jl")     # 本地开发
# 或 Pkg.add("PhyExp")                    # 发布后
```

依赖：`Unitful`、`Measurements`、`Plots`（绘图）、`Statistics`（标准库）。

---

## 核心 API 速查（与 Python 版对照）

| 功能 | Python 版（中文） | **Julia 版（英文）** | Julia 版（中文别名） |
|---|---|---|---|
| 构造带不确定度 | `ufloat(5.2, 0.1)` | **`5.2 ± 0.1`** | — |
| 带单位 | `Q_(..., "m")` | **`(... ± ...)u"m"`** | — |
| A 类不确定度 | `求A类不确定度(x)` | **`ua(x)`** | `求A类不确定度` |
| B 类不确定度 | `仪器误差限转B类不确定度(Δ)` | **`ub(Δ)`** | `仪器误差限转B类不确定度` |
| 不确定度合成 | `不确定度合成(a, b)` | **`combine(a, b)`** | `不确定度合成` |
| 多次测量综合 | `多次测量结果(x)` | **`mean_measurement(x; ub=...)`** | `多次测量结果` |
| 一元线性回归 | `一元线性回归(x, y)` | **`fit_line(x, y)`** | `一元线性回归` |
| York 回归 | — | **`fit_line_york(x, y)`** | — |
| 修约显示 | `修约(m)` | **`roundsig(m; str=true)`** | `修约` |
| 提取标称值 | `提取标称值(m)` | **`nominal(m)`** | `提取标称值` |
| 提取不确定度 | `提取不确定度(m)` | **`uncertval(m)`** | `提取不确定度` |
| 相对误差 | `相对误差(测, 真)` | **`relerr(测, 真)`** | `相对误差` |
| 异常值检验 | — | **`outliers(x; method=:grubbs)`** | — |
| 剔除异常值 | — | **`reject(x)`** | — |
| Markdown 表格 | `table(data)` | **`mdtable(data)`** | — |
| LaTeX 表格 | — | **`latextable(data)`** | — |
| 回归拟合图 | `绘制回归图(x, y)` | **`plotfit(x, y)`** | — |

英文短名面向日常使用；中文别名与 Python 版 API 一一对应，便于迁移。

---

## 端到端示例

```julia
using PhyExp
using Unitful, Measurements, Statistics, Plots

# 1. 带不确定度的物理量（原生语法 + 单位）
t = (2.12 ± 0.05)u"s"
d = (6.00 ± 0.02)u"mm"
v = d / t                                 # 速度，误差自动传播
roundsig(v; str=true)                     # "(2.83 ± 0.07) mm s⁻¹"

# 2. 多次测量：自动 A 类 + B 类 + 合成
t_list = [2.12, 2.15, 2.09, 2.13, 2.11]
mean_measurement(t_list; ub=0.02)         # 2.12 ± 0.022 s

# 3. 异常值检验
outliers([10.1, 10.2, 10.0, 15.7]; method=:grubbs)   # [false,false,false,true]

# 4. 一元线性回归（结果自带不确定度）
U = [45.0, 50.0, 55.0, 60.0, 65.0, 70.0]u"V"
I = ([2.35, 2.78, 3.21, 3.64, 4.08, 4.50] .± 0.05)u"μA"
fit = fit_line(U, I)
slope(fit)            # (0.0862 ± 0.0024) μA V⁻¹
intercept(fit)        # (-1.53 ± 0.14) μA
correlation(fit)      # 1.0
println(fit)          # 一键打印完整拟合报告

# 5. 表格与绘图
mdtable(Dict("栅极电压 U" => U, "板极电流 I" => I);
        title="弗兰克-赫兹实验数据")
p = plotfit(U, I; xlabel="栅极电压 U", ylabel="板极电流 I")
savefig(p, "fit.png")
```

完整运行：`julia --project=. demo/demo.jl`（生成 `demo/fit_demo.png`）。

---

## 设计要点

### 1. 多重分派：一种操作，多种输入

```julia
ua([1.0, 2.0, 3.0, 4.0])                  # 纯数值
ua([1.0, 2.0, 3.0, 4.0]u"V")             # 整体带单位
ua([1.0u"V", 2.0u"V", 3.0u"V", 4.0u"V"])# 元素带单位
ua([1.0±0.1, 2.0±0.2, 3.0±0.3, 4.0±0.4])# 测量值数组
```
同一个 `ua` 函数自动适配。不需要写 `isinstance` 分支。

### 2. 误差自动传播

```julia
v = d / t            # Quantity{Measurement{...}}，误差自动算
sin(t * u"Hz")       # 任意数学函数自动传播（Measurements.jl 提供）
```

### 3. 单位与不确定度原生组合

```julia
x = (5.2 ± 0.1)u"m" / (3.0 ± 0.05)u"s"^2   # 156 ± 5.5 m s⁻²
```

### 4. 修约规范符合中国大学物理实验教学惯例

- 不确定度首位 1/2 → 取 2 位有效；其他 → 取 1 位
- 测量值末位与不确定度末位对齐
- 数值过小（>4 位小数）或过大（> 百位）自动改用 `(a ± b)×10ᵏ` 科学计数法

```julia
roundsig(9.8429 ± 0.0512; str=true)    # "9.84 ± 0.05"
roundsig(1.234 ± 0.012; str=true)      # "1.234 ± 0.012"
roundsig(1.234e-6 ± 0.051e-6; str=true)# "(1.234 ± 0.051)×10⁻⁶"
```

### 5. 回归结果带不确定度

Python 版 `一元线性回归(x, y)` 只返回裸 `a, b`；Julia 版 `fit_line(x, y)` 返回 `LinearFitResult`，其中 `slope`、`intercept` 都是 `Measurement`：

```julia
fit = fit_line(x, y)
slope(fit)                 # b ± σb
intercept(fit)             # a ± σa
intercept(fit) + slope(fit) * x_new   # 自动传播预测值的不确定度
```

### 6. York 回归：x、y 都带不确定度时

```julia
fit_y = fit_line_york(x .± σx, y .± σy)
```

York et al. (2004) 迭代算法，是大学物理实验数据处理文献中推荐的标准方法。

---

## 模块结构

```
PhyExp.jl/
├── Project.toml           # 依赖声明
├── src/
│   ├── PhyExp.jl          # 主模块：导入、include、导出
│   ├── uncert.jl          # A/B 类不确定度、合成、异常值（格拉布斯 + 3σ）
│   ├── rounding.jl        # nominal/uncertval/roundsig/relerr/mean_measurement
│   ├── regression.jl      # fit_line（加权最小二乘）+ fit_line_york（York）
│   ├── report.jl          # mdtable / latextable
│   └── plotting.jl        # plotfit / errorbarplot（Plots.jl）
├── test/
│   └── runtests.jl        # 104 个测试，覆盖所有 API
└── demo/
    ├── demo.jl            # 端到端演示
    └── fit_demo.png       # 演示生成的拟合图
```

---

## 依赖

| 包 | 用途 |
|---|---|
| [`Unitful`](https://github.com/PainterQubits/Unitful.jl) | 物理单位系统 |
| [`Measurements`](https://github.com/JuliaPhysics/Measurements.jl) | 不确定度与误差传播 |
| [`Plots`](https://github.com/JuliaPlots/Plots.jl) | 绘图（GR 后端） |
| `Statistics` | 标准统计（标准差、相关系数） |
| `Test`（仅测试） | 单元测试 |

---

## 与 Python 版的不兼容点（重构的取舍）

- Python 版的 `Q_` 自定义子类、ureg 冲突处理 → Julia 版不需要（生态原生支持）
- Python 版的"输入转换装饰器" → Julia 版用多重分派
- Python 版异常值检验未实现 → Julia 版已实现
- Python 版只有 Markdown → Julia 版新增 LaTeX
- Python 版中文函数名为主 → Julia 版**英文为主 + 中文别名**（兼容调用）
- 包管理：Python `pip install -e .` → Julia `Pkg.develop(path=...)`

---

## 测试与演示

```bash
cd PhyExp.jl
julia --project=. -e 'using Pkg; Pkg.test()'      # 104 个测试
julia --project=. demo/demo.jl                     # 端到端演示
```

---

## 许可

MIT（与 Python 版 phyexp 一致）。