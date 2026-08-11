module PhyExp

using Unitful
using Measurements
using Statistics
import Plots

import Unitful: ustrip, unit, Quantity, @u_str
import Measurements: ±, measurement, value, uncertainty

export @u_str, ±, Quantity, unit, ustrip, measurement, value, uncertainty

include("rounding.jl")
include("uncert.jl")
include("regression.jl")
include("report.jl")
include("plotting.jl")

# ---- 中文别名（与 Python 版 phyexp API 一一对应，方便迁移）----
const 求A类不确定度 = ua
const 仪器误差限转B类不确定度 = ub
const 不确定度合成 = combine
const 一元线性回归 = fit_line
const 提取标称值 = nominal
const 提取不确定度 = uncertval
const 修约 = roundsig
const 相对误差 = relerr
const 多次测量结果 = mean_measurement
const 斜率 = slope
const 截距 = intercept
const 相关系数 = correlation

export ua, ub, combine, zscore, outliers, reject, grubbs
export roundsig, nominal, uncertval, relerr, mean_measurement
export fit_line, fit_line_york, slope, intercept, correlation, LinearFitResult
export mdtable, latextable
export plotfit, errorbarplot
export 求A类不确定度, 仪器误差限转B类不确定度, 不确定度合成, 一元线性回归,
       提取标称值, 提取不确定度, 修约, 相对误差, 多次测量结果, 斜率, 截距, 相关系数

end
