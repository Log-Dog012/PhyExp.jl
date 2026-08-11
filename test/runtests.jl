using PhyExp
using Unitful
using Measurements
using Statistics
using Test

@testset "PhyExp.jl" begin

@testset "A 类不确定度 ua" begin
    # 手算：x=[1,2,3,4]，s=√(5/3)，u_A = s/√4
    x = [1.0, 2.0, 3.0, 4.0]
    expected = std(x) / 2
    @test ua(x) ≈ expected
    @test ua([1.0, 1.0, 1.0, 1.0]) ≈ 0.0

    # 带单位（整体 Quantity 数组）
    q = [1.0, 2.0, 3.0, 4.0]u"V"
    uq = ua(q)
    @test ustrip(uq) ≈ expected
    @test unit(uq) == u"V"

    # 元素带单位
    q2 = [1.0u"m", 2.0u"m", 3.0u"m", 4.0u"m"]
    @test ustrip(ua(q2)) ≈ expected

    # 元素为 Measurement：只取标称值
    xm = [1.0 ± 0.1, 2.0 ± 0.2, 3.0 ± 0.3, 4.0 ± 0.4]
    @test ua(xm) ≈ expected

    # 错误情况
    @test_throws ArgumentError ua([1.0])
    @test_throws ArgumentError ua(3.14)
end

@testset "B 类不确定度 ub" begin
    @test ub(0.1) ≈ 0.1 / √3
    @test ub(0.1; dist=:normal) ≈ 0.1 / 3
    @test ub(0.1; K=2) ≈ 0.05
    @test ub(0.1; dist="均匀") ≈ 0.1 / √3
    @test ub(0.1; dist="正态") ≈ 0.1 / 3
    @test_throws ArgumentError ub(0.1; dist=:triangular)
    @test ub(0.5u"mm") |> ustrip ≈ 0.5 / √3
end

@testset "不确定度合成 combine" begin
    @test combine(0.3, 0.4) ≈ 0.5
    @test combine(0.3, 0.4, 1.2) ≈ sqrt(0.3^2 + 0.4^2 + 1.2^2)
    # 输入为 Measurement：自动提取不确定度
    @test combine(1.0 ± 0.3, 2.0 ± 0.4) ≈ 0.5
    # 带单位
    @test ustrip(combine(0.3u"V", 0.4u"V")) ≈ 0.5
    # 空合成（0 分量）
    @test combine(0.3, 0.0) ≈ 0.3
end

@testset "有效数字修约 roundsig" begin
    # 首位 > 2 取 1 位
    @test roundsig(9.8429 ± 0.0512; str=true) == "9.84 ± 0.05"
    # 首位 1/2 取 2 位
    @test roundsig(1.234 ± 0.012; str=true) == "1.234 ± 0.012"
    @test roundsig(2.3456 ± 0.0213; str=true) == "2.346 ± 0.021"
    # 修约后仍为 Measurement
    r = roundsig(9.8429 ± 0.0512)
    @test r isa Measurement
    @test value(r) == 9.84
    @test uncertainty(r) == 0.05
    # 不确定度为 0
    @test roundsig(3.14 ± 0.0; str=true) == "3.14"
    # 带单位
    @test roundsig((9.8429 ± 0.0512)u"m"; str=true) == "(9.84 ± 0.05) m"
    # 科学计数法（极小值，首位>2 取 1 位）
    @test roundsig((1.234e-6 ± 0.0512e-6); str=true) == "(1.23 ± 0.05)×10⁻⁶"
    # 科学计数法（极小值，首位=1 取 2 位）
    @test roundsig((1.234e-6 ± 0.012e-6); str=true) == "(1.234 ± 0.012)×10⁻⁶"
    # 末位百位：直接显示
    @test roundsig(12345.678 ± 1200.0; str=true) == "12300 ± 1200"
    # 科学计数法（极大值）
    @test roundsig(123456.78 ± 12000.0; str=true) == "(1.23 ± 0.12)×10⁵"
end

@testset "测量值提取 nominal / uncertval" begin
    @test nominal(3.14) == 3.14
    @test nominal(1.5 ± 0.2) == 1.5
    @test nominal((1.5 ± 0.2)u"m") == 1.5u"m"
    @test nominal(1.5u"m") == 1.5

    @test uncertval(3.14) == 0.0
    @test uncertval(1.5 ± 0.2) == 0.2
    @test uncertval((1.5 ± 0.2)u"m") == 0.2u"m"
    @test uncertval(1.5u"m") == 0.0u"m"
end

@testset "相对误差 relerr" begin
    @test relerr(1.02, 1.00) ≈ 0.02
    @test relerr(1.02, 1.00) isa Real
    # 相对不确定度
    @test relerr(100.0 ± 2.0) ≈ 0.02
    @test relerr((100.0 ± 2.0)u"m") ≈ 0.02
end

@testset "多次测量 mean_measurement" begin
    x = [1.01, 1.02, 1.00, 1.01]
    m = mean_measurement(x)
    @test value(m) ≈ mean(x)
    @test uncertainty(m) ≈ combine(ua(x), 0.0)
    m2 = mean_measurement([1.0, 2.0, 3.0, 4.0]; ub=0.1)
    @test uncertainty(m2) ≈ combine(ua([1.0, 2.0, 3.0, 4.0]), 0.1)
    # 带单位
    mq = mean_measurement([1.0, 2.0, 3.0, 4.0]u"V"; ub=0.1u"V")
    @test unit(mq) == u"V"
    @test ustrip(uncertval(mq)) ≈ combine(ua([1.0, 2.0, 3.0, 4.0]), 0.1)
end

@testset "线性回归 fit_line" begin
    # 精确线性关系：无不确定度数据
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = 2.0 .* x .+ 1.0
    r = fit_line(x, y)
    @test nominal(intercept(r)) ≈ 1.0 atol=1e-10
    @test nominal(slope(r)) ≈ 2.0 atol=1e-10
    @test correlation(r) ≈ 1.0 atol=1e-10

    # 带不确定度数据：y 有误差
    y2 = (2.0 .* x .+ 1.0) .± 0.1
    r2 = fit_line(x, y2)
    @test nominal(intercept(r2)) ≈ 1.0 atol=1e-8
    @test nominal(slope(r2)) ≈ 2.0 atol=1e-8
    @test uncertainty(intercept(r2)) > 0
    @test uncertainty(slope(r2)) > 0

    # 带单位：x 单位 s，y 单位 m → 斜率单位 m/s
    xq = [1.0, 2.0, 3.0, 4.0]u"s"
    yq = ([1.5, 3.5, 5.5, 7.5])u"m"
    r3 = fit_line(xq, yq)
    @test unit(slope(r3)) == u"m/s"
    @test unit(intercept(r3)) == u"m"
    @test ustrip(nominal(slope(r3))) ≈ 2.0 atol=1e-8
    @test ustrip(nominal(intercept(r3))) ≈ -0.5 atol=1e-8

    # 加权效果：不确定度小的点权重更大（数据有噪声时加权结果应偏向低误差点）
    x4 = [1.0, 2.0, 3.0, 4.0]
    y4 = [2.9, 5.1, 7.2, 9.8]                 # 大致 y = 2.3x + 0.6
    σ4 = [0.01, 0.5, 0.5, 0.5]                # 第一个点很准
    r4 = fit_line(x4, y4 .± σ4)
    # 第一个点 (1, 2.9) 权重极大 → 直线应接近过该点
    @test abs(nominal(intercept(r4)) + nominal(slope(r4)) - 2.9) < 0.05

    # 长度不匹配报错
    @test_throws DimensionMismatch fit_line([1.0, 2.0], [1.0, 2.0, 3.0])
end

@testset "York 正交回归 fit_line_york" begin
    x = collect(1.0:10.0)
    y = 3.0 .+ 2.0 .* x
    xm = x .± 0.1
    ym = y .± 0.2
    r = fit_line_york(xm, ym)
    @test nominal(intercept(r)) ≈ 3.0 atol=1e-6
    @test nominal(slope(r)) ≈ 2.0 atol=1e-6
    @test uncertainty(intercept(r)) > 0
    @test uncertainty(slope(r)) > 0

    # 全无不确定度 → 退化 fit_line
    r2 = fit_line_york(x, y)
    @test nominal(slope(r2)) ≈ 2.0 atol=1e-10

    # 带单位（元素带单位）
    xq = (x .± 0.1) .* u"s"
    yq = (y .± 0.2) .* u"m"
    r3 = fit_line_york(xq, yq)
    @test unit(slope(r3)) == u"m/s"
    @test ustrip(nominal(slope(r3))) ≈ 2.0 atol=1e-6
    @test ustrip(nominal(intercept(r3))) ≈ 3.0 atol=1e-6
end

@testset "异常值检验" begin
    x = [10.1, 10.2, 10.0, 10.3, 10.1, 15.7, 10.2]
    mask = outliers(x)
    @test count(mask) == 1
    @test x[findall(mask)][1] == 15.7
    @test reject(x) == [10.1, 10.2, 10.0, 10.3, 10.1, 10.2]
    @test grubbs(x)
    # 3σ 原则（大样本才可能检出：单点最大 |z| = (n-1)/√n）
    xbig = vcat(fill(10.0, 30), 30.0)
    mask3 = outliers(xbig; method=:zscore)
    @test sum(mask3) == 1
    # 小样本（n<11）3σ 检不出任何点（统计特性）
    @test !any(outliers([10.1, 10.2, 10.0, 10.3, 10.1, 15.7, 10.2]; method=:zscore))
    # zscore 数值正确
    zs = zscore([1.0, 2.0, 3.0])
    @test zs[1] ≈ (1.0 - 2.0) / std([1.0, 2.0, 3.0])
    # 无离群值的数据
    @test !any(outliers([10.0, 10.1, 10.0, 10.1]))
    # 带单位
    maskq = outliers([10.1, 10.2, 10.0, 10.3, 10.1, 15.7, 10.2]u"V")
    @test count(maskq) == 1
end

@testset "Markdown 表格" begin
    data = Dict("电压 U" => [1.0, 2.0, 3.0]u"V",
                "电流 I" => [0.5, 1.1, 1.5]u"mA")
    md = mdtable(data)
    # Dict 无序：只断言各列名（含单位）都在表头中出现
    @test occursin("电压 U / V", md)
    @test occursin("电流 I / mA", md)
    @test occursin("| --- |", md)
    @test occursin("| 1 |", md)               # 序号列
    # 测量值自动修约
    md2 = mdtable(Dict("x" => [9.8429 ± 0.0512]))
    @test occursin("9.84 ± 0.05", md2)
    # 无序号列
    @test !occursin("Index", mdtable(data; index=false))
    # 表题/表注/标签
    md3 = mdtable(data; title="测量数据", caption="三次测量", label="tab:data")
    @test occursin("表：测量数据", md3)
    @test occursin("*三次测量*", md3)
    @test occursin("tab:data", md3)
    # 长度不一致报错
    @test_throws DimensionMismatch mdtable(Dict("a" => [1, 2], "b" => [1, 2, 3]))
end

@testset "LaTeX 表格" begin
    data = Dict("电压 U" => [1.0, 2.0]u"V", "电流 I" => [0.5, 1.1]u"mA")
    lt = latextable(data)
    @test occursin("\\begin{tabular}", lt)
    @test occursin("U~(\\mathrm{V})", lt)
    # 带不确定度的数据 → \pm 输出
    lt2 = latextable(Dict("x" => [9.8429 ± 0.0512, 1.234 ± 0.012]))
    @test occursin("\\pm", lt2)
    @test occursin("9.84 \\pm 0.05", lt2)
end

@testset "中文别名（兼容 Python 版 API）" begin
    x = [1.0, 2.0, 3.0, 4.0]
    @test 求A类不确定度(x) ≈ ua(x)
    @test 仪器误差限转B类不确定度(0.1) ≈ ub(0.1)
    @test 不确定度合成(0.3, 0.4) ≈ 0.5
    @test 修约(9.8429 ± 0.0512; str=true) == "9.84 ± 0.05"
    r = 一元线性回归(x, 2.0 .* x .+ 1.0)
    @test nominal(斜率(r)) ≈ 2.0 atol=1e-10
    @test 提取标称值(1.5 ± 0.2) == 1.5
    @test 提取不确定度(1.5 ± 0.2) == 0.2
    @test 相对误差(1.02, 1.00) ≈ 0.02
    m = 多次测量结果([1.0, 2.0, 3.0, 4.0])
    @test value(m) ≈ 2.5
end

@testset "线性拟合结果打印" begin
    r = fit_line([1.0, 2.0, 3.0], [2.0, 4.0, 6.0])
    s = sprint(show, r)
    @test occursin("一元线性回归结果", s)
    @test occursin("斜率 b = 2.0", s)
end

end
