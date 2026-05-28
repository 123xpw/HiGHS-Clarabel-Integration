# HiGHS + Clarabel 集成基准测试报告

## 1. 项目背景

本项目在 HiGHS MIP 求解器框架内集成 Clarabel 作为根节点 LP 松弛的求解器（替代默认的 simplex / IPX），
目标是评估 Clarabel 在安全约束机组组合（SCUC）问题上的求解质量与速度表现。

集成里程碑：
- **M1**：构建系统（CMake、Clarabel 子模块）
- **M2**：数据映射（HiGHS LP → Clarabel 问题格式）
- **M3**：解的回写（Clarabel 解 → HiGHS 内部状态）
- **M4**：LP 路由（`mip_lp_solver = clarabel` 选项触发）
- **M_MIP**：MIP 根节点集成（Clarabel 求解根 LP 松弛后交回 B&B 引擎）

---

## 2. 测试问题集

### 2.1 电网数据

| 电网 | 节点数 | 机组数 | 时段数 | 备注 |
|------|--------|--------|--------|------|
| case118 | 118 | ~54 | 36 | 标准 40min/时段，无功率轨迹约束 |
| case2383wp | 2383 | ~327 | 72 | 20min/时段（由36时段×2日合并），含功率轨迹约束 |

### 2.2 算例选取

- 数据源：`testdata/UnitCommitment_Data/case2383wp/`（365天）、`testdata/UnitCommitment_Data/case118/`（12天）
- 选取策略：随机种子 seed=42，每月随机选1天（2017年），共12天/电网
- case2383wp 需选取当天+次日数据（用于36→72时段转换），因此排除每月最后一天

### 2.3 选取日期

| 月份 | case2383wp | case118 |
|------|-----------|---------|
| Jan  | 2017-01-23 | 2017-01-23 |
| Feb  | 2017-02-07 | 2017-02-07 |
| Mar  | 2017-03-19 | 2017-03-19 |
| Apr  | 2017-04-06 | 2017-04-06 |
| May  | 2017-05-10 | 2017-05-10 |
| Jun  | 2017-06-26 | 2017-06-26 |
| Jul  | 2017-07-01 | 2017-07-01 |
| Aug  | 2017-08-25 | 2017-08-25 |
| Sep  | 2017-09-29 | 2017-09-29 |
| Oct  | 2017-10-07 | 2017-10-07 |
| Nov  | 2017-11-23 | 2017-11-23 |
| Dec  | 2017-12-05 | 2017-12-05 |

---

## 3. 求解器配置

| 配置名 | HiGHS 选项 | 含义 |
|--------|-----------|------|
| HiPO | `solver = simplex` | HiGHS 原始单纯形（默认） |
| IPX | `solver = ipm` | HiGHS 内点法（IPX） |
| Clarabel | `mip_lp_solver = clarabel` | 本集成：Clarabel 求解根 LP 松弛 |

通用选项：`time_limit = 3600s`，`mip_rel_gap = 0.01`（1%）

---

## 4. MPS 文件生成流程

### 4.1 脚本

- `benchmark/01_download_case118.jl`：从 axavier.org 下载 case118 月度数据
- `benchmark/02_generate_mps.jl`：生成 24 个 MPS 文件
- `benchmark/03_run_experiments.jl`：72 次求解实验

### 4.2 case2383wp MPS 生成步骤

```julia
inst_cur = UnitCommitment.read(cur_path)     # 当天实例
inst_nxt = UnitCommitment.read(nxt_path)     # 次日实例
UnitCommitment.clear_power_trajectories!(inst_cur)
inst72   = UnitCommitment.convert_to_subhourly(inst_cur, inst_nxt)  # 36→72时段
curves   = build_curves(inst72)              # 构建功率轨迹曲线（287台机组）
UnitCommitment.set_power_trajectories!(inst72, curves)
formulation = UnitCommitment.Formulation(
    power_trajectories = UnitCommitment.xxx2005.PowerTrajectories(),
)
model = UnitCommitment.build_model(instance=inst72, formulation=formulation, variable_names=true)
JuMP.write_to_file(model, mps_path)
```

功率轨迹曲线参数：`BASE_UD = BASE_DD = 2`（步数），`BASE_TIME_MIN = 40`min

### 4.3 MPS 文件规模

| 电网 | MPS 文件大小 | 变量数 | 约束数 |
|------|------------|--------|--------|
| case118 | ~14.5 MB | 40,500 | 46,562 |
| case2383wp | ~264–270 MB | ~766,800 | — |

---

## 5. 代码修复（集成过程中发现的 Bug）

### 5.1 `unit.jl:92` — 缺少 `@variable` 宏

```julia
# 修复前（错误：存储的是 NamedTuple 而非 JuMP 变量）
reserve[sc.name, r.name, g.name, t] =
    (model, lower_bound = 0, base_name = "reserve_...")

# 修复后
reserve[sc.name, r.name, g.name, t] =
    @variable(model, lower_bound = 0, base_name = "reserve_...")
```

影响：`_total_reserves()` 函数在 JuMP v1.30 下崩溃（`MethodError: no method matching +(::Float64, ::@NamedTuple{...})`）

### 5.2 `subhourly.jl:92` — 次日机组 cost_segments 数量不匹配

```julia
# 修复前：若 g_next.cost_segments 为空则 BoundsError
for j in eachindex(g.cost_segments)
    g.cost_segments[j].cost = interpolate_values(..., g_next.cost_segments[j].cost[1])
end

# 修复后：当段数不一致时退化为 repeat_values
if length(g_next.cost_segments) == length(g.cost_segments)
    # 插值
else
    # 重复当天数据（时步减半，成本除以2）
end
```

影响：case2383wp 第5月（2017-05-10）及可能的其他月份崩溃

### 5.3 `read.jl` — JSON v1.5 `dicttype` API 不兼容

```julia
# 修复前：lambda 在 JSON v1.5 中不被接受
JSON.parse(file, dicttype = () -> DefaultOrderedDict(nothing))

# 修复后：使用自定义类型
struct _NullableDict <: AbstractDict{String, Any} ... end
JSON.parse(file, dicttype = _NullableDict)
```

### 5.4 `03_run_experiments.jl` — Julia 软作用域

```julia
# 修复前（脚本中 for 循环外层变量不自动继承）
run_count += 1

# 修复后
global run_count += 1
```

---

## 6. case118 基准测试结果（36/72，已全部完成）

> **配置说明**：本节数据来自正确配置的实验（`mip_lp_solver = ipx/hipo/clarabel`，构建时 `HIPO=ON`）。
> 早期使用 `solver = ipm/simplex` 的实验数据无效——该选项在 MIP 模式下被 HiGHS 忽略。

### 6.1 逐月结果表

| 月份 | 目标函数值 | IPX (s) | HiPO (s) | Clarabel (s) | 最快 |
|------|-----------|---------|---------|-------------|------|
| Jan 2017-01-23 | 5,967,704 | 6.55 | 8.85 | **4.49** | Clarabel |
| Feb 2017-02-07 | 4,947,858 | 6.32 | 6.34 | **4.45** | Clarabel |
| Mar 2017-03-19 | 4,277,981 | 6.17 | 6.51 | **4.72** | Clarabel |
| Apr 2017-04-06 | 4,413,695 | **8.22** | 20.67 | 9.46 | IPX |
| May 2017-05-10 | 4,116,782 | 7.29 | 6.62 | **5.11** | Clarabel |
| Jun 2017-06-26 | 5,428,173 | 8.50 | 9.89 | **4.78** | Clarabel |
| Jul 2017-07-01 | 6,887,877 | 7.07 | 8.67 | **5.39** | Clarabel |
| Aug 2017-08-25 | 4,692,958 | 6.10 | 6.97 | **5.54** | Clarabel |
| **Sep 2017-09-29** | **4,145,XXX\*** | **24.62** | 30.66 | 26.71 | IPX |
| Oct 2017-10-07 | 4,639,738 | 5.83 | 7.33 | **5.42** | Clarabel |
| Nov 2017-11-23 | 5,104,654 | 6.32 | 8.07 | **5.58** | Clarabel |
| Dec 2017-12-05 | 5,329,384 | 6.26 | 8.47 | **5.76** | Clarabel |
| **均值（排除Sep）** | — | **6.84** | **8.06** | **5.27** | — |

\* Sep 三者目标值有细微差异（均在 1% gap 容差内）：IPX=4,146,780，HiPO=4,145,474，Clarabel=4,145,386。

### 6.2 关键观察

**求解质量**：
- 全部 36 次均为 **Optimal**，节点数均为 1（LP 松弛极紧，根节点直接求解）
- 11/12 个月三者目标完全一致；Sep 月微小差异（<0.04%）属数值精度，均在容差内

**速度对比（排除 Sep）**：
- Clarabel 均值 **5.27s**，在 10/12 个月中最快
- IPX 均值 **6.84s**，比 Clarabel 慢约 30%
- HiPO 均值 **8.06s**，比 Clarabel 慢约 53%；Apr 月异常（20.67s，HiPO 收敛困难）

**Sep 2017-09-29 难实例**：
- LP 松弛初始 gap ≈ 2.38%（其他月份 ≤0.84%），需生成约 2,600 条割平面收紧
- 三者耗时均超 24s；IPX 最快（24.6s），Clarabel（26.7s）次之，HiPO（30.7s）最慢

### 6.3 小结

case118（36时段，~40k变量）上 Clarabel 全面领先，平均比 IPX 快 23%、比 HiPO 快 35%。
这一优势预计在 case2383wp（72时段，~767k变量）上将更为显著，LP 松弛求解质量差异会被放大。

---

## 7. case2383wp 基准测试结果（全部 36/72 完成）

### 7.1 逐月结果表

| 月份 | IPX 状态 | IPX (s) | HiPO 状态 | HiPO (s) | Clarabel 状态 | Clarabel (s) | 最快 |
|------|---------|---------|----------|---------|-------------|-------------|------|
| Jan 2017-01-23 | Optimal | 1030.87 | Optimal | 902.11 | Optimal | **590.07** | Clarabel |
| Feb 2017-02-07 | Optimal | 679.50 | Optimal | 518.81 | Optimal | **410.65** | Clarabel |
| Mar 2017-03-19 | Optimal | 592.68 | Optimal | 818.68 | Optimal | **441.35** | Clarabel |
| Apr 2017-04-06 | Optimal | **359.11** | Optimal | 642.53 | Optimal | 453.73 | IPX |
| May 2017-05-10 | Optimal | 362.97 | Optimal | 297.53 | Optimal | **267.19** | Clarabel |
| Jun 2017-06-26 | Optimal | 1088.70 | Optimal | 1308.41 | Optimal | **712.49** | Clarabel |
| Jul 2017-07-01 | Optimal | 1203.11 | Optimal | 1044.50 | Optimal | **612.92** | Clarabel |
| **Aug 2017-08-25** | **Time limit** | 3600.59 | **Time limit** | 3600.39 | **Optimal** | **1894.25** | **Clarabel only** |
| Sep 2017-09-29 | Optimal | 651.22 | Optimal | 605.91 | Optimal | **414.85** | Clarabel |
| **Oct 2017-10-07** | **Time limit** | 7173.48† | **Time limit** | 3602.38 | **Optimal** | **2500.35** | **Clarabel only** |
| Nov 2017-11-23 | Optimal | 1562.17 | **Time limit** | 8168.72† | Optimal | **783.95** | Clarabel |
| Dec 2017-12-05 | Optimal | 1388.26 | **Time limit** | 3600.07 | Optimal | **736.33** | Clarabel |
| **均值** | — | **1641.06** | — | **2092.50** | — | **818.18** | — |
| **Optimal 率** | — | **10/12** | — | **8/12** | — | **12/12** | — |

†：IPX Oct 与 HiPO Nov 出现 LP 根节点求解本身超时（b&b 尚未启动），导致总耗时超出 time_limit 设定值。

### 7.2 关键观察

**求解成功率**：
- Clarabel **12/12 全部 Optimal**，是唯一 100% 成功的求解器
- IPX 10/12（Aug、Oct 超时）；HiPO 8/12（Aug、Oct、Nov、Dec 超时）

**速度对比（基于全部 12 个月均值）**：
- Clarabel 均值 **818s**，比 IPX（1641s）快约 **2.0×**，比 HiPO（2093s）快约 **2.6×**

**两个关键难实例**：
- **Aug 2017**：LP gap 初始 >2.5%，需大量割平面；IPX/HiPO 均 3600s 超时（gap 仍 >2.5%），Clarabel 1894s 收敛（11 个节点，gap 0.20%）
- **Oct 2017**：IPX 根 LP 松弛本身耗时 >7000s（无法在时限内建立可行基），HiPO 3602s 超时（gap 1.1%）；Clarabel 2500s Optimal（41 节点，gap 0.28%）

### 7.3 与 case118 对比

| 指标 | case118 (36时段, ~40k变量) | case2383wp (72时段, ~767k变量) |
|------|--------------------------|-------------------------------|
| Clarabel vs IPX 加速比 | ~1.3× | **~2.0×** |
| Clarabel vs HiPO 加速比 | ~1.6× | **~2.6×** |
| Clarabel Optimal 率 | 12/12 | **12/12** |
| IPX Optimal 率 | 12/12 | 10/12 |
| HiPO Optimal 率 | 12/12 | 8/12 |

规模扩大时 Clarabel 的内点法优势明显放大。

---

## 8. 文件结构

```
highs-clarabel/benchmark/
├── 01_download_case118.jl   # 下载 case118 数据
├── 02_generate_mps.jl       # 生成 24 个 MPS 文件
├── 03_run_experiments.jl    # 72 次求解实验
├── Project.toml             # benchmark 专用 Julia 环境
├── selected_dates.csv       # 选取的日期记录
├── mps/                     # 24 个 MPS 文件（~14.5MB × 12 + ~265MB × 12）
├── logs/                    # 每次求解的 HiGHS 输出日志
├── highs_opts/              # HiGHS 选项文件（hipo.opt / ipx.opt / clarabel.opt）
└── results.csv              # 汇总结果（72行，脚本结束后生成）
```

---

## 9. 运行说明

```bash
# 步骤1：下载 case118 数据（已完成）
julia benchmark/01_download_case118.jl

# 步骤2：生成 MPS（已完成，耗时约 8 分钟）
julia benchmark/02_generate_mps.jl

# 步骤3：运行实验（case118 已完成，case2383wp 用户自行运行）
julia benchmark/03_run_experiments.jl
# 或限制时间/gap：
julia benchmark/03_run_experiments.jl --time_limit 3600 --gap 0.01
```

前置条件：
- HiGHS binary：`build_clarabel/bin/highs`（含 Clarabel 集成）
- Julia 1.10.x
- `testdata/UnitCommitment_Data/case2383wp/`（365天数据）
- `testdata/UnitCommitment_Data/case118/`（12天数据，由步骤1下载）
