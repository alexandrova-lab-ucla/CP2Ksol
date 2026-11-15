# 关键字参考

[English](keywords.md) · [中文](keywords.zh.md) · [返回 README](../README.zh.md)

本页列出本实现新增或修改的全部输入关键字。默认值按输入所用的单位给出，别名写在括号里。默认值和取值范围以 `src/input_cp2k_dft.F`、`src/input_cp2k_scf.F`、`src/input_cp2k_print_dft.F`、`src/start/input_cp2k.F` 为准，检查条件在 `src/cp_control_utils.F`、`src/scf_control_types.F`、`src/qs_sccs.F`。`&SCCS &DEVEL` 是开发者专用的关键字，保持默认即可。

### `&DFT &SCCS`

#### 模型与求解器

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `METHOD` | `ANDREUSSI` | 介电函数的平滑方式。新增两个取值：`ISODENSITY`、`SOLVENT_EXCLUDED` |
| `SOLVER` | `JACOBI`；`ISODENSITY` 和 `SOLVENT_EXCLUDED` 未给定时为 `NEWTON_CG` | `JACOBI`：对极化电荷的阻尼不动点迭代，按每步极化电荷的变化判断收敛（`EPS_SCCS`）。`NEWTON_CG`：对真实残差做 Newton 迭代，内层 CG，线搜索（`SOLVER_TOLERANCE`、`SOLVER_MAX_ITER`、`SCF_TOL_FACTOR`）。`NEWTON_CG` 要求 `RELATIVE_PERMITTIVITY` > 1，且 `METHOD` 不是 `SAA_ANDREUSSI` |
| `SOLVER_GRID_LEVEL` | 1 | 求解器所在的 MGRID 层；1 是辅助基组网格。n > 1 时把溶质密度限制到第 n 层，再把溶剂化势延拓回来；真空 Hartree 能仍在细网格上算。只用于 `NEWTON_CG`。`SOLVENT_EXCLUDED` 在 n > 1 时中止，`ISODENSITY` 警告收敛很差 |
| `SOLVER_STRIDE` | 1 | 每隔 n 个被激活的 SCF 步求解一次，其间电势冻结。必须 ≥ 1 |
| `PSI_FORM` | F；未给定时，`ISODENSITY` 为 T，`SOLVENT_EXCLUDED` 在 `NONLINEAR_DIELECTRIC F`、`NEWTON_CG`、且没有 `&EXCLUSION` 时为 T | ψ = √ε φ 离散化（https://doi.org/10.1063/1.4865107；https://doi.org/10.1063/1.4939125）：算符 s(−∇²+q)s，q = ∇²(√ε)/√ε，定义离散的 −∇·(ε∇)，内层 CG 每次迭代 2 次 FFT。只用于 `NEWTON_CG`；`NONLINEAR_DIELECTRIC` 必须关闭。`PSI_FORM F` 选用普通算符 |
| `SOLVER_TOLERANCE` | 1.0E-6 | 冷启动时、以及 `SCF_TOL_FACTOR` 为 0 时，求解器的 RMS 残差容差 [e/bohr³] |
| `SOLVER_MAX_ITER` | 50 | Newton 迭代次数上限 |
| `SCF_TOL_FACTOR` | 0.1；未给定时，线性介电且没有非线性离子为 0.3 | 从第二个被激活的 SCF 步起，求解器在 RMS ≤ max(`SOLVER_TOLERANCE`，因子 × 入口残差的 RMS) 时停止。0 关闭这一耦合。显式给定的值始终优先 |
| `SOLVER_ENERGY_TOL` | 1.0E-10 Ha | 线性介电 / 线性化离子路径：最近 20 次 CG 迭代使凹的自由能泛函总共增加不到这个值时，CG 也停止。0 = 只用 RMS 判据。必须 ≥ 0 |
| `SOLVER_MAX_FAILURES` | 0 | 连续这么多个 SCF 步的求解都没有收敛就中止；0 = 只警告。电势冻结的步既不计数也不清零 |
| `GAMMA`（`SURFACE_TENSION`） | 0 mN/m | CP2K 关键字，见[表面张力](#表面张力) |

#### 继承自 CP2K 的关键字

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `RELATIVE_PERMITTIVITY`（`DIELECTRIC_CONSTANT`、`EPSILON_RELATIVE`、`EPSILON_SOLVENT`） | 80 | 溶剂的体相相对介电常数 ε_b。298 K 的水用 78.4 |
| `EPS_SCCS`（`EPS_ITER`、`TAU_POL`） | 1.0E-6 | `JACOBI` 迭代中极化密度的收敛容差 |
| `EPS_SCF` | 0.5 | SCF 收敛到这个阈值后才开始 SCCS 循环 |
| `MAX_ITER` | 100 | `JACOBI` 求解器的迭代次数上限 |
| `MIXING`（`ETA`） | 0.6 | `JACOBI` 迭代的阻尼 |
| `ALPHA` | 0 mN/m | 排斥项 αS |
| `BETA` | 0 GPa | 色散项 βV |
| `DELTA_RHO` | 2.0E-5 | 计算空腔量子表面用的数值增量 |
| `DERIVATIVE_METHOD` | `FFT` | 实空间网格上的数值导数：`FFT`、`CD3`、`CD5`、`CD7` |

`METHOD` 的取值 `ANDREUSSI`、`FATTEBERT-GYGI`、`SAA_ANDREUSSI` 及其分节来自 CP2K 本身，见 CP2K 手册。

#### 电解质

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `DEBYE_LENGTH`（`LAMBDA_D`、`DEBYE_SCREENING_LENGTH`） | 0 Å = 关闭 | 线性化 Poisson–Boltzmann 的屏蔽长度。298 K 的 1 M 1:1 水溶液：3.04 Å。需要 `POISSON_SOLVER PERIODIC`、`RELATIVE_PERMITTIVITY` > 1，且 `METHOD` 不是 `SAA_ANDREUSSI` |
| `ION_CONCENTRATION` | 0 = 关闭 | 体相浓度，单位 mol/l；大于 0 时取代 `DEBYE_LENGTH`，后者由 λ_D = √(ε_b k_B T/(8π c z²)) 算出 |
| `ION_VALENCE` | 1 | 对称 z:z 电解质的价数 z（用于 `NONLINEAR_IONS` 和 `ION_CONCENTRATION`） |
| `TEMPERATURE` | 298 K | Boltzmann 因子（`NONLINEAR_IONS`）和由 `ION_CONCENTRATION` 算 λ_D 时用的电解质温度 |
| `NONLINEAR_IONS`（`NLPB`、`SIZE_MODIFIED_PB`） | F | 用带尺寸的非线性离子代替线性化。需要 `DEBYE_LENGTH` > 0 |
| `ION_STERIC_RADIUS`（`ION_RADIUS`） | 0 Å | 晶格气体的水合离子半径；0 给出点离子（普通 sinh 型非线性 PB） |
| `ION_PACKING_DIAMETER` | −1 = 2^(5/6) × `ION_STERIC_RADIUS` | 密堆积离子直径 d_ion；n_max = 1/d_ion³。θ_b ≥ 1 时输入中止 |

#### 非线性介电

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `NONLINEAR_DIELECTRIC`（`NLDIEL`、`DIELECTRIC_SATURATION`） | F；`SOLVENT_EXCLUDED` 未给定时为 T | 介电的转动饱和。需要 `SOLVER NEWTON_CG` 或 `NONLINEAR_IONS`。不能与 `PSI_FORM` 同用 |
| `EPS_INF` | 1.78 | 高频相对介电常数 |
| `SOLVENT_DENSITY` | 0.0335 Å⁻³（4.9649E-3 bohr⁻³） | 溶剂分子数密度 |
| `SOLVENT_DIPOLE` | 0.9449 e·bohr | 溶剂分子偶极矩 |
| `SMOOTHING_WIDTH` | 0.125 Å | 对极化电荷和束缚电荷做 Gauss 平滑的宽度；0 = 关闭 |

#### 赝核空腔密度

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `ADD_CORE_DENSITY`（`PSEUDO_CORE`） | T | 空腔密度 = 价电子密度 + 赝核密度。`ISODENSITY` 在 `F` 或 `CORE_DENSITY_SCALE` ≤ 0 时中止 |
| `CORE_DENSITY_SCALE`（`PSEUDO_CORE_SCALE`） | 1 | 缩放 Z − Z_eff |
| `CORE_DENSITY_ZMIN`（`PSEUDO_CORE_ZMIN`） | 1 | 原子序数更小的种类不加核 |
| `CORE_DENSITY_WIDTH`（`PSEUDO_CORE_SIGMA`） | 0 = GTH 核宽度 | 核密度的 Gauss 宽度 σ |
| `CORE_DENSITY_WIDTH_SCALE`（`PSEUDO_CORE_SIGMA_SCALE`） | 1 | 乘在 GTH 核宽度上 |
| `CORE_DENSITY_CHARGE`（`PSEUDO_CORE_CHARGE`） | −1 = `CORE_DENSITY_SCALE` × (Z − Z_eff) | 核 Gauss 里每个原子的电子数 |

#### Poisson–Boltzmann 求解的提速选项

除 `NEWTON_CG_FORCING` 外都默认开启；关掉它们就走普通代码路径。它们改变 CG 的路径，所以能量会在 SCF 停止判据的量级上有变动。

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `COARSE_DIRECT_SUM` | T | ψ 形式 CG 的两级预条件（粗空间：低平面波 s·exp(iGr)）：投影和展开用局部网格块上的可分离相位和代替三维 FFT，每次 CG 迭代 4 次 FFT 而不是 6 次。相位约定在第一次求解时检查，检查失败就改用 FFT 路线 |
| `LAP_CARRY` | T | ψ 形式：CG 之后以及线搜索试探点处的 ∇²(sφ)，由 CG 已经算过的项更新，每个残差省 2 次 FFT |
| `LINEAR_SINGLE_CG` | T | ψ 形式、线性介电、无离子：一次 CG 直接到 SCF 耦合的容差。`F` 则每个 Newton 步把残差降到十分之一 |
| `HISTORY_GUESS` | 2 | 线性介电、无离子：同一次 SCF 运行里前几次电势增量的个数（0–3），它们张成初始猜测；只有把入口残差降到原来的 0.9 以下时才使用。0 = 关闭 |
| `LINEAR_CG_FLOOR` | 0.8 | 线性介电、无离子：CG 在 RMS ≤ max(0.1 × 入口 RMS，值 × 容差) 时停止。范围 0–1；0 = 在 0.1 × 入口 RMS 处停止 |
| `COARSE_SPACE_IONS` | T | ψ 形式、线性化 Debye 离子：使用两级预条件；每次 CG 迭代 4 次 FFT 而不是 3 次，离子可达区域只占晶胞一小部分时迭代次数约减半。只有 N_e 能在 SCF 期间变化时才起作用（`&SCF &GCE` 里 `MAX_OUTER_UPDATES` > 0、`UPDATE_EVERY SCF_STEP`、方案 1 或 `EARLY_UPDATE_DELTA` > 0），条件与 `SCF_TOL_JUMP` 相同；固定 N_e 或非 GCE 运行始终用普通预条件 |
| `SCF_TOL_JUMP` | 100 | SCF 耦合的容差跟随上次重置以来入口残差的最小值。某个入口残差超过最小值的这个倍数，就重启最小值。0 = 关闭，否则 > 1。只在有 Debye 离子、且 `&SCF &GCE` 能在 SCF 期间改变 N_e 时起作用（`MAX_OUTER_UPDATES` > 0、`UPDATE_EVERY SCF_STEP`、方案 `MU_MIXING` 或 `EARLY_UPDATE_DELTA` > 0）；无离子或固定 N_e 的运行等同于 0 |
| `COLD_START_FROM_LAST` | T | 一次 SCF 运行的第一次求解从上一次运行结束时的电势出发；容差仍是 `SOLVER_TOLERANCE` |
| `NEWTON_CG_FORCING` | 0 = 关闭 | 非线性路径：每个 Newton 步的 CG 在 RMS ≤ max(0.1 × 入口 RMS，值 × 容差) 时停止。范围 0–1。约快 15 %，`EPS_SCF` 为 1E-6 时能量最多变动约 1E-5 Ha |

#### 诊断

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `PRINT_POL_DIAGNOSTICS` | F | 打印极化能的各项和每次 CG 求解的一行；出现极化能为正的警告时，也会打印这些能量项 |
| `SOLVER_DEBUG` | F | 在每次求解的第一个 Newton 步上，对 Jacobian 和自由能做有限差分检查 |

### `&SCCS &ISODENSITY`

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `RHO_CAVITY`（`NC`、`RHO_CUTOFF`） | 0.0030 Å⁻³（4.44554E-4 bohr⁻³） | 空腔边界的密度等值；在 GPW/GTH 密度上对 11 个实验溶剂化自由能拟合得到 |
| `SWITCH_WIDTH`（`SIGMA`） | 0.6 | ln(n/`RHO_CAVITY`) 轴上 erfc 开关的宽度 |

### `&SCCS &SOLVENT_EXCLUDED`

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `RHO_CAVITY`（`NC`） | 0.0158 Å⁻³（2.34132E-3 bohr⁻³） | 最内层表面的密度等值 |
| `SWITCH_WIDTH`（`SIGMA`） | 0.6 | erfc 开关的宽度，每个表面都用 |
| `PROBE_KERNEL_WIDTH`（`KERNEL_A`） | 0.125 Å | 探针壳核的衰减长度 a；壳的厚度为 a/`SWITCH_WIDTH` |
| `SOLVENT_RADIUS` | 1.4 Å | 溶剂分子中心与 van der Waals 表面的最近距离 |
| `ION_RADIUS` | 4.0 Å | 电解质离子的最近距离；这个表面与介电表面之间有溶剂但没有离子 |
| `DIELECTRIC_RADIUS` | 1.0 Å（SE 用 0.7 Å） | 介电屏蔽起始表面向内的位移；决定电极处无屏蔽间隙的宽度。线性介电配大于 0 的值会警告 SCF 可能塌缩 |
| `CAVITATION_RADIUS` | 0 Å | 承载空化项的表面向内的位移 |

### `&SCCS &REGION`

在含显式水的晶胞里，对隐式溶剂做平面限制。介电空腔和离子空腔沿晶胞轴 `AXIS` 乘以 m(z) = ½[1 + erf((z − LOWER_EDGE)/(WIDTH√2))] × ½[1 + erf((UPPER_EDGE − z)/(WIDTH√2))]，并对周期像求和。`LOWER_EDGE`、`UPPER_EDGE` 和 `WIDTH` 是沿该轴晶格矢量的长度，所以这个矢量必须与另两个垂直。分节出现就生效。每次 SCCS 调用都报告落在区域外的空腔体积，施加掩模之后报告窗口内、整个晶胞内和 `PROBE_ATOMS` 周围球内的读数。空化表面不加掩模。`SAA_ANDREUSSI` 不可用。

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `AXIS` | 3 | 垂直于显式/隐式边界的晶胞轴（1、2 或 3） |
| `LOWER_EDGE` | 无 | 区域下沿，范围 [0, L)。没有 `&FLEXIBLE` 时必须给；两个边要么都给要么都不给 |
| `UPPER_EDGE` | 无 | 区域上沿，> `LOWER_EDGE`；可以超过 L，使区域穿过晶胞边界；`UPPER_EDGE` − `LOWER_EDGE` < L − 6 `WIDTH` |
| `WIDTH` | 1.0 Å | erf 边缘的宽度；必须 > 0。膜用 0.5 Å |
| `APPLY` | T | `F` 不改变模型，只报告泄漏 |
| `PROBE_ATOMS` | 无 | 另外在这些原子周围取读数（从 1 起编号）；只用于报告 |
| `PROBE_RADIUS` | 3.0 Å | 这些球的半径；必须 > 0 |

#### `&REGION &FLEXIBLE`

几何优化或 MD 中，区域边界跟着显式溶剂走。列出的原子定义粗粒化密度 ρ̄(r) = Σᵢ wᵢ g(r − Rᵢ)，g 是宽度为 `COARSE_GRAINING_LENGTH` 的归一化 Gauss，对所有周期像求和。令 t = erfc⁻¹(2ρ̄/`DENSITY`)，掩模为 m(r) = ½ erfc[(`OFFSET` − √2·`COARSE_GRAINING_LENGTH`·t)/(√2·`WIDTH`)]，其中 `WIDTH` 是 `&REGION` 的关键字。显式溶剂内部 m 为 0（分子之间的空隙也算在内），远离它时为 1。对密度为 `DENSITY` 的平坦均匀水膜，这个掩模等于 `&REGION` 的 erf 平面，位置在膜的 Gibbs 分界面之上 `OFFSET` 处，所以平面设置可以用同一个 `WIDTH` 直接沿用。掩模乘在 S_diel 和 S_ion 上，列出的原子一动就重建，并对列出的原子给出解析力。同时给了 `LOWER_EDGE` 和 `UPPER_EDGE` 时，平面掩模再乘上去。只用于 `METHOD SOLVENT_EXCLUDED`。ρ̄ 在倒空间构造（每次重建一次 FFT，算力再一次），所以代价不随 `COARSE_GRAINING_LENGTH` 增长。

膜里的水氧和离子按权重 1 列出。不要列氢，也不要列电极原子。在 Pt(111) 上三层水双层加一个 K⁺ 的体系中（88 个原子；8 个静态帧和 100 步 MD），这样列出能把隐式溶剂挡在电极与水膜之间的间隙之外：那里的 max(ε − 1) 保持在 1.0e-4 以下。把电极原子也列进去，权重取 `DENSITY` 除以金属的原子密度，会让平板裸露一侧的边界向外移动（在那个测试里，从最低一层 Pt 之下 2.71 Å 移到 3.50 Å），并使带电平板的功函数降低 0.29 eV。`ATOM_LIST` 可以把列表限制在所列种类的部分原子上。

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `KIND <name> <weight>` | 无 | 可重复。膜列 `Ow 1.0` 和 `K 1.0`。`&KIND` 里的种类名及其在 ρ̄ 中的权重。至少一行；未知名字会中止 |
| `COARSE_GRAINING_LENGTH` | 2.4 Å | 粗粒化的 Gauss 宽度；必须大于 `&REGION` 的 `WIDTH` |
| `DENSITY` | 0.0334 Å⁻³ | 显式溶剂的数密度，按单位体积的加权原子数计（液态水） |
| `OFFSET` | 2.4 Å | 边界（m = ½）到显式溶剂 Gibbs 分界面的距离 |
| `ATOM_LIST` | 所列种类的全部原子 | 把 ρ̄ 限制在这些原子上（从 1 起编号）；选不出原子的列表会中止 |
| `FROM_DENSITY` | F | 改由电子密度构造场：ρ̄ = g ∗ (1 − S_solv)，S_solv 是溶剂空腔函数，`DENSITY` = 1。掩模随之进入 Kohn–Sham 势，不列原子，`OFFSET` 从 S_solv = ½ 的表面算起，默认 0 |

每次重建打印 `SCCS| flexible region: atoms n  int rho_bar dV x  sum of weights y  int(1-m) dV [bohr^3] z`；高斯在网格上分辨得好时 x 和 y 一致，z 是对隐式溶剂关闭的体积。

### `&SCCS &EXCLUSION`

对显式溶剂原子（会动的显式水、K⁺ 或 Na⁺）周围的隐式介电和离子做跟随原子的排除。只用于 `METHOD SOLVENT_EXCLUDED`，其他方法会中止。M(r) = ∏ₐ[1 − ½ erfc((|r−Rₐ|−rₐ)/(√2·`WIDTH`))] 乘在 S_diel 和 S_ion 上（`&REGION` 的掩模再乘上去）。每次 SCCS 调用都按当前位置重建掩模，并对被排除的原子给出解析力。空化表面不加掩模。没有默认半径。分节参数：`&EXCLUSION on`。

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `KIND <name> <r_diel> <r_ion>` | 无 | 可重复。203 原子 Cu/水/K⁺ 膜用 `Ow 1.80 2.40`、`Hw 1.10 1.70`、`K 2.00 2.60`。`&KIND` 里的种类名，两个半径单位都是 Å。至少一行；没有对应行的种类不排除；未知或重复的名字会中止 |
| `WIDTH` | 0.25 Å | erfc 边缘的宽度；必须 > 0。膜用 0.25 Å |
| `ATOM_LIST` | 所列种类的全部原子 | 把排除限制在这些原子上（从 1 起编号） |
| `PRINT_SUMMARY` | F | 每次 SCCS 调用打印被排除的原子数以及 ∫(1 − M_diel) dV、∫(1 − M_ion) dV |

### 表面张力

`&SCCS GAMMA` 保持 CP2K 的默认值 0。算分子的溶剂化自由能时，`ISODENSITY` 设 `GAMMA [mN/m] 18.47`，`SOLVENT_EXCLUDED` 在 `DIELECTRIC_RADIUS` 0.7 Å 时设 `GAMMA [mN/m] 11.09`（分别是 1.1531 和 0.6923 meV/Å²，1 meV/Å² = 16.0218 mN/m），这是用本文档的空腔设置在烷烃上拟合的。拟合时从 `GAMMA [mN/m] 1.0` 的运行里读出空化项再按比例缩放；在输入里设这个值还会把空化势反馈进 SCF，这个二阶差别没有单独测量。平板、吸附能、零电荷电位和电容保持为 0：`GAMMA [mN/m] 14.50` 时，这一项只让 Pt(111) 上 H 的吸附能移动了 −1.8 meV。

### `&DFT &SCF &GCE`

CP2K 自带的 `TARGET_WORKFUNCTION` 和 `MIXING_COEF` 属于方案 `MU_MIXING`。下面是本代码新增的关键字。电极的对电荷要么是 `PLANAR_COUNTER_CHARGE`，要么是 SCCS 电解质（`DEBYE_LENGTH` > 0）；两者同时给会中止。恒电势运行需要 `&DIAGONALIZATION`、`&SMEAR METHOD FERMI_DIRAC` 和 `ADDED_MOS`。不支持轨道变换。

#### 方案与目标

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `UPDATE_SCHEME` | `MU_MIXING` | `MU_MIXING`：Fermi 能级向目标功函数混合。`CAPACITANCE`：N_e ← N_e + C (WF − WF_target)，C 是割线更新的电容，固定 N_e 下用 Fermi–Dirac 占据。`DYNAMICS`：电子数动力学，见下 |
| `UPDATE_EVERY` | `SCF_STEP` | `CAPACITANCE`：每个 SCF 步更新 N_e，或每个收敛的 SCF 更新一次（`SCF_CYCLE`，用于 MD 和 `GEO_OPT`）。`SCF_CYCLE` 是受支持的模式，`DYNAMICS` 和 `ENFORCE_TARGET` 只能用它 |
| `TARGET_POTENTIAL` | 未设置 | 电极电位，V vs SHE，按 eV 输入；目标功函数设为 `PHI_SHE` + 该值。`CAPACITANCE` 和 `DYNAMICS` 需要 `TARGET_POTENTIAL`、`TARGET_WORKFUNCTION`、`TARGET_CHEMICAL_POTENTIAL` 之一，固定 N_e 的块除外（`MAX_OUTER_UPDATES` 0 配 `UPDATE_EVERY SCF_CYCLE`） |
| `PHI_SHE` | 4.600 eV | 配合 `TARGET_POTENTIAL` 使用的 SHE 绝对电位。4.600 eV 是本实现针对九个 Cu、Ag、Au 晶面实验 pzc 标定的值（拟合值 4.602 eV） |
| `TARGET_CHEMICAL_POTENTIAL` | 未设置 | 电子绝对化学势 [eV]，μ_e = −WF；覆盖 `TARGET_WORKFUNCTION` 和 `TARGET_POTENTIAL` |
| `EXCESS_ELECTRONS` | 0 | 开始时加到 `CHARGE` 态上的分数过剩电子；配 `MAX_OUTER_UPDATES 0` 得到分数电子数下的恒电荷运行 |
| `SOLVENT_ONLY_FOR_CHARGE` | F | 隐式溶剂只用来确定 N_e；更新步上 N_e 一固定，就关掉溶剂，SCF 在真空里重新收敛。没有更新的步（`UPDATE_INTERVAL`）上溶剂保持关闭，所以能量和力不含溶剂 |

#### 更新循环（`CAPACITANCE`，`UPDATE_EVERY SCF_CYCLE`）

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `MAX_OUTER_UPDATES` | 20 | 每个更新步里 N_e 更新次数的上限，每次更新都在内层 SCF 收敛之后进行。`MAX_SCF` 计的是上次更新以来的内层步数（第一个 N_e 则从 SCF 开始算），所以每个更新周期都有完整的 `MAX_SCF` |
| `UPDATE_INTERVAL` | 1（MD 推荐 5，在速度和精度之间取得平衡） | 只在第 1、n、2n、3n、… 次力计算时更新，初始那次算第 1 次。n > 1 时，一次更新计算里的第一次更新从上次采纳的电容出发，忽略之前离子步的割线点（它含有中间各步的几何漂移）；割线在这次计算内部重建。只做一次更新的计算（`MAX_OUTER_UPDATES` 1、`CONVERGENCE_MODE` 2，或模式 1 的第一步之后）使用更早那次更新计算的点对，新电容限制在上次采纳值的 0.5 到 2 倍之间 |
| `CONVERGENCE_MODE` | 3 | 1：只在第一个离子步把功函数收敛，之后每个更新步更新一次；2：每个更新步始终更新一次；3：每个更新步都收敛（最多 `MAX_OUTER_UPDATES` 次）。`MAX_OUTER_UPDATES 0` 在任何模式下都保持 N_e 固定。`DYNAMICS` 忽略它 |
| `FERMI_TOL` | 0.01 eV | \|WF − WF_target\| 的容差，用于更新循环、`CONVERGENCE_CHECK` 和 `ENFORCE_TARGET` |
| `CONVERGENCE_CHECK` | F | SCF 收敛还要求 \|WF − WF_target\| < `FERMI_TOL` |
| `FINAL_TOL_FRAC` | 1 | 内层 SCF 收敛后，只要 \|WF − WF_target\| > `FERMI_TOL` × 该值，就继续更新 N_e。范围 (0, 1] |
| `CAPACITANCE_INIT` / `_MIN` / `_MAX` | 1 / 0.05 / 2 e/V | 割线电容的初值和上下限；都为正，MIN ≤ MAX |
| `NE_STEP_MAX` | 0.5 e | 每次更新 N_e 的最大变化 |
| `FIRST_UPDATE_SIZE` | 0 = 关闭 | 第一次更新的固定大小 [e]，符号取自功函数误差；0 则按 `CAPACITANCE_INIT` 定步长 |
| `FERMI_UPDATE_TOL` | 0.02 eV | 只有功函数自上次更新以来移动超过这个值，才刷新割线电容 |
| `CAP_REFRESH_MIN` | 0 eV | 只有功函数移动超过 max(`FERMI_UPDATE_TOL`，该值)，才刷新割线 |
| `MIN_STEPS_AFTER_UPDATE` | 12 | 一次更新之后必须再走完的内层 SCF 步数，之后才能宣布 SCF 收敛。`PIN_MIXER_G0` 关闭且用 Broyden 混合时，会带警告提高到 `NBROYDEN` + 2 |
| `UPDATE_SCF_DELTA` | 0.01 | 只用于 `UPDATE_EVERY SCF_STEP`：只有上一步的 SCF 密度变化低于这个值时才更新 N_e |
| `PIN_MIXER_G0` | F；`DYNAMICS` 未给定时为 T | 同时把混合器输入的 G = 0 系数设为 N_e/V，使新 N_e 立刻体现在密度里。只用于自旋限制运行；UKS 会带警告关闭 |
| `RESPONSE_GUESS` | F | 在 N_e 更新处，进入下一个内层 SCF 的密度是 ρ + dN_e·∂ρ/∂N_e，∂ρ/∂N_e 取同一次力计算里最近两次更新所得密度的差商（保留给之后的力计算）。混合器输入设为同一个密度。在两个这样的状态都有之前，更新只移动 G = 0 系数。只用于自旋限制运行。只有一次力计算里做两次或更多次 N_e 更新时才起作用：在 323 原子的 Cu(111)/水膜里，它把第一次力计算的 SCF 步数从 163 降到 138，与提前更新一起用时降到 74（只用提前更新为 122）。每次力计算只更新一次的较小测试晶胞没有变化，在它们上强制开启提前更新还会多花 SCF 步。每次使用都打印 `GCE\| response guess: density moved along d rho/d N_e by x e` |
| `KEEP_MIXER_BELOW` | 0；`DYNAMICS` 未给定时为 0.05 e | \|dN_e\| [e] 低于这个值的 N_e 更新保留 Broyden 历史，不再重启 |
| `EARLY_UPDATE_DELTA` | 0 = 关闭 | 大于 0 时，内层 SCF 收敛之前也会更新 N_e：密度变化低于该值，并且在当前 N_e 下已走过 `EARLY_UPDATE_MIN_STEPS` 步。最终的 N_e 仍在收敛的密度上判断 |
| `EARLY_UPDATE_MIN_STEPS` | 3 | 允许提前更新之前，在新 N_e 下要走的步数 |
| `EARLY_UPDATE_MAX` | 2 | 每次 SCF 运行提前更新的最大次数；不计入 `MAX_OUTER_UPDATES` |
| `EARLY_UPDATE_WF_STABLE` | 0 eV | 大于 0 时，提前更新还要求相邻两个内层步之间功函数的变化小于该值 |
| `EARLY_UPDATE_SECANT` | F | 第一次提前更新记录它的 (N_e, WF) 点对，使下一次更新有割线电容，之后的提前更新再刷新它。`F`：提前更新不动割线历史，留给第一次收敛的更新 |
| `CARRY_GAIN` | 0 = 关闭 | 每个离子步收敛之后，N_e 为下一步移动 该值 × C × (WF − WF_target)。0.5 到 1 比较合适。不能与 `DYNAMICS` 同用 |
| `ENFORCE_TARGET` | F | 一个更新步用完了全部 `MAX_OUTER_UPDATES` 次更新，仍偏离目标超过 `FERMI_TOL` 时，带着原因停止运行（目标落在能隙里、步长太小等）。在 `GEO_OPT` 里每次力计算都是完整的更新步，所以最后一次在 `FERMI_TOL` 以内。`UPDATE_EVERY SCF_STEP` 或 `MAX_OUTER_UPDATES` 0 时无效 |

#### Restart 状态

每个离子步之后，代码把这些值写进运行中输入的 `&GCE` 分节，所以每个 MD 和 `GEO_OPT` 的 restart 文件都带着它们。不需要手动编辑。

| 关键字 | 默认值 | 内容 |
|---|---|---|
| `RESTART_VALID` | F | `RESTART_*` 的值有效（`CAPACITANCE` 和 `DYNAMICS` 共用） |
| `RESTART_N_E`、`RESTART_N_INPUT`、`RESTART_CAPACITANCE`、`RESTART_ION_STEP` | 0 | 最近一次计算的 N_e [e]、输入电荷态的 N_e [e]、电容、已完成的力计算次数（两种方案共用） |
| `RESTART_N_LAST`、`RESTART_WF_LAST`、`RESTART_SECANT_READY`、`RESTART_LAST_WF`、`RESTART_STEPS_SINCE_UPDATE` | 0 / F | `CAPACITANCE`：上一个割线点 (N_e, WF)、割线标志、上一个离子步的功函数、自上次更新以来的内层步数 |
| `RESTART_P_N`、`RESTART_DU_LAST`、`RESTART_NEXT_LEN`、`RESTART_OPEN_PENDING`、`RESTART_XI`、`RESTART_P_XI` | 0 / 1 / F | `DYNAMICS`：N_e 的动量、上一次的 U − U_target、距离下一次电位测量的 MD 步数、是否该做开头半步、链的位置和动量 |

`&EXT_RESTART` 有关键字 `RESTART_GCE`（和该分节其他 `RESTART_*` 关键字一样跟随 `RESTART_DEFAULT`，默认 T），它把 restart 文件里显式的 `RESTART_*` 关键字拷进新输入。

### `&DFT &SCF &GCE UPDATE_SCHEME DYNAMICS`

只用于 `RUN_TYPE MD`，`UPDATE_EVERY SCF_CYCLE`，固定的 MD 时间步。电子数 n 是动力学变量（质量 M_n，动量 P_n），由 U − U_T 驱动，每个离子步一次 SCF，第一次力计算之后 SCF 内部不再更新 N_e。第一次计算用常规的割线循环把 n 收敛到目标（`MAX_OUTER_UPDATES`；0 则从 `EXCESS_ELECTRONS` 出发）：

```
dn/dt = P_n/M_n      dP_n/dt = (U − U_T) − P_n P_ξ1/Q_1      (chain: dξ_j/dt = P_ξj/Q_j, dP_ξj/dt = P_ξ(j−1)²/Q_(j−1) − kT − P_ξj P_ξ(j+1)/Q_(j+1))
M_n = 1/(C ω²),  ω = 2π/POTENTIOSTAT_PERIOD,  Q_j = kT/ω_th²
```

链用 https://doi.org/10.1080/00268979600100761 的对称 Suzuki–Yoshida 分解。MD 时间步和温度从 `MOTION&MD` 读取。

| 关键字 | 默认值 | 含义 |
|---|---|---|
| `POTENTIOSTAT_PERIOD` | 必填，> 0 fs | n 自由振荡的周期 2π/ω；决定 M_n。取值要远大于 10 个 MD 时间步（低于 10 个更新周期时代码会警告） |
| `POTENTIOSTAT_CAPACITANCE` | 0 = 第一次计算测得的割线值 | 晶胞的 C = dN_e/dU [e/V]；没有可用值时运行停止。生产运行请用静态 N_e 扫描测量 |
| `POTENTIOSTAT_CHAIN_LENGTH` | 0 | 0：恒能量（K_n + Ω 守恒）；1–16：在 `POTENTIOSTAT_TEMPERATURE` 下的 Nosé–Hoover 链（做统计用 3 个及以上；单个链节对谐振模式的抽样不遍历） |
| `POTENTIOSTAT_TEMPERATURE` | `MOTION&MD TEMPERATURE` | 链的 kT |
| `POTENTIOSTAT_THERMOSTAT_PERIOD` | `POTENTIOSTAT_PERIOD` | Q_j = kT/ω_th² |
| `POTENTIOSTAT_INITIAL_RATE` | 0 e/fs | 初始 dn/dt；续算运行恢复上一次运行的状态，不用这个值 |

用了 `SOLVENT_ONLY_FOR_CHARGE` 时，电位在开着溶剂的步上测量；在这些步上，`GCE| trace` 行、`RESTART_LAST_WF` 和 `CARRY_GAIN` 用开溶剂时的功函数（不开溶剂的步没有这个值，报的是不含溶剂的功函数，它不是电极电位）；`UPDATE_INTERVAL` n 决定间隔，kick 以 n 步的冲量施加。守恒量是 (MD 的 `Cons Qty`) + `Hextra`，Ω = E + WF_T (n − n_input)。

### 恒电势运行的续算

用 MD 驱动写出的 `.restart` 文件作输入（`cp2k -i <project>-1.restart`），或者用固定的输入加 `&EXT_RESTART`。第一次力计算以相同的 N_e 和相同的更新步相位重复上一次运行的最后一次计算，不做额外的 N_e 更新，之后的步与不间断的运行一致。如果上一次运行的最后一次计算结束时偏离目标（用完了全部 `MAX_OUTER_UPDATES`），重复的那次计算会再做一次更新。恢复的 N_e 只是起始值：更新步和 `CARRY_GAIN` 拿功函数与新输入的 `TARGET_POTENTIAL` 比较，所以换了目标电位续算时，N_e 会在 `UPDATE_INTERVAL` 次力计算之内移到新目标。`RESTART_GCE F` 跳过电子数状态，N_e 重新从 `CHARGE` 和 `EXCESS_ELECTRONS` 开始。`CHARGE` 或原子数不同的续算会带着信息中止。用 `DYNAMICS` 时，`POTENTIOSTAT_CHAIN_LENGTH` 必须等于 restart 文件里存的长度。

### 自旋

用 `UPDATE_SCHEME CAPACITANCE` 或 `DYNAMICS` 时，`UKS` 的两种自旋共用一个电子数和一个 Fermi 能级。`FIXED_MAGNETIC_MOMENT` 配恒电势运行会中止。在共用 Fermi 能级的情况下，`MULTIPLICITY` 只设初始占据（有警告说明）：从三重态出发的运行可能收敛到非磁态，同一三重态出发的固定 N_e 运行可能到达能量更低的高自旋解。方案 `MU_MIXING` 配 `LSD` 保留按自旋的占据调用。

### `&DFT &PRINT &SCCS`

`IONIC_CHARGE_DENSITY`（debug 打印级别）写出 Poisson–Boltzmann 电解质（`DEBYE_LENGTH` > 0）离子电荷密度的 cube 文件，带 `STRIDE` 和 `APPEND`。

## 原名对照

下列关键字改了名。原名在所有输入和 restart 文件里仍然可用；代码两个名字都读，写出时用新名。

| 分节 | 原名 | 新名 |
|---|---|---|
| `&SCCS` | `ION_CHARGE` | `ION_VALENCE` |
| `&SCCS` | `EPS_NLPB` | `SOLVER_TOLERANCE` |
| `&SCCS` | `MAX_ITER_NLPB` | `SOLVER_MAX_ITER` |
| `&SCCS` | `NLPB_ENERGY_TOL` | `SOLVER_ENERGY_TOL` |
| `&SCCS` | `MAX_NLPB_FAILURES` | `SOLVER_MAX_FAILURES` |
| `&SCCS` | `DEBUG_NLPB` | `SOLVER_DEBUG` |
| `&SOLVENT_EXCLUDED` | `RADIUS_SOLVENT` | `SOLVENT_RADIUS` |
| `&SOLVENT_EXCLUDED` | `RADIUS_ION` | `ION_RADIUS` |
| `&SOLVENT_EXCLUDED` | `RADIUS_DIELECTRIC` | `DIELECTRIC_RADIUS` |
| `&SOLVENT_EXCLUDED` | `RADIUS_CAVITATION` | `CAVITATION_RADIUS` |
| `&REGION` | `Z_LOW` | `LOWER_EDGE` |
| `&REGION` | `Z_HIGH` | `UPPER_EDGE` |
| `&REGION &FLEXIBLE` | `SIGMA` | `COARSE_GRAINING_LENGTH` |
| `&GCE` | `NELECTRON_EXCESS` | `EXCESS_ELECTRONS` |

直接写在 `&SCCS` 下的 `ION_RADIUS` 仍是 `ION_STERIC_RADIUS`（带尺寸离子的晶格半径）的别名。在 `&SOLVENT_EXCLUDED` 里，它是离子到溶质表面的最近距离。

## 符号

符号、物理量、输入关键字和单位。没有输入关键字的量不列。

| 符号 | 物理量 | 关键字 | 单位 |
|---|---|---|---|
| **`&SCCS`** | | | |
| ε_b | 体相相对介电常数 | `RELATIVE_PERMITTIVITY` | – |
| ε_∞ | 溶剂的光学介电常数 | `EPS_INF` | – |
| λ_D | Debye 长度 | `DEBYE_LENGTH` | Å |
| c | 离子浓度 | `ION_CONCENTRATION` | M |
| z_ion | 离子价数（对称 z_ion:z_ion 电解质） | `ION_VALENCE` | – |
| T | 电解质温度 | `TEMPERATURE` | K |
| R_h | 晶格气体的水合离子半径 | `ION_STERIC_RADIUS` | Å |
| ψ | Poisson 算符的对称化形式 | `PSI_FORM` | – |
| τ | 空化表面张力 | `GAMMA` | mN/m |
| **`&SCCS ISODENSITY and SOLVENT_EXCLUDED`** | | | |
| n_c | 溶质表面的电子密度等值 | `RHO_CAVITY` | Å⁻³ |
| σ_c | 溶质表面在 ln(n/n_c) 轴上的宽度 | `SWITCH_WIDTH` | – |
| **`&SCCS SOLVENT_EXCLUDED`** | | | |
| λ_p | 探针壳核的衰减长度 | `PROBE_KERNEL_WIDTH` | Å |
| R_solv | 溶剂分子中心到溶质表面的最近距离 | `SOLVENT_RADIUS` | Å |
| R_ion | 电解质离子到溶质表面的最近距离 | `ION_RADIUS` | Å |
| R_diel | 介电屏蔽从溶剂中心表面向内延伸的距离 | `DIELECTRIC_RADIUS` | Å |
| R_cav | 承载空化项的表面向内的位移 | `CAVITATION_RADIUS` | Å |
| **`&SCCS EXCLUSION`** | | | |
| w_ex | 排除球的边缘宽度 | `WIDTH` | Å |
| **`&SCCS REGION`** | | | |
| z_low, z_high | 隐式区域的下沿和上沿 | `LOWER_EDGE`, `UPPER_EDGE` | Å |
| w | 区域的边缘宽度 | `WIDTH` | Å |
| **`&SCCS REGION FLEXIBLE`** | | | |
| ℓ | 粗粒化溶剂密度的 Gauss 宽度 | `COARSE_GRAINING_LENGTH` | Å |
| ρ_0 | 显式溶剂的数密度 | `DENSITY` | Å⁻³ |
| d_off | 边界到 Gibbs 分界面的距离 | `OFFSET` | Å |
| **`&DFT SCF GCE`** | | | |
| Φ_SHE | SHE 的绝对电位 | `PHI_SHE` | eV |
| U_T | 目标电极电位 vs SHE | `TARGET_POTENTIAL` | eV（按 V 用） |
| μ_e | 目标电子化学势 | `TARGET_CHEMICAL_POTENTIAL` | eV |
| N_exc | 相对中性晶胞的过剩电子 | `EXCESS_ELECTRONS` | e |
| δ_W | \|W − W_T\| 的容差 | `FERMI_TOL` | eV |
| C_0 | 初始微分电容 | `CAPACITANCE_INIT` | e/V |
| C_min, C_max | 割线电容的下限和上限 | `CAPACITANCE_MIN`, `CAPACITANCE_MAX` | e/V |
| ΔN_max | 每次更新 N_e 的最大变化 | `NE_STEP_MAX` | e |
| k_upd | 两次外层更新之间的离子步数 | `UPDATE_INTERVAL` | 步 |
| n_out | 每个更新步里 N_e 外层更新的最大次数 | `MAX_OUTER_UPDATES` | – |
| g_ff | 带到下一个离子步的电子数修正的前馈增益 | `CARRY_GAIN` | – |
| τ_pot | 恒电势器周期 | `POTENTIOSTAT_PERIOD` | fs |
| **`PLANAR_COUNTER_CHARGE (CP2K)`** | | | |
| d_pl | 对电荷平面到晶胞边界的距离 | `DIST_EDGE` | Å |
| σ_pl | 平面对电荷的 Gauss 宽度 | `GAU_C` | Å |
