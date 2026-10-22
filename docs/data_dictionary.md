# 数据、输出与结果管理

模拟结果保存在 `results/<批次>/`，原始数据、衍生汇总和图片分别组织。模拟输出不能作为实测校准数据。

```text
settings.json                       实际参数、种子设计和源码哈希
 data/status.json                   整批运行状态
 data/run_status.csv                各处理、重复及终止代
 data/<case>/repXX/
   metadata.json / status.json       配置、版本、随机流、文件哈希、终止原因
   populations.csv                   各物种逐代人口行（N_start/N_next）
   plants.csv                        个体植株供给、共享与受损状态
   networks.csv                      物种对互作边与实际转移（含 candidate_exposure）
   events.csv                        单虫遭遇、候选、选择与取食事件
   parents.csv                       个体存活、繁殖成本与期望后代
   offspring.csv                     已抽子代及亲缘
   diagnostics.csv                   固定祖先参考面板
   genetic_summary.csv               逐代种群性状/等位基因汇总（record_type 区分）
   genetic_loci.csv                  位点定义、功能模块、目标性状与突变参数
   allele_catalog.csv                等位基因目录
   mutations.csv                     突变事件记录
   genotypes.npz / phenotypes.npz    基因型和表型
   snapshots/gXXXXXX.npz            逐代个体状态
 data/summary/                      衍生汇总表、作图数据及 <脚本名>_metadata.json 来源元信息
 figures/                           PNG/PDF 结果图
```

测试与数据检查报告放 `tests/reports/`；分析代码放 `evo_sim/models/plant_insect/analysis/`；模型说明放 `docs/`。
`metadata.json` 的文件哈希检查仍保留。历史代码、基线清单及源码副本不再保存。

## 关键字段

| 文件/字段 | 含义与边界 |
| --- | --- |
| snapshots: D / I / A | 植物防御、昆虫抗性、气味偏好中心；各行是个体，各列是功能轴 |
| plant_voc_independent_trait | 可遗传独立气味性状；ρ=1 时不影响 V_emit |
| V_emit / V_cue | 实际释放性状与用于选择的线索；打乱处理只改后者 |
| *_ids / *_species | 类群内跨代唯一 ID / 从 0 开始的固定物种标签 |
| events: host_row / host_id | 当前数组行 / 植物稳定 ID；拒绝时均为 -1 |
| q | 潜在摄食比例；无寄主为 NA，不应与接受后 q=0 混同 |
| d / x / g | 潜在需求 dmax·q / 实际组织转移 / x÷dmax |
| plants: R / T / gamma / X / h | 单株供给 / 总需求 / 共享缩放 / 实际移除量 / X÷R |
| parents: S_i / F_i / W_i | 存活概率 / 相对繁殖成本因子 / S_i·F_i，W 不是后代数 |
| phi_i / lambda_s / μ | 密度因子 / 基础成功招募尺度 / 存活者 Poisson 均值 λFφ |
| expected_offspring / offspring_variance | 条件于当代状态的期望与方差，含存活随机性 |
| offspring_count | 已抽到的整数子代数；超过计算上限时也保留完整计数 |
| parent_genotype_id | 历史字段名，实际为 P/H 前缀加亲代个体 ID，不是基因型编码 |
| *_mutated_* | 该位置被选中突变；不保证数值改变，例如突变步长为零 |
| N_start / N_next | 本代取食前人数 / 已提交的下一代人数；未提交时 N_next 为 NA |
| extinct_at_generation | 首次灭绝代；零为吸收态，初始零种为第 0 代 |
| N_accept | 接受该寄主物种的虫数，即使 q=0 也可增加 |
| Q_actual | 物种对实际转移组织量 Σx；不是接受次数 |
| L_per_insect | Q_actual÷该虫种当代人数，虫种灭绝时为 NA |
| snapshots: plant_cells / insect_cells | 每个体的栅格行列坐标；仅空间配置存在 |
| events: search_row / search_col | 本代昆虫搜索起点；第 0 代由独立初始分布产生 |
| events: candidate_plant_ids / candidate_count | 实际抽到的候选植株 ID JSON 列表及数量；可与快照核对局部遭遇 |
| events: reason / host_assignment | `no_candidate`、`rejected`、`accepted`；当前自由初态与后续均为 `local_choice` |
| networks: candidate_exposure | 该虫种所有个体的候选集中，该植物种出现的次数 |
| diagnostics.csv | 固定祖先参考面板上的偏好与无限资源利用率，不是当代实际遭遇 |

个体身份使用 `(kind, individual_id)`；当前数组行号不可跨代追踪。所有物种槽位保留，包括初始为零或已灭绝物种。NA 表示不可定义/未抽样/未提交，不能全局填零。

## 状态和单位

`completed` 和 `all_species_extinct` 是正常终止。`invalid_parameter` 是输入不合法；`arithmetic_error` 是无法表示的算术；`computational_limit` 是资源/整数采样范围限制，不是生物灭绝。计算失败或提前停止必须单独报告，不能当作存活到最后一代或从统计中无声删除。

R、dmax、d、x、X、Q_actual 使用同一取食窗口和干组织质量尺度。g、h、q、S、F、φ 无量纲。α 在固定样地以人数建模时为每株系数；把克换成毫克不应同时缩放 α。λ 是无相关成本/拥挤时每个存活亲代的预期成功招募数，并非原始产卵量。

审计重建取食与人口账本、边计数、父子关系、局部候选集、坐标和哈希；它不等于实证校准。长期运行的累计内存限制、专化判据和适用生命周期见科学审核；当前初态方法见 [host_network.md](host_network.md)。
