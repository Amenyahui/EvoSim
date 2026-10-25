# evoSim 参数配置指南

模型输入由 **YAML 通用参数 + CSV 数据表**组成。[`model_parameters.yaml`](../evo_sim/models/plant_insect/configs/model_parameters.yaml) 保留模型规则、运行参数和输入路径；物种、位点、等位效应、初始频率等逐条记录放在 [`inputs/`](../evo_sim/models/plant_insect/inputs/) 的 CSV 中。实验文件通过 `extends` 继承公共值，Python 读取表格、计算和检查约束，不用隐藏的模型默认值补齐缺项。

## 1. 配置文件分工

| 文件（相对于 `evo_sim/models/plant_insect/`） | 内容 |
| --- | --- |
| `configs/model_parameters.yaml` | 通用遗传、生态、空间、交配、观测、随机种子、安全上限、基础实验设置和 CSV 路径 |
| `inputs/species.csv` | 初始物种、数量、资源/需求和基础繁殖率；当前 3 行数据 |
| `inputs/competition.csv` | 植物竞争系数及来源；当前 4 行数据 |
| `inputs/loci.csv` | 位点、连锁组、作用坐标和突变参数；当前 16 行数据 |
| `inputs/alleles.csv` | 创始等位基因及效应；当前 48 行数据 |
| `inputs/founder_frequencies.csv` | 各物种各位点的初始等位频率；当前 48 行数据 |
| `inputs/reference_traits.csv` | 标准化观测用的固定参考植物性状；当前 2 行数据 |
| `configs/sexual_qtl_development.yaml` | 继承公共配置，运行 5 代 |
| `configs/sexual_qtl_smoke.yaml` | 继承开发配置，运行 1 代 |
| `configs/sexual_qtl_200gen.yaml` | 继承开发配置，运行 200 代、每处理 3 次重复 |
| `configs/initial_conditions.yaml` | 7 组初始数量对照、0 代、4 次重复；供初始条件诊断使用 |
| `config.yaml` | 调度流程、状态路径、运行设备等框架配置；单独使用它不能构成完整的模型实验 |

公共 YAML 从 908 行缩减到 139 行，数据转入六张表，所有数值沿用原演示场景，没有依据种群变化重新调参。等位效应系数仍为 `0.5`，遗传表型和 VOC 截断范围仍为 `[0, 1]`，突变效应均值仍全为 `0`。数值属于理论示例，不是经验标定。

## 2. CSV 列与编辑方法

所有表都是普通逗号分隔 CSV，一格一个值，没有 JSON 列表。可用文本编辑器或 Excel 修改，保存为 UTF-8 CSV；读取支持带 BOM 的 UTF-8。必须保留表头和所需列，不要把空白数值当作 0。

| 文件 | 列 | 编辑要点 |
| --- | --- | --- |
| `species.csv` | `species_id, organism, initial_size, resource_capacity, max_feeding_demand, offspring_lambda, source_kind, source` | `organism` 为 `plant` 或 `insect`；植物填资源量、留空需求，昆虫填需求、留空资源量；数量为非负整数。`source_kind/source` 记录基础繁殖率来源，资源来源仍见 YAML 中 `parameter_provenance.resources` |
| `competition.csv` | `recipient, competitor, alpha, source_kind, source` | 一行表示竞争者对受影响植物物种的系数；必须覆盖全部植物物种对，包括同种竞争，`alpha≥0` |
| `loci.csv` | `locus_id, organism, functional_module, chromosome_id, position_cM, mutation_rate`，以及各坐标的 `target_j, mutation_mean_j, mutation_sd_j` | `target_j` 为 0/1，表示是否作用于坐标 j；非目标坐标的均值与标准差必须为 0；突变率在 [0,1]、标准差非负 |
| `alleles.csv` | `allele_id, locus_id`，以及各坐标的 `effect_j` | 一行一个等位基因；效应可正可负，非目标坐标必须为 0；ID 按行从 0 连续递增，`locus_id` 必须存在 |
| `founder_frequencies.csv` | `species_id, locus_id, allele_id, frequency` | 一行一条初始频率；每个物种必须覆盖所属类群的所有位点，每个物种—位点组合的频率和为 1；同一组合内不能重复等位 ID |
| `reference_traits.csv` | `species_id`，以及 `D_1...D_n, Z_1...Z_K` | 每个植物物种一行、值在 [0,1]；这是固定的观测参考面板，不直接给创始个体赋表型 |

`P1, P2, …` 和 `H1, H2, …` 仍是预先定义的植物、昆虫物种槽位，编号分别从 1 连续排列；不是程序自动识别的物种或演化生成的物种。物种表的行顺序不改变编号，数量为零的物种也保留记录。更换物种集合时，要同步更新竞争表、初始频率和植物参考表。

遗传 CSV 中的坐标 j 从 **0** 到 **n+K−1**，其中 n 为 YAML 的 `n`，K 为 B 的行数。植物按 D 后 Z 排列，昆虫按 I 后 A 排列。当前 `n=2, K=2` 时：

| 遗传坐标 | 0 | 1 | 2 | 3 |
| --- | --- | --- | --- | --- |
| 植物 | D₁ | D₂ | Z₁ | Z₂ |
| 昆虫 | I₁ | I₂ | A₁ | A₂ |

例如 `target_0=1`、其余 target 列为 0，表示该位点只作用于植物 D₁ 或昆虫 I₁；`effect_0` 是该坐标上的等位效应。`functional_module` 是功能注释，实际作用坐标由 target 列决定。改变 n 或 B 的行数时，需同步调整遗传 CSV 的坐标列、背景向量和参考表列。

常见修改入口：

- 改 P1/P2/H1 的初始数量、资源或基础繁殖率：修改 `species.csv` 对应行。
- 改突变概率、效应均值/标准差或连锁位置：修改 `loci.csv` 对应行。
- 改等位基因功能效应：修改 `alleles.csv` 的 effect 列。
- 改物种的初始遗传组成：修改 `founder_frequencies.csv`，保持每物种每位点频率和为 1。

为了精确复现，保留位点、等位记录及同一物种—位点内频率记录的顺序；这些顺序参与当前抽样和数据排列。重复键、缺列、非有限数字、未知引用或非法频率会被拒绝，不自动补齐或归一化错误输入。

## 3. 展开后的参数及其意义

读取后，CSV 会转换为现有 `parameters` 结构；下表列的是**展开后的配置键**，不表示这些数据仍须手写在 YAML 中。初始数量、资源和繁殖率来自 species 表，竞争矩阵来自 competition 表，遗传记录来自三张遗传表，固定观测面板来自 reference_traits 表；其余通用参数保留在 YAML 中。

| 类别 | 展开后的配置键 | 含义或约束 |
| --- | --- | --- |
| 初始种群 | `plant_sizes`, `insect_sizes` | 每个物种的初始数量，非负整数；保留数量为零的物种位置 |
| 昆虫总数 | `insect_pop` | 默认明确写 `null`，总数由 `insect_sizes` 派生；若填数字必须等于其和 |
| 性状维度 | `n`, `B`, `weights`, `mapping` | `n` 为 D/I 维数，B 的行数决定 Z/A/VOC 维数；B 必须显式提供，权重非负且和为 1。`mapping` 为保留的旧标签，实际映射始终由 B 决定，不再自动生成矩阵 |
| 防御与解毒成本 | `plant_cost`, `insect_cost`, `voc_cost` | 非负的性状成本系数 |
| 生存与生活史 | `sP0`, `sI0`, `unmatched_insect_fitness`, `insect_life_history` | `sP0` 在 [0,1]；当前自主选择—取食—繁殖模型要求后三项分别为 0、0、`self_feeding_pre_reproductive`；不包含母体产卵选择—幼虫发育过程 |
| 寄主选择 | `candidate_size`, `olfaction_mode`, `shuffle_voc`, `choice_mode` | 候选数量、气味选择模式、气味打乱处理、是否允许拒绝所有寄主 |
| 气味效用 | `sigma`, `b`, `tau`, `outside_utility` | 距离核宽度 σ>0、截距、判别强度 τ≥0、外部选项效用 |
| VOC 映射 | `rho_voc`, `voc_bounds` | 防御与独立气味性状的混合权重，以及 VOC 截断区间 |
| 取食资源 | `resource_capacity`, `max_feeding_demand` | 各植物单代资源量、各昆虫最大需求，必须正且有限 |
| 资源量纲 | `resource_unit`, `demand_unit`, `feeding_window`, `demand_window` | 双方必须采用相同可食干物质量纲和时间窗口 |
| 人口过程 | `plant_offspring_lambda`, `insect_offspring_lambda`, `plant_competition_alpha` | 各物种基础母方繁殖率和植物邻域竞争矩阵 |
| 模型协议 | `model_version`, `life_history`, `demography` | 明确声明现已实现的非重叠世代、Bernoulli 存活和 Poisson 出生等过程 |
| 场景与来源 | `illustrative_only`, `study_site`, `parameter_provenance` | 是否为理论示例、研究场景、人口和生态参数来源 |
| 空间 | `spatial.grid_shape`, `boundary` | 网格形状；当前仅支持周期边界 |
| 空间作用范围 | `spatial.search_radius`, `competition_radius`, `plant_dispersal_radius`, `insect_dispersal_radius` | 搜寻、竞争和两类群扩散半径，单位为格 |
| 初始位置 | `spatial.initial_plant_cells`, `initial_insect_cells` | `null` 表示按现有均匀位置规则抽样；也可逐个体指定坐标 |
| 遗传协议 | `genetics.mode`, `schema_version`, `ploidy`, `inheritance`, `recombination_model` | 显式二倍体、有性双亲、Haldane–Poisson 重组协议；当前不支持任意倍性 |
| 遗传映射 | `genetics.background`, `phenotype_map.effect_scale`, `phenotype_map.bounds` | 两类群背景向量、所有等位效应总和的非负尺度系数、遗传表型截断区间 |
| 染色体与位点 | `genetics.loci[]` 中的 `locus_id`, `organism`, `functional_module`, `chromosome_id`, `position_cM`, `targets` | 位点身份、类群、功能注释、连锁组、遗传位置和作用坐标；位点数量由列表长度决定 |
| 突变 | 每个位点的 `mutation_rate`, `mutation_mean`, `mutation_sd`，以及 `genetics.mutation_enabled` | 逐拷贝逐位点事件概率，以及效应增量的均值/标准差向量；未声明的目标坐标必须均为 0 |
| 创始遗传组成 | `genetics.alleles`, `founder_frequencies` | 等位 ID、所属位点、功能效应，以及每物种每位点的初始频率；每个位点频率和为 1 |
| 遗传来源 | `genetics.provenance` | 遗传映射、图谱、突变、创始组成及计算限制的依据 |
| 交配规则 | `sexual_reproduction.plant_mating_system`, `plant_selfing`, `insect_mating_system`, `father_sampling` | 当前为植物严格异交、昆虫雌雄异体、合格父本均匀有放回抽样 |
| 交配半径与性比 | `sexual_reproduction.pollen_radius`, `insect_mate_radius`, `insect_founder_female_probability`, `insect_offspring_female_probability` | 交配空间范围、创始和子代雌性概率；另有 `provenance` 记录假设来源 |
| 观测 | `record_observations`, `observation_reference_traits`, `record_every` | 是否记录标准化观测、固定参考植物 D/Z 面板；完整账本目前要求逐代记录，即 `record_every: 1` |
| 计算保护 | `max_individuals_safety`, `genetics.max_alleles_safety`, `genetics.max_genotype_values_safety` | 个体、等位目录和累计基因型存储上限，超限停止而不截断出生数 |
| 随机复现 | `seed` | 主随机种子；第 r 次重复使用 `seed + r - 1`，子流编号为固定复现协议 |
| 实验时间 | 顶层 `simulation_metadata.num_steps_per_episode` | 模拟代数；转入 Config 的 `generations`，不另设第二份模型默认值 |
| 批量实验 | 顶层 `batch.purpose`, `replicates`, `cases[].name`, `cases[].overrides` | 实验目的、重复次数、处理名称及参数差异 |

归一化性状范围是模型的量纲约定。两个可配置截断区间均需满足 `0 ≤ lower < upper ≤ 1`；它们允许在当前归一化尺度内调整截断，不支持把实测浓度直接放进现有归一化公式。

`effect_scale` 默认 0.5 是二倍体效应尺度约定，不是显性系数。突变均值默认 0 表示无方向偏置；修改为非零将引入定向突变假设，需同步说明依据，不能把这种方向性解释成自然选择的结果。

## 4. YAML 引用与实验覆盖

公共 YAML 使用顶层 `input_files` 引用六张表：

```yaml
input_files:
  species: ../inputs/species.csv
  competition: ../inputs/competition.csv
  loci: ../inputs/loci.csv
  alleles: ../inputs/alleles.csv
  founder_frequencies: ../inputs/founder_frequencies.csv
  reference_traits: ../inputs/reference_traits.csv
```

路径相对于**声明该路径的 YAML 所在目录**解析，与运行命令时的工作目录无关。独立实验可以复制需要改动的 CSV，再只覆盖对应路径；未列出的输入继续继承。例如把另一份频率表保存为 `inputs/founders_alternative.csv`，在 `configs/` 下的新实验中写：

```yaml
extends: model_parameters.yaml
input_files:
  founder_frequencies: ../inputs/founders_alternative.csv
```

通用参数和实验差异仍可在 YAML 中覆盖：

```yaml
extends: model_parameters.yaml
simulation_metadata:
  num_steps_per_episode: 20
parameters:
  seed: 20261025
  spatial:
    search_radius: 3
  genetics:
    phenotype_map:
      effect_scale: 0.5
    mutation_enabled:
      plant: false
      insect: true
batch:
  replicates: 2
  cases:
  - name: plant_mutation_off
    overrides: {}
```

假设文件名是 `my_experiment.yaml`，从项目根目录运行：

```powershell
python -m evo_sim --config evo_sim/models/plant_insect/configs/my_experiment.yaml
```

读取顺序为：父配置 → 本层 CSV 整表替换 → 本层显式 YAML 参数 → 后续子配置；运行时最后应用 `batch.cases[].overrides`。因此显式 YAML 参数可以覆盖刚导入的 CSV 值，修改 CSV 未生效时应先检查实验里是否仍有对应覆盖。

1. 更换 CSV 路径会替换对应整张表的目标数据，不把新旧表逐行合并。引用一张更小的频率表不会保留旧表中删掉的记录，但仍须通过完整性检查。
2. YAML 参数字典递归合并，只覆盖指定键；例如关闭突变不会覆盖染色体、位点、等位效应或来源说明。
3. YAML 列表整体替换，不按编号或 ID 合并；标量和 `null` 显式替换原值。直接编辑 CSV 的一个单元格不需要重写整张表。
4. CSV 路径只写在 YAML 顶层 `input_files`，不写进 `parameters` 或 `batch.cases[].overrides`。批次处理继续使用内联参数差异；需要换一套 CSV 时使用另一份实验 YAML。
5. 必填参数缺失、未知参数或非法取值会报错，不从 Python 补默认值。

批次 `settings.json` 保存展开后、处理覆盖前的配置和各处理差异；每次运行的 `metadata.json.config` 保存实际 Config、重复种子及解析后的矩阵等。通过配置文件运行时，`input_sources` 还保存输入表的路径和 SHA-256 指纹。输出已包含展开后的实际数据，独立审计无需重新读取原 CSV；输入路径用于追溯，不代替参数快照。修改 YAML 或 CSV 后应使用新输出目录，已有目录不会覆盖。

初始条件诊断也读取 YAML，可使用默认的 `initial_conditions.yaml` 或通过 `--config` 指定另一份零代实验：

```powershell
python -m evo_sim.models.plant_insect.analysis.diagnose_initial_conditions --out results/initial_conditions_new
```

## 5. 哪些内容仍由代码定义

方程形式、资源守恒、概率合法范围、单位换算（100 cM = 1 Morgan）、数组索引、性别编码和浮点误差检查属于数学定义或实现约束。突变效应均值已在 loci CSV 中显式给出，等位效应尺度和性状截断保留在 YAML；没有在代码中补一套隐藏默认值。

YAML 中列出的协议项会校验可支持的机制。把 `ploidy: 2` 改为 4、把 `immigration: false` 改为 true，并不会自动获得多倍体或迁入模型；当前会明确拒绝。这次修改统一参数来源，没有实现新的物种形成或互作网络机制。

完整内联 YAML 和 YAML + CSV 均须遵守当前 `model_version: 0.1.2` 与 `insect_life_history: self_feeding_pre_reproductive` 契约，不提供旧配置自动迁移。默认取食窗口为 `one_generation_pre_reproductive_feeding`。历史结果没有修改；重建历史结果应使用生成它们时的代码和配置。

## 6. v0.1.2 对照的可比性

`configs/v012*.yaml` 与 `inputs/v012/` 提供小规模已有适配、连锁、无新突变、完整物种重标记、较少预设形成及多物种场景。各配置只覆盖必要的参数或 CSV，不另建配置体系。实际场景及检查命令见 [因果核对报告](v0.1.2_causal_check_report.md)。

- 连锁比较只改图谱；不改等位效应、创始频率、位点行序和突变输入。有限创始样本仍可能存在关联，需用相同位点对检查。
- 无新突变不代表无演化；位点概率是每传递 QTL 拷贝的概率，不是核苷酸率。
- 较少预设的形成场景改变初始变异分布，因此不是与已有适配基线仅差一个参数的处理。须报告性状支持范围、p/q/J 分布及有限样本偏离。
- 完整重标记同步处理种群、资源、繁殖、频率、参考面板、竞争矩阵和来源映射；同一 seed 不保证标签交换前后的轨迹逐点相同。
- `rho_voc` 是表型映射耦合，不能用作连锁开关。

所有新增数值仍是探索性输入。`theoretical_placeholder` 中需文字区分归一化尺度约定和具体探索参数；文献采用过某个数值不使它成为本项目物种的经验估计。
