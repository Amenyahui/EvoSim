# evoSim v0.1.7 — 参数与输入指南

本文以当前 `slim-qtl-4`、`shared-chemical-1` 和 `expression-tuning-2` 为准；模型的公式和阶段顺序见 [模型说明](model.md)，简化及证据范围见 [假设边界](assumptions.md)。以下数值是配置中实际使用的**理论示例**，不是已校准的自然系统参数。每一代是一轮模型生活史，不预设对应的日历时长。

## 配置入口与输入

以下路径相对于 `evo_sim/models/plant_insect/`。父 YAML 先载入，子层指定的 CSV 按表替换，YAML 显式参数再覆盖；路径相对于声明路径的 YAML，映射递归合并，列表整体替换。

| 入口 | 作用和默认初始条件 |
| --- | --- |
| [configs/v017_modules.yaml](../evo_sim/models/plant_insect/configs/v017_modules.yaml) | v0.1.7 公共配置入口，继承 `model_parameters.yaml`，以 `v017_no_flowering/` 的遗传和参考表替换旧表；继承公共物种和竞争表。P1、P2 各 200 株，H1 为 40 只；植物 `lambda=2`，昆虫潜在卵 `lambda=4`；植物同池 `alpha=0.02`、跨池 `alpha=0.005`。两植物池创始等位频率相同。这是理论默认，不自动构成正式实验设计。 |
| [configs/v017_single_ancestor.yaml](../evo_sim/models/plant_insect/configs/v017_single_ancestor.yaml) | 继承上表并替换物种、竞争和创始频率表，供单祖先兼容池机制检查：P1 40 株、H1 40 只，`lambda` 分别为 2 和 4，植物 `alpha=0.05`。运行轮数覆盖为 3。它的竞争系数和初始数量不能当成公共配置默认值。 |
| [configs/v017_200gen_400insects.yaml](../evo_sim/models/plant_insect/configs/v017_200gen_400insects.yaml) | 用户安排的另一种设计，继承公共配置，覆盖为 200 轮、P1/P2 各 200 株、H1 初始 400 只及 3 个重复；不改变功能参数。 |

公共配置从 [inputs/species.csv](../evo_sim/models/plant_insect/inputs/species.csv) 读取兼容池、初始数量、资源、需求与繁殖率，从 [inputs/competition.csv](../evo_sim/models/plant_insect/inputs/competition.csv) 读取完整植物竞争矩阵。`species_id` 的 P/H 编号约束双亲兼容性，既不是只读谱系标签，也不自动对应不同参考化学型。

`inputs/v017_no_flowering/` 包含当前位点、等位、创始频率及参考表；[reference_traits.csv](../evo_sim/models/plant_insect/inputs/v017_no_flowering/reference_traits.csv) 的 `species_id` 是参考面板 ID，默认 `chem_A=(0.25,0.65,0.25)`、`chem_B=(0.65,0.25,0.65)`，可与兼容池数量不同，只供反事实诊断，不参与真实交配或取食。`inputs/v017/` 是历史含花期草稿，不能冒充当前三轴目录；`model_parameters.yaml` 本身仍为 v0.1.6 / `slim-qtl-2`，不能只改版本号混用遗传表。

## 可遗传坐标和遗传输入

| 类群 | 原发坐标顺序，从 0 开始 | 公共配置的位点数 |
| --- | --- | --- |
| 植物 | `x_N, x_S, x_V`：防御独有、共享、挥发独有的化学功能轴 | 6 QTL，三轴各 2 |
| 昆虫 | `met_E1, met_E2, met_mu1, met_mu2, sense_E1, sense_E2, sense_mu1, sense_mu2`：处理和感觉各自的表达量与调谐 | 16 QTL，八轴各 2 |

坐标是无量纲模型尺度；默认各轴背景为 0，加性拷贝效应乘 `genetics.phenotype_map.effect_scale=0.5` 后加背景，按 `bounds=[0,1]` 截断。表达及调谐是可演化个体性状，下面列出的核宽、类型特征和价性等是运行中固定参数。派生的 `I` 与 `A` 不是独立 QTL 性状，亦不保证分别落在 `[0,1]`。

| 输入或键 | 当前含义、单位和默认值 |
| --- | --- |
| `loci.csv` | 每轴两个位点位于同一抽象连锁群的 `0`、`10 cM`；各轴的连锁群独立。`target_j` 指定受作用的轴，`mutation_mean_j`、`mutation_sd_j` 给出该轴效应增量的分布。共 6 植物与 16 昆虫位点。 |
| `alleles.csv`、`founder_frequencies.csv` | 创始等位效应及每兼容池、每位点的创始等位频率；默认每位点两种等位各 0.5。植物 P1/P2 的频率相同；单祖先表只保留 P1/H1。等位 ID 是账本标识，不是 DNA 序列。 |
| `mutation_rate`、`mutation_mean_j`、`mutation_sd_j` | 默认每位点每拷贝每代突变概率 0.01，目标轴均值 0、标准差 0.02 个效应单位；是理论设定，不是每核苷酸突变率。默认植物 12、昆虫 32 份位点拷贝进入突变过程；均以 μ 表示时，至少一次突变概率分别为 `1−(1−μ)^12`、`1−(1−μ)^32`。 |
| `mutation_enabled.plant/insect` | 分别控制新突变；公共配置均为 `true`，批次的 `mutation_off` 处理将两者设为 `false`。关闭新突变不关闭站立变异、选择或重组。 |
| `ploidy`、`inheritance`、`recombination_model` | 分别为 `2`、`biparental_meiosis`、`slim_haldane_marker_map`；遗传执行器为 SLiM，默认命令 `slim`，`EVOSIM_SLIM` 可覆盖路径。 |

遗传 CSV 为 UTF-8，公共向量列编号 0–7：植物仅填 0–2，3–7 必须为空；昆虫八列均填写。位点目标、效应和突变向量须与本类群宽度一致，不能把植物空列补成第四至第八轴。

## 化学与昆虫功能模块

这些参数在一次运行内固定；`x`、模块表达 `E` 和调谐 `mu` 随遗传而变化。所有化学量和特征均为归一化模型单位，并非测得的浓度、流量或分子结构距离。

| 键 | 公共默认值及合法范围 | 作用阶段和含义 |
| --- | --- | --- |
| `mapping`、`n`、`eta_shared` | `shared_chemical_modules`、`2`、`1.0`；`eta_shared` 有限且在 `[0,1]` | 植物输出固定为 `D=(x_N,x_S)`、`V=(eta_shared*x_S,x_V)`；`n` 是两条防御通道，`eta_shared` 无量纲，1 是归一化参考，不是 100% 排放率。 |
| `functional_modules.chemical_features` | `[0.0,1.0]`；两个有限数，各在 `[0,1]` | 两种化学类型的固定特征坐标，与各类型的 `D`/`V` 数量分开；供处理与感觉核读取。 |
| `processing_width` | `0.25`；有限且大于 0 | 高斯处理核的固定宽度，单位为化学特征坐标；与可演化的 `met_mu` 共同决定幼虫处理分配。 |
| `processing_idle` | `0.1`；有限且大于 0 | 无量纲闲置通道；以 `T/(idle+Σ全目录T)` 分配各处理模块的表达容量，允许未分配部分并限制总容量。 |
| `sensory_width` | `0.25`；有限且大于 0 | 高斯感觉核的固定宽度，单位为化学特征坐标；与可演化的 `sense_mu` 共同决定母体气味激活。 |
| `sensory_half` | `0.2`；有限且大于 0 | 感觉激活量的半饱和常数，单位与模型中的 VOC 加权激活相同；表达量作为饱和后的增益。 |
| `sensory_valence` | `[1.0,-1.0]`；两个有限数 | 两个感觉模块固定的行为读出系数，分别为正和负；不是可演化受体类别或真实神经网络。 |
| `plant_cost`、`voc_cost` | `0.125`、`0.05`；非负有限 | 植物母源产量乘子分别按 `x_N+x_S`、`x_V` 收线性成本；系数为对应归一化性状单位的倒数，`x_S` 的基本成本只计一次。 |
| `insect_cost`、`sensory_cost` | `0.125`、`0.05`；非负有限 | 昆虫母源潜在卵数分别按 `met_E1+met_E2`、`sense_E1+sense_E2` 收线性成本；系数为相应表达单位的倒数，不对派生 `I` 再收费。 |

昆虫母源乘子为 `max(0,1−insect_cost*Σmet_E−sensory_cost*Σsense_E)`；植物乘子为 `max(0,1−plant_cost*(x_N+x_S)−voc_cost*x_V)`。它们不筛选父本，亦不是能量守恒模型。默认 `B`、`rho_voc`、`voc_bounds` 均为 `null`，属于旧映射的退役键，不能把它们替换成当前 `eta_shared` 的同义项。

处理与感觉读出按索引共用 `chemical_features[k]`，但 `D_k` 与 `V_k` 不因此代表同一种真实化合物；共享性状 `x_S` 在上述映射中分别贡献 `D_2` 与 `V_1`。

## 产卵选择与干预

| 键 | 公共默认值及含义 |
| --- | --- |
| `candidate_size` | `8`，正整数；每枚潜在卵从搜索半径内植物无放回抽取至多 8 株候选。 |
| `olfaction_mode` | `template`；当前协议中该名称表示使用模块效用，**不是**旧版模板距离。`neutral` 时所有候选具有相同基线效用。 |
| `b`、`tau`、`outside_utility` | 分别为 `1.0`、`1.0`、`0.0`，均为模型效用尺度；候选效用乘非负 `tau`，拒绝选项使用 `outside_utility`。`tau=0` 时候选之间无效用差异，拒绝仍由自身效用决定。 |
| `choice_mode` | `outside`：可选择拒绝全部候选；`forced`：只要候选非空就选其中一株。无候选时均无寄主。 |
| `shuffle_voc` | `false`；若启用，每轮打乱分配给植物的产卵气味线索，不改变其防御、植物成本、真实资源或成虫配偶池。 |
| `sigma`、`weights` | `0.2`、`[0.5,0.5]`；分别仍要求正值与非负且和为 1 的两维向量，以维持输入验证和旧协议兼容，**不参与当前模块效用计算**。 |

母体以自身成年位置和感觉模块为每枚潜在卵选择；子代自身 `I` 之后作用于幼虫取食。当前配置既无 `flowering_kappa` 或植物花期状态，也无 `adult_settlement_enabled` 开关；含花期的 `slim-qtl-3` 和含交配前定居的 `expression-tuning-1` 按原身份保留供历史读取或审计，不能当作当前运行合同。

## 资源、存活、繁殖和空间

| 键 | 公共默认值及作用 |
| --- | --- |
| `resource_capacity` | 每株 P1/P2 均为 `1.0` 个 `resource_unit`，可为零；该单株量供所有本轮幼虫共同消费。 |
| `max_feeding_demand` | 每只 H1 幼虫 `1.0` 个 `demand_unit`，须大于零；防御后需求 `d=dmax*q` 与同株幼虫按比例分配，实际摄入 `x` 对应的收益为 `g=x/dmax`，分母不是 `d`。 |
| `resource_unit`、`demand_unit` | 均为 `normalized_edible_dry_biomass`；输入必须使用相同受支持单位，也可均为克或毫克可食干生物量，不能混用。默认量没有经验质量标定。 |
| `feeding_window`、`demand_window` | 均为 `one_generation_larval_feeding`，资源及需求必须是相同的幼虫取食窗口。 |
| `sP0`、`sI0`、`unmatched_insect_fitness` | 分别为 `0.5`、`0`、`0`；`sP0` 是全食害时的植物存活基线，植物按食害一次抽样；幼虫存活为 `g`，拒绝产卵不生成幼虫。现有成虫不再次抽幼虫存活。 |
| `plant_offspring_lambda` | 公共表中 P1/P2 均为 `2.0`，单位为每合格母株每繁殖轮的有效植物后代基准；实际 Poisson 均值乘植物母源成本因子与局部密度因子。母株须存活且有合格父本。 |
| `insect_offspring_lambda` | H1 为 `4.0`，单位为每合格雌成虫每繁殖轮的潜在卵基准；实际 Poisson 均值乘该母体表达成本因子，卵之后可能被拒绝或幼虫死亡。不能将潜在卵等同成年补充。 |
| `plant_competition_alpha` | 公共表按受竞争种为行、竞争种为列；同池 `0.02`、跨池 `0.005`，单位为每个局部植株的倒数。局部数量含母株，`phi=1/(1+Σalpha*N)` 折减母源产量，不是落点着生概率。 |
| `sexual_reproduction` | 严格异交、不自交的两性植物及雌雄异体昆虫；父本从局部同池候选均匀有放回抽取。`insect_founder_female_probability` 和 `insect_offspring_female_probability` 均为 `0.5`，作用于各自的独立性比抽样。 |
| `grid_shape`、`boundary` | `[10,10]`、`periodic`；空间距离为周期格上的切比雪夫距离，格无经验长度。`initial_plant_cells` / `initial_insect_cells` 默认为 `null`，表示创始位置随机分配。 |
| `search_radius`、`competition_radius` | 分别为 2、1 格；前者用于母体逐卵找候选植株，后者计算母株局部竞争。 |
| `pollen_radius`、`insect_mate_radius` | 分别为 1、2 格；确定植物存活父本及成年昆虫同池异性候选，均不使用未来卵寄主筛父本。 |
| `plant_dispersal_radius`、`insect_dispersal_radius` | 均为 1 格；植物子代从母株格扩散，成功幼虫从出生寄主格扩散至下一轮成年位置。 |

`species.csv` 的初始数和上述 `lambda`、资源值可被明确的子配置覆盖；`insect_pop: null` 表示由各昆虫池初始数求和。资源为零时实际消费与幼虫收益为零；植物食害比例在零资源时按零定义。半径是非负整数格数，网格和局部性不是对实测距离的承诺。

## 运行身份、观察和输出

公共入口继承 `simulation_metadata.num_steps_per_episode=5`、`record_every=1`、`seed=20261021` 和 1 个批次重复；批次包含默认突变与 `mutation_off` 两种处理。模型标识、展开配置、输入及源码指纹应按各次运行的元数据读取，不以文档描述覆盖历史结果。`genetics.slim.timeout_seconds=120` 是计算时间上限；`max_individuals_safety=5000`、`max_alleles_safety=100000`、`max_genotype_values_safety=4000000` 为计算保护，触发停止而非生态承载量或人为截断出生。

[analysis/v017_lineages.yaml](../evo_sim/models/plant_insect/analysis/v017_lineages.yaml) 只控制已保存结果的观察，不进入模拟规则或随机过程：

| 观察键 | 默认值、单位与含义 |
| --- | --- |
| `distance_threshold` | `0.25`，无量纲，范围 `[0,1]`；同兼容池无相位等位拷贝不匹配距离不超过阈值的个体相连，连通分量可能因桥接而合并。 |
| `minimum_cluster_size` | `2` 个个体，整数且至少为 2；更小的分量标为 `unassigned`，不代表灭绝。 |
| `persistence_generations` | `3` 代，正整数；候选谱系持续标记的观察窗口，不能证明自然物种形成。 |

这些阈值、参考面板和汇总网络都不改变交配、取食、存活或繁殖。旧 [v014_analysis.yaml](../evo_sim/models/plant_insect/analysis/v014_analysis.yaml) 是旧版汇总的观察设计，不能把其中试验资源及边阈值当作 v0.1.7 生物参数；新生活史请使用 [输出说明](results.md) 中对应的阶段汇总。

开发、冒烟、机制检查和审计批次写入 `codex/validation/`；批量默认路径为 `codex/validation/v<version>/run<N>`。仅用户明确安排的正式科研实验显式使用 `--out results/...`。v0.1.6 旧输入合同、配置和历史核定完整保留于 [参数归档](history/parameters_v0.1.6.md)，旧数据、源码指纹与运行身份不随文档迁移而改写。
