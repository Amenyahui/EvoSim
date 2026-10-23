# evoSim v0.1.1 — Architecture：Agent、Rule 与 Step 的实现规范

> **适用版本**：evoSim v0.1.1；遗传协议 `sexual-qtl-1`。  
> **文件属性**：唯一目标架构 / Agent 实现与审计契约。下面描述 v0.1.1 **必须实现的最终架构**，不代表附件源码已完成全部实现。  
> **研究类型**：植物—植食性昆虫的理论个体基础共同演化模型（IBM）；不使用真实基因序列。  
> **相互引用**：`model.md` 给出全部公式（G、E、D 系列）；`assumptions.md` 解释假设（A 系列）；`results.md` 规定输出和图件（O 系列）。

## 1. 总体设计与核心约束

evoSim 采用“**Agent 状态 → Observation → Policy → Transition → 世代结算**”组织一个多物种生态遗传仿真器。`Registry` 负责登记可调用规则，`Runner` 决定执行顺序，`Controller` 执行所声明的规则；植物与昆虫 Agent 的遗传和生态状态由统一 `state` 管理。此组织方式参考模块化 Agent 系统；显式遗传、nonWF 式人口学自主控制参考 SLiM，IBM 可重复描述参考 ODD 协议 [R01–R03]。**这些文献支持建模范式，不意味着代码已采用 SLiM 引擎或 PyTorch。**

全局硬约束：

1. 只存在 `sexual-qtl-1`：二倍体、双亲有性繁殖、孟德尔分离、连锁群、Haldane 重组、配子突变、加性 QTL。不得加入克隆或无遗传的可执行分支。
2. 生态互作机制：先局部感知/寄主选择，再有限资源竞争和取食，之后一次存活抽签、配偶搜索、母本产子、双亲遗传和扩散。
3. 同一完整世代内，所有 Agent 在相同起始群体下取食；本代出生的子代不参与本代的取食、竞争、交配。
4. 严格区分 **生存至繁殖阶段**、**存在配偶**、**母方出生数** 和 **父方成功提供配子数**；结果均独立记录。
5. 原子提交：构造完整下一代成功后再替换群体；内存/容量/数值错误不得产生半代群体或孤儿等位基因。

## 2. Agent：定义、状态和生命周期

**Agent** 是具有独立身份、位置、基因型和性状，并能参与互作与繁殖的植物或昆虫个体。`P1、P2、H1` 是固定物种标签，**物种不是 Agent，也不会在 v0.1.1 中动态生成**。

| 类别 | 固定属性 | 动态状态 | 繁殖角色和退出条件 |
|---|---|---|---|
| `PlantAgent`（向量化容器 `PlantPopulation`） | `id`, `species` | `cell`, `genotype[2,L_P]`, `D[2]`, `Z[2]`, `V_emit[2]`（后者可计算缓存）、本代损伤与资源账本 | 理论两性花，可为种子母本和其他植株的花粉父本；必须异交；一代结束全部退出 |
| `InsectAgent`（`InsectPopulation`） | `id`, `species`, `sex`（稳定、出生时抽样） | `search_cell`, 本代 `host_id/host_cell`, `genotype[2,L_H]`, `I[2]`, `A[2]`, `q,d,x,g` | 雌虫为母本、雄虫为父本；两者均需取食并通过生存检验；一代结束全部退出 |

遗传 `genotype` 保存**等位基因 ID**，不是 DNA 字符串。每个个体两份拷贝的语义在出生时固定为 **copy 0=母源、copy 1=父源**。每个物种使用其 `locus_id` 与独立等位目录（或强制全局唯一 ID），同一行的 `id/species/cell/sex/genotype/traits` 必须严格对齐。`sex` 不由 QTL 决定，仅作为本版繁殖角色变量。

**环境 Agent？** 栅格是 `Environment` 的空间单元而非可繁殖 Agent；每格可存在多株植物。植物资源归属于**具体个体**，不能把同格多株的资源提前合并为一个“格子库存”。

## 3. Rule：规则类别、输入与输出

**Rule** 是对输入状态进行确定性计算或一次随机抽样的最小单元。禁止让单个规则私自同时执行多阶段更新；随机规则必须显式接收 RNG；声明其依赖及写入字段。

| Rule 组 | 规则/职责 | 读入 | 写出及对应公式 | 实现边界 |
|---|---|---|---|---|
| 初始化 `Initialization` | `InitializeCommunity` | 配置、seed、起始计数/位点目录 | 个体、物种、性别、图谱、等位目录、初始基因型 | 只接受新协议；所有输出校验 |
| 遗传 `Genetics` | `FormGamete` | 亲本基因型、连锁群、有序 `position_cM`、RNG | 1套单倍体等位 ID；G2–G5 | 各连锁群独立起始，组内按 Poisson 奇偶交换 |
| 遗传 `Genetics` | `MutateGamete` | 已形成配子、位点级 `mu, sigma`、等位目录 | 新等位 ID、突变事件；G7–G8 | 仅新等位变异，旧等位对象不可被原位修改 |
| 遗传 `Genetics` | `Fertilize` | 两个有资质的同种亲本及其配子 | `child_genotype[2,L]`、双亲 ID；G6 | 不能缺少任何父源/母源拷贝 |
| 遗传 `Genetics` | `MapPhenotype` | 子代基因型、等位效应目录 | `D/Z/I/A`、原始表型、裁剪标记；G9 | 纯确定性函数，不直接突变表型 |
| 观察 `Observation` | `ObserveCues` | 植物 D/Z、位置、B、rho | `V_emit, V_cue` 及候选集合；E1、S1 | 不改变任何 Agent |
| 策略 `Policy` | `ChooseHosts` | 候选植物 V、昆虫 A、选择参数 | 每昆虫 `host_id` 或无寄主；E3–E4 | 所有昆虫基于本代固定候选集选择 |
| 结算 `Transition` | `SettleFeeding` | 所有寄主选择、植物资源、D/I | `q,d,T,gamma,x,g,h,S,F`；E5–E11 | 共享资源按比例分配，严格守恒 |
| 人口学 `Transition` | `SampleSurvival` | 所有植物/昆虫 S、RNG | `survived_to_mating`；D1 | 全体一次抽签，后续不得重抽 |
| 交配 `Rule` | `FindMates` | 存活者、性别/物种、空间、F | 合格花粉/雄性父本候选集合；D3–D4 | 物种、角色、空间与 F>0 同时成立 |
| 人口学 `Rule` | `SampleMaternalBirths` | 生存母本、候选池、`lambda,F,phi` | 子代数、条件均值；D5–D9 | 无配偶严格零；雄虫不计母方出生数 |
| 结算 `Transition` | `SexualInheritanceAndDispersal` | `births, parents, gametes, traits, dispersal_rng` | 完整下一代、谱系和传递日志；G2–G9、D10 | 全部成功才整体提交 |
| 输出 `Audit/IO` | `Snapshot/Record/Validate` | 世代前后状态与事件 | `*.csv/*.npz/metadata.json`；见 results.md | 日志和正文指标语义一致 |

Rule 表内的名称为**v0.1.1 建议稳定接口名**，Agent 可合理拆分类，只要输入、事件顺序和科学语义保持一致。禁止把不同阶段合并而造成存活重复抽样或父本来源缺失。

## 4. Step：完整世代的唯一执行顺序

### Step 0 — 创始种群构建（只执行一次）

读取 `sexual-qtl-1` 配置，验证 `L=8`/类群、每个位点的连锁群与 `position_cM`、等位频率、性别参数、空间范围、所有理论参数与来源标签。初始化个体位置；为两个亲源拷贝独立抽始祖等位基因；计算性状。初始化种群中两份拷贝只是“创始基因型的同源副本”，不可声称其历史父母在模拟中出现。

### Step 1 — 感知与寄主选择

- 植物按 E1 输出 VOC；候选寄主通过周期边界距离 S1 构建。
- 每昆虫从本代候选池抽取最多 `candidate_size` 株，并按 E3/E4 的概率选择实际寄主或外部选项。
- 使用“同时选择、统一结算”：前一只昆虫选中寄主不立即使后者从其他昆虫候选池消失。
- 昆虫无寄主时 `host_id=-1`，其 `q,d,x,g=0` 且无法提供合格交配位置。

### Step 2 — 取食与共享资源结算

计算对应 `D/I` 的利用系数 `q` 和潜在需求 `d`，汇总同一植物全部需求 `T`，按 E5–E8 分配实际消耗 `x`，满足 `sum_{h:host(h)=p}x_h ≤ R_p`。计算每植物食害 `h_p`、每昆虫资源收益 `g_h`，再计算 `S^P,S^H,F^P,F^H,phi`。此处 `S` 是之后的抽样概率，不是直接保留亲本进下一代。

### Step 3 — 全体存活检验

对本代每一名植物/昆虫独立执行一次 `Bernoulli(S_i)`，存储布尔向量 `R_i`。**未存活个体不能出现在任一父母候选池中**；不得在后面的 `SampleMaternalBirths` 再抽存活。

### Step 4 — 交配可获得性

- 植物：所有存活植株可以尝试当种子母本；花粉父本必须是**另一株存活的同种植物**且 `F>0`、位置落在 `pollen_radius` 内；严格禁止自交。
- 昆虫：只有存活雌虫可以尝试产子；雄虫必须存活、同物种、`F>0`，且其**本代实际寄主位置**处在雌虫 `insect_mate_radius` 内。
- 周期栅格 `d_infty` 由 S1 定义；候选按 ID 去重，父本可以跨多个母本反复作为父本（不设置父本容量）。
- 本阶段建立母本 `C_i`；没有合格候选的母本标记 `no_compatible_mate`。

### Step 5 — 母本出生数与双亲选择

- 对 `eligible_mother && survived && len(C_i)>0` 的母本执行**一次** `Poisson(mu_i)`；其他 `B_i=0`。
- 对每个实际出生子代，从母本的本代 `C_i` **有放回、等概率**选择一名父本；植物母方与父方永远不为同一 ID。
- 记录 `maternal_id, paternal_id, parent_species, maternal_cell, paternal_cell, maternal_offspring_count, paternal_sired_count`；不要由记录倒推出没有参与的父本。

### Step 6 — 遗传、重组、突变、受精、子代性状

对每一名子代：两个亲本**各自独立**生成配子；每条连锁群等概率选择起始同源拷贝，沿位点之间的遗传距离执行 `K~Poisson(d)` 和奇偶交换；每个传递等位以 `mu_l` 判断突变并创建新 ID；两个配子按母源/父源组成子代 `[2,L]`；据此重算加性性状 `D/Z/I/A`。必要时记录 `raw_trait` 和被截断比例。**绝不能直接在子代表型数组上进行额外随机变异。**

### Step 7 — 位置、性别与全世代原子替换

- 植物种子从**母本所在格**按 `plant_dispersal_radius` 扩散。
- 昆虫子代从**雌虫本代所选寄主所在格**按 `insect_dispersal_radius` 扩散。
- 昆虫子代性别 `female~Bernoulli(p_f)`；不从父母 QTL 推导性别。
- 完成全部新生个体、等位目录、谱系与日志**临时对象**校验后一次性替换全体 Agent；本代所有亲本退出。`N_{t+1}=sum_i B_i`。
- 所有物种总个体数为零时进入吸收性灭绝状态；允许植物存续、昆虫局部或全部灭绝，不自动补充移民。

### Step 8 — 观察与终止

每代记录 `populations, events, parents, offspring, gamete_transmissions, mutations, allele_catalog, genetic_summary, snapshots` 等（字段见 `results.md`）。**终点世代仅记录生态与统计，不再产生“多出来的一代”子代或突变**。必须记录 `update_applied` 以区分统计终点与实际状态更新。

## 5. 文件/组件的建议职责边界

| 目标文件 | 必须承担的职责 | 不应承担的职责 |
|---|---|---|
| `evo_sim/core/registry.py` | 注册和取得 Initialization/Observation/Policy/Transition | 隐式科学方程、随机采样 |
| `evo_sim/core/runner.py` | 正确调度阶段、终止、RNG/配置边界 | 直接定义遗传/取食公式 |
| `evo_sim/core/controller.py` | 规则依赖传递、全代提交保护 | 静默覆盖有冲突的状态字段 |
| `evo_sim/models/plant_insect/helpers/agents.py` | Agent 容器、个体 ID、行对齐和快照 | 生成假双亲谱系 |
| `.../helpers/genetics.py` | 位点目录、配子、重组、突变、受精、性状与审计 | 寄主选择、交配空间资格 |
| `.../helpers/spatial.py` | 周期距离、候选池、竞争邻域、花粉/交配邻域与扩散 | 性状突变、随机重选幸存者 |
| `.../helpers/ecology.py`、`voc.py`、`traits.py` | VOC、选择、取食、成本、共享资源与生存概率 | 繁殖遗传分离 |
| `.../helpers/mating.py`（建议新增） | 双亲角色资格、存活父本候选、父本抽样 | 在配偶不存在时制造后代 |
| `.../helpers/demography.py` | 母本条件繁殖均值、Poisson 出生与矩审计 | 有性遗传/变异本体 |
| `.../substeps/*.py` | 组织纯规则为符合 Step 1–8 的实际事务 | 内嵌第二套科学公式 |
| `.../helpers/config.py`、`configs/*.yaml` | 唯一 schema 和完整参数校验 | 容忍缺失遗传图谱/交配字段 |
| `.../audit*.py`、`.../analysis/*.py` | 谱系/守恒/统计检查和研究指标汇总 | 通过修改结果来使模型通过检查 |

项目只保留 `sexual-qtl-1` 有性遗传、双亲出生及其对应配置和执行入口；不实现其他遗传协议的运行路径。

## 6. 状态转换契约与数据验证

### 6.1 初始化/读取

- `config.genetics.schema_version` 必须等于 `sexual-qtl-1`，`ploidy=2`，`inheritance=sexual`；无位点/无图谱/无性别字段均拒绝加载。
- `chromosome_id`、`position_cM`、`locus_id` 的顺序和唯一性必须明确；**不同连锁群之间独立分配，同组仅按递增遗传位置交换**。
- 每个 `allele_id` 有且仅有一个所属 `locus_id`；等位效应数组维度必须与所属物种性状长度一致；父源谱系 `parent_allele_id` 能溯源。
- 同物种同代个体和父方/母方 ID 必须唯一可追踪；植物配偶 ID≠母方 ID；跨物种永不配对。

### 6.2 终止与原子性

- 运行安全阈值 `max_individuals_safety`、基因型内存上限只是**计算安全约束**，不得当作生态承载量使用。超限应失败并保持完整已提交的上一代状态。
- 建议使用 staged `next_agents`, `new_alleles`, `new_mutations`, `new_events`，预检后一次性提交；`Controller` 不应在一个繁殖步骤中先更新部分个体、随后写失败日志。
- `state.generation` 标记本代已经提交的完整状态；输出文件另记 `evaluated_generation` 与 `transition_to_generation`，避免终止时越界。

### 6.3 随机数流

以主 seed 和固定 `generation/taxon/process` 标签派生互不干扰的 RNG（例如 SeedSequence 派生）；分别用于创始基因型、创始位置/性别、寄主候选、寄主选择、存活、配偶父本、Poisson 出生、母源/父源配子、突变、出生性别、扩散。明确**固定配置和 seed 的逐记录可复现性**；不同代码版本之间不强求字节级相同。

## 7. 代码实现验收清单（必须有单元与集成测试）

1. `[G1]` 每个个体有 `[2,L]`、所属位点与等位 ID 一致，创始等位频率有效。
2. `[G2–G5]` 杂合位点分离 `1:1`、两个不同连锁群 `1:1:1:1`；纯合永不分离出不存在的等位。
3. `[G3–G5]` 同组 10 cM 的重组分数约 `0.0906346`，0 cM 的**纯函数极限测试**不交换；较远位点趋近 0.5。
4. `[G3–G5]` 交换次数 `K=0,1,2,3` 分别无/有/无/有表型可观测位点切换；不能把“发生2次”错记为“2个重组型”。
5. `[G6–G9]` 母父来源可验证、每子代两个配子、突变产生新 ID、无直接表型突变。
6. `[D1–D8]` 植物严格异交、昆虫仅雌虫母方产仔、同种空间配偶可获得性、无配偶=零出生；父本可贡献多个子代。
7. `[E5–E8]` 植物资源守恒，具有同一个寄主的多昆虫按比例共享资源；取食顺序不改变最终分配。
8. `[S3、D10]` 母源扩散起点、周期边界、多植物同格、空种群都无越界错误。
9. 完整随机重放、原子失败、终点不额外出生、输出数据字典和公式字段一致。
10. 实际运行短程 smoke 及结果审计，报告真实退出码与运行日志；不根据种群数量方向调整科学参数。

## 8. 公式—规则—结果三向索引

| 机制 | model.md | architecture.md | results.md |
|---|---|---|---|
| 遗传 / 连锁重组 / 突变 / 表型 | G1–G9 | Step 0、6 | `gamete_transmissions`, `mutations`, `genetic_summary` |
| 寄主气味 / 选择 | E1–E4 | Step 1 | `host_choices`, `events`, `host_network` |
| 食物分配 / 成本 / 生存 | E5–E11 | Step 2–3 | `feeding`, `plants`, `parents`, `populations` |
| 配偶 / 繁殖 / 期望方差 | D1–D10 | Step 4–7 | `parents`, `offspring`, `mating_summary`, `populations` |
| 空间 / 扩散 / 代替换 | S1–S3、D10 | Step 1、4、7 | `spatial`, `snapshots` |

## 9. 理论依据的使用纪律

-  **遗传学/方法依据**：孟德尔分离与有性遗传、Haldane 连锁映射、SLiM nonWF 正向仿真 [R02–R04]。
-  **生态/行为机制线索**：植物化学性状的遗传基础、气味对寄主选择的作用、昆虫解毒的底物差异 [R08–R11]。
-  **项目特定数学函数**：VOC 混合、气味距离核、最低维度限制解毒、资源比例分配、繁殖成本及交配均匀抽样。文献不能自动验证其函数形状或系数；具体见 `model.md` 和 `assumptions.md`。
-  **模型报告规范**：遵循 ODD 的目的、实体/尺度、顺序、设计概念、初始化、输入、子模型和可重现性描述 [R01]。

---

## 参考文献与理论出处（按文内编号）

**[R01]** Grimm V, Railsback SF, Vincenot CE, et al. 2020. The ODD protocol for describing agent-based and other simulation models: A second update to improve clarity, replication, and structural realism. *Journal of Artificial Societies and Social Simulation*, 23(2):7. https://doi.org/10.18564/jasss.4259  
**适用范围**：IBM 的 ODD 模型描述、方法透明度和可复现性；不是具体生态方程的来源。

**[R02]** Haller BC, Messer PW. 2023. SLiM 4: Multispecies eco-evolutionary modeling. *The American Naturalist*, 201(5):E127–E139. https://doi.org/10.1086/723601 （[开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10793872/)）  
**适用范围**：多物种显式遗传、QTL、生态互作、nonWF 风格的前向个体模型。

**[R03]** Haller BC, Messer PW. 2019. SLiM 3: Forward genetic simulations beyond the Wright–Fisher model. *Molecular Biology and Evolution*, 36(3):632–637. https://doi.org/10.1093/molbev/msy228  
**适用范围**：个体出生/死亡/适合度的非 Wright–Fisher 正向模拟。evoSim 不是 SLiM 软件的封装。

**[R04]** Haldane JBS. 1919. The combination of linkage values, and the calculation of distances between loci of linked factors. *Journal of Genetics*, 8:299–309. [Haldane 著作目录](https://jbshaldane.org/bib/JBS-Haldane-Bibliography.html)  
**适用范围**：无交换干涉的重组映射 `r=(1−e^(−2d))/2`；适用条件必须说明。

**[R05]** Chu TT, Kristensen PS, Jensen J. 2024. Simulation of functional additive and non-additive genetic effects using statistical estimates from quantitative genetic models. *Heredity*, 133:33–42. https://doi.org/10.1038/s41437-024-00690-5 。同年勘误：https://doi.org/10.1038/s41437-024-00725-x 。  
**适用范围**：个体功能遗传效应与群体统计遗传效应的区分；不等于分子酶功能模拟。

**[R06]** Barrett SCH, Harder LD. 2017. The ecology of mating and its evolutionary consequences in seed plants. *Annual Review of Ecology, Evolution, and Systematics*, 48:135–157. https://doi.org/10.1146/annurev-ecolsys-110316-023021  
**适用范围**：种子植物交配系统多样性和花粉受限的重要性；不证明严格异交适用于所有植物。

**[R07]** McFadden D. 1974. Conditional logit analysis of qualitative choice behavior. In: Zarembka P, ed. *Frontiers in Econometrics*. Academic Press, pp. 105–142. [作者提供的章节](https://eml.berkeley.edu/reprints/mcfadden/zarembka.pdf)  
**适用范围**：多项 Logit / Softmax 离散选择的数学基础；不证明昆虫行为遵从该具体效用核。

**[R08]** Kliebenstein DJ, Kroymann J, Brown P, et al. 2001. Genetic control of natural variation in Arabidopsis glucosinolate accumulation. *Plant Physiology*, 126:811–825. https://doi.org/10.1104/pp.126.2.811  
**适用范围**：植物防御化学谱具有遗传基础；不提供本项目两个 D 轴的校准值。

**[R09]** Köllner TG, Schnee C, Gershenzon J, Degenhardt J. 2004. The variability of sesquiterpenes emitted from two *Zea mays* cultivars is controlled by allelic variation of two terpene synthase genes encoding stereoselective multiple product enzymes. *The Plant Cell*, 16:1115–1131. https://doi.org/10.1105/tpc.019877  
**适用范围**：等位变异改变多组分 VOC 产物谱；不验证本研究的线性混合系数。

**[R10]** Auer TO, Khallaf MA, Silbering AF, et al. 2020. Olfactory receptor and circuit evolution promote host specialization. *Nature*, 579:402–408. https://doi.org/10.1038/s41586-020-2073-7  
**适用范围**：昆虫嗅觉受体及神经环路参与寄主相关行为；不意味着 A 就是某个 OR 的酶学参数。

**[R11]** Li X, Baudry J, Berenbaum MR, Schuler MA. 2004. Structural and functional divergence of insect CYP6B proteins: From specialist to generalist cytochrome P450. *PNAS*, 101:2939–2944. https://doi.org/10.1073/pnas.0308691101  
**适用范围**：昆虫不同解毒酶具有不同的底物广度与效率；不提供本模型的 `q=min(I/D)` 公式。

**[R12]** Chesson P. 2000. Mechanisms of maintenance of species diversity. *Annual Review of Ecology and Systematics*, 31:343–366. https://doi.org/10.1146/annurev.ecolsys.31.1.343  
**适用范围**：稳定化共存机制和种间/种内竞争区别；不证明本项目竞争矩阵。

**[R13]** Beverton RJH, Holt SJ. 1957. *On the Dynamics of Exploited Fish Populations*. Fisheries Investigations, Series II, Vol. 19. [Beverton–Holt 函数与原始文献说明](https://search.r-project.org/CRAN/refmans/LeMaRns/html/rec_BH.html)  
**适用范围**：负密度依赖的有理函数建模传统；本项目局部个体竞争是从该理论形状借鉴的变体。

**证据声明**：引用文献只对相应的“遗传规律、方法形式或经验机制”负责。`μ=0.01`、`σ_mut=0.02`、`λ=2`、`α=0.02/0.005`、化学效用系数、配偶半径、严格异交与 Poisson 出生的特定参数等，凡未直接从独立数据拟合，都必须标记 `theoretical_placeholder` 或 `implementation_choice`。这些假设不能因为列出参考文献而被写成实测事实。
