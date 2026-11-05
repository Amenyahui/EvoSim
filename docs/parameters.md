# 参数、输入与文献依据

本文集中维护配置入口、参数含义、单位、取值和证据。公式见[模型](model.md)，当前任务见[项目入口](../README.md#current-work)，全部字段见[227条参数清单](parameter_inventory.csv)。不同自然系统分开列；机制证据不等于当前数字已获经验标定。

## 配置入口与输入

以下路径相对于 `evo_sim/models/plant_insect/`。父 YAML 先载入，子层指定的 CSV 按表替换，YAML 显式参数再覆盖；路径相对于声明路径的 YAML，映射递归合并，列表整体替换。

| 入口 | 作用和默认初始条件 |
| --- | --- |
| [configs/current.yaml](../evo_sim/models/plant_insect/configs/current.yaml) | 多年生工作入口：继承单祖先设置，P1/H1各40，K=10、dmax=1；K=10/30/100为共享理论来源的短程容量对照，不是实测值。关闭新突变。 |
| [configs/current_base.yaml](../evo_sim/models/plant_insect/configs/current_base.yaml) | 模型0.1.9公共基础：P1、P2各40株、H1为40只，10×10格、每格1个共享建立位点；植物跨代留存。K=1，植物`lambda=2`/繁殖脉冲、昆虫潜在卵`lambda=4`/昆虫代。两植物池创始频率相同，处理表达参考见G1。 |
| [configs/current_single_ancestor.yaml](../evo_sim/models/plant_insect/configs/current_single_ancestor.yaml) | 继承公共基础并替换物种与创始频率表：P1/H1各40，K=1、3轮。它不改变建立位点规则；短于5轮的运行不能验收植物繁殖脉冲。 |
| [configs/current_200gen.yaml](../evo_sim/models/plant_insect/configs/current_200gen.yaml) | 继承当前多虫工作入口并延长为200轮、4重复；新参考下的诊断设计，不是先前两池400昆虫历史实验的等价重跑。 |
| [configs/demographic_sensitivity.yaml](../evo_sim/models/plant_insect/configs/demographic_sensitivity.yaml) | 单祖先K=1基线与最大叶量2/4、卵均值8、搜索/交配3及扩散2的理论敏感性对照；不改变多年生人口学规则。 |
| [configs/initial_conditions.yaml](../evo_sim/models/plant_insect/configs/initial_conditions.yaml) | 零转移的有限位点初态诊断：总P20/P60/P100（两池各10/30/50），与H60/H6/H0组合形成7种预定情景。旧P180/P600超出100位点，已退役；该更正不调长期工作入口。 |

公共配置从[inputs/current/species.csv](../evo_sim/models/plant_insect/inputs/current/species.csv)读取兼容池、初始数量、最大资源K、需求与繁殖率。旧competition表、alpha和母源密度φ退役；空位由所有兼容池共享，不同时叠加旧竞争矩阵。`species_id`的P/H编号约束双亲兼容性，不是观察谱系或参考化学型。

`inputs/current/` 包含当前位点、等位、创始频率及参考表；[reference_traits.csv](../evo_sim/models/plant_insect/inputs/current/reference_traits.csv) 的 `species_id` 是参考面板ID，可与兼容池数量不同，只供反事实诊断，不参与真实交配或取食。旧协议配置、饱和输入及源码保存原文、哈希和路径映射后归档到 `codex/history/`，不由当前入口混用，也不回写历史结果身份。

## 输入格式与运行模式

- 遗传 CSV 使用 UTF-8。`loci.csv` 的 `target_j` 指目标轴，`mutation_mean_j`、`mutation_sd_j` 为效应增量分布；`alleles.csv` 和 `founder_frequencies.csv` 分别保存效应与创始频率。公共向量列为0–7：植物仅填0–2，3–7须为空；昆虫八列均填，不能把空列补为额外植物性状。
- `ploidy=2`、`inheritance=biparental_meiosis`、`recombination_model=slim_haldane_marker_map`；默认引擎命令 `slim`，`EVOSIM_SLIM` 可覆盖。`mutation_enabled.plant/insect` 分别控制新突变，工作配置关闭两者；公共配置保留开启及关闭对照，站立变异与重组仍在。
- `olfaction_mode=template` 在当前协议表示模块效用，`neutral` 令候选效用相同；`choice_mode=outside` 允许拒绝，`forced` 在候选非空时必选，无候选均不产卵。`tau≥0`，候选效用乘tau，拒绝项不乘；`shuffle_voc` 仅置换逐卵线索。
- `shuffle_voc` 使用独立阶段8随机流，避免置换次数本身推进选择顺序流；处理产生行为差异后仍不承诺整条轨迹完全配对。默认false，旧批次沿用其保存的旧随机协议。
- 核宽、闲置项、半饱和常数、dmax及最大资源K须为正，当前叶量存量B允许零；危害率、恢复率和成本非负，eta在[0,1]，新株资源比例在[0,1]。繁殖间隔和每格位点数为正整数。`sigma/weights`仍为不参与模块效用的接口字段，配置中的退役空键`B`不是植物叶量存量。模型不接受花期或交配前定居开关。
- `species.csv` 的数量、资源和lambda可由子配置覆盖，`insect_pop=null` 表示初始昆虫数求和。`initial_plant_cells`、`initial_insect_cells` 默认为null、即随机位置；半径为非负整数格数。资源与需求须同单位、同窗口，可同为归一化干生物量或受支持的克/毫克可食干生物量，不能混用。
- 公共入口为5代、`record_every=1`、种子20261021、1重复；工作入口的时长、种子及重复由各自YAML覆盖。`genetics.slim.timeout_seconds=120`、`max_individuals_safety=5000`、`max_alleles_safety=100000`、`max_genotype_values_safety=4000000` 是计算保护，不是生态承载量；停止状态单独报告。
- [analysis/lineages.yaml](../evo_sim/models/plant_insect/analysis/lineages.yaml) 控制只读观察；旧 v014 分析设计不能用于新生活史。观察阈值的含义和限制见下文O1，不反馈模拟行为。

## 分组参数与证据

本文件对应[逐字段清单](parameter_inventory.csv)的evidence ID。清单覆盖六个当前展开配置的227条叶字段、输入表列和观察设置，列表保留全部元素；不是227个独立生物学参数。`current_working`含容量情景覆盖，`inputs.*`为原CSV，YAML覆盖后的有效量以`parameters.*`为准。各条保留证据组及当前值支持层级。早期定向核查为2026-11-02，新增多年生机制核查为2026-11-04；不是系统综述。不同物种、地点、处理和阶段独立解释，不拼成一个经验系统。**机制有实证支持不等于默认数字已校准。**

清单的`current_public/current_single/current_sensitivity/current_working/current_200gen/current_initial`逐列对应公共基础、单祖先、人口敏感性、工作入口、200代延长与零代初态诊断；`source_file`标明各展开配置和实际输入表。CSV列保留原表值，有效参数列记录YAML覆盖后的值。观察设置来自`analysis/lineages.yaml`，不反馈模拟。2026-11-04的处理表达修订延续到多年生协议；旧批次保留原输入，不能据当前清单反推旧实验初态。

<a id="r1"></a>

### R1 叶量存量、容量与需求

| 参数 | 值、单位与实施决定 |
|---|---|
| `resource_capacity` K | 每株最大归一化可食叶量，严格大于零；多虫参考K=10、30/100为容量对照，公共基础K=1。是上界，不再每轮重赋资源。 |
| 叶量存量B | 单位与K相同，范围[0,K]；创始株B=K，旧株食后扣减并恢复，新株初始B=.25K。分配核的R读取当轮B。 |
| `max_feeding_demand` dmax | 每幼虫完整取食窗口1个相同单位；保持归一化。实际分配受利用和同株竞争影响，不能把需求当成必然摄入。 |
| `resource_unit`、`demand_unit`、窗口 | 相同的 `normalized_edible_dry_biomass` 与 `one_generation_larval_feeding`；可另以同单位真实干重输入，但默认不是克。见N1。 |

Sarfraz等[1]全文方法与图3：温室4周龄油菜每株接10只初龄小菜蛾，5株/处理，5种施肥组约35–76%化蛹；**10为初始接虫量而非资源容量**，图估读而非容量测定。Wang与Lu[2]全文方法、图1：25±1°C第四龄自由取食白菜/甘蓝分别4.995±0.245、4.426±0.189 cm²/虫（均值±SE）；不是全龄干质量。两项支持同株多虫与取食有限的可能，不能由此推得R/dmax=10或最大值。

### R2 幼虫招募与资源限制

取食需求为 dmax×q，同株按资源约束分配，实际摄入 x；幼虫存活概率 `g=x/dmax`，`sI0=0` 和 `unmatched_insect_fitness=0`，零摄入无招募，满摄入本环节存活1。Wang与Lu[2]的第四龄限食影响化蛹与羽化，有限范围近线性不能验证全龄S=g及两端点。保留该闭合作对照；现有数值不是测得的背景死亡率，未来同实验分阶段拟合充分供食背景生存及摄入响应，避免重复扣死亡。

同dmax时每株招募**期望**为`min(Σq,B_start/dmax)`，不是实际数的硬上限；异质dmax时依需求组成变化。昆虫不另设干扰密度因子，不能把资源闭合下的密度响应直接解释为寄主局部适应。植物本轮死亡不附加幼虫死亡。B_start=0时g和当轮消费比例h均为零，但叶量亏损为1，植物仍承担背景和亏损危害；旧零容量下S_P=1的约定不属于本协议。

### P1 多年生存活、恢复和独立繁殖时间

以下均为`researcher_selected_theoretical`参考，不是自然系统数值拟合。时间单位固定为一个昆虫世代Δ=1，没有换算成天或年。

| 参数 | 当前值、量纲及明确含义 |
| --- | --- |
| `background_mortality_rate` μ₀ | `−ln(.99)=.01005033585350145`/昆虫代；满叶无虫时单轮存活.99，连续5轮背景存活`.99⁵≈.95099`。 |
| `defoliation_mortality_rate` μ_D | `ln(2)/5=.13862943611198905`/昆虫代；若完整5轮始终全缺叶，则相对背景存活再乘.5，不是单次剪叶死亡率。 |
| `resource_recovery_rate` r | `ln(2)=.6931471805599453`/昆虫代；幸存株叶量亏损每轮减半，不保证受食株能恢复满叶。 |
| `reproduction_interval` k_P | 5昆虫轮/植物繁殖脉冲；仅`(t+1)%5=0`发生提案，新株下一脉冲才可繁殖；同步钟不是实测花期。 |
| `newborn_resource_fraction` f_new | .25，无量纲；新建立株B=.25K，未显式模拟真实种子萌发或株高。 |
| `plant_offspring_lambda` λ_P | 2有效建立提案/合格母株/植物繁殖脉冲；均值λ_P F_m，不乘旧φ。数值虽沿用2，时间分母和终点已变。 |

死亡使用`defoliation=1−Bpost/K`，而`h=X/Bstart`只描述当轮现有资源的消费比例；即使当轮没有新消费，先前未恢复的亏损仍影响死亡。恢复为`Bnext=K−(K−Bpost)exp(−rΔ)`，允许零可食叶量的幸存株由未显式追踪的根茎/储备恢复；没有储碳状态，不能宣称已模拟碳饥饿阈值。繁殖不另扣叶量，也不按叶量再乘母源损伤因子；这是非闭合能量近似。

`defoliation`在此表示潜在冠层叶量亏损，不全是虫致损伤。新株以.25K建立，即使无虫也先有δ=.75，承担额外建立压力；无恢复对照持续承担亏损。模型没有成熟体型、株高或幼株储备状态，这项风险是明确的简化，不是实测幼株死亡率。

**机制证据分系统解释。** Stam等[27]在野生多年生甘蓝上跨两生长季跟踪节肢动物、植物表现与第二年结籽，支持同株跨季和繁殖时间区别；早季引虫处理未显著改变后续群落，不能概括成所有前期取食都有长期效应。Gols等[28]比较两季甘蓝化学动态，并以菜粉蝶进行16天取食实验；它支持时间和诱导反应的重要性，不标定固定D/V、恢复r或5轮繁殖钟。

Heichel与Turner[29]对红栎/红枫人工去叶后观察再展叶，支持存活个体可重建叶量，不验证指数恢复或每轮恢复一半。Barker Plotkin等[30]在舞毒蛾—栎树系统发现去叶、非结构性碳储备和后续死亡相关；叶量和储碳不是同一状态，树木碳阈值不得移为甘蓝参数或本模型μ_D。各系统未被拼成一套实测多年生寄主参数。早期拟南芥剪叶及白菜花数文献[3,4]也不支持退役sP0的数值。

### P2 共享建立位点与空位抽签

`microsites_per_cell=m=1`个建立位点/格，10×10格最多100株；这是生态空间假设，不是计算保护。所有池共享，幸存旧株占原位，死亡释放空位。每格从到达提案中均匀无放回抽最多空位数的新株；零提案仍为空，不强填、不移民、不补种，不给各池保留席位。α、competition输入和φ全部退役，避免在母源均值和建立阶段重复限制密度。

Metcalf等[31]对两个**单次结实后死亡**的多年生蓟类群体研究，将密度依赖放在微位点建立；这支持建立阶段可以重要，不证明本项目重复繁殖生活史或m=1。Peng等[32]是空间群落**理论模拟**，按到达空位的种子及竞争能力抽签；本项目用均匀提案竞争且允许保留空位，不据其推自然寄主的位点密度。旧年度繁殖和竞争文献[5–7]保留原系统背景，不能再标为现行α的数值依据。

### D1 昆虫繁殖、性别与交配

`insect_offspring_lambda=4` 为每有局部合格雄虫的雌虫、每模型轮潜在卵的Poisson基准；不等于成年招募。创始及子代雌概率各0.5；潜在卵经寄主接受与幼虫生存，父本对每枚成功卵从同池局部雄虫均匀有放回抽取，无精子耗竭。Yang等[8]全文表1：23±1°C、萝卜条件100卵中68成虫（30雌38雄），30雌产卵173.20±11.40卵/雌；初代豌豆17成虫、9雌，57.33±19.55卵/雌。羽化雌比0.441是条件于生存的比例，终生卵数不是潜在卵或招募的直接参数。Song等[9]全文表1：26±2°C，1日龄配对组交配100%、产卵雌占已交配雌85%，而6日龄对应50%与47.37%；多次相遇的实验不等于有候选即逐卵抽父本。原值保留为理论设定，待相同生命周期、机会和分母下拟合，检验泊松方差并避免重复乘死亡或成本。

### D2 生活史与繁殖结构

植物严格异交、同池且非自身，旧株跨昆虫代存活并可重复繁殖；昆虫雌雄异体、非重叠代。没有移民或种子库，未建立提案当轮消失。既有兼容池不是新物种。Postma与Ågren[10]的拟南芥种子库和Yamamori[11]的油菜异交研究是不同系统，只说明被省略的过程可能重要；它们不验证本模型无种子库、硬花粉半径或严格异交。当前先验证多年生人口学，不新增寄主配对限制、杂交死亡或双寄主实验。

<a id="g1"></a>

### G1 创始遗传与表型映射

植物3轴各2位点、昆虫8轴各2位点，二倍体每位点两种创始等位，频率各0.5；各轴背景0，表型为 `clip(background + .5×拷贝效应和, 0, 1)`。创始同位点两个等位效应相差0.05，即相对中点±0.025；四份独立拷贝经scale=.5后每轴站立标准差=.025、范围为中心±.05。突变SD=.02是另一参数。坐标范围、scale及效应存在重参数化，0和1不是实测浓度。Fisher[23]原始理论论文支持孟德尔因素与加性方差分析的理论框架，不规定此处的位点数或频率。Henniges-Janssen等[12]全文分析：小菜蛾豌豆适应的家系推断2–5个位点，支持复杂遗传的可能，不证明每轴恰两个位点或频率0.5。

旧输入的 `E^I` 两中心均为`.8`，`D_N∈[.30,.40]`、`D_S∈[.50,.60]`、`μ^I_1∈[.05,.15]`、`μ^I_2∈[.85,.95]`，使全部15,625个相关表型组合`q=1`，未截断的最小`I/D≈1.11352046`。两防御轴在这个无新突变的可达范围内无法经利用率形成选择差异，升高R也不解除截断。该结论不意味着VOC选择和表达成本同时消失，也不对有新突变的未来状态作恒饱和保证。

当前参考固定原植物分布、`μ^I`中心`(.1,.9)`、核宽`.25`、闲置`.1`及站立方差，仅反解处理表达中心。令名义中心的处理分配矩阵

\[
P^I=\begin{pmatrix}
.9009088028932118&.001496911573913638\\
.001496911573913638&.9009088028932120
\end{pmatrix},\quad
(P^I)^\mathsf T\bar E^I=(.35,.55)^\mathsf T.
\]

解为 `E^I`中心`(.38748329034947154,.6098508195430548)`，约`.3875/.6099`；每个相关位点两等位效应取`E_r/2±.025`，详见[当前等位表](../evo_sim/models/plant_insect/inputs/current/alleles.csv)。计算精度记录设计方程的解，不代表测量精度。这使名义中心`I=D`；由于核非线性，不声称随机创始个体的平均I恰等于D中心。它是跨越截断阈值的理论机制参考，不是自然祖先估计、经验校准或按不灭绝目标拟合；感觉、植物、调谐、映射scale和突变率均不因此更改。

这是整体理论基线的修订，**并非只改变q的演化实验**：按未更改的`insect_cost=.125`，处理表达中心总量从1.6降为约`.997334`，母源产出乘子在其他轴固定时提高约`.0753332`。因此新旧人口轨迹差异不能全部归于防御反馈；创始q对照只诊断固定核和防御条件下的利用率响应，未分离繁殖成本效应。本轮不为抵消该代价变化而追加生态参数补偿。

枚举六个相关轴各五级的`5⁶=15,625`个表型组合（不是全部QTL基因型），得到两通道`I/D`范围约`.75197–1.32962`和`.83096–1.19960`，`q∈[.75197,1]`。按`q<1−1e−12`计，等权组合覆盖`11828/15625=75.6992%`；按创始独立四拷贝的`Binomial(4,.5)`权重计算，随机昆虫—植物配对的未饱和概率约`75.1226%`，平均q约`.948645`。两通道分别增大D `1e−4`可降低q的加权覆盖约`41.8398%`和`35.4719%`；阈值附近两方向探针的覆盖可以重叠。这些是功能空间及创始频率加权的覆盖，均不是已选择寄主上实际幼虫比例。启动诊断另按该次实际抽到的创始个体报告，后续资源拥挤、选择和遗传仍需实际账本检验。

### G2 新突变

每出生个体每位点每拷贝突变概率 `.01`，目标效应增量均值0、SD `.02`（未乘映射scale）；植物12份、昆虫32份拷贝的至少一次突变概率约11.36%和27.50%。Ossowski等[13]原文结果：拟南芥5条单籽传代30代，替换率约7×10⁻⁹ **每碱基位点每代**；不知功能靶大小及效应分布不能转换成抽象QTL事件概率。保留原值和关闭新突变对照，后续报告每代有效事件供给、效应方差及截断；关突变不关闭站立变异和重组。

### G3 连锁图与重组

同轴位点0、10 cM，不同轴的连锁群独立，双亲均依SLiM和Haldane图重组。Haldane无干扰映射是项目采用的经典理论假设，Haldane 1919 本次仅核对原始文献书目、未核全文；在10 cM给重组率 `(1−e⁻⁰·²)/2≈0.0906`；50 cM也不是自由重组，不同独立连锁群才为0.5。Jiggins等[14]全文方法、结果：Heliconius一F2家系73只用于连锁分析，描述鳞翅目雌性无交叉；该类群性别差异与本模型双亲相同图不同。10 cM非实测本模型性状图距；保留结构，研究具体鳞翅目时需另立性别特异遗传机制。

### C1 化学映射

`n=2` 防御和挥发目录；植物三轴映射 D=(x_N,x_S)、V=(eta_shared×x_S,x_V)，`eta_shared=1` 是归一化比例，不是分子浓度、排放率或实测耦合系数。`chemical_features=[0,1]` 为抽象类型坐标，两通道对应关系不保证同一化合物。保持此项目结构，解释共享轴只作为理论耦合；若映射自然系统须指定浓度、化学空间与时间窗口并独立验证。

Frey等[25]玉米研究连接诱导吲哚释放与防御相关代谢的酶学分支，提供化学功能联系的实例；相关同源酶及分支不等于本模型同一坐标同时线性控制两种输出，不能据此标定n=2或eta=1。

### C2 处理和感觉模块

`processing_width=.25`、`sensory_width=.25` 是化学坐标高斯核宽，`processing_idle=.1` 是分配分母的闲置项，`sensory_half=.2` 为项目VOC加权活化的半饱和常数，`sensory_valence=[1,−1]` 为固定读出符号。宽度对应半高距离约0.294坐标单位，而非浓度EC50。Liu等[15]摘要与论文讨论中小菜蛾Or35/Or49对49种气味筛选、3种主要异硫氰酸酯相关的受体/产卵反应；Liu等[16]摘要示另一受体对三种芳香物及驱避有关。它们支持特定化学反应和行为方向，不支持本模型的具体核宽、闲置、饱和值及价性。保持当前值为理论选择，需联合剂量、化学结构与行为实验辨识。

价性输入只要求两个有限实数；同正、同负或零可定义合法理论对照，不能因默认`[1,−1]`而强制异号。默认`b=1`且感觉表达不超过1时`u≥0`，故默认`tau=1,outside=0`在n候选下拒绝概率至多`1/(n+1)`；相对回避仍存在。判断候选是否低于拒绝项须比较`tau*u`与`outside_utility`，详见[模型§3](model.md#3-产卵线索与选择)。

### C3 表达与繁殖成本

植物 `plant_cost=.125`、`voc_cost=.05`，昆虫 `insect_cost=.125`、`sensory_cost=.05`；线性乘原发表达量并以零截断母源产量，不对父源资格收费。Baldwin[17]野外745对烟草诱导处理，未受攻击诱导株种子更少，但诱导多条通路，不能分离到某一性状的每单位线性代价。Zhang等[18]摘要：无杀虫剂时抗性小菜蛾雌株平均少产卵约20%，且成本随温度变化；ace1靶标抗性不等同解毒或感觉表达。保持系数为待校准假设，不能将两研究结果直接代入。

### B1 寄主候选、产卵效用

`candidate_size=8`，每潜在卵从半径内无放回抽至多8株；`olfaction_mode=template` 在本协议指模块效用；`choice_mode=outside`，`b=1`、`tau=1`、`outside_utility=0` 决定候选与拒绝的softmax相对效用。零候选拒绝；候选间及接受概率取决实际数量，b、outside、tau不能仅从已产卵份额分别识别。Little等[19]全文方法：温室2–4小时，*P. brassicae* 平均每株1–2块、每块10–40卵；*P. rapae* 每株约20枚单卵，涉及多虫而非每雌代卵数，不证明候选8或独立逐卵。条件logit数学形式可追溯至McFadden[26]，主审核对了作者机构扫描件题名与书目，未逐页复核推导；这属于概率选择理论，不是昆虫行为实测。当前数值是行为机制假设，需逐母轨迹、候选及拒绝数据拟合。`shuffle_voc=false` 是保留真实线索分配的干预设计开关，不是需测量的自然概率。

### S1 空间、位置与初始数量

10×10周期网格、Chebyshev距离，无实测米制；搜索2、昆虫交配2、昆虫扩散1、植物花粉1、植物提案扩散1格。成年植物原位留存，创始位置在有限建立位点中无放回抽取。公共P1/P2各40株、H1 40只，单祖先P1/H1各40；新公共80株适配100个共享位点，不兼容保留旧400株再自动裁剪。初始密度是研究设计，不是实测。Mo等[20]摘要：小菜蛾作物田标记重捕的平均位移估计14–35 m；Bolton等[21]摘要：2.83 ha田地雄虫重捕率7.93%、距离中位14 m，田界95 m；Ishii等[22]摘要：下风向信息素有效空间估计0.6–0.9 m。各测不同过程且有陷阱与边界限制，不是逐代格位移或硬交配半径。保留格值为模型尺度；如拟合自然系统，先定义格长、时间和边界并采同阶段轨迹。

大半径下周期取模可使目的格的位移重数不同，每轴最多相差1；这符合整数位移均匀的定义，不是三角分布错误。不以此无依据缩小合法半径，也不将目的格重采样为均匀。

### O1 观察阈值与参考面板

只读操作谱系的无相位QTL拷贝不匹配阈值0.25、最小群2、持续3代；参考植物面板 `chem_A=(.25,.65,.25)`、`chem_B=(.65,.25,.65)`，只作功能响应比较，不参与配对或实际取食。阈值、单链接桥接、样本量及窗口影响结果，不构成物种证据；保留原设并做预先声明的阈值稳定性比较，未知保留，不能反馈兼容池。

### N1 单位与归一化

`resource_unit`和`demand_unit`同单位，取食和需求窗口同为一代；性状、化学类型位置、表达和效用均无量纲，格与代没有自然长度或天数。数值1、0及坐标端点是尺度约定，不是生物观测。换到真实单位需同源测量并校准映射，不能混用鲜重、干重、叶面积或单株与单虫分母。

### T1 运行与数据协议

版本/输入路径、来源元数据、种子、重复、运行代数、记录间隔、最大个体/等位基因/基因型值保护、SLiM路径与超时、状态、子步骤与输出ID为计算和实验设计，不是生物实测参数。计算上限触发中止而不是灭绝；不改历史源码哈希。清单对应项按实际配置记录，只有明确完成的实验才能解释结果；验证输出使用 `codex/validation/`，正式实验须明确指定 `results/`。

### L1 退役字段和旧时间尺度

`sigma=.2`、`weights=[.5,.5]` 为旧接口输入校验，在当前模块效用不生效；`B`、`rho_voc`、`voc_bounds` 为退役空键。不得当成活跃参数标定，也不得把历史协议输出冒充当前模型。旧sP0、alpha、competition输入、competition_radius和母源φ不进入0.1.9；植物λ仍记为2但已改为每脉冲的建立提案，不能按旧每轮新生语义解释。

## 原始研究及核验层级

1. Sarfraz, Dosdall & Keddie (2009), *Bottom-up effects of host plant nutritional quality on Plutella xylostella (Lepidoptera: Plutellidae) and top-down effects of herbivore attack on plant compensatory ability*. 方法、图3及表2：[DOI](https://doi.org/10.14411/eje.2009.073)。
2. Wang & Lu (2025), *Starvation biology of Plutella xylostella from a post-harvest crop sanitation perspective*. 全文方法及图1：[DOI](https://doi.org/10.1127/entomologia/2025/2883)。
3. Akiyama & Ågren (2012), *Magnitude and Timing of Leaf Damage Affect Seed Production in a Natural Population of Arabidopsis thaliana (Brassicaceae)*. 全文图2、表2：[DOI](https://doi.org/10.1371/journal.pone.0030015)。
4. Brosset 等 (2023), *Airborne cues accelerate flowering and promote photosynthesis in Brassica rapa*. 全文方法与讨论：[DOI](https://doi.org/10.1111/1365-2745.14046)。
5. Stouffer (2022), *A critical examination of models of annual-plant population dynamics and density-dependent fecundity*. 理论推导全文§2.1：[DOI](https://doi.org/10.1111/2041-210X.13965)。
6. Hauser 等 (2003), *Frequency-dependent fitness of hybrids between oilseed rape (Brassica napus) and weedy B. rapa (Brassicaceae)*. 全文表6：[DOI](https://doi.org/10.3732/ajb.90.4.571)。
7. Muñoz-Parra 等 (2017), *Plant–plant interactions influence developmental phase transitions, grain productivity and root system architecture in Arabidopsis via auxin and PFT1/MED25 signalling*. 全文方法及图3：[DOI](https://doi.org/10.1111/pce.12993)。
8. Yang 等 (2021), *Fitness comparison of Plutella xylostella on original and marginal hosts using age-stage, two-sex life tables*. 全文方法与表1：[DOI](https://doi.org/10.1002/ece3.7804)。
9. Song 等 (2014), *Analysis of the Mating and Reproductive Traits of Plutella xylostella (Lepidoptera: Plutellidae)*. 全文方法与表1：[DOI](https://doi.org/10.1093/jisesa/ieu129)。
10. Postma & Ågren (2022), *Effects of primary seed dormancy on lifetime fitness of Arabidopsis thaliana in the field*. 全文方法：[DOI](https://doi.org/10.1093/aob/mcac010)。
11. Yamamori (2011), *Outcrossability of Brassica napus L. and B. rapa L. in an Experimental Field*. 机构摘要及原文：[DOI](https://doi.org/10.6090/jarq.45.173)。
12. Henniges-Janssen 等 (2011), *Complex inheritance of larval adaptation in Plutella xylostella to a novel host plant*. 全文连锁分析、讨论：[DOI](https://doi.org/10.1038/hdy.2011.27)。
13. Ossowski 等 (2010), *The Rate and Molecular Spectrum of Spontaneous Mutations in Arabidopsis thaliana*. 原文摘要及结果：[DOI](https://doi.org/10.1126/science.1180677)。
14. Jiggins 等 (2005), *A Genetic Linkage Map of the Mimetic Butterfly Heliconius melpomene*. 全文方法与结果：[DOI](https://doi.org/10.1534/genetics.104.034686)。
15. Liu 等 (2020), *The Molecular Basis of Host Selection in a Crucifer-Specialized Moth*. 摘要及出版方讨论：[DOI](https://doi.org/10.1016/j.cub.2020.08.047)。
16. Liu 等 (2022), *Odorant Receptor PxylOR11 Mediates Repellency of Plutella xylostella to Aromatic Volatiles*. 摘要：[DOI](https://doi.org/10.3389/fphys.2022.938555)。
17. Baldwin (1998), *Jasmonate-induced responses are costly but benefit plants under attack in native populations*. 全文摘要与讨论：[DOI](https://doi.org/10.1073/pnas.95.14.8113)。
18. Zhang 等 (2015), *Temperature-sensitive fitness cost of insecticide resistance in Chinese populations of the diamondback moth Plutella xylostella*. 摘要：[DOI](https://doi.org/10.1111/mec.13133)。
19. Little 等 (2007), *Oviposition by Pierid Butterflies Triggers Defense Responses in Arabidopsis*. 全文方法：[DOI](https://doi.org/10.1104/pp.106.090837)。
20. Mo, Baker, Keller & Roush (2003), *Local dispersal of the diamondback moth (Plutella xylostella (L.)) (Lepidoptera: Plutellidae)*. 摘要：[DOI](https://doi.org/10.1603/0046-225X-32.1.71)。
21. Bolton 等 (2019), *Response to a Synthetic Pheromone Source by OX4319L, a Self-Limiting Diamondback Moth (Lepidoptera: Plutellidae) Strain, and Field Dispersal Characteristics of its Progenitor Strain*. 摘要：[DOI](https://doi.org/10.1093/jee/toz056)。
22. Ishii 等 (1981), *Active Space of Sex Pheromone of the Diamondback Moth, Plutella xylostella (L.) (Lepidoptera: Yponomeutidae)*. 英文摘要：[DOI](https://doi.org/10.1303/jjaez.25.71)。

23. Fisher (1918), *The Correlation between Relatives on the Supposition of Mendelian Inheritance*. 原始PDF理论框架核验（主审）：[DOI](https://doi.org/10.1017/S0080456800012163)。
24. Haldane (1919), *The combination of linkage values, and the calculation of distances between the loci of linked factors*. 本次仅核原始书目，未核全文：[FlyBase记录](https://flybase.org/reports/FBrf0000954)。
25. Frey等 (2000), *An herbivore elicitor activates the gene for indole emission in maize*. 摘要及原文酶学分支说明：[DOI](https://doi.org/10.1073/pnas.260499897)。
26. McFadden (1973，作者机构所列年份), *Conditional Logit Analysis of Qualitative Choice Behavior*, in *Frontiers in Econometrics*, 105–142；常见书目亦列1974。核对作者机构标题与扫描件，未逐页复核推导：[作者机构扫描本](https://eml.berkeley.edu/reprints/mcfadden/zarembka.pdf)。


27. Stam, Kos, Dicke & Poelman (2019), *Cross-seasonal legacy effects of arthropod community on plant fitness in perennial plants*. 原始研究全文方法、结果及限制：[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6774310/)，[DOI](https://doi.org/10.1111/1365-2745.13231)。
28. Gols 等 (2018), *Seasonal and herbivore-induced dynamics of foliar glucosinolates in wild cabbage (Brassica oleracea)*. 原始研究全文，两季田间与16天温室实验：[DOI](https://doi.org/10.1007/s00049-018-0258-4)。
29. Heichel & Turner (1983), *CO₂ assimilation of primary and regrowth foliage of red maple (Acer rubrum L.) and red oak (Quercus rubra L.): response to defoliation*. 原始论文摘要核验，未核全文：[PubMed](https://pubmed.ncbi.nlm.nih.gov/28310150/)，[DOI](https://doi.org/10.1007/BF00379555)。
30. Barker Plotkin 等 (2021), *Defoliated trees die below a critical threshold of stored carbon*. 原始研究全文方法与死亡/储备结果；不将木本碳阈值标定到本模型：[DOI](https://doi.org/10.1111/1365-2435.13891)。
31. Metcalf 等 (2008), *Evolution of flowering decisions in a stochastic, density-dependent environment*. 原始研究全文Population Biology及Structured Dynamics节；物种为单次结实多年生蓟，不是重复结实甘蓝：[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC2492471/)，[DOI](https://doi.org/10.1073/pnas.0800777105)。
32. Peng, Zhou & Zhang (2012), *Dispersal and recruitment limitation contribute differently to community assembly*. 原始理论研究空位竞争段；是模拟机制而非田间位点数测量：[DOI](https://doi.org/10.1093/jpe/rtr041)。


## 补充系统与来源可靠性

以下为早期定向核查中未被上表覆盖的独立记录，保留其系统、分母和限制；本次仅迁移文字，未重新检索文献。

- **Golizadeh 等 2009，小菜蛾—甘蓝**：[DOI 10.1016/j.aspen.2009.05.002](https://doi.org/10.1016/j.aspen.2009.05.002)。核验作者全文方法与表1、3：25°C条件下报告卵到成虫累计存活80.95%、生命表平均世代时间18.22±0.33天；每温度检查200卵，幼虫发育另取50只，繁殖用15对。累计存活含多个阶段，不能直接乘入已表示相同死亡的g；生命表T也不是所有个体发育时长。
- **Green 等 2015，菜粉蝶—野生甘蓝**：[DOI 10.1371/journal.pone.0136884](https://doi.org/10.1371/journal.pone.0136884)。核验全文方法与表3：每株最初放2只二龄幼虫，死亡个体持续补入；幼虫24日龄时各叶色组的每株累计受食面积约180.1±23.95至253.4±32.52 cm²。它是补虫条件下每株累计面积，不是单虫全龄干质量，不能与小菜蛾数据拼接成资源比值。
- **Saeed 等 2019，小菜蛾温度生命表**：[DOI 10.1016/j.sjbs.2018.08.026](https://doi.org/10.1016/j.sjbs.2018.08.026)。全文表与摘要存在内部疑点：20°C时摘要R₀/GRR与表2恰差两倍；表1卵到雌成虫时间与首次产卵时间、阶段之和不一致。未澄清前不用于标定。文献不是本模型数值的自动授权。
