# 参数、输入与文献依据

本文集中维护配置入口、参数含义、单位、取值和证据。公式见[模型](model.md)，当前任务见[项目入口](../README.md#current-work)，全部字段见[237条参数清单](parameter_inventory.csv)。不同自然系统分开列；机制证据不等于当前数字已获经验标定。

## 配置入口与输入

以下路径相对于 `evo_sim/models/plant_insect/`。父 YAML 先载入，子层指定的 CSV 按表替换，YAML 显式参数再覆盖；路径相对于声明路径的 YAML，映射递归合并，列表整体替换。

| 入口 | 作用和默认初始条件 |
| --- | --- |
| [configs/current.yaml](../evo_sim/models/plant_insect/configs/current.yaml) | 当前多虫理论工作入口：继承单祖先设置，P1/H1各40，R=10、dmax=1；R=10/30/100为共享理论来源的短程诊断点，不是实测容量。关闭新突变。 |
| [configs/current_base.yaml](../evo_sim/models/plant_insect/configs/current_base.yaml) | 当前唯一模型合同的公共基础。P1、P2各200株，H1为40只；R=1，植物`lambda=2`、昆虫潜在卵`lambda=4`；同池`alpha=.02`、跨池`.005`。两植物池创始频率相同，新处理表达参考见G1。含开启/关闭新突变对照，不自动构成正式实验设计。 |
| [configs/current_single_ancestor.yaml](../evo_sim/models/plant_insect/configs/current_single_ancestor.yaml) | 继承公共基础并替换物种、竞争和创始频率表：P1/H1各40，R=1、植物`alpha=.05`，3轮。单池与两池初始数量及竞争设计不同。 |
| [configs/current_200gen.yaml](../evo_sim/models/plant_insect/configs/current_200gen.yaml) | 继承当前多虫工作入口并延长为200轮、4重复；新参考下的诊断设计，不是先前两池400昆虫历史实验的等价重跑。 |
| [configs/demographic_sensitivity.yaml](../evo_sim/models/plant_insect/configs/demographic_sensitivity.yaml) | 单祖先R=1基线与资源2/4、卵均值8、搜索/交配3及扩散2的理论敏感性对照；同样使用当前创始参考。 |

公共配置从 [inputs/current/species.csv](../evo_sim/models/plant_insect/inputs/current/species.csv) 读取兼容池、初始数量、资源、需求与繁殖率，从 [inputs/current/competition.csv](../evo_sim/models/plant_insect/inputs/current/competition.csv) 读取完整植物竞争矩阵。`species_id` 的 P/H 编号约束双亲兼容性，既不是只读谱系标签，也不自动对应不同参考化学型。

`inputs/current/` 包含当前位点、等位、创始频率及参考表；[reference_traits.csv](../evo_sim/models/plant_insect/inputs/current/reference_traits.csv) 的 `species_id` 是参考面板ID，可与兼容池数量不同，只供反事实诊断，不参与真实交配或取食。旧协议配置、饱和输入及源码保存原文、哈希和路径映射后归档到 `codex/history/`，不由当前入口混用，也不回写历史结果身份。

## 输入格式与运行模式

- 遗传 CSV 使用 UTF-8。`loci.csv` 的 `target_j` 指目标轴，`mutation_mean_j`、`mutation_sd_j` 为效应增量分布；`alleles.csv` 和 `founder_frequencies.csv` 分别保存效应与创始频率。公共向量列为0–7：植物仅填0–2，3–7须为空；昆虫八列均填，不能把空列补为额外植物性状。
- `ploidy=2`、`inheritance=biparental_meiosis`、`recombination_model=slim_haldane_marker_map`；默认引擎命令 `slim`，`EVOSIM_SLIM` 可覆盖。`mutation_enabled.plant/insect` 分别控制新突变，工作配置关闭两者；公共配置保留开启及关闭对照，站立变异与重组仍在。
- `olfaction_mode=template` 在当前协议表示模块效用，`neutral` 令候选效用相同；`choice_mode=outside` 允许拒绝，`forced` 在候选非空时必选，无候选均不产卵。`tau≥0`，候选效用乘tau，拒绝项不乘；`shuffle_voc` 仅置换逐卵线索。
- `shuffle_voc` 使用独立阶段8随机流，避免置换次数本身推进选择顺序流；处理产生行为差异后仍不承诺整条轨迹完全配对。默认false，旧批次沿用其保存的旧随机协议。
- 核宽、闲置项、半饱和常数及 dmax 须为正；资源、成本非负，eta在[0,1]。`sigma` 仍要求正值，`weights` 非负且和为1，但二者不进入当前模块效用。`B/rho_voc/voc_bounds` 退役为空，不与eta混用。模型不接受花期或交配前定居开关。
- `species.csv` 的数量、资源和lambda可由子配置覆盖，`insect_pop=null` 表示初始昆虫数求和。`initial_plant_cells`、`initial_insect_cells` 默认为null、即随机位置；半径为非负整数格数。资源与需求须同单位、同窗口，可同为归一化干生物量或受支持的克/毫克可食干生物量，不能混用。
- 公共入口为5代、`record_every=1`、种子20261021、1重复；工作入口的时长、种子及重复由各自YAML覆盖。`genetics.slim.timeout_seconds=120`、`max_individuals_safety=5000`、`max_alleles_safety=100000`、`max_genotype_values_safety=4000000` 是计算保护，不是生态承载量；停止状态单独报告。
- [analysis/lineages.yaml](../evo_sim/models/plant_insect/analysis/lineages.yaml) 控制只读观察；旧 v014 分析设计不能用于新生活史。观察阈值的含义和限制见下文O1，不反馈模拟行为。

## 分组参数与证据

本文件对应 [逐字段清单](parameter_inventory.csv) 的 evidence ID，记录公开、单祖先、人口敏感性及新多虫工作配置的参数含义及决定。清单覆盖237条配置叶字段、输入表列和观察设置；列表保留全部元素，不能把237条解释为237个独立生物学参数。`current_working` 记录新入口及三个情景的实际覆盖值；`inputs.*` 记录原始CSV，YAML覆盖后的有效量以 `parameters.*` 为准。每条均有证据组、处置决定和当前数值是否获直接标定的状态。检索日期2026-11-02，主审整理2026-11-03，系定向研究而非系统综述。项目是通用理论模型，文献中的不同物种、地点、处理和生活史阶段独立陈列，不拼成同一个经验系统。**机制有实证支持不等于默认数字已实测校准。** 公式与阶段见[模型](model.md)。

清单的 `current_public/current_single/current_sensitivity/current_working` 分别对应上表公共基础、单祖先、人口敏感性、多虫入口；`source_file`逐列标明当前来源，CSV列保留原始表值，已加载的参数列记录有效值。2026-11-04的处理表达修订适用这些当前入口，旧批次使用保存的旧输入，不能据当前清单反推旧实验初态。

<a id="r1"></a>

### R1 资源与需求

| 参数 | 值、单位与实施决定 |
|---|---|
| `resource_capacity` R | 旧公共配置每株1个归一化可食干生物量；新多虫情景R=10，30/100为检查点，不是实测。每株每幼虫取食窗口重赋；R与dmax必须同单位。 |
| `max_feeding_demand` dmax | 每幼虫完整取食窗口1个相同单位；保持归一化。实际分配受利用和同株竞争影响，不能把需求当成必然摄入。 |
| `resource_unit`、`demand_unit`、窗口 | 相同的 `normalized_edible_dry_biomass` 与 `one_generation_larval_feeding`；可另以同单位真实干重输入，但默认不是克。见N1。 |

Sarfraz等[1]全文方法与图3：温室4周龄油菜每株接10只初龄小菜蛾，5株/处理，5种施肥组约35–76%化蛹；**10为初始接虫量而非资源容量**，图估读而非容量测定。Wang与Lu[2]全文方法、图1：25±1°C第四龄自由取食白菜/甘蓝分别4.995±0.245、4.426±0.189 cm²/虫（均值±SE）；不是全龄干质量。两项支持同株多虫与取食有限的可能，不能由此推得R/dmax=10或最大值。

### R2 幼虫招募与资源限制

取食需求为 dmax×q，同株按资源约束分配，实际摄入 x；幼虫存活概率 `g=x/dmax`，`sI0=0` 和 `unmatched_insect_fitness=0`，零摄入无招募，满摄入本环节存活1。Wang与Lu[2]的第四龄限食影响化蛹与羽化，有限范围近线性不能验证全龄S=g及两端点。保留该闭合作对照；现有数值不是测得的背景死亡率，未来同实验分阶段拟合充分供食背景生存及摄入响应，避免重复扣死亡。

同dmax时每株招募**期望**为 `min(Σq,R/dmax)`，不是实际招募数的硬上限；异质dmax时由需求组成决定。昆虫不另设干扰密度因子，资源共享本身已造成竞争；不能把该资源闭合下的密度响应直接解释为寄主局部适应。植物本轮死亡不附加幼虫死亡；`R=0`定义`h=0,S_P=1,g=0`，这些时间窗口与端点简化详见[模型§4](model.md#4-防御资源与存活)。

### P1 植物繁殖、损伤与更新

公共植物潜在有效后代基准 `plant_offspring_lambda=2`，`sP0=.5` 是完全食害时单轮生存，生存函数为 `1−(1−sP0)h`；h为已耗资源比例（R=0按0）。植物存活且有局部合格父本后，母源有效后代按Poisson抽样并受成本与密度影响；该λ已压缩种子到建成，**不是种子数**。Akiyama与Ågren[3]全文图2及表2：瑞典自然拟南芥50%人工剪叶在莲座、抽薹期种子数分别下降约60%、22%，开花期未见降低；剪叶面积、结籽和成株生存不同，不能推出sP0。Brosset等[4]全文方法与讨论：21日龄白菜每株8只小菜蛾5天，平均叶面积损失约18%，花数未显著改变；花数不是有效后代。保留旧规则，待相同损伤尺度及阶段数据联合拟合，不能以这些数替换λ或sP0。

### P2 植物密度竞争

公共竞争矩阵同池α=.02、异池=.005，单祖先同池=.05，单位为每局部植株倒数；`competition_radius=1`格，计入焦点株，产量乘 `1/(1+ΣαN)`。Stouffer[5]全文§2.1式(4)提供年度植物倒数型密度繁殖的理论联系，指出焦点计数及λ与竞争强度关系，未验证本模型局部系数。Hauser等[6]全文表6：油菜与野生芸薹20–90株/m²，两个方向相对替代系数8.71（95%可信区间5.65–18.53）及0.71（0.33–1.15）；每处理仅一个地块，且可能含花粉效应，不能换为模型α。Muñoz-Parra等[7]的拟南芥1–32株/培养皿实验显示密度影响产量，但种子数非建成数。保留矩阵为假设，未来在明确格长的同系统密度×组成实验联合估计λ与α。

当前标签 `plant_density_stage=maternal_output_mean` 明确密度折减的是母源产量均值，不是子代落点存活；旧标签 `offspring_recruitment` 只按旧合同追溯，不作为当前别名。

### D1 昆虫繁殖、性别与交配

`insect_offspring_lambda=4` 为每有局部合格雄虫的雌虫、每模型轮潜在卵的Poisson基准；不等于成年招募。创始及子代雌概率各0.5；潜在卵经寄主接受与幼虫生存，父本对每枚成功卵从同池局部雄虫均匀有放回抽取，无精子耗竭。Yang等[8]全文表1：23±1°C、萝卜条件100卵中68成虫（30雌38雄），30雌产卵173.20±11.40卵/雌；初代豌豆17成虫、9雌，57.33±19.55卵/雌。羽化雌比0.441是条件于生存的比例，终生卵数不是潜在卵或招募的直接参数。Song等[9]全文表1：26±2°C，1日龄配对组交配100%、产卵雌占已交配雌85%，而6日龄对应50%与47.37%；多次相遇的实验不等于有候选即逐卵抽父本。原值保留为理论设定，待相同生命周期、机会和分母下拟合，检验泊松方差并避免重复乘死亡或成本。

### D2 生活史与繁殖结构

当前植物为严格异交、同兼容池且非自身；昆虫雌雄异体，非重叠代，无显式成虫跨代、移民或种子库。子代按既有双亲遗传，既有兼容池不是新物种。Postma与Ågren[10]野外研究在不同年份跟踪拟南芥种子休眠、萌发和成活，说明种子库可重要；拟南芥多自交，不能验证本模型严格异交。Yamamori[11]田间标记油菜后代的异交比例在0.25、1、5、10、30、60 m分别4.1%、1.4%、0.4%、0.15%、0.09%、0.01%；这些是后代比例而非花粉核，油菜可自交，不支持本模型硬半径和强制异交。保留为有条件的机制简化，研究越冬、迁入或混合交配须另行定义过程与数据。

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

10×10周期网格、Chebyshev距离，无实测米制；搜索2、昆虫交配2、昆虫扩散1、植物花粉1、植物扩散1、竞争1格。扩散为每轴整数均匀抽样；初始位置默认随机。公共P1/P2各200株、H1 40只，单祖先P1/H1各40；初始密度随格与数量共同变化，是研究设计而非实测。Mo等[20]摘要：小菜蛾作物田标记重捕的平均位移估计14–35 m；Bolton等[21]摘要：2.83 ha田地雄虫重捕率7.93%、距离中位14 m，田界95 m；Ishii等[22]摘要：下风向信息素有效空间估计0.6–0.9 m。各测不同过程且有陷阱与边界限制，不是逐代格位移或硬交配半径。保留格值为模型尺度；如拟合自然系统，先定义格长、时间和边界并采同阶段轨迹。

大半径下周期取模可使目的格的位移重数不同，每轴最多相差1；这符合整数位移均匀的定义，不是三角分布错误。不以此无依据缩小合法半径，也不将目的格重采样为均匀。

### O1 观察阈值与参考面板

只读操作谱系的无相位QTL拷贝不匹配阈值0.25、最小群2、持续3代；参考植物面板 `chem_A=(.25,.65,.25)`、`chem_B=(.65,.25,.65)`，只作功能响应比较，不参与配对或实际取食。阈值、单链接桥接、样本量及窗口影响结果，不构成物种证据；保留原设并做预先声明的阈值稳定性比较，未知保留，不能反馈兼容池。

### N1 单位与归一化

`resource_unit`和`demand_unit`同单位，取食和需求窗口同为一代；性状、化学类型位置、表达和效用均无量纲，格与代没有自然长度或天数。数值1、0及坐标端点是尺度约定，不是生物观测。换到真实单位需同源测量并校准映射，不能混用鲜重、干重、叶面积或单株与单虫分母。

### T1 运行与数据协议

版本/输入路径、来源元数据、种子、重复、运行代数、记录间隔、最大个体/等位基因/基因型值保护、SLiM路径与超时、状态、子步骤与输出ID为计算和实验设计，不是生物实测参数。计算上限触发中止而不是灭绝；不改历史源码哈希。清单对应项按实际配置记录，只有明确完成的实验才能解释结果；验证输出使用 `codex/validation/`，正式实验须明确指定 `results/`。

### L1 退役字段

`sigma=.2`、`weights=[.5,.5]` 为旧接口输入校验，在当前模块效用不生效；`B`、`rho_voc`、`voc_bounds` 为退役空键。不得当成活跃参数标定，也不得把历史协议输出冒充当前模型。

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


## 补充系统与来源可靠性

以下为早期定向核查中未被上表覆盖的独立记录，保留其系统、分母和限制；本次仅迁移文字，未重新检索文献。

- **Golizadeh 等 2009，小菜蛾—甘蓝**：[DOI 10.1016/j.aspen.2009.05.002](https://doi.org/10.1016/j.aspen.2009.05.002)。核验作者全文方法与表1、3：25°C条件下报告卵到成虫累计存活80.95%、生命表平均世代时间18.22±0.33天；每温度检查200卵，幼虫发育另取50只，繁殖用15对。累计存活含多个阶段，不能直接乘入已表示相同死亡的g；生命表T也不是所有个体发育时长。
- **Green 等 2015，菜粉蝶—野生甘蓝**：[DOI 10.1371/journal.pone.0136884](https://doi.org/10.1371/journal.pone.0136884)。核验全文方法与表3：每株最初放2只二龄幼虫，死亡个体持续补入；幼虫24日龄时各叶色组的每株累计受食面积约180.1±23.95至253.4±32.52 cm²。它是补虫条件下每株累计面积，不是单虫全龄干质量，不能与小菜蛾数据拼接成资源比值。
- **Saeed 等 2019，小菜蛾温度生命表**：[DOI 10.1016/j.sjbs.2018.08.026](https://doi.org/10.1016/j.sjbs.2018.08.026)。全文表与摘要存在内部疑点：20°C时摘要R₀/GRR与表2恰差两倍；表1卵到雌成虫时间与首次产卵时间、阶段之和不一致。未澄清前不用于标定。文献不是本模型数值的自动授权。
