# 参数、输入与文献依据

本文集中维护配置入口、参数含义、单位、取值和证据。公式见[模型](model.md)，当前任务见[项目入口](../README.md#current-work)，全部字段见[292条参数清单](parameter_inventory.csv)。不同自然系统分开列；机制证据不等于当前数字已获经验标定。

## 配置入口与输入

[configs/current.yaml](../evo_sim/models/plant_insect/configs/current.yaml) 是默认入口：P1–P10十种初始化学类型各20株，H1初始200只；10×10周期网格每格3个建立位点，共300位点。默认5轮、两种有虫和无虫情景各4重复；`current_200gen.yaml`继承该设计并延长至200轮，配置不代表已运行。类型标签不是植物生殖隔离。数值是理论探索，不按共存结果拟合；展开参数与来源见[逐字段清单](parameter_inventory.csv)。旧配置及历史结果须由[冻结源码](../codex/history/before_simplification_0.1.11_20261106T034842Z/README.md)解释。

## 资源与植物人口学

令K为每株最大资源、B为可食存量、X为本轮消费。创始株B=K；新株B=.25K。`B_post=B−X`，亏损`δ=1−B_post/K`；无同龄无食害参照。存活概率`exp(−μ₀−μDδ)`，幸存株恢复为`K−(K−B_post)exp(−r)`。每轮幸存株提案均值`λP max(0,1−C_P) B_post/K`，其中λP=2；它现在是每轮而非旧版每5轮的基准，数值不等价。参数未按共存结果调节。每格3个位点只限制建立，不代表共享光水养分，也不构成全株碳收支。没有年龄状态、年度、花期、繁殖脉冲或显式花粉机制。植物父本从全局其他幸存植物中随机选择，不按距离或初始来源标签限制。

空位提案按植物扩散核到达目标格并抽签，未建立者消失；没有移民或种子库。K和需求为归一化可食干物质量，恢复率、死亡率、繁殖率和空间单位未针对自然系统标定。恢复和空位建立的机制可参考下方文献[27–32]，这些研究不能直接标定本模型数值。

## 遗传、昆虫与输入

植物和昆虫均通过SLiM双亲遗传与重组；植物后代的标签表示母系初始来源，性状由基因型重新生成。昆虫保持SLiM双亲有性遗传、局部异性交配和非重叠世代；出生幼虫依据资源竞争及自身处理能力存活。化学输入、QTL、表型映射和昆虫机制的含义如下；参考面板只用于反事实诊断，不决定真实配对或取食。

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

解为 `E^I`中心`(.38748329034947154,.6098508195430548)`，约`.3875/.6099`；每个相关位点两等位效应取`E_r/2±.025`，P1及昆虫参考详见[当前等位表](../evo_sim/models/plant_insect/inputs/current/ten_plants_alleles.csv)。计算精度记录设计方程的解，不代表测量精度。这使名义中心`I=D`；由于核非线性，不声称随机创始个体的平均I恰等于D中心。它是跨越截断阈值的理论机制参考，不是自然祖先估计、经验校准或按不灭绝目标拟合；感觉、植物、调谐、映射scale和突变率均不因此更改。

这是整体理论基线的修订，**并非只改变q的演化实验**：按未更改的`insect_cost=.125`，处理表达中心总量从1.6降为约`.997334`，母源产出乘子在其他轴固定时提高约`.0753332`。因此新旧人口轨迹差异不能全部归于防御反馈；创始q对照只诊断固定核和防御条件下的利用率响应，未分离繁殖成本效应。本轮不为抵消该代价变化而追加生态参数补偿。

枚举六个相关轴各五级的`5⁶=15,625`个表型组合（不是全部QTL基因型），得到两通道`I/D`范围约`.75197–1.32962`和`.83096–1.19960`，`q∈[.75197,1]`。按`q<1−1e−12`计，等权组合覆盖`11828/15625=75.6992%`；按创始独立四拷贝的`Binomial(4,.5)`权重计算，随机昆虫—植物配对的未饱和概率约`75.1226%`，平均q约`.948645`。两通道分别增大D `1e−4`可降低q的加权覆盖约`41.8398%`和`35.4719%`；阈值附近两方向探针的覆盖可以重叠。这些是功能空间及创始频率加权的覆盖，均不是已选择寄主上实际幼虫比例。启动诊断另按该次实际抽到的创始个体报告，后续资源拥挤、选择和遗传仍需实际账本检验。

### G2 新突变与重组

植物和昆虫均有标记重组，可按配置启用新突变；突变率定义为每出生后代每个位点每份传递拷贝的概率。公共基础配置双方开，十类型及默认入口双方关闭，数值不因简化而调整。现有突变率为每拷贝每位点每出生后代0.01；植物六个位点的12份传递拷贝至少发生一次突变的概率为 `1−0.99^12≈0.113615`，以各拷贝独立且突变启用为条件，不是实际性状演化速率。关闭突变没有新等位变异，但仍有选择和重组；进化速率由结果估计，不设额外速度参数。抽象QTL事件率不是核苷酸突变率；Ossowski等[13]的测量不能直接转换为本模型概率。植物和昆虫同轴标记之间的Haldane重组是理论设定，不是实测连锁图。

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

默认10×10周期网格，每格3个共享建立位点，初始十种植物各20株、H1 200只。种子扩散及昆虫搜寻、交配、扩散半径是抽象格尺度，不是实测米制；具体展开值见配置和清单。不同物种的移动观测[20–22]不能直接转换为本模型逐代位移或硬半径。

### O1 观察阈值与参考面板

只读操作谱系的无相位QTL拷贝不匹配阈值0.25、最小群2、持续3代；参考植物面板由当前输入表指定，只作功能响应比较，不参与配对或实际取食。阈值、单链接桥接、样本量及窗口影响结果，不构成物种证据；保留原设并做预先声明的阈值稳定性比较，未知保留，不能反馈兼容池。

### N1 单位与归一化

`resource_unit`和`demand_unit`同单位，取食和需求窗口同为一代；性状、化学类型位置、表达和效用均无量纲，格与代没有自然长度或天数。数值1、0及坐标端点是尺度约定，不是生物观测。换到真实单位需同源测量并校准映射，不能混用鲜重、干重、叶面积或单株与单虫分母。

### T1 运行与数据协议

版本/输入路径、来源元数据、种子、重复、运行代数、记录间隔、最大个体/等位基因/基因型值保护、SLiM路径与超时、状态、子步骤与输出ID为计算和实验设计，不是生物实测参数。计算上限触发中止而不是灭绝；不改历史源码哈希。清单对应项按实际配置记录，只有明确完成的实验才能解释结果；验证输出使用 `codex/validation/`，正式实验须明确指定 `results/`。

### L1 退役字段和旧时间尺度

`sigma=.2`、`weights=[.5,.5]` 为旧接口输入校验，在当前模块效用不生效；`B`、`rho_voc`、`voc_bounds` 为退役空键。不得当成活跃参数标定，也不得把历史协议输出冒充当前模型。旧sP0、alpha、competition输入、competition_radius和母源φ不进入当前协议；植物λ仍记为2但现为每轮建立提案的基准，不能与旧版每5轮的数值视作等价。

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
