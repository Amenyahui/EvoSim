> 历史文档：以下“当前”“本版”等词仅指 v0.1.6；链接已调整到归档位置。

# v0.1.6 — 保留输入合同与核定配置（原§1–7）

更新：2026-10-28。当前化学/遗传契约为 shared-chemical-1 / slim-qtl-2。YAML保存通用参数、CSV路径和场景覆盖，CSV保存逐条物种/遗传输入。新规则与历史四轴配置不能仅凭相同版本号混用。

## 1. 文件分工与数据形状

| 文件（相对evo_sim/models/plant_insect） | 当前契约 |
| --- | --- |
| configs/model_parameters.yaml | 公共参数、模型0.1.6、遗传协议2、CSV路径 |
| config.yaml | 既有调度；植物增加chemical_traits状态属性 |
| inputs/species.csv | 物种、数量、单株R或单虫dmax、基础λ及来源；P/H编号各自连续 |
| inputs/competition.csv | 完整植物物种对的alpha及来源 |
| inputs/loci.csv | 14条：6植物、8昆虫；默认每原发轴2QTL，图谱/突变/目标列 |
| inputs/alleles.csv | 40条创始等位，ID从0连续；效应可正可负 |
| inputs/founder_frequencies.csv | 40条；每物种每位点和为1，涵盖所属类群所有位点 |
| inputs/reference_traits.csv | 每植物物种一行species_id,x_N,x_S,x_V |
| inputs/v012/*.csv | 同步迁移的已有场景输入；名称表示历史来源，不保证与旧版生态等价 |
| analysis/v014_analysis.yaml | 只读分析的时期、共同资源和边持续窗口设计 |

CSV为UTF-8（可含BOM），一格一个值，不使用JSON向量。species/competition/founder表原列名继续使用；loci/alleles的共同表有0–3四组向量列，因为昆虫仍4维。**植物只填写0–2，第3坐标的所有格必须为空，连0也不能填；昆虫四格全部必填。** 读入的植物向量长度确为3，不把空格补零当第四轴。

| 类群 | 坐标0 | 坐标1 | 坐标2 | 坐标3 |
| --- | --- | --- | --- | --- |
| 植物 | x_N | x_S | x_V | 不存在，CSV空格 |
| 昆虫 | I_N | I_S | A_S | A_V |

loci列为locus_id,organism,functional_module,chromosome_id,position_cM,mutation_rate，以及target_j、mutation_mean_j、mutation_sd_j。target取0/1；目标索引和效应/突变向量必须与本类群维数一致，非目标效应/增量均为零。alleles为allele_id,locus_id,effect_j。引用未知位点、重复键、缺列、非有限值或错误频率均拒绝。

## 2. 新化学参数与退役项

| 键 | 意义、单位与默认值 | 有效域/边界 |
| --- | --- | --- |
| mapping | 显式shared_chemical_modules | 不与identity/mixed旧映射混用 |
| eta_shared | 共享x_S到V_S的固定有效映射系数；默认1为标度约定 | 数值且有限，[0,1]；不接受bool或缺失，不是实测挥发率 |
| n | 两条D/I通道，当前固定2 | 本三模块构型不支持任意改维数 |
| B、rho_voc、voc_bounds | 新配置中明确null，保留键仅为历史反序列化 | 不能偷换成eta、不能覆盖为旧数值后静默参与新规则 |
| weights | 两条V/A距离权重，默认[0.5,0.5] | 非负有限、和为1 |
| plant_cost=c_D | 对x_N+x_S的母源产出惩罚，默认0.125 | 非负有限，归一化化学性状单位的倒数；共享模块只收一次 |
| voc_cost=c_V | 对x_V的母源产出惩罚，默认0.05 | 非负有限；不对V均值计费，也不另收共享释放成本 |
| insect_cost=c_I | 对成年母亲ΣI的潜在卵产出惩罚，默认0.125 | 非负有限；昆虫函数及角色解释保留 |

当前化学能力和线索均为一个幼虫取食窗口内的归一化功能尺度，窗口写入chemical_observation_window。真实储量、排放流量及释放耗竭未模拟。Fm是无量纲母源乘子，非总能量余额，也不是父源授精能力。

eta=0时共享模块仍提供防御，但无对应挥发线索；eta=1是参考映射。改变eta不直接改变母源成本。成本系数可大于1，Fm非负截断保证[0,1]；不能为了避免期望外轨迹而事后调节。所有新示例仍非实测校准。

## 3. 遗传、人口与空间参数

| 类别 | 键及意义 |
| --- | --- |
| 遗传身份 | genetics.schema_version=slim-qtl-2；ploidy=2、biparental_meiosis、slim_haldane_marker_map保留 |
| 原发背景 | genetics.background.plant长度3、insect长度4，默认零 |
| 表型映射 | effect_scale默认0.5、bounds默认[0,1]；前者为拷贝效应尺度，不是显性系数 |
| 突变 | loci CSV逐拷贝概率及均值/SD；默认0.01/0/0.02为理论示例，不是核苷酸率 |
| 创始组成 | allele/founder CSV；初始频率必须和为1，不按运行结果修改 |
| 引擎 | genetics.slim.executable指向SLiM5.2，EVOSIM_SLIM可覆盖；timeout_seconds是计算时限 |
| 数量 | plant_sizes/insect_sizes来自species CSV；insect_pop=null时由各物种数派生 |
| 资源 | resource_capacity≥0、max_feeding_demand>0；零资源时消费、食害比例与幼虫存活为零，共享resource_unit/demand_unit及feeding_window/demand_window |
| 选择 | candidate_size、sigma、b、tau、outside_utility、olfaction_mode、choice_mode、shuffle_voc用于每枚潜在卵的母体选择；不读取子代I |
| 存活 | sP0∈[0,1]；幼虫S=g，sI0及unmatched_insect_fitness强制0；现有成虫不再抽幼虫存活 |
| 母源繁殖 | 植物offspring_lambda>0为每合格母株每轮有效后代（ν=λPFmφ）；昆虫insect_offspring_lambda>0为每合格雌成虫每轮潜在卵（ν=λegg Fm(I_m)） |
| 密度 | 植物完整alpha矩阵≥0，按母株邻域计算φ；offspring_recruitment为历史标签，不是落点招募 |
| 空间 | 周期grid_shape，搜寻/竞争/母源扩散/花粉与昆虫交配半径分别指定，单位为格 |
| 生殖角色 | 两性植物严格异交；昆虫雌雄异体；父本均匀有放回，不以Fm筛选或加权 |
| 性比 | 创始与子代雌性概率分别在[0,1]，不强制数量相等 |
| 计算保护 | max_individuals_safety、max_alleles_safety、max_genotype_values_safety仅停止保护，不截断出生 |
| 复现/窗口 | seed、replicates及simulation_metadata.num_steps_per_episode；当前逐代完整记录；reproduction_round与adult_cohort分开 |

## 4. 配置继承及历史迁移

继承顺序为父YAML→本层CSV整表替换→本层显式参数→子层及处理覆盖。路径相对声明该路径的YAML解析；字典递归合并、列表整体替换。展开配置、输入CSV路径/哈希、完整遗传架构和新坐标顺序写入元数据。新运行要求模型0.1.6及协议2；旧协议只读审计/分析仍可使用，重建旧模拟使用原代码/配置。

v0.1.5 已显式迁移x_N←旧D1、x_S←旧D2、x_V←旧Z2，删除旧Z1两位点和相关等位，重新连续编号等位ID；默认植物6位点，昆虫8位点。保留其余位点效应、频率和成本系数。参考面板同样迁移；不能自动把已有快照Z转为新x。祖先V因此改变，旧“P1适配”标签不再是既定事实。

input_sources和协议/配置指纹区分新旧数据。eta_shared的Python缺省None仅用于读取缺少该键的历史配置；新协议中省略会验证失败，不是隐藏的生物默认值。历史config_payload去除该空键以保持原配置哈希。

## 5. 依据和观察参数

三模块分类、线性映射、等权c_D和eta=1为理论/标度选择；综述及理论先例支持的范围见 [model.md §7](model_v0.1.6.md)。分析设计YAML仍只用于观察，早期/后期、共同R*/d*、边阈值和连续窗口均不反馈模拟。未提交贡献、无后代的频率及全零q宽度等统计保留缺失语义，见 [results.md](../results.md)。

## 6. v0.1.6 生活史与阶段参数合同

| 键/输入 | 当前含义与限制 |
| --- | --- |
| model_version | 0.1.6；新运行仍为 slim-qtl-2 / shared-chemical-1 / SLiM 5.2 |
| insect_life_history | adult_oviposition_larval_feeding |
| feeding_window、demand_window | 均为 one_generation_larval_feeding，同一幼虫窗口 |
| resource_unit、demand_unit | 默认均为 normalized_edible_dry_biomass；允许配置同一受支持干物质量单位；R≥0、dmax>0 |
| insect_offspring_lambda | 由 species.csv 昆虫行 offspring_lambda 读取或显式覆盖；每合格雌成虫每轮潜在卵，单位 eggs/female/round |
| plant_offspring_lambda | 每合格母株每轮有效植物后代；仍受 Fm 和局部 phi，单位 effective_offspring/maternal_plant/round |
| insect_mate_radius | 当前成年坐标的同种异性候选半径，默认2格；不以本轮寄主或Fm限制雄虫 |
| search_radius、candidate_size | 每枚潜在卵以母体坐标搜索；默认2格、8株候选，不是幼虫移动范围 |
| insect_dispersal_radius | 默认1格，成功幼虫从自身出生寄主扩散到下一轮成年位置 |
| 创始出生寄主 | 当前未提供已知寄主初始化输入；统一 unknown/-1，不用初始位置推断 |

lambda 沿用的旧示例数值没有重新校准；新增拒绝和幼虫存活后，不能与旧有效补充数视作等价。Fm(I_m) 仅使用母亲自身抗性，不重复乘她前轮资源收益。每卵独立、固定幼虫寄主、线性 S=g、均匀父本和各半来源贡献均是操作性假设；不能为保证人口、专化或网络方向而事后调参。

时间/寄主/事件字段是追踪状态而非新选择参数。出生背景群或亲子标签不改变交配兼容性，寄主源→目的矩阵不是有效基因流率。历史观察 YAML 继续用于原版本只读分析；依赖旧成年摄食/单阶段贡献定义的入口拒绝新生活史，参见 [results.md](../results.md)。

## 7. 显式初始条件与核定配置

| 配置 | 用途及覆盖 |
| --- | --- |
| v016_ancestral_matched.yaml | 已有功能匹配的理论起点；使用ancestral_matched_alleles.csv，仅A_S两位点的四条等位效应各加0.2，期望A_S增加0.4；I、频率、零背景及lambda=4保持原值 |
| v016_recognition_mismatch.yaml | 原始识别不匹配情景；显式读取原alleles.csv，利用已偏向P1，不是完全未适应祖先 |
| v016_reference_validation.yaml | 匹配起点、50代、种子20261021–20261024、关闭双方新突变；P1单独与双寄主各4重复；只写codex/validation |

前两者继承短程开发配置的代数与突变开/关处理，正式实验必须另行明确实验窗口和设计。匹配起点不修改历史默认、出生寄主标签或参考诊断面板。lambda=4仍为潜在卵的理论标度，无实测鳞翅目卵数依据；理想化参考替代能力不构成生活史校准。核定结果见 [报告](../plan/v0.1.6_reference_qualification_2026-10-28.md)。
