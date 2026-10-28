# evoSim v0.1.5 — 当前参数与输入指南

更新：2026-10-28。当前化学/遗传契约为 shared-chemical-1 / slim-qtl-2。YAML保存通用参数、CSV路径和场景覆盖，CSV保存逐条物种/遗传输入。新规则与历史四轴配置不能仅凭相同版本号混用。

## 1. 文件分工与数据形状

| 文件（相对evo_sim/models/plant_insect） | 当前契约 |
| --- | --- |
| configs/model_parameters.yaml | 公共参数、模型0.1.5、遗传协议2、CSV路径 |
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
| insect_cost=c_I | 对ΣI的母源产出惩罚，默认0.125 | 非负有限；昆虫函数及角色解释保留 |

当前化学能力和线索均为一个繁殖前取食世代内的归一化功能尺度，窗口写入chemical_observation_window。真实储量、排放流量及释放耗竭未模拟。Fm是无量纲母源乘子，非总能量余额，也不是父源授精能力。

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
| 资源 | resource_capacity、max_feeding_demand均正，共享resource_unit/demand_unit及feeding_window/demand_window |
| 选择 | candidate_size、sigma、b、tau、outside_utility、olfaction_mode、choice_mode、shuffle_voc均保持原意义 |
| 存活 | sP0∈[0,1]，当前sI0及unmatched_insect_fitness强制0；自主取食生活史 |
| 母源繁殖 | 各物种offspring_lambda>0，单位为每合格母本每模型代有效后代；ν=λFmφ |
| 密度 | 植物完整alpha矩阵≥0，按母株邻域计算φ；offspring_recruitment为历史标签，不是落点招募 |
| 空间 | 周期grid_shape，搜寻/竞争/母源扩散/花粉与昆虫交配半径分别指定，单位为格 |
| 生殖角色 | 两性植物严格异交；昆虫雌雄异体；父本均匀有放回，不以Fm筛选或加权 |
| 性比 | 创始与子代雌性概率分别在[0,1]，不强制数量相等 |
| 计算保护 | max_individuals_safety、max_alleles_safety、max_genotype_values_safety仅停止保护，不截断出生 |
| 复现/窗口 | seed、replicates及simulation_metadata.num_steps_per_episode；当前逐代完整记录 |

## 4. 配置继承及历史迁移

继承顺序为父YAML→本层CSV整表替换→本层显式参数→子层及处理覆盖。路径相对声明该路径的YAML解析；字典递归合并、列表整体替换。展开配置、输入CSV路径/哈希、完整遗传架构和新坐标顺序写入元数据。新运行要求模型0.1.5及协议2；旧协议只读审计/分析仍可使用，重建旧模拟使用原代码/配置。

本轮显式迁移x_N←旧D1、x_S←旧D2、x_V←旧Z2，删除旧Z1两位点和相关等位，重新连续编号等位ID；默认植物6位点，昆虫8位点。保留其余位点效应、频率和成本系数。参考面板同样迁移；不能自动把已有快照Z转为新x。祖先V因此改变，旧“P1适配”标签不再是既定事实。

input_sources和协议/配置指纹区分新旧数据。eta_shared的Python缺省None仅用于读取缺少该键的历史配置；新协议中省略会验证失败，不是隐藏的生物默认值。历史config_payload去除该空键以保持原配置哈希。

## 5. 依据和观察参数

三模块分类、线性映射、等权c_D和eta=1为理论/标度选择；综述及理论先例支持的范围见 [model.md §7](model.md)。分析设计YAML仍只用于观察，早期/后期、共同R*/d*、边阈值和连续窗口均不反馈模拟。未提交贡献、无后代的频率及全零q宽度等统计保留缺失语义，见 [results.md](results.md)。
