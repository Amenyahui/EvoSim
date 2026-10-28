# evoSim v0.1.6 — 当前输出与解释边界

更新：2026-10-28。当前生活史 adult_oviposition_larval_feeding，化学 shared-chemical-1，遗传 slim-qtl-2，执行器 SLiM 5.2。本文是阶段字典；验收状态见 [实施记录](plan/v0.1.6_implementation.md)，不包含演化结果分析。

## 1. 目录、身份与时间

新批次默认 results/results_v0.1.6/runN/，显式目录必须不存在。批次有 settings.json、data/run_status.csv，逐运行在 data/<case>/repNN/，审计在 audits/，派生汇总在 data/summary/。计算停止与生物灭绝分开。

主键含 case/replicate；个体再含 kind/individual_id，两类群 ID 不保证全局唯一。reproduction_round=t 是 Adult(t) 的繁殖轮，adult_cohort=t+1 是本轮子代成功补充的成年队列；offspring.generation=t+1 指子代队列。events/larvae.generation=t 指幼虫事件轮，不能仅按同名 generation 合并成年和幼虫记录。

## 2. 原始文件

| 文件 | 实际阶段及语义 |
| --- | --- |
| metadata.json、status.json | 展开配置、输入/源码/架构哈希、生活史/协议/SLiM身份、停止状态及文件哈希 |
| snapshots/gNNNNNN.npz | 该代 P_t/Adult(t) 的原发/派生性状、基因型、位置、性别、出生寄主与队列；终点也有快照 |
| populations.csv | 每次实际更新的代初数量、潜在卵、出生数、下一成年数量、首次物种灭绝代与更新标识 |
| oviposition.csv | 全部潜在卵的母体选择/拒绝及候选；成功事件才有 offspring_id |
| offspring.csv | 全部实际出生（含随后死亡幼虫）的真实双亲、基因型、出生寄主、事件映射、队列及存活结果；不是仅幸存者 |
| larvae.csv、events.csv | 实际出生幼虫的自身 I、q、需求、消费、收益、植物关联和一次存活结果；同一幼虫事件的两种输出入口 |
| plants.csv | P_t 受本轮幼虫食害的资源、需求、消费、损伤、成本及空间/密度 |
| parents.csv | Adult(t) 或母/父株候选、母源率、潜在卵、出生与补充两套双亲计数 |
| host_groups.csv | 每代成年队列按出生背景的数量及 present/absent，含 unknown 和固定所有寄主背景 |
| host_flows.csv | 出生/补充两个阶段的亲本出生来源→子代出生目的贡献矩阵，完整含 unknown |
| networks.csv | 固定物种层本轮幼虫消费、出生及成功补充关联；无遗传贡献同义替换 |
| lagged_networks.csv | 当前真实双亲按其自己的出生背景归属本轮出生/补充贡献，含 unknown |
| genetic_loci.csv、allele_catalog.csv、mutations.csv | QTL、类群特异效应长度与全部出生者的突变，不限于补充者 |
| gamete_transmissions.csv | 全部出生者双亲位点来源、同源标记切换和突变关联；不是完整树序列 |
| genetic_summary.csv、haplotype_summary.csv | 成年/植物代初状态的原发性状、等位/单倍型统计 |
| allele_contributions.csv | 按 stage=birth/recruitment 分别比较阶段子代及双亲贡献/传递分解 |
| diagnostics.csv | 成年 A/I 基因型的反事实参考诊断；不是成年当轮真实取食 |
| genotypes.npz、phenotypes.npz | 逐代成年/植物状态汇集；死亡幼虫基因型见 offspring/位点传递而非成年快照 |

终点不追加下一轮 parents、oviposition、larvae、events、plants 或 networks 事件。终点没有相应事件行表示未执行，不能解释为该轮真实摄食/补充恰为零。初始诊断的幼虫损伤、繁殖期望也不能由成年快照伪造。

## 3. 事件映射、背景与位置

oviposition.event_id=`round:mother_id:egg_index`，在一次运行内唯一；offspring_id 在拒绝时为空。offspring/events/larvae 用 oviposition_event_id 与 egg_index 回链。卵被拒绝时 host_id=-1、plant_species=unknown 只表示没有选择，不是存在未知植物上的幼虫。

offspring 的 birth_host_id/birth_host_species 是实际产卵植物；maternal_birth_host_species/paternal_birth_host_species 是双亲各自出生背景，不能互换。创始成年 birth_host_id/species 内部为-1，文本背景为 unknown。任一亲本背景未知时 cross_background_pair 为 NA，pair_background_status=unknown；两者已知时才分别 same/different 与 false/true，不能把 unknown 当同背景。

昆虫 offspring_cell/cell_row/col 表示出生寄主格；成年快照位置为幸存后扩散位置，两者可不同。幼虫出生至成年保持同一 ID。植物 offspring 仍从母株格扩散，其 survived_to_adult=True 是两类群通用记录字段，表示没有额外后代筛选，不将植物称为成虫。

host_groups 是每代出生背景分类，无遗传独立性判据，空群可随后重新占据。populations.extinct_at_generation 是固定物种首次灭绝代；物种灭绝为吸收态，背景暂时 absent 不是物种灭绝或谱系终止。全部 offspring 的双亲构成 DAG，不生成候选独立谱系标签。

寄主背景群视图由 host_groups、larvae/events、host_flows 和 lagged_networks 按轮次、物种与背景连接得到；没有另行保存的候选谱系网络或独立背景网络 CSV。

## 4. 母源率、出生与补充贡献

parents 中 F_i 是母源 Fm；昆虫 S_i=1、survived_to_mating_stage=True 表示现有成虫已经合格，不是重新测得其幼虫存活率。昆虫 conditional_mu 是 lambda_egg*Fm(I_m)，expected_maternal_offspring / maternal_offspring_variance 表示给定成年状态和父本机会的**潜在卵**矩，不含产卵拒绝或幼虫竞争。potential_egg_count 是实际潜在卵数；雄虫母源矩 NA、实际母源计数零。植物对应字段仍是含其存活、配偶与密度的有效植物后代矩。

| parents 字段 | 阶段终点 |
| --- | --- |
| offspring_count_as_mother、offspring_count_sired | 实际出生子代的母/父计数 |
| realized_genetic_contribution | 上述两计数之和/2，出生阶段 |
| recruited_count_as_mother、recruited_count_sired | 幸存补充子代的母/父计数 |
| realized_recruitment_contribution | 上述两计数之和/2，补充阶段 |

同种每轮两套账分别 ΣBm=ΣBf=Σw=N_birth / N_next；幼虫死亡时两套数不同。植物没有额外幼苗筛选时两套相同。F_m=0 仍可父源供源，不把 SFm 再乘实际贡献。parents.W_i=S_i*Fm 是阶段代理，昆虫现有成虫为 Fm；larvae/events.W_insect 本版为 g（幼虫资源存活代理），不再是旧成体摄食模型 g*Fm，也不是完整适合度。

populations.N_survived_to_reproduce 对昆虫是当前成年数量，不是本轮幼虫补充数；N_offspring 是出生幼虫数，N_next 是成功补充成年数，N_potential_eggs 为潜在卵数。zero_fitness_count 是历史代理字段，不表示所有遗传功能丧失。未提交或未执行的实际贡献缺失；实际已提交但没有后代时计数为零，零总数的归一化频率仍未定义。

allele_contributions.stage 必须参与主键。每阶段 p_offspring 及 N_offspring 指该阶段全部出生者或成功成年者，parent_weighted 采用同阶段 w；不能用出生权重分解成年频率。位点传递可用于实际等位来源核验，父母各1/2只是角色来源份额。

## 5. 来源矩阵和三种网络边

host_flows 的 reproduction_round、adult_cohort、insect_species、stage、source_host_species、destination_host_species 确定矩阵单元。每子代把母/父各1/2从亲本出生背景投向子代出生背景；stage=birth 全矩阵合计出生数，recruitment 合计成年数，unknown 不得删去。P1/P1→P2 与 P1/P2→P3 都在该账本内。它不是有效基因流率、迁移率或成种证据。

networks.Q_actual 是同物种幼虫到植物的真实消费；N_accept/birth_count 是实际出生幼虫数，successful_recruits 是其中成功成年数；edge_semantics=larval_resource_and_recruitment。resource_share 分母是该昆虫物种本轮全部实际消费，零消费时缺失。L_per_insect 本版分母是**该寄主上的出生幼虫数**，不再是旧成年物种总人数。

lagged_networks.plant_species 指**亲本出生寄主背景**，不是本轮卵目的地。maternal_offspring/paternal_offspring 是按该来源亲本归类的出生母/父计数；genetic_contribution 为出生角色各半贡献，recruitment_contribution 为幸存补充角色各半贡献；edge_semantics=parent_birth_host_lagged_contribution。unknown 单列。某个子代成年补充只构成补充关联，只有它后来实际成为双亲才有自己的滞后繁殖贡献。

实际资源量、补充数和滞后来源计数不可混作一个边权；背景标签和网络不改变生态参数或交配资格。

## 6. 三模块、参考诊断与历史兼容

plant_chemical_traits 每行是 x_N/x_S/x_V，N_P×3；D、V_emit 及反事实线索为 N_P×2，I/A 各 N_H×2。原发统计 trait 为 plant_chemical_traits_1/2/3 与 I_1/2/A_1/2；D/V 是确定性输出，不另算独立 QTL。等位目录效应向量植物3维/昆虫4维。

元数据保存 chemical_schema、genetics_version、plant_trait_order、plant_defense_order、odor_trait_order、chemical_observation_window、eta_shared 及生活史。plants.defense_cost=cD(x_N+x_S)、voc_cost=cV*x_V；eta 改变线索但不直接改变成本。resolved_B=null，旧 B/rho/Z 不参与新规则。

standardized 的成年 A/I 参考 p/q 联合谱可作基因型潜在功能诊断；A 的母体选择作用和 I 的幼虫作用属于不同生命阶段，p*q 不代表该成虫真实当轮摄食或繁殖适合度。参考 p 不含外部拒绝，参考 q 无资源竞争；逐个体宽度与群体平均宽度不同，全零 q 宽度未定义。

新批次先 audit，再运行 analysis.summarize_lifecycle；它合并以下原始账本到 data/summary/，添加 case/replicate，不新增网络指标：lifecycle_counts（populations）、host_background_flows（host_flows）、host_background_groups（host_groups）、lifecycle_networks（networks）、parent_host_contributions（lagged_networks）。

旧 summarize_qtl_run、summarize_coevolution 和依赖 Batch 的 plot_diagnostics **明确拒绝 v0.1.6**，避免沿用成年自主取食/单阶段贡献假设；历史结果的只读路径继续保留，需匹配原协议和审计。历史200代分析或旧测试数不能当本版验证，不改旧元数据/源码哈希，不重命名旧 Z 或旧贡献字段来冒充新运行。
