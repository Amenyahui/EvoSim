# 输出、审计与绘图

适用协议与科学边界见[模型](model.md)，输入和数值见[参数](parameters.md)；这里定义当前输出字段和分析命令。当前代码只执行和审计模型0.1.7 / slim-qtl-4 / expression-tuning-2。旧数据保持原协议、原始文件和指纹，需用对应冻结源码读取和审计。

## 当前模块与阶段语义

| 文件/字段 | 定义 |
| --- | --- |
| metadata.version、software_version | 前者是实际模型协议0.1.7，后者是运行软件版本0.1.8；两者分别记录 |
| snapshots/phenotypes中的insect_module_traits | N×8原发感觉/处理表达与调谐；顺序由元数据insect_module_order保存 |
| plant_chemical_traits | N×3原发x_N/x_S/x_V；当前协议不保存花期/theta状态 |
| I / A | I为两类底物的派生有效容量；A为固定目录上的无符号感觉表达profile，实际选择还经过饱和和行为价性，不是旧偏好模板 |
| parents.csv新增字段 | mating_row/col、candidate_mate_ids/probabilities、两项昆虫母源表达成本；植物存活、局部同池非自身候选的概率均匀 |
| genetic_summary | 植物三化学轴与昆虫八模块轴的raw/clip及遗传统计；当前没有wrapped_fraction统计 |
| reference_plant_species | 新协议中是显式参考化学型ID，可与祖先兼容池数量不同；不是实际寄主事件 |

当前无adult_settlement.csv及metadata.settlement_semantics。成虫从出生寄主扩散后直接按成年位置建立局部兼容异性池，随后母体逐卵择主；parents的交配位置/候选字段仍用于审计该顺序，并非定居专属输出。历史expression-tuning-1中的定居账本按原身份保留，不补入新批次，也不据此回填新配偶池。

`summarize_lineages --input <batch> --design <yaml>`仅观察已审计批次，输出`data/summary/lineages/`。阈值、最小规模及持续窗口不进入模拟：

- memberships/clusters：同祖先池的无相位QTL拷贝距离图分群；小群标为unassigned，不能据此说已灭绝。
- transitions/lineage_events/populations：补充子代双亲贡献、多对多分裂/合并及单元终止；祖先池灭绝与观察未检出分别记录。一个桥接基因型即可造成单链接合并，不能直接当作真实基因库融合。
- mating：实际出生/补充子代加权的跨观察群双亲数及明确分母；无可分类双亲时比例缺失。
- mating_opportunities：每个合格母体等权的保存候选池概率，不与子代加权的实际交配比例混为同一量。
- networks：`maternal_source`记录母源群t→植物群t的幼虫消费/出生/补充；`recruit_cluster`记录成年补充群t+1→出生寄主植物群t，仅记录补充，消费字段缺失。
- network_balances/metadata：两套网络与原始资源/出生/补充总账、观察设计/代码/输入指纹和限制。消费包括死亡幼虫；unknown/unassigned仍参与总量。

这些节点是操作性候选谱系。持续窗口只使用当前及过去轮次，不能将未来存续倒填为早已确证物种。植物化学性状仍可演化，但当前没有额外植物适应性成种机制；没有输出经验意义的物种形成率或新物种数。统一绘图入口支持当前v0.1.7的A–H图及敏感性模式，生命周期和谱系数据汇总均可用。历史slim-qtl-3文件中的花期/回绕字段仅按原保存协议解释，不补入当前新输出。

`plot_diagnostics`沿用A人口、B性状、C寄主利用/繁殖阶段、D空间、E遗传的布局和PNG/PDF/SVG导出。B使用植物x_S/x_V、x_N/x_S，以及昆虫met E₁/E₂、met μ₁/μ₂、sense E₁/E₂、sense μ₁/μ₂六对原发轴，均相对各运行第0代中心化；派生I/A不当作旧识别模板。C以幼虫实际消费计算寄主份额；没有H1灭绝时下排改为最后20代成年丰度和末轮成年→交配雌虫→潜在卵→幼虫出生→成年补充，不能把期末删失当灭绝。A追加终点成年快照，E保留位点平均H_e，不填充灭绝后的性状或遗传值。

当前分析读取`audits/integrity_audit.json`。绘图前核对实际使用原始CSV/快照与metadata.output_sha256，绘图元数据另存当前绘图/汇总代码和输入指纹，保留原始模拟源码哈希。历史修正批次的专属读取逻辑随原源码冻结，不通过当前入口改写其身份。

## 1. 目录、身份与时间

开发与验收批次按软件版本默认写入codex/validation/v0.1.8/runN/，显式目录必须不存在。测试、冒烟和临时验证结果不写入 results/；正式科研实验需显式指定 --out results/...。批次有 settings.json、data/run_status.csv，逐运行在 data/<case>/repNN/，审计在 audits/，派生汇总在 data/summary/。计算停止与生物灭绝分开。

主键含 case/replicate；个体再含 kind/individual_id，两类群 ID 不保证全局唯一。reproduction_round=t 是 Adult(t) 的繁殖轮，adult_cohort=t+1 是拟补充子代队列，即使幼虫随后死亡仍保留该标签。generation 的阶段定义如下，不能只按同名 generation 合并：

| 表 | generation | reproduction_round | adult_cohort |
| --- | --- | --- | --- |
| offspring、mutations、gamete_transmissions | t+1，出生子代队列 | t | t+1 |
| larvae、events、networks | t，实际幼虫资源轮 | t | t+1 |
| parents | t，亲本繁殖轮 | t | 不适用 |
| snapshots、genetic_summary、host_groups | t，代初成年/植物状态 | 不适用 | 昆虫快照另保存各个体所属队列 |

新输出以metadata.stage_key_schema标识显式阶段键；独立审计核对出生表的generation=adult_cohort=reproduction_round+1。幼虫与昆虫子代表用(case, replicate, reproduction_round, adult_cohort, insect_id/offspring_id)逐个体连接；连接位点传递时另含kind、role和locus_id，不能只按轮次连接而制造笛卡尔积。旧CSV不回填或重命名。

批次CLI只在全部运行以completed或all_species_extinct结束时退出0；backend_error、computational_limit、arithmetic_error和invalid_parameter导致非零退出并保留失败账本。完整性审计仍可验证失败记录，因此`integrity_status=passed`不代表模拟成功。报告另列execution_success、successful_runs、zero_transition_runs和acceptance_passed；运行验收须明确成功状态及每个运行的最小已提交转移数，例如：

```powershell
python -m evo_sim.models.plant_insect.audit --input codex/validation/multilarva_example --require-success --min-transitions 1
```

这只检验至少一轮成功提交；固定期限完成还须核对请求代数和停止状态，不能据此宣称达到200代。零代初态诊断可以完整成功且转移数为0；用diagnose_initial_conditions读取，生命周期诊断绘图明确拒绝含零转移运行的批次。敏感性汇总保留0代观测并区分已知灭绝、存续和未知。

max_genotype_values_safety限制累计保留的全部出生者基因型值与逐代快照基因型值，metadata.stored_genotype_values记录同一口径。死亡幼虫仍有出生账本；幸存者又出现在成年快照，两个表示均需计入，不能以峰值活体量代替累计存储量。每次提交前预留下一快照，初始和终点也受保护；超限整轮停止，不裁剪人口。gamete_transmissions每行对应一个出生QTL拷贝，其行数受同一出生预算约束；该预算是计数保护，不是CSV字节数或进程峰值内存上限。

shuffle_voc的线索置换使用独立phase=8随机流，宿主选择顺序仍使用phase=1；关闭时不额外抽样。同线索条件下切换不再消耗选择顺序流，但状态分化后不保证不同处理的全部事件仍逐个配对。

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

host_groups 是每代出生背景分类，无遗传独立性判据，空群可随后重新占据。populations.extinct_at_generation 是固定物种（v0.1.7祖先兼容池）首次灭绝代；该灭绝为吸收态，背景暂时 absent 不是物种灭绝或谱系终止。全部 offspring 的双亲构成 DAG，原始记录不写入候选谱系资格；v0.1.7操作群在只读分析中另行识别。

寄主背景群视图由 host_groups、larvae/events、host_flows 和 lagged_networks 按轮次、物种与背景连接得到；上述原始文件不另造背景网络CSV。v0.1.7候选谱系网络另写在批次data/summary/lineages/，不替换原始寄主背景账本。

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

## 6. 当前模块、参考诊断与历史兼容

当前植物原发性状为 N_P×3，派生 D、V 各 N_P×2；昆虫模块性状为 N_H×8，派生 I、A 各 N_H×2。原发性状统计按模块坐标解释，等位效应向量植物3维、昆虫8维；QTL和顺序以元数据与输入表为准。旧版 I/A 四维映射只用于其原协议。

元数据保存 chemical_schema、genetics_version、功能协议、性状顺序、观察窗口和生活史。plants.defense_cost=cD(x_N+x_S)、voc_cost=cV*x_V；eta 改变线索但不直接改变成本，旧 B/rho/Z 不参与当前规则。

成年参考 p/q 是基因型反事实诊断；感觉影响母体选择，处理影响幼虫表现。参考 p 不含外部拒绝，参考 q 无资源竞争，p*q 不是成虫本轮真实取食或完整繁殖适合度；全零 q 的宽度未定义。

新批次先 audit，再运行 analysis.summarize_lifecycle；它将 populations、host_flows、host_groups、networks、lagged_networks 合并到 data/summary/，分别写为 lifecycle_counts、host_background_flows、host_background_groups、lifecycle_networks、parent_host_contributions，并添加 case/replicate。summarize_qtl_run复用当前成年/植物快照的遗传和性状统计；旧summarize_coevolution入口已冻结，不在当前流程。

plot_diagnostics 通过阶段适配使用终点成年快照、幼虫资源及成年→合格雌虫→潜在卵→出生幼虫→补充成年账本。A–E要求mutation_off/qtl_mutating对照案例；默认多资源情景使用--sensitivity。单寄主池没有P2时P2_resource_fraction为NA；已配置P2但实际消费为零且总消费为正时才为0。原始CSV的空表与非空表使用相同字段及顺序，未发生或不适用的字段保留NA。历史结果保持原字段、源码哈希和协议，不把旧200代结果或旧测试数当成本版验证。

## 匹配机会与固定参考代 QTL 来源

这些分析表可通过统一绘图入口生成中文 F–H 图（PNG/PDF/SVG）：

以下为命令示例，应替换为实际存在且已审计的批次路径。先生成谱系及可选祖源汇总，再绘图：

```powershell
python -m evo_sim.models.plant_insect.analysis.summarize_lineages --input codex/validation/v0.1.8/lineage_check --design evo_sim/models/plant_insect/analysis/lineages.yaml --ancestry-reference-generation 1
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/v0.1.8/lineage_check --lineages
```

先对该批次运行 `summarize_lineages`，需祖源图时显式指定 `--ancestry-reference-generation`。默认图片位于批次 `figures/lineages/`，可用 `--out` 指定验证目录。F 为匹配交配机会和风险集覆盖，G 为全体成年个体的参考代 QTL 来源份额，H 为母本来源的消费网络和出生/招募计数；只绘制已有数据，不增加模拟。NA 保留为不可估计，不按零隔离绘制。不同运行按母本权重、QTL 拷贝数或事件数汇总；运行内谱系标签不能解释为跨运行的同一谱系。绘图源表和输入指纹保存在 `data/summary/lineages/plot_*.csv` 与 `plot_metadata.json`，原 A–E 绘图命令不变。

`summarize_lineages` 的 `matched_opportunities.csv` 分出生寄主与候选群、兼容池、母本群及父本目标群；`__all_cross__` 表示所有异群父本。局部与同池全局计算使用同一完整分类母本集和相同权重，列出三种 `weight_type`、总数/权重、无配偶、未知候选边、排除与匹配覆盖。全局比例为零或权重和为零时 RI 缺失；负值允许。潜在卵数由局部有配偶者抽取，成功卵权重条件化于实际产卵；两者不能外推所有雌虫。

追加 `--ancestry-reference-generation <代数>` 后，`qtl_ancestry.csv` 以该代寄主或候选群为固定来源，逐条校验并追踪每个位点实际传递的父母拷贝，缺失、重复或不一致的传递记录会报错。`reference` 是 k=0 基线；`birth` 目标为出生寄主（包含死亡幼虫），`recruitment` 目标为当代候选群及其当前操作 lineage ID。分母 `total_qtl_copies` 包含全部来源及 unknown/unassigned，`fraction` 为对应来源份额，`unknown_fraction` 为未知份额。空总体产生 `empty_population` 状态行且比例缺失；不能将后续来源份额直接称为新增基因流。网络仍仅有 maternal_source 和 recruit_cluster，不生成物种确证节点。

## 人口参数敏感性分析

20代人口敏感性配置为 `configs/demographic_sensitivity.yaml`：8种情景各8重复，种子20261102–20261109；比较基线、R=2/4、潜在卵基准8、搜索3、交配3、扩散2和R=2×搜索3。默认`configs/current.yaml`另比较R=10/30/100、每种4重复和5代。两类配置都是理论诊断，配置存在不等于已完成实验。按[当前计划](../README.md#current-work)，运行并审计批次后可绘图（路径须换成实际批次）：

```powershell
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/v0.1.8/run2 --sensitivity
```

该命令在 `data/summary/` 输出 `sensitivity_runs.csv`、`sensitivity_trajectories.csv`、`sensitivity_survival.csv` 和 `sensitivity_diagnostics.csv`，并保存审计与绘图输入指纹；中文轨迹、终点和阶段诊断图默认写入 `figures/sensitivity/`，各有 PNG、PDF、SVG 三种格式。每个情景内的重复是独立运行，跨情景相同种子用于配对，不能把不同情景当成独立重复。

固定期限存续汇总分别列已知存续、已知灭绝与未知；运行失败或未完成后状态未知时报告全体存续比例上下界，仅全部可判定时报告 Wilson 95% 区间。终点个体数分布只包含完成的运行；已经确认灭绝的轨迹保留零，未完成的存活轨迹截止最后观测代。阶段比率的分母为零时为 NA。运行配置、种子或源码指纹与批次记录不符时视为未审计、状态未知，不当作灭绝。
