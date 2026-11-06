# 输出、审计与绘图

当前协议为模型0.1.11 / `slim-qtl-5` / `expression-tuning-2` / `resource-coupled-plants-larval-recruitment-1`。模型公式见[模型](model.md)，默认值见[参数](parameters.md)。本版实施验收见[当前记录](#current-validation)；下面的0.1.9和0.1.10记录均为历史结果，不验证本版。历史结果保持原文件、元数据及源码哈希，使用[冻结包](../codex/history/before_simplification_0.1.11_20261106T034842Z/README.md)读取。

## 目录与审计

开发、测试和验收默认写入 `codex/validation/v0.1.11/run<N>/`，显式输出目录必须不存在。仅用户明确安排的正式科研实验写入 `results/`。运行后以相同目录执行：

```powershell
python -m evo_sim --config evo_sim/models/plant_insect/configs/current.yaml --out codex/validation/resource_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/resource_example --require-success --min-transitions 1
```

完整性审计通过不等于达到请求期限或共同存续；计算保护停止不等于灭绝。派生汇总保存于批次 `data/summary/`，原始数据不回写。

## 阶段与输出解释

每轮成年昆虫交配、逐卵选择植物，出生幼虫取食及存活；植物按食后资源抽存活，幸存母株每轮产生提案并恢复，空位建立的新株加入下一状态。`plants.csv`记录轮初资源、消费、亏损、存活及恢复；`plant_proposals.csv`记录提案与建立；`offspring.csv`中植物仅包含建立新株，旧株留存另计。终点无下一轮事件表示未执行，不表示真实零事件。

`B_post=B_start−X`，`defoliation=1−B_post/K`；不使用同龄无食害参照。新株资源为0.25K。恢复和新株资源是显式输入，死亡株残余移出；资源账不代表全株碳守恒。每格位点限制建立，不表示共享光水养分。

植物与昆虫均由SLiM进行双亲遗传，母父贡献各半；`gamete_transmissions`记录双方实际位点传递、重组及启用时的突变。植物父本为全局其他幸存株，初始类型不限制交配；标签表示母系初始来源，不能替代个体化学性状。出生与成年补充的贡献分别计数。运行元数据记录协议和遗传语义。植物类型是初始化学标签，不是生殖隔离或物种确证。

按植物标签汇总 `allele_contributions` 时，`p_parent` 是母系来源标签的参照群；`p_parent_weighted` 包含实际双亲，可以跨标签。因此 `delta_parent_contribution` 包含父源混合和繁殖贡献，不能单独解释为纯自然选择；没有跨标签交配时才适用原来的同标签解释。关闭突变时没有新等位变异，但选择和重组仍可改变群体。十种初始来源不维持为十个离散化学型，后代连续性状由重组及启用时突变重新生成。

`summarize_lifecycle`只读汇总阶段账本，`summarize_qtl_run`汇总性状与遗传，`summarize_lineages`及`plot_diagnostics`提供只读谱系与诊断；具体字段和当前批次适用性以实际运行元数据和审计为准。参考面板、候选谱系和网络不改变模拟行为；固定期限存续不证明长期共存或共同演化。

## 批次身份、阶段键与原始表

批次根目录保存 `settings.json`、`data/run_status.csv`；各运行位于 `data/<case>/repNN/`，审计位于 `audits/`。`case/replicate` 标识运行，两类群个体ID不保证全局唯一，连接个体时还需 `kind/individual_id`。`reproduction_round=t` 是本轮繁殖，`adult_cohort=t+1` 是子代队列；不能仅按同名generation连接不同阶段：

| 表 | generation | reproduction_round | adult_cohort |
| --- | --- | --- | --- |
| offspring、mutations、gamete_transmissions | t+1，出生队列 | t | t+1 |
| larvae、events、networks | t，幼虫资源轮 | t | t+1 |
| parents | t，亲本繁殖轮 | t | 不适用 |
| snapshots、genetic_summary、host_groups | t，轮初状态 | 不适用 | 快照保存昆虫队列 |

幼虫与昆虫子代以 `(case, replicate, reproduction_round, adult_cohort, insect_id/offspring_id)` 连接；位点传递另用 `kind`、`role`、`copy` 和 `locus_id`。`metadata.stage_key_schema`记录阶段键。`oviposition.event_id=round:mother_id:egg_index`，拒绝卵无offspring ID；`offspring/events/larvae`可通过`oviposition_event_id`和`egg_index`回链。植物提案`event_id=round:mother_id:propagule_index`，仅建立提案有位点及新生ID。

| 文件 | 含义 |
| --- | --- |
| metadata.json、status.json | 展开配置、版本、协议、源码/输入/输出哈希、SLiM身份和停止状态 |
| snapshots/gNNNNNN.npz、genotypes.npz、phenotypes.npz | 轮初及终点的植物和成年昆虫状态、位置、基因型、原发和派生性状；植物N×3化学坐标、昆虫N×8模块坐标 |
| populations.csv | 轮初及下一状态数量；植物N_retained为幸存旧株，N_propagules为建立前提案，N_offspring为新建株；昆虫N_offspring为出生幼虫、N_next为成年补充 |
| parents.csv、offspring.csv | 亲本资格、候选池、母源率和出生/补充贡献；昆虫记录全部出生幼虫，植物只记录已建立新株，不把留存旧株重记出生 |
| oviposition.csv、larvae.csv、events.csv | 潜在卵选择与拒绝；实际幼虫的出生寄主、需求、摄入、资源利用和一次存活 |
| plants.csv、plant_proposals.csv | 轮初植物资源、消费、亏损、存活、恢复及位点；每个扩散提案的目标和建立结果 |
| genetic_loci.csv、allele_catalog.csv、mutations.csv | QTL、等位效应与全部出生者突变，包括未招募的幼虫 |
| gamete_transmissions.csv | 每位点每份出生拷贝的双亲来源、标记切换与突变关联，不是完整树序列 |
| genetic_summary.csv、haplotype_summary.csv、allele_contributions.csv | 轮初遗传统计、单倍型统计及按出生/补充分阶段的贡献分解 |
| host_groups.csv、host_flows.csv、networks.csv、lagged_networks.csv | 昆虫出生寄主背景、双亲来源、实际消费/出生/补充及亲本出生背景贡献；unknown单列 |
| diagnostics.csv | 成年基因型对参考植物的反事实p/q，不是当轮真实取食或完整适合度 |

终点没有下一轮事件行表示未执行，不是该轮摄食或繁殖为零。`birth_host_id/species`是昆虫子代实际出生寄主，`maternal_birth_host_species`和`paternal_birth_host_species`分别是双亲出生背景；unknown不能当同背景。昆虫出生格和扩散后成年格可不同，植物旧株保留原ID和位置。

## 植物资源、父本及贡献

`plants.R=B_start`，`resource_capacity=K`，`resource_after_feeding=B_start−X`，`defoliation=1−B_post/K`；当轮消费比例`h=X/B_start`在B_start=0时定义为0。`survival_probability=exp(−μ₀−μD defoliation)`，`survived`为本轮唯一植物存活抽样；幸存旧株恢复至`K−(K−B_post)exp(−r)`，死株恢复与下一资源为NA。新株初始资源为0.25K；没有同龄无食害参照资源。资源账检查轮初=消费+食后，下一状态加上幸存恢复及新株资源并移出死亡残余；它不是全株碳收支。

`parents.F_i=max(0,1−defense_cost−voc_cost)×B_post/K`是植物母源产出乘子，`conditional_mu=λP F_i`。每轮在其他幸存植物中全局均匀抽取父本，不按距离或初始类型限制，孤株不能异交；零母源乘数的幸存植物仍可作父本。`mate_probability`是根据轮初潜在父本存活率计算的至少一位父本概率；`expected_propagules`和`propagule_variance`是存活及配偶机会条件下建立前提案的无条件矩，不是已建立后代数的矩。`propagule_count`是实际提案数，空位抽签后才产生新株。

`offspring_count_as_mother`和`offspring_count_sired`分别统计出生母/父次数，`realized_genetic_contribution`是二者之和的一半；`recruited_count_as_mother`、`recruited_count_sired`及`realized_recruitment_contribution`对应补充阶段。植物两阶段都只计已建立新株，旧株留存另计；昆虫出生含死亡幼虫，补充只计幸存成虫。全体植物母、父计数各等于已建立株数，但按母系来源标签分组时父本可以来自其他标签，不能要求每标签母父计数相等。零后代的归一化频率为NA。

`gamete_transmissions.copy=0/1`分别对应母源/父源，`role`为maternal/paternal；`parent_id`指实际亲本，`transmitted_homolog_before_mutation`是该位点亲本同源拷贝，`switch_before_locus`标记相邻位点同源来源变化，`inherited_allele_id`及`child_allele_id`区分突变前后。`crossover_count`不可由标记切换数重建；逐位点记录用于核验实际来源，不等同全基因组祖源。

## 汇总和绘图

先对已有批次审计；以下路径是示例，须替换为实际已运行的批次：

```powershell
python -m evo_sim.models.plant_insect.audit --input codex/validation/resource_example --require-success --min-transitions 1
python -m evo_sim.models.plant_insect.analysis.summarize_lifecycle --input codex/validation/resource_example
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input codex/validation/resource_example
python -m evo_sim.models.plant_insect.analysis.summarize_lineages --input codex/validation/resource_example --design evo_sim/models/plant_insect/analysis/lineages.yaml --ancestry-reference-generation 1
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/resource_example --lineages
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/resource_example --sensitivity
```

`summarize_lifecycle`写 `data/summary/` 下的`lifecycle_counts`、`plant_resource_cycles`、`plant_establishment_proposals`、`host_background_flows/groups`、`lifecycle_networks`和`parent_host_contributions`，保留case/replicate；`summarize_qtl_run`写种群、性状、遗传时间序列及运行汇总，只作描述。谱系汇总写`data/summary/lineages/`：候选群和跨代事件、交配及机会、网络和账本核对均为只读；指定参考代后另输出逐QTL拷贝祖源。unknown/unassigned保留，观察群不是确证物种，参考代来源份额不是新增基因流。

A–E默认绘图命令为 `python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input <batch>`，分别展示人口、性状、寄主利用与繁殖、空间、遗传；它要求完整同期限、至少4轮且案例名为`mutation_off`和`qtl_mutating`，空间示例重复可用`--spatial-reps`指定。图只显示P1、P2及H1，不能代表十来源完整互作；默认`current.yaml`是`ten_types`/`ten_types_no_insects`，不满足A–E对照案例要求。F–H使用`--lineages`，先生成谱系汇总并在需要祖源图时指定参考代，分别显示匹配机会、参考代QTL来源及母本来源互作网络。`--sensitivity`为多情景轨迹、终点、阶段与植物资源图，并在`data/summary/`写敏感性表；图可用`--out`指定目录，输出PNG/PDF/SVG。零分母和不可估计量保留NA，不以零代替。分类标签是初始母系来源而非保持不变的化学型；任何P1/P2子集图都不能推断其余八类互作。

完整性审计可验证失败账本：`integrity_status=passed`不代表运行成功。`execution_success`和`acceptance_passed`分别反映执行状态及额外验收条件；`--require-success --min-transitions 1`要求成功状态且至少一次已提交转移，并不证明达到固定期限。计算保护停止后的状态为未知，不当作灭绝。绘图和汇总检查输入指纹；不回写原始metadata。

<a id="current-validation"></a>

## 0.1.11 实施验收

[最终测试日志](../codex/validation/v0.1.11/tests_release_final.log)109组assert测试通过，包含全局跨初始标签双亲遗传、重组与突变、资源耦合繁殖、留存基因型不变及篡改拒绝；Ruff检查和格式检查覆盖74个Python文件。植物的`W_plant`诊断与亲本`W_i`均为`S×max(0,1−C)×B_post/K`，不把恢复后的资源用于本轮母源产出。

[最终默认批次](../codex/validation/v0.1.11/run3/data/run_status.csv)为有虫/无虫各4重复、每次5轮；[严格独立审计](../codex/validation/v0.1.11/run3/audits/integrity_audit.json)通过8次运行、40次转移、7744条幼虫记录和11081条植物记录，资源超用及相对收支误差均为0，`source_matches_current=true`。[机制核对](../codex/validation/v0.1.11/run3/validation_mechanisms.json)确认留存植株基因型逐轮不变、繁殖乘子与食后资源一致，并观察到1094个跨初始来源标签的已建立植物后代。[终点表](../codex/validation/v0.1.11/run3/validation_endpoints.csv)保留全部重复；这些短程结果不验证长期共存或互相适应，也未运行200轮正式实验。

[植物资源图](../codex/validation/v0.1.11/run3/figures/sensitivity/sensitivity_resource_plants.png)及其他敏感性诊断已生成PNG/PDF/SVG并检查显示。每坐标的[性状均值逐轮变化](../codex/validation/v0.1.11/run3/data/summary/trait_change_per_round.csv)按同类群相邻快照均值之差计算，单位为归一化性状/轮；[定义和输入指纹](../codex/validation/v0.1.11/run3/data/summary/trait_change_metadata.json)随表保存。变化含选择、传递、漂变和留存组成变化，不当作相互适应证明；不把3维植物和8维昆虫合成一个未经标准化的速度排名。默认该对照关闭新突变，突变启用机制由独立遗传测试检查。

[安装包核对](../codex/validation/v0.1.11/package_release/verification.json)检查了全部21项YAML、CSV和SLiM资源，并在隔离安装目录用默认配置完成零轮初始化。简化前104项源码、配置、输入及测试已冻结并核对SHA256，历史结果与元数据未改写。所有本次产物均位于`codex/validation/`，属于开发验证。

<a id="ten-plant-types-results"></a>

## 历史0.1.10十类植物—昆虫互作

以下均为历史0.1.10，与现行全局植物异交和每轮繁殖规则不同。十类历史情景的原配置和运行入口保存在冻结包中，须用冻结源码读取，不能以当前同名配置重跑并视为原协议。

正式科研实验须用户明确安排并显式使用`--out results/...`。计算保护停止或未完成的后续状态为未知，不填为灭绝。已有A–E图针对P1/P2，不可视为十类完整互作结果。

互作按重复和轮次读取`networks.csv`的`plant_species`、`birth_count`、`Q_actual`、`successful_recruits`，分别表示各类植物上的出生幼虫、实际消费量和招募成虫数；用`populations.csv`及终点快照检查各植物池数量，保持全部零、缺失和停止状态。卵选择与拒绝读取`oviposition.csv`；拒绝卵没有被选中的植物类型，不能强行归类。候选出现次数可由该文件的`candidate_plant_ids`与同轮`plants.csv`中的ID/类型匹配，分母为0时为NA。出生或消费占比受丰度、空间相遇、拒绝与竞争共同影响，不能单称为遗传偏好；招募/出生是幼虫条件存活，不是完整成虫适合度。

**本次开发验证。** [30代批次](../codex/validation/v0.1.10/run4/data/run_status.csv)的有虫与无虫各4重复均完成30代；[独立审计](../codex/validation/v0.1.10/run4/audits/integrity_audit.json)覆盖8次运行、240次转移和122560条幼虫事件，资源超用及相对守恒误差均为0。108组测试通过，新配置和四份CSV已检查进入安装包。没有运行正式科研实验。

有虫4重复在第30代的植物总数依次为268、255、265、261株，昆虫为844、825、751、796只；无虫对照4次植物总数均达300株的位点上限。有虫重复1的P10在第27代归零：第24轮繁殖脉冲时只剩1株且没有同池配偶。这个例子不证明所有类型或更长时间都能存续。

[互作概览图](../codex/validation/v0.1.10/ten_types_overview_20261105/interaction_overview.png)和[有虫／无虫种群轨迹](../codex/validation/v0.1.10/ten_types_final_30gen_20261105/population_timeseries.png)保留各重复；[逐类型、轮次和重复的互作表](../codex/validation/v0.1.10/ten_types_final_30gen_20261105/interactions.csv)及[输入指纹](../codex/validation/v0.1.10/ten_types_final_30gen_20261105/metadata.json)随图保存。图中P2和P6的名义初始坐标分别为(.2,.3,.2)和(.6,.3,.2)，只有第一防御轴中心不同；P2上的幼虫补充/出生为74.8%–83.1%，P6为53.7%–59.0%（各4重复的范围）。这只是探索性观测，植物丰度、空间相遇和同株竞争仍可能不同，不是独立控制的因果效应或实测物种结论。

本次只读绘图脚本位于本地验证目录，可对已完成且已审计的同合同批次使用；输出目录须不存在：

```powershell
python codex/validation/ten_plants_visualization_20261105/plot_interactions.py --input codex/validation/v0.1.10/run4 --out codex/validation/v0.1.10/ten_types_new_plots
python codex/validation/ten_plants_visualization_20261105/overview.py --input codex/validation/v0.1.10/ten_types_new_plots
```

<a id="feedback-validation"></a>

## 历史0.1.10 同龄参照验证记录

[完整测试日志](../codex/validation/feedback_implementation_0.1.10/full_tests_after_guards.log)107组通过；Ruff检查与格式检查覆盖72个Python文件。最终wheel已核对包含当前YAML、CSV及SLiM脚本，并在隔离安装目录以默认配置完成零代创始快照检查。所有输出属于开发验证，不是正式科研实验。

[20代实施验收](../codex/validation/v0.1.10/run1/audits/integrity_audit.json)覆盖10个预设情景、每情景2个配对种子，共20次运行和400次转移，全部完成并通过严格审计；无虫对照损伤为0。[敏感性汇总](../codex/validation/v0.1.10/run1/data/summary/sensitivity_runs.csv)仅用于检查20代的阶段与分类，不能当作200代敏感性结果。

[K10基线200代审计](../codex/validation/v0.1.10/run2/audits/integrity_audit.json)显示4次均完成、800次转移通过、最大资源超用和相对守恒误差为0。各运行在200代的植物/昆虫数分别为20/67、28/121、10/36、26/54；[固定期限汇总](../codex/validation/v0.1.10/run2/data/summary/sensitivity_runs.csv)4次均为`coexistence_at_horizon`。这是有限期限人口学结果，不保证长期稳定，也不证明适应或共同演化。预设10情景×4重复×200代完整敏感性设计尚未执行。

为避免重复进行昂贵的逐条遗传重放，上述200代汇总复用已经完成的独立批次审计，并在汇总前重新核验每个运行的全部已记录输出SHA256；具体审计来源和哈希记录在[汇总元数据](../codex/validation/v0.1.10/run2/data/summary/sensitivity_audit.json)。两批模拟结束后修复了当前审计器和敏感性读取器对旧协议失败/空批次的拒绝边界，[源码身份核对](../codex/validation/feedback_implementation_0.1.10/execution_identity_check.json)证实仅这两个读取文件与运行时指纹不同；当时原文及SHA256保存在[读取器副本](../codex/validation/reader_guard_before_fix/SHA256SUMS.txt)。执行生态、遗传和生命周期文件均匹配原指纹。200代批次审计报告的`source_matches_current=false`据实保留；短程批次审计时该字段为true，随后同样发生了上述两个读取器变化。原始metadata均未回写。

<a id="perennial-validation"></a>

## 历史0.1.9多年生协议验证记录

以下是历史0.1.9开发验收，均位于`codex/validation/`，不是正式科研实验；须用[冻结读取包](../codex/history/before_feedback_0.1.10_20261105T095448Z/READ_FROZEN_SOURCE.txt)读取，不能用当前协议重解释旧defoliation。[最终测试记录](../codex/validation/perennial_design/full_tests_final.log)101组通过，71个Python文件的Ruff检查及[最终安装包隔离校验](../codex/validation/perennial_design/package_verification_final.json)通过；包包含60项文件、43个模块和6个配置，wheel SHA256为`f25a9a6ccec8bde9feed45ce5ca0f6748ead66a96c30bc7cdcdadb792918127b`。[短程对照](../codex/validation/v0.1.9/population_checks/validation_endpoints.csv)共6次10代、严格审计60次转移，资源超用和相对守恒误差均为0。[无虫植物对照](../codex/validation/v0.1.9/plants_only_200/data/run_status.csv)两次均完成200代、期末96/97株，严格审计共400次转移通过；这不是昆虫与植物共同存续的证据。

[主200代请求的终点表](../codex/validation/v0.1.9/run1/validation_endpoints.csv)包含全部12个重复：K=10四次昆虫在64/96/166/87代灭绝，前三次植物继续到200代（期末1/96/1株），第四次87代双方全灭；K=30四次在112/119/157/83代全灭；K=100在37/14/13/71代计算截断，37/71为累计基因型保存量限制、14/13为潜在卵数限制。后四次200代结局未知，不并入已知灭绝；不存在已验证的双方200代共同存续重复。

[主批次完整审计](../codex/validation/v0.1.9/run1/audits/integrity_audit.json)已核对全部12次运行、1293次转移、394312条幼虫事件和396590条新生记录，`integrity_status=passed`，资源超用及相对守恒误差均为0。但4次计算截断使`execution_success=false`；另有5次全灭终止，全部运行至少200次转移的要求未满足，故`acceptance_passed=false`。审计中的8次成功终止包含全灭，不表示8次完成200代，更不表示共同存续。

审计读者补上了生产分配器原已存在的舍入预算纠偏重放，未放宽整体容差或改变生态计算。`source_matches_current=false`只来自`audit_lifecycle.py`读取器修复；[执行身份核对](../codex/validation/perennial_design/execution_identity_check.json)确认20次运行的执行源码等于[冻结执行包](../codex/history/perennial_run_0.1.9_20261105/READ_FROZEN_SOURCE.txt)，生态、遗传和空间核未改变。原运行metadata、输出及源码哈希未回写，冻结执行包和最终读取器包的身份分别保留；不得通过改哈希掩盖差异。

只读生成的[200代人口轨迹](../codex/validation/v0.1.9/run1/population_200gen.png)与[短程/无虫对照图](../codex/validation/v0.1.9/population_checks/population_checks.png)保留实际停止代；计算截断后的曲线不补零、不外推。当时的结论是0.1.9实现验收通过、长期共存人口学基线未合格，不能将这些轨迹当成成种实验支持；后续设计见[当前计划](../README.md#current-work)。
