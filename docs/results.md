# 输出、审计与绘图

当前协议为模型0.1.12 / slim-qtl-6 / expression-tuning-contact-3 / dose-tolerance-lifecycle-2。
原始文件名及个体/事件/阶段键保留，但下述列的科学意义已经修改；
旧0.1.11结果须使用冻结源码解释，不能以新读者重算。
[模型](model.md)、[参数](parameters.md)和[科学修改设计](ecology_v0.1.12_design.md)给出定义。

## 运行、审计和汇总

```powershell
python -m evo_sim --config evo_sim/models/plant_insect/configs/current.yaml
python -m evo_sim.models.plant_insect.audit --input <batch> --require-success --min-transitions 1
python -m evo_sim.models.plant_insect.analysis.summarize_lifecycle --input <batch> --sensitivity
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input <batch>
```

开发输出默认创建全新的 `codex/validation/v0.1.12/run<N>/`，不覆盖已有路径。
正式科研实验仅在用户明确安排时写results/。原始记录由审计独立重建后用于分析，
汇总位于data/summary，不修改metadata或原始文件。完成执行、达到期限、有限期限存续和
稳定共存是不同判断；计算保护中止之后保持unknown，不能填零当作灭绝。

## 新生态记录

| 文件与字段 | v0.1.12意义 |
| --- | --- |
| snapshots.plant_chemical_traits | 五列PSM1、PSM2、Z、VOC1、VOC2 |
| snapshots.insect_module_traits | 十列，末两列是[0,1]接触坐标；签名权重2*raw-1 |
| oviposition.contact_plant_ids | 每卵实际接触植物的有序JSON列表，最多候选数次 |
| contact_probabilities / contact_accepted | 对应接触接受概率及Bernoulli结果，拒绝后继续剩余候选 |
| larvae/events.d、x、g | 实际食物需求dmax、分配量、x/dmax；I不参与需求计算 |
| dose_1/2 | x*PSM*dose_scale，抽象负荷 |
| processing_capacity_1/2 | I*processing_scale，摄入后处理能力 |
| processed_1/2、residual_1/2 | min(dose,capacity)、dose-processed |
| antinutrient_1/2 | 本版等于全部摄入剂量，独立于处理能力 |
| nutrition_survival / toxicity_survival / antinutrition_survival | 三项无量纲存活分量 |
| W_insect | 三项存活分量的乘积，不再等于g |
| q / reference_q / q_unlimited | 足量食物、无竞争时的化学存活测定，不再是min(I/D)或需求乘子 |
| plants.R、X、resource_after_feeding | 取食前资源、实际消费、取食后资源 |
| defoliation | X/R，R=0时0；不再是1-post/K |
| tolerance_Z、recovery_rate | 遗传耐受坐标及基础恢复率+Z系数*Z |
| resource_after_recovery / resource_next | 幸存旧株恢复后资源，死亡株为NA |
| W_plant、parents.F_i | 基于恢复后资源的植物繁殖代理/母源乘子 |

outside嗅觉外部选项可在接触前终止搜索；forced仍执行接触接受。
contact_plant_ids与候选列表不同，不把未接触植物写成接触拒绝。
拒绝卵无出生子代ID，不计算其虚构幼虫表现。

阶段键维持：亲本/幼虫资源轮t，出生遗传记录generation=t+1，
reproduction_round=t、adult_cohort=t+1。连接事件需case/replicate及个体类群/ID，
两类群ID可能相同。植物offspring只包含建立成功新株；留存旧株另计。

参考p为等机会固定面板中第一次嗅觉接近概率乘接触接受概率，
不是完整重试后的最终产卵分布；不强制母虫偏好与后代表现匹配。
谱系、亲本贡献、预设兼容池与生殖隔离分析规则不变，标签不等于新物种。

## 诊断输出与解释

`sensitivity_runs/trajectories/diagnostics/survival.csv`保留已有结果分类。
幼虫预期招募为sum(W_insect)，实际为存活计数。无幼虫时比例为NA。
资源耗尽率仅在R>0植物中计算，受攻击耗尽率再限制herbivore_load>0；分母单列，
耗尽容差1e-12*R。未启动、审计失败和计算停止保持区别。
plant_start_deficit_mean加plant_consumption_capacity_fraction_mean等于
plant_post_resource_deficit_mean；它不再等于实际defoliation。

`trait_variance_partition.csv`以各类型个体数加权，将原发坐标及已记录D/I的总体方差
拆成类型内和类型间方差，ddof=0、单位为相应性状单位平方。
不把这些描述性方差当作加性遗传方差、遗传力或跨物种进化速度。
输入哈希和解释写入summarize_qtl_run_metadata.json。

A–E绘图入口仍针对P1/P2及mutation_off/qtl_mutating对照；B图改用独立PSM/VOC轴、
加入接触坐标，不能把第三植物轴Z画成VOC。F–H谱系分析与--sensitivity入口保留。
新protocol不支持直接读取0.1.11历史输出。

## 本版验证

2026-11-09 已执行：

- `python -m ruff check evo_sim tests example.py`：通过。
- `python -m ruff format --check evo_sim tests example.py`：114个文件通过。
- `python -m tests`：[完整日志](../codex/validation/v0.1.12/tests_complete.log)116组通过。
  新增的 `tests/scientific/test_ecology_v012.py` 覆盖八项科学验收、接触重试边界、
  固定种子重复和重新散列后的记录篡改拒绝；实际调用SLiM验证新增坐标遗传、重组与突变。
- [最终小规模核验](../codex/validation/v0.1.12/run2/verification.json)：
  基线、关闭突变、无虫和全接触拒绝四情景，各2重复、3轮，共8次运行全部完成。
  [独立审计](../codex/validation/v0.1.12/run2/audits/integrity_audit.json)
  通过24次转移、32个状态、104条幼虫事件、490条植物事件和690条新生记录，
  核对118次突变；资源超用及相对收支误差均为0，`source_matches_current=true`。
  重复运行的基因型与快照一致，接触次数有界、处理能力不改变食物需求，
  分析前后原始记录哈希不变。
- [隔离安装包核对](../codex/validation/v0.1.12/package_final/verification.json)：
  wheel中的99项Python、YAML、CSV与SLiM资源逐项匹配最终源码；
  隔离安装后版本及5/10维配置校验通过。SLiM后端与核心脚本的SHA256与基线相同。

以上均为开发实现核验，未执行本版200/300轮实验、完整参数不确定性分析或经验校准。
新增系数仍是理论假设；抗营养负荷本版等于摄入剂量，不假定处理模块消除该效应。
不根据短程存续声称长期共存、共同适应或物种形成。

<a id="current-validation"></a>

## 历史0.1.11 实施验收

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
