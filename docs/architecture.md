# evoSim v0.1.7 / v0.1.6 — 实现架构

更新：2026-10-29。SLiM 5.2不变，v0.1.7当前执行slim-qtl-4 / expression-tuning-2，v0.1.6保留slim-qtl-2执行/回读。新协议植物3轴/6QTL、昆虫8轴/16QTL；植物保存三化学性状，昆虫保存模块性状，I/A为派生诊断。`helpers/modules.py`负责功能核和母源表达成本，生命周期恢复原成年空间交配；植物保留原存活、局部同池严格异交与均匀父本抽样。当前没有交配前定居或花期参数、theta状态、wrapped_fraction。旧slim-qtl-3花期草稿和expression-tuning-1定居版本仅历史回读，不允许当前执行。原草稿CSV/原始运行及source/metadata哈希不变。`audit_v016.py`按实际协议独立重建账本，新增字段也须通过独立审计。

`analysis/summarize_lineages.py`只读已审计输出，在模拟外计算操作分群、亲子来源、分裂/合并/终止与时标明确的网络；这些标签不进入父本候选或SLiM。`Batch.open`默认拒绝新生活史，阶段适配入口显式允许且统一绝对路径。实际验收见 [v0.1.7基线修订记录](plan/v0.1.7_baseline_revision_2026-10-29.md)。下面的职责表保留v0.1.6基础；新协议覆盖其固定3/4维、直接I/A与固定参考面板解释。

## 1. 职责与执行路径

| 位置（相对 evo_sim/models/plant_insect） | 职责 |
| --- | --- |
| configs/*.yaml、inputs/*.csv | 生活史、版本、参数、物种/QTL及参考面板；数值为理论示例 |
| helpers/config.py、inputs.py | 配置展开、输入/域/协议检查，保留历史读取所需身份 |
| helpers/agents.py | 植物三模块与派生防御；昆虫成年 I/A、性别、出生寄主、队列及基因型 |
| helpers/chemistry.py | 同一 x→D/V 的确定性映射；历史映射独立按协议读 |
| helpers/voc.py、ecology.py、demography.py、spatial.py | 复用选择/拒绝、瓶颈、资源、成本、空间候选与植物人口学计算 |
| substeps/lifecycle.py | AdvanceGeneration 完成成年交配→逐卵择主→SLiM幼虫→资源/存活→植物更新→成年补充的耦合转换 |
| helpers/genetics.py、slim_backend.py、slim_reproduction.slim | 按协议3/4维或3/8维加性 QTL、原生 addCrossed、标记重组、项目突变及位点传递 |
| helpers/contributions.py | 阶段双亲计数、等位贡献分解；本版生命周期分别调用出生/补充两次 |
| batch.py、simulator.py | 调度、批次、状态、文件及指纹；历史结果不覆盖 |
| audit.py、audit_genetics.py、audit_v016.py | 按版本分派、独立核对遗传/资源/阶段来源/事件映射 |
| analysis/summarize_lifecycle.py | 新阶段原始账本只读汇总；不增加宏观网络指标 |
| analysis/common.py | 历史 Batch 入口明确拒绝 v0.1.6，防止把成年当作本轮摄食者 |

现行 config.yaml 以一个耦合 transition 调用 advance_generation；不再依赖旧 ecology→demography→reproduction 三阶段时序。保留的旧源码/历史观察路径不代表当前运行仍沿用旧生物学流程。

## 2. 状态与身份

植物 x 为 N_P×3，D/V 为 N_P×2，默认基因型 N_P×2×6；当前v0.1.7昆虫原发模块为N_H×8、I/A各N_H×2，默认基因型N_H×2×16（v0.1.6原发I/A共4维、8位点）。当前输入来自`inputs/v017_no_flowering/`，共同CSV植物只填写向量0–2，3–7留空；位点数由CSV给定，不用空列补出额外植物性状。SLiM每类群当前效应矩阵分别3/8维，不强行同宽；每子代12/32份拷贝进入QTL突变核。历史slim-qtl-3的4/8维只在原结果回读中保留。

昆虫新增 birth_host_id、birth_host_species、adult_cohort。内部 unknown=-1，输出标签 unknown；创始未知不能按位置补造。reproduction_round=t 记录事件，offspring.generation 和 adult_cohort=t+1 记录队列。个体 ID 类群内唯一，子代出生至成年不换 ID，跨类群引用含 kind。

oviposition.event_id 使用 round:mother_id:egg_index，成功事件链接 offspring_id；offspring/events/larvae 的 oviposition_event_id 回链。全部出生者有双亲和存活结果，成年快照只保留幸存者；因此 DAG 不会因幼虫死亡丢失出生关系。双亲 ID 与位点传递不等于完整树序列。

## 3. 原子更新与两个贡献终点

Adult(t) 已完成上一轮幼虫存活和从出生寄主的扩散。直接按成年位置、同祖先兼容池、异性、交配半径形成候选，不以本轮产卵寄主作为资格，也不再抽存活。母源处理/感觉表达成本决定潜在卵，每卵从母体成年位置独立择主/拒绝，成功卵逐子代从既定池均匀抽父并生成自身基因型。产卵选择不移动成虫或重建配偶池，未来寄主不能决定此前交配。shuffle_voc只影响逐卵线索。幼虫以自身I固定出生寄主、同批比例资源分配、S=g只一次，幸存后从出生寄主格扩散。植物受食后抽存活，合格母株需有局部同池非自身存活父本，以Poisson(lambda_P Fm phi)产出、逐后代均匀抽父本，经SLiM从母株格扩散，非重叠替换；不增加物候模块。

当前已删除交配前定居函数、开关、metadata.settlement_semantics及专属adult_settlement表，带旧开关配置明确拒绝。parents.mating_row/col与candidate_mate_ids/probabilities仍保存实际位置和候选池，兼用于植物及昆虫审计/观察。只读历史审计按expression-tuning-1重建旧定居，不修改旧数据。

临时转换内同时构建全部出生状态和幸存成年状态；出生、补充分别累计双亲各1/2贡献和等位分解。资源、遗传、存活、空间及所有输出验证后整轮提交。失败不提交半轮实际账本，保护上限停止而不截断人口。终点只 observe 快照与参考诊断，不调用下一轮幼虫转换。

## 4. 事件与观察分层

oviposition 保存全部潜在卵选择，larvae/events 保存实际出生幼虫资源与存活。host_flows 在 birth/recruitment 两阶段按亲本出生源→子代出生目的含 unknown 计账。networks 是祖先兼容池层幼虫资源/成功补充，lagged_networks 是亲本出生背景的实际滞后贡献，host_groups 是可消失再现的出生背景群。全部出生双亲构成 DAG；原始模拟不生成谱系资格或新物种，v0.1.7只读分析另识别操作候选群。

standardized 继续计算成年 A/I 的反事实参考谱；不是成年个体本轮取食，也不回写真实幼虫事件。旧 summarize_qtl_run/summarize_coevolution 命令不接受新生活史；plot_diagnostics 由 prepare_diagnostics 适配终点成年、幼虫资源与新繁殖阶段，复用已有 A–E 图函数和阶段无关的遗传统计。详细列定义见 [results.md](results.md)。

元数据保存实际模型版本、生活史、化学/遗传协议、坐标顺序、观察窗口、展开配置/输入/源码哈希和SLiM身份。历史D/Z、slim-qtl-1、含花期slim-qtl-3草稿和旧阶段含义均按原保存协议读取；不能重写旧指纹或仅凭版本号断言生态等价。植物化学共同演化仍在当前执行路径中，不宣称本版另有植物适应性成种机制。
