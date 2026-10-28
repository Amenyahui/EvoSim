# evoSim v0.1.6 — 当前实现架构

更新：2026-10-28。shared-chemical-1 / slim-qtl-2 / SLiM 5.2；实际检查状态见 [实施记录](plan/v0.1.6_implementation.md)。

## 1. 职责与执行路径

| 位置（相对 evo_sim/models/plant_insect） | 职责 |
| --- | --- |
| configs/*.yaml、inputs/*.csv | 生活史、版本、参数、物种/QTL及参考面板；数值为理论示例 |
| helpers/config.py、inputs.py | 配置展开、输入/域/协议检查，保留历史读取所需身份 |
| helpers/agents.py | 植物三模块与派生防御；昆虫成年 I/A、性别、出生寄主、队列及基因型 |
| helpers/chemistry.py | 同一 x→D/V 的确定性映射；历史映射独立按协议读 |
| helpers/voc.py、ecology.py、demography.py、spatial.py | 复用选择/拒绝、瓶颈、资源、成本、空间候选与植物人口学计算 |
| substeps/lifecycle.py | AdvanceGeneration 完成成年交配→逐卵择主→SLiM幼虫→资源/存活→植物更新→成年补充的耦合转换 |
| helpers/genetics.py、slim_backend.py、slim_reproduction.slim | 3/4维加性 QTL、原生 addCrossed、标记重组、项目突变及位点传递 |
| helpers/contributions.py | 阶段双亲计数、等位贡献分解；本版生命周期分别调用出生/补充两次 |
| batch.py、simulator.py | 调度、批次、状态、文件及指纹；历史结果不覆盖 |
| audit.py、audit_genetics.py、audit_v016.py | 按版本分派、独立核对遗传/资源/阶段来源/事件映射 |
| analysis/summarize_lifecycle.py | 新阶段原始账本只读汇总；不增加宏观网络指标 |
| analysis/common.py | 历史 Batch 入口明确拒绝 v0.1.6，防止把成年当作本轮摄食者 |

现行 config.yaml 以一个耦合 transition 调用 advance_generation；不再依赖旧 ecology→demography→reproduction 三阶段时序。保留的旧源码/历史观察路径不代表当前运行仍沿用旧生物学流程。

## 2. 状态与身份

植物 x 为 N_P×3，D/V 为 N_P×2，默认基因型 N_P×2×6；昆虫 I/A 各 N_H×2、默认基因型 N_H×2×8。位点数由 CSV 给定，当前化学构型固定三原发/两功能轴。SLiM 每类群效应矩阵分别3/4维，全局目录异宽不能强行当同宽；另类群未使用 ID 的零占位不构成植物第四性状。

昆虫新增 birth_host_id、birth_host_species、adult_cohort。内部 unknown=-1，输出标签 unknown；创始未知不能按位置补造。reproduction_round=t 记录事件，offspring.generation 和 adult_cohort=t+1 记录队列。个体 ID 类群内唯一，子代出生至成年不换 ID，跨类群引用含 kind。

oviposition.event_id 使用 round:mother_id:egg_index，成功事件链接 offspring_id；offspring/events/larvae 的 oviposition_event_id 回链。全部出生者有双亲和存活结果，成年快照只保留幸存者；因此 DAG 不会因幼虫死亡丢失出生关系。双亲 ID 与位点传递不等于完整树序列。

## 3. 原子更新与两个贡献终点

Adult(t) 已完成上一轮幼虫存活和扩散，候选只按当前成年位置、同种、雌雄、交配半径；不需本轮寄主、不再抽存活。母 Fm 决定潜在卵，每卵独立择主/拒绝，成功卵逐子代选父并生成自身基因型。幼虫固定出生寄主、同批比例资源分配、S=g 只一次，幸存后从出生寄主格扩散。植物受该批食害后按已有成本/存活/密度/异交生成后代。

临时转换内同时构建全部出生状态和幸存成年状态；出生、补充分别累计双亲各1/2贡献和等位分解。资源、遗传、存活、空间及所有输出验证后整轮提交。失败不提交半轮实际账本，保护上限停止而不截断人口。终点只 observe 快照与参考诊断，不调用下一轮幼虫转换。

## 4. 事件与观察分层

oviposition 保存全部潜在卵选择，larvae/events 保存实际出生幼虫资源与存活。host_flows 在 birth/recruitment 两阶段按亲本出生源→子代出生目的含 unknown 计账。networks 是固定物种层幼虫资源/成功补充，lagged_networks 是亲本出生背景的实际滞后贡献，host_groups 是可消失再现的出生背景群。全部出生双亲构成 DAG，不生成候选独立谱系或自动新物种。

standardized 继续计算成年 A/I 的反事实参考谱；不是成年个体本轮取食，也不回写真实幼虫事件。旧 summarize_qtl_run/summarize_coevolution 及依赖旧 Batch 的绘图不接受新生活史；新入口只合并明确阶段的账本。详细列定义见 [results.md](results.md)。

元数据保存0.1.6身份、生活史、化学/遗传协议、坐标顺序、观察窗口、展开配置/输入/源码哈希和 SLiM 身份。历史 D/Z、slim-qtl-1 和旧0.1.5阶段含义仍按原保存协议读取；不能重写旧指纹或仅凭名称断言生态等价。
