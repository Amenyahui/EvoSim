# evoSim v0.1.5 — 当前输出与解释边界

更新：2026-10-28。当前数据为 shared-chemical-1 / slim-qtl-2。本文描述实际文件与新旧身份；本轮不进行演化结果分析。

## 1. 目录和原始账本

新批次默认 results/results_v0.1.5/runN/，显式目标目录必须不存在。批次根目录有settings.json，data/run_status.csv，逐运行data/<case>/repNN/，审计在audits/，派生汇总在data/summary/。计算停止与生物灭绝分开记录。

| 文件 | 内容及关键语义 |
| --- | --- |
| metadata.json、status.json | 完整展开配置、输入/源码/架构哈希、协议与SLiM身份、实际停止状态 |
| populations.csv | 每物种每代起始、生存、出生、下一代及更新状态 |
| events.csv、plants.csv | 昆虫候选/选择/摄入与植物损伤、资源、成本及代理 |
| parents.csv | S、Fm的旧存储名F_i、W_i=SFm、母源均值/矩、角色资格、母源和父源计数 |
| offspring.csv | 有效后代的真实双亲、物种、位置与性别 |
| snapshots/gNNNNNN.npz | 原发x、派生D/V、I/A、基因型、位置、性别及启用时的参考反应 |
| genotypes.npz、phenotypes.npz | 逐代基因型/表型汇集，供独立审计；新表型包含D和plant_chemical_traits，不含独立Z |
| genetic_loci.csv、allele_catalog.csv、mutations.csv | 位点、类群特异长度的等位效应及实际突变 |
| gamete_transmissions.csv | 双亲来源、同源拷贝与标记切换；完整交换次数缺失 |
| genetic_summary.csv、haplotype_summary.csv | 原发性状与遗传多样性/单倍型统计 |
| allele_contributions.csv | 已提交的等位频率及贡献/传递分解 |
| networks.csv、diagnostics.csv | 实际资源、寄主相关双亲贡献，固定参考面板的联合p/q统计 |

主键含case/replicate，个体再含kind/individual_id，逐代观察再含generation。同名0.1.5不代表数据结构相同，必须同时检查协议。

## 2. 三模块状态、统计及元数据

新快照plant_chemical_traits每行依次为x_N,x_S,x_V，形状N_P×3；D及V_emit/V_cue形状N_P×2。I/A各N_H×2。D/V是同一原发遗传状态的确定性输出，不能解释为独立遗传变异源。新状态不保存plant_voc_independent_trait。

genetic_summary的植物trait名为plant_chemical_traits_1/2/3，依次对应x_N/x_S/x_V；昆虫仍为I_1/I_2/A_1/A_2。均值、方差、raw_mean、clipped_fraction针对原发轴，不把派生D重算为额外QTL轴。等位目录effect JSON向量植物3维、昆虫4维；不能再把整个目录强制转为同宽矩阵。

元数据必须记录：chemical_schema=shared-chemical-1，genetics_version=slim-qtl-2，plant_trait_order=[x_N,x_S,x_V]，plant_defense_order=[D_N,D_S]，odor_trait_order=[V_S,V_V]，chemical_observation_window及eta_shared的配置值。resolved_B为null，不能称仍使用某个B映射。

plants.csv中的defense_cost=c_D(x_N+x_S)、voc_cost=c_V*x_V；后者是历史字段名，当前指独立挥发模块投入，不是V的排放量收费。eta变化可以改变V但不直接改变成本。Fm=0只使母源有效产出为零，供源资格仍可成立。

## 3. 代理和实际双亲贡献

F_i/fp/fi存储母源Fm；W_i/W_plant/W_insect存储SFm，不代表完整适合度。zero_fitness_count只统计代理零值，不代表所有生殖功能丧失。

conditional_mu=λFmφ是逐个体共有函数的代数诊断；实际只有合格母源角色进入抽样，雄虫不是产卵者。母源期望对全部植物/雌虫定义，无寄主雌虫为零，雄虫为NA。eligible_as_mother不保证Fm正、存在父本或产生子代；eligible_as_donor不保证被任何母本选择。

已提交时Bm/Bf分别为实际母源/父源计数，w=(Bm+Bf)/2；同物种ΣBm=ΣBf=Σw=N_next。Fm=0、SFm=0仍可Bf>0，这与母源限定定义一致。未提交终代/失败更新的实际贡献是NA，不是零。无后代时计数可为零，但归一化频率/贡献未定义。不能再把SFm乘到计数上。

## 4. 保留的观察与绘图入口

固定参考面板现为每物种一组三模块x，生成D/V后计算条件参考p及无竞争q；p不含外部拒绝概率。individual p×q先逐虫相乘再平均，不能用两个平均值的乘积替代。逐个体有效寄主宽度仍与群体平均反应宽度不同，q全零宽度未定义且保留绝对利用值。

summarize_coevolution仍输出cross_assay、contribution_differentials、host_relations、host_edges、host_breadth CSV和分析元数据。新植物贡献差异按三个原发轴计算；历史数据按其D/Z规则读取。单株单虫跨时间试验中的SFm不含配偶/父源和密度，因此不能等同完整繁殖成功。关系首次出现、持续窗口与因果寄主转换不自动等价。

plot_diagnostics仍是唯一绘图入口；新图使用共享/独立挥发模块x_S/x_V及防御x_N/x_S的正确坐标，历史图仍识别旧性状名。本轮不运行长期模拟或产出趋势分析图。

## 5. 历史兼容与本轮验收

历史slim-qtl-1结果保留原等位、D/Z、成本与指纹；不能重命名Z或重算V后冒充新模拟。独立审计分别核对协议和元数据。新/旧数据即使同为0.1.5，必须区分遗传和化学schema。

本轮仅执行解析边界、短程资源/繁殖账本及真实SLiM遗传检查；正式测试结果记录在 [实施核对](plan/v0.1.5_chemical_shared_function.md)。历史v0.1.4的40组测试和200代重分析不能当作本轮新化学模型验证。没有重新运行200代或按演化方向调参。

当前验收（2026-10-28）：48组测试全部通过，静态及格式检查通过（66个Python文件）；历史slim-qtl-1两代记录只读审计通过。完整记录和未决假设见 [实施核对](plan/v0.1.5_chemical_shared_function.md)。本轮没有新增长期模拟或结果分析。
