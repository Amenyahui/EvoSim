# evoSim v0.1.4 — 输出内容与解释边界

更新：2026-10-27。当前遗传协议 `slim-qtl-1`。本文只把实际写出的文件列为输出契约；旧版拟议表名不再当作现有文件。

## 1. 结果目录

默认新批次写入 `results/results_v0.1.4/runN/`；已有目录不会覆盖。历史 v0.1.3 结果保持原路径和元数据。

```text
<batch>/
  settings.json
  data/
    run_status.csv
    <case>/repNN/
      metadata.json
      status.json
      populations.csv
      events.csv
      plants.csv
      parents.csv
      offspring.csv
      networks.csv
      diagnostics.csv
      allele_contributions.csv
      genetic_summary.csv
      genetic_loci.csv
      allele_catalog.csv
      mutations.csv
      gamete_transmissions.csv
      haplotype_summary.csv
      genotypes.npz
      phenotypes.npz
      snapshots/g000000.npz ...
    summary/                   # 分析命令生成
```

审计与绘图采用既有批次工具的目录规则。因果检查目录没有常规批次 `settings.json`，不能直接当作普通批次汇总。

## 2. 原始账本及主键

| 文件 | 内容与关键边界 |
| --- | --- |
| `populations.csv` | 每物种每代起始数量、生存、子代、下一代、更新状态与灭绝；安全停止单独记录 |
| `events.csv` | 昆虫候选、选择/拒绝、寄主、q、需求、实际摄入 x/g 和资源状态；不是独立的 host_choices/feeding 两张表 |
| `plants.csv` | 植物损伤 h、资源转移、D/V 成本及生态表现 |
| `parents.csv` | S、F、母源条件均值/期望、配偶机会、实际母源和父源贡献、是否提交 |
| `offspring.csv` | 每个后代的真实双亲身份、类群和物种 |
| `gamete_transmissions.csv` | 逐后代拷贝的来源、等位传递、标记间来源切换；`crossover_count` 缺失 |
| `allele_catalog.csv`、`mutations.csv` | 等位来源、效应以及实际突变事件 |
| `genetic_summary.csv` | trait 与 allele 两类记录、频率、多样性和表型统计 |
| `allele_contributions.csv` | 已提交代际更新的等位频率贡献与传递分解 |
| `networks.csv` | H×P 的接受次数、摄入资源、寄主相关母源/父源/遗传贡献；不是植物后代网络 |
| `diagnostics.csv` | 固定祖先植物面板上的种群平均 p、q、逐虫 p×q 及方差/协方差 |
| `snapshots/*.npz` | 按 ID 对齐的双方性状、基因型、性别、位置，以及启用时的个体 reference_p/reference_q |

跨运行主键包括 `(case, replicate)`；个体还需 `(kind, individual_id)`，同一个体的观察需 generation。`realized_genetic_contribution=(母源子代数+父源子代数)/2`。母源期望对雄虫缺失，对无寄主雌虫为零。实际贡献只在更新提交后有定义；终点未繁殖时为 NA。

## 3. v0.1.4 新增功能分析

先完成该批次的常规审计，再运行：

```powershell
python -m evo_sim.models.plant_insect.analysis.summarize_coevolution --input <batch>
```

命令仅接收 `--input`；以下文件写入 `<batch>/data/summary/`。重复执行覆盖这些派生汇总，原始数据不变。分析定义见 [参数](parameters.md)，数学定义见 [模型](model.md)。

| 文件 | 行单位与主要指标 |
| --- | --- |
| `cross_assay.csv` | case/replicate × 植物代 × 昆虫代 × P×H；双方样本数、pair_count、available，acceptance、q、resource、insect_survival、plant_damage、plant_survival、insect_SF、plant_SF |
| `contribution_differentials.csv` | 每代/物种/性状坐标；均值、方差、实际平均贡献、贡献加权亲代均值、realized_differential、下一代均值及 transmission_residual |
| `host_relations.csv` | 每代 P×H；连接原始网络、固定面板诊断、双方丰度及昆虫实际接受率/无候选率。接受率是 H 整体量，不能对 P 行求和 |
| `host_edges.csv` | 每代 P×H × resource/contribution 层；边权、是否有观测、是否超过阈值、连续代数、是否满足窗口、截至本代的首次正边时间 |
| `host_breadth.csv` | 每代 H × reference_p/reference_q；有响应个体数、逐个体有效寄主数的均值和方差 |
| `coevolution_metadata.json` | 分析设计及哈希、各运行实际端点和量纲、输入哈希、当前分析源码哈希、原审计和运行状态文件哈希 |

默认早期为 0，后期为最后实际记录代；若两者相同，仅有一个唯一组合。未完成运行使用实际可用终点并保留状态，不解释为完成预定实验。指定不存在的代数会报错，不自动换代。

交叉试验是独立的单株—单虫接受试验，采用该运行的行为参数、固有 V 及统一指定的 R/dmax；不执行空间遇见、气味随机打乱、资源竞争或交配。SF 仅为生态表现代理，不是完整适合度或实际繁殖成功。不同处理若行为参数不同，其比较还包含这些参数差异。

## 4. 专化、形成、转换和维持的操作含义

固定面板 p 是条件于选择面板内一个寄主的相对概率，不包含外部拒绝。q 是无竞争利用上限。`host_breadth` 对每只虫的 p，或归一化 q，计算 exp(Shannon entropy)，再平均；它不等于先平均反应再算宽度。全零 q 个体的宽度未定义，单列正响应个体数，避免把不能取食称作专化。归一化 q 的宽度不表示绝对利用能力，须与 q 一起解释。

资源边默认 Q_actual>0；贡献边默认 genetic_contribution>0。默认连续三个**相邻且有观测的代**才满足持续窗口。阈值与三代窗口是可修改的操作定义，不是文献证实的生物常数。缺代打断连续计数；未提交贡献既不是有边也不是无边。首次正边仅是首次观测利用/贡献，不自动宣称关系已经形成。

`host_relations` 联合展示功能反应、丰度和实际利用，供识别变化的多种来源；本版不自动把占优寄主变化归类为遗传性寄主转换。丰度、原寄主灭绝、空间机会和功能变化可共同影响实际网络，不能用其中一项单独归因。

有限窗口的边持续性不等于扰动后恢复能力或永久稳定。既有 `analyze_macro` 可输出资源网络周转等指标；2×2 场景适合机制检查，不能支持复杂网络结构的一般结论。

## 5. 图件与验证层级

既有 `plot_diagnostics` 是唯一绘图入口，按其要求的两个处理和完整批次输出 A–E 图；新增分析当前只输出数据。无需为了得到预期曲线改动生态参数。

解析和集成检查验证计算规则，而不是证明共同演化。历史 200 代数据可用新分析重读，但须明确“v0.1.3 模拟、v0.1.4 分析”，不可重命名为重新完成的 v0.1.4 演化实验。历史原始文件与源码身份不改写。

双方性状变化、贡献差异或跨时间表现差异均不足以独立排除漂变及其他选择路径。关闭新突变仍保留既有变异的选择、漂变和重组；当前没有冻结功能分布的四组演化控制。物种标签固定，不能从网络分群声称出现新物种。

## 6. 理论出处与适用范围

- Price (1970), [10.1038/227520a0](https://doi.org/10.1038/227520a0)：贡献协方差与传递项的区分；本模型具体双亲计数是独立定义。
- Forister (2005), [10.1111/j.0014-3820.2005.tb01051.x](https://doi.org/10.1111/j.0014-3820.2005.tb01051.x)：偏好与表现可独立遗传；不能视为自主取食生活史的直接验证。
- Blüthgen (2010), [10.1016/j.baae.2010.01.001](https://doi.org/10.1016/j.baae.2010.01.001)：网络指标应结合丰度和机会解释。
- Maliet et al. (2020), [10.1111/ele.13592](https://doi.org/10.1111/ele.13592)：由个体生态—演化过程产生二分网络的理论先例，不保证本项目必然产生某种拓扑。

## 7. 本轮实际验证记录（2026-10-27）

完整断言测试 40 组通过，格式和静态检查通过。新增检查覆盖单寄主接受概率、资源上限、双方功能差异、个体联合分布、贡献分解恒等式、缺代/未提交边及原始数据被修改后的拒绝行为。

使用 v0.1.3 已保存的 `results/results_v0.1.3/slim_200gen_20261026` 六次 200 代运行完成新分析，核对 1,242 份实际输入哈希。新表行数：cross_assay 48，contribution_differentials 12,960，host_relations 2,412，host_edges 4,824，host_breadth 2,412。全部六次的分析终点均为 200；12 行交叉组合因昆虫灭绝而不可用，60 行未提交终代的性状贡献差异缺失。

上述是“v0.1.3 模拟、v0.1.4 分析”的验证记录；本轮没有重新运行 200 代演化，也没有据此判定共同演化或寄主转换成立。派生文件及元数据位于原批次 `data/summary/`；这些本地结果遵循项目规则不纳入 Git。
