# evoSim v0.1.2 — Results：输出内容、指标定义、图件与结果解释

> **版本**：当前文件契约为 `evoSim v0.1.3 / slim-qtl-1`。
> **文件性质**：输出规范与分析计划。O01–O10、F1–F12 保留原始规划身份，不表示所有规划表和图均已实现。
> **目的**：清楚说明“模型将产生什么原始数据、汇总统计和论文图件”，并对接 `architecture.md` 的 Agent/Step 和 `model.md` 的公式。

## 当前实现与本规范的对应（2026-10-25）

默认批次入口输出到 `results/results_v0.1.2/runN/`，不覆盖已有目录。因果检查入口由 `--out` 指定目录，输出相同逐运行文件及 `causal_checks.json/md`，不生成常规批次的 `settings.json`；各运行已经单独审计，该目录不直接用于常规批次汇总命令。

```text
results/results_v0.1.2/runN/
├── settings.json
├── data/
│   ├── run_status.csv、status.json
│   ├── <case>/rep01/
│   │   ├── metadata.json、status.json
│   │   ├── populations.csv、plants.csv、events.csv、parents.csv
│   │   ├── offspring.csv、gamete_transmissions.csv、mutations.csv
│   │   ├── genetic_loci.csv、allele_catalog.csv、genetic_summary.csv
│   │   ├── haplotype_summary.csv、networks.csv、diagnostics.csv
│   │   ├── allele_contributions.csv
│   │   └── snapshots/、genotypes.npz、phenotypes.npz
│   └── summary/                 # 汇总命令生成
├── audits/integrity_audit.json   # 审计命令生成
└── figures/                     # 绘图命令生成
```

当前个体主键为 `(case, replicate, kind, individual_id)`，处理和重复由目录提供。`events.csv` 合并选择与取食事件，性别、位置和性状保存在快照。重组信息在 `gamete_transmissions.csv`；未单独输出 `insects/host_choices/feeding/mating_summary/recombination_intervals/traits_summary/spatial_summary.csv`。验证范围见 [本版检查报告](v0.1.2_causal_check_report.md)。

### v0.1.2 新增字段及缺失语义

| 文件 | 字段 | 定义 |
| --- | --- | --- |
| `parents.csv` | `is_maternal_role` | 植物及雌虫为真；与有无寄主、生存、配偶分开 |
| 同上 | `realized_genetic_contribution` | 母源、父源后代数之和除以 2；未提交时缺失 |
| `parents.csv`、`networks.csv` | `contribution_committed` | 完成代际更新才为真；终代、失败更新不能伪造贡献 |
| `allele_contributions.csv` | `p_parent`、`p_parent_weighted`、`p_offspring` | 每物种/位点的亲代频率、按实际双亲贡献加权频率、实际后代频率；记录亲子等位基因并集 |
| 同上 | `delta_parent_contribution`、`delta_transmission` | 加权频率减亲代频率、后代频率减加权频率；不是纯选择与纯突变分解 |
| `networks.csv` | `resource_share` | 昆虫物种内实际资源份额；总量为零时缺失 |
| 同上 | `maternal_offspring`、`paternal_offspring`、`genetic_contribution` | 按母本、父本各自当代寄主归属；合并两者各乘 1/2 |
| 同上 | `maternal_share`、`paternal_share`、`genetic_contribution_share` | 各自在昆虫物种内归一化；相应总量零时缺失 |
| `diagnostics.csv` | `joint_pq`、`cov_pq`、`variance_p`、`variance_q`、`variance_pq` | 个体先相乘后平均；协方差/方差用总体分母 N；p 不含拒绝，q 不含资源竞争 |
| `snapshots/*.npz` | `reference_p`、`reference_q`、`reference_plant_species` | 矩阵行与 `insect_ids` 对齐，列与参考物种标签对齐；关闭观测时不保存 |

`expected_maternal_offspring` 对所有雌虫和植物有定义：没有寄主的雌虫记零，雄虫记缺失。全体个体的母源期望均值把非母源角色按零计入；全体雌虫均值则不包含雄虫，须报告分母。`W_i/W_insect/W_plant` 仍为 `S×F` 的兼容字段，仅是表现代理量。

没有后代时实际贡献可以为零，但下一代等位频率及分解缺失；在非空后代中消失的等位基因频率才是零。未提交的更新不输出频率分解行。空账本仍保存列头。新增观测不抽随机数，不改变生态和遗传抽样顺序。

以下保留原 v0.1.1 分析规划；尚未实现的独立表和图不属于本版完成声明。

## 1. 结果分类：四级证据，不可混写

| 级别 | 作用 | 何时可以报告 | 典型内容 |
|---|---|---|---|
| **A. 实现正确性验证** | 证明程序按既定公式执行 | 完成 v0.1.1 单元测试/冒烟运行后 | 孟德尔分离、Haldane r(d)、资源守恒、父本/母本一致性、无配偶零出生 |
| **B. 人口学与模型诊断** | 了解规则在随机运行中如何工作、定位异常 | 实际执行相应配置后 | 种群规模、雌雄率、配偶可获得性、剪裁率、生态/遗传数值边界 |
| **C. 生态遗传描述** | 描述演化过程中发生了什么，不宣称机制因果 | 预定重复数和分析条件后 | 等位频率、性状分布、寄主资源流、群落网络与空间格局 |
| **D. 机制假说评价** | 检验遗传与生态机制如何影响宏观结果 | 模型整体建成、事先明确假说、独立校准与设计后 | 寄主专化、寄主转换、群落稳定、网络结构的对照实验 |

**本版不产生 `speciation_events` 或真实分子酶活输出**，因为物种标签是固定的、QTL 无真实序列意义。图片标题必须与可支持的分析级别一致。[R01–R03,R08–R11]

## 2. 结果目录与文件契约（v0.1.1 唯一规范）

建议每次完整运行使用：

```text
results/results_v<version>/
└── run<N>/
    ├── metadata.json
    ├── status.json
    ├── tables/
    │   ├── populations.csv
    │   ├── plants.csv
    │   ├── insects.csv
    │   ├── host_choices.csv
    │   ├── feeding.csv
    │   ├── parents.csv
    │   ├── mating_summary.csv
    │   ├── offspring.csv
    │   ├── gamete_transmissions.csv
    │   ├── recombination_intervals.csv
    │   ├── mutations.csv
    │   ├── allele_catalog.csv
    │   ├── genetic_summary.csv
    │   ├── traits_summary.csv
    │   ├── spatial_summary.csv
    │   └── host_network.csv
    ├── snapshots/g000000.npz ...
    ├── audits/
    │   ├── integrity.json
    │   ├── pedigree.json
    │   ├── resource_conservation.json
    │   └── tests.log
    └── figures/           # 完成审计后从 tables/snapshots 派生，不是原始真相
```

每一行至少含 `run_id, replicate, generation, species/kind`（适用时）；所有 individual ID 在**同一个 run 内唯一且不可复用**，不允许只用数组行号当长期个体 ID。建议 `metadata` 记录完整 YAML、随机种子、软件/科学版本、代码提交哈希、配置哈希、公式版本、参数 `source_kind`。`status.json` 至少有 `completed/failed/terminated`, `last_committed_generation`, `failed_phase`, `exception`, `max_individuals_safety_reached`。

### O01. `populations.csv`：种群规模与灭绝

每行是一个物种 × 一个世代。至少包含：

`generation, kind, species, N_start, N_female, N_male, N_survived_to_mating, N_eligible_mothers, N_mothers_with_mate, N_no_compatible_mate, N_eligible_fathers, N_offspring, N_next, extinction_generation, update_applied`。

- 植物没有雌雄性别，`N_female/N_male` 为 **NA**；可以同时有种子母本与花粉父本贡献。
- 雄虫没有母方生殖计数；`N_offspring=Σ当代母方 B_i`，不能把父源数量再相加造成重复。
- `N_next` 仅在执行了 `t→t+1` 更新的记录上有定义；终点生态评估使用 `update_applied=false`，`N_next=NA`，不能编造下一代。
- `extinction_generation` 以首次提交 `N_{s,t}=0` 的代为准；不存在的物种初始为零时须说明初始缺失。

**对应公式**：D1–D10；**图件**：F1、F2。

### O02. `plants.csv`：植物个体生理与繁殖成本

每行一个植物 × 世代：

`plant_id, species, cell_row, cell_col, D1, D2, Z1, Z2, V1, V2, resource_capacity_R, total_demand_T, gamma, consumption_X, herbivory_h, defense_cost, voc_cost, fecundity_factor_FP, survival_probability_SP, local_density_phi, survived_to_mating, mate_count, born_as_mother, sired_as_father`。

诊断 `FP` 是生殖能力系数，不等于在本代一定产生后代；植物可作为花粉父本给其他母本贡献遗传来源。**对应公式**：E1、E6–E9、E11、S2、D3–D7；**图件**：F3、F5。

### O03. `insects.csv`：昆虫个体寄主选择与收益

每行一个昆虫 × 世代：

`insect_id, species, sex, search_cell, chosen_host_id, host_cell, I1, I2, A1, A2, candidate_count, host_choice_reason, q, demand_d, allocated_x, feeding_gain_g, FI, SH, survived_to_mating, eligible_as_mother, eligible_as_father, candidate_mate_count, born_as_mother, sired_as_father`。

无寄主 `chosen_host_id=NA`、`allocated_x=g=0`，因本版 `SH=g` 不可通过生存检验；不能用预期选择分数替代真实摄食收益。**对应公式**：E2–E11、D1–D7；**图件**：F1、F3、F6。

### O04. `host_choices.csv` 与 `feeding.csv`：机会—选择—真实利用三阶段

- `host_choices.csv`：`generation, insect_id, candidate_plant_ids, candidate_species, candidate_probabilities, outside_probability, chosen_host_id, was_rejected`；每次选择概率之和为1。
- `feeding.csv`：`generation, insect_id, plant_id, host_plant_species, q, d, gamma, x, g, plant_X, plant_R, plant_h`；所有 `x` 不得超过同一植物的 `R`。

**重要**：`选择`与`成功利用`不是同一条生态边。以后画网络时必须声明以候选暴露、接受选择还是实际摄食量 `x` 定义边，不能互相替换。**对应公式**：E2–E8；**图件**：F4。

### O05. `parents.csv` 与 `mating_summary.csv`：交配与母方产子

`parents.csv` 每行一个当代个体 × 世代，包含：

`id, kind, species, sex_or_role, cell, S_i, F_i, phi_i, survived_to_mating, can_be_seed_or_female_parent, can_be_pollen_or_male_parent, potential_father_count, surviving_compatible_father_count, has_mate, conditional_nu, mate_probability_q, maternal_expected_B, maternal_variance_B, offspring_count_as_mother, offspring_count_sired_as_father, failed_mating_reason`。

`mating_summary.csv` 按物种汇总母本数、无配偶数、配偶成功比例、父本贡献分布、同格交配份额、配偶距离均值等。`maternal_expected_B/maternal_variance_B` 按 D8–D9 的**给定交配前状态**定义；雄虫和不具母方角色者为 NA。已发生存活与择偶之后的 `conditional_nu` 应另有记录，不能混用。

**对应公式**：D1–D9；**图件**：F2、F3。可检查配偶缺失对小种群演化的影响，但不能单凭某次灭绝声称出现了普适 Allee 效应。

### O06. `offspring.csv`：真实双亲谱系

每个实际出生的子代一行：

`generation, kind, species, offspring_id, maternal_id, paternal_id, maternal_species, paternal_species, maternal_cell, paternal_cell, offspring_cell, offspring_sex, maternal_gamete_id, paternal_gamete_id, maternal_mutation_count, paternal_mutation_count`。

- `generation` 为**子代所属代 `t+1`**，母方和父方记录在 `t`，必须明确区分。`offspring_id` 不复用，父源贡献 `offspring_count_sired` 用真实子代记录汇总。
- 植物子代 `offspring_sex=NA`，昆虫子代按 `p_female` 抽取稳定性别。
- 父母必须存活至繁殖、同种且符合性别/异交限制；不能记录并不存在的虚构父本。

**对应公式**：G6、D3–D7、S3；**图件**：F7（谱系示例）。

### O07. `gamete_transmissions.csv` 与 `recombination_intervals.csv`：配子与重组

`gamete_transmissions.csv` 每行一个子代 × 一个角色 × 一个 QTL：

`generation, offspring_id, role, parent_id, gamete_id, chromosome_id, locus_id, position_cM, source_homolog, inherited_allele_id, child_allele_id, mutation_event_id`。

`recombination_intervals.csv` 每行是一次配子形成时同一连锁群的一个相邻 QTL 间隔：

`generation, offspring_id, role, chromosome_id, left_locus, right_locus, distance_cM, K_crossovers, parity, homolog_before, homolog_after`。

**注意区分**：`K_crossovers` 是抽到的交换次数；`parity` 是两个标记之间来源是否切换；它们不是同一变量。当前压缩算法只模拟标记区间内的交换次数，不输出从未生成过的实际物理断点坐标。**对应公式**：G2–G6；**图件**：F8（重组率测试）。

### O08. `mutations.csv` 与 `allele_catalog.csv`：突变事件和等位来源

`mutations.csv`：`event_id, generation, offspring_id, role, locus_id, inherited_allele_id, new_allele_id, effect_delta_1..4`。

`allele_catalog.csv`：`allele_id, organism, locus_id, chromosome_id, parent_allele_id, birth_generation, effect_1..4, is_founder`。

每次突变产生新 ID，不覆盖既有等位；所有 `new_allele_id` 必须能找到源等位基因的父节点；突变**是 QTL 效应变化**，不能写成核酸点突变的实测次数。**对应公式**：G7–G8；**图件**：F9。

### O09. `genetic_summary.csv` 与 `traits_summary.csv`：遗传和性状演化

`genetic_summary.csv` 按世代 × 物种 × QTL 输出：

`N, locus_id, n_alleles, allele_frequencies, observed_heterozygosity_Ho, gene_diversity_He, major_allele_frequency, founder_allele_fraction, new_mutation_count`。

定义：对于具有 `N>0` 个二倍体个体的某 QTL，记某等位基因的两个拷贝总计数为 `n_a`：

\[
p_a=\frac{n_a}{2N},\qquad H_E=1-\sum_a p_a^2,\qquad
H_O=\frac{1}{N}\sum_{i=1}^{N}\mathbf 1(a_{i,0}\neq a_{i,1}).\tag{O1}
\]

`N=0` 时 H、p 应记 NA；不能为了绘图填成0。存在多个 QTL 的连锁群才统计真实单倍型频率。`H_E` 是基因多样性而不是自动证明 Hardy–Weinberg 平衡的期望。

`traits_summary.csv` 输出各 D/Z/V/I/A 坐标的 `raw_mean, raw_sd, clipped_mean, clipped_sd, q05, median, q95, clipped_low_fraction, clipped_high_fraction`。必须报告裁剪比例，否则加性效应可能隐藏在0/1平台。**对应公式**：G1/G9；**图件**：F9、F10。

### O10. `spatial_summary.csv`、`host_network.csv`：空间与互作网络

`spatial_summary.csv` 每一代每格记录 `n_P1,n_P2,...,n_H1_female,n_H1_male,...,n_total_plants, plant_species_richness, mixed_plant_cell, n_accepted_hosts`。对每次报告的群落图可以重建“空格、只有P1、只有P2、两者共存”四类，而不依赖手工读图片。

`host_network.csv` 每代每种昆虫 `h` × 每种植物 `p` 记录：

\[
E^{choice}_{hpt}=\sum_{i:\,species(i)=h}\mathbf 1[species(H_i)=p],\qquad
W^{feed}_{hpt}=\sum_{i:\,species(i)=h,\ species(H_i)=p}x_i.\tag{O2}
\]

另记 `candidate_exposure`, `n_accepted`, `n_consumers`, `total_demand`, `W_feed`, `W_gain`。如网络权重用 `W_feed`，单位为**归一化可食生物量**；`E_choice` 是接受计数，二者不得画成同一种边。

将昆虫 `h` 的实际利用份额定义为：

\[
\pi_{hpt}=\frac{W^{feed}_{hpt}}{\sum_p W^{feed}_{hpt}},\qquad
N_{eff,h,t}=\exp\left(-\sum_{p:\pi>0}\pi_{hpt}\ln\pi_{hpt}\right).\tag{O3}
\]

`N_eff` 是摄食权重计算的有效寄主物种数；若该昆虫物种当代**没有实际摄食**，则 `π` 和 `N_eff=NA`，不是“专食程度=0”。当前 P1/P2/H1 的网络规模有限，嵌套性和模块度可能定义不稳；只有在未来多物种网络有足够节点和非零边后才将相应指标设为正式检验结果。

**对应公式**：S1–S3、E4–E8；**图件**：F4、F5、F11。

## 3. 建议论文图件：哪些真正值得画？

| 图号 | 推荐图名 | 横/纵轴与颜色 | 主要数据来源 | 科学解释边界 |
|---|---|---|---|---|
| **F1** | 群落数量与存续 | x=世代；y=各物种 N；颜色=物种 | O01 | 数量变化不自动代表共进化机制得证 |
| **F2** | 生存—交配—繁殖漏斗 | x=世代/物种；y=存活母本、成功获配偶母本、出生数、无配偶率 | O01/O05 | 分解繁殖失败原因；不是假说显著性检验 |
| **F3** | 植物与昆虫资源/适合度 | x=世代；y=食害 `h`、虫获食 `g`、生存率、母本 `ν` | O02/O03/O05 | 区分潜在需求和实际摄食 |
| **F4** | 寄主选择与实际摄食双网络 | 行=昆虫种、列=植物种；颜色=接受数或 `W_feed` | O04/O10 | 选择网络与摄食网络分开面板、共同注明单位 |
| **F5** | 多世代空间共现图 | 指定世代 `0、25、50、100...`；植物四类格颜色，昆虫数黑色圆点 | O10/snapshots | 染色规则预定义；热图面积不等于真实地理面积 |
| **F6** | 气味匹配、选择和利用的分离 | x=气味距离/偏好分数；y=选择概率、q、g 等 | O03/O04 | 避免把相似气味等同于高适口性 |
| **F7** | 真实父母—子代谱系 | 连线=母源/父源传递，附位点示例 | O06/O07 | 显示分离与双亲来源，不是系统发育树 |
| **F8** | 遗传图谱与重组验收 | x=cM；y=观察 r 和 Haldane 理论曲线 | O07 | 属于代码验证，不是生态科学结果 |
| **F9** | 等位频率轨迹与变异累积 | x=世代；y=主要等位频率、突变数、H_O/H_E | O08/O09 | 新等位不是实测 SNP，注意漂变与选择混合 |
| **F10** | D/Z/I/A 与裁剪诊断 | x=世代；y=性状分布及剪裁率 | O09 | 裁剪饱和可能隐藏潜在基因型差异 |
| **F11** | 寄主利用广度 | x=世代/预定义情景；y=`N_eff`、不同宿主资源份额 | O10 | 无摄食 NA，广度受寄主资源可利用性限制 |
| **F12** | 机制参数—结果关系（未来） | x=预注册的 X 类情景；y=主要 Y 指标，重复/不确定性 | 所有审计后汇总表 | **须待模型建成、预先定义假说后**才属于正式机制检验 |

绘图输入参数分两类：**(i) 运行自变量 X**，如 `rho_voc、chromosome_id/position_cM、competition_radius、植物/昆虫扩散半径、初始等位频率`；**(ii) 输出指标 Y**，如网络权重、存活率、性状方差。这两类不能与 `lambda/alpha/mutation_rate` 等待校准参数混为一谈；详情见 `model.md §6`。

## 4. v0.1.1 必须通过的分析与审计指标

### 4.1 遗传与数值正确性：不依赖生态运行趋势

| 检验 | 理论结果 / 条件 | 合格证据 |
|---|---|---|
| 杂合单 QTL 分离 | `Aa` 形成配子 A:a=1:1 | 固定 seed、大样本频率和容差；纯合 AA 必传 A |
| 两个独立连锁群 | 双杂合 `Aa/Bb` 配子 AB/Ab/aB/ab 各1/4 | 不应只出现亲本型单倍型 |
| Haldane 10 cM | `r=(1-exp(-0.2))/2≈0.0906346` | 重组型统计符合理论与抽样误差 |
| 0 cM 区间纯函数 | `K=0`，不切换 | 正式图谱要求位置唯一，0 cM 在纯函数接口验收 |
| 突变0与1的边界 | `μ=0` 无新变异，`μ=1` 每传递拷贝每位点均生新等位 | 所有等位谱系准确、旧等位不变 |
| 生存与交配 | `R=0` 不得当父母；无合格父本 `B=0` | `parents+offspring` 精确审计 |
| 母本出生数 | 合格母本 `E B=Var B=ν`；混合模型按 D8–D9 | 条件/无条件统计分别核验 |
| 资源守恒 | `Σ_i x_i = min(T_p,R_p)` 数值容差内 | `0≤X_p≤R_p`，无顺序特权 |
| 世代替换 | `N_{t+1}=Σ_i B_i`，终点不出额外后代 | 统计及谱系行数完全一致 |
| 重放/原子性 | 同 seed+config 字节稳定；超限没有半代 | 哈希、失败路径和状态检查 |

这些**是目标理论计算应满足的验收值，不是已经观察到的模拟实验结果**。文献 [R02–R04] 提供遗传模拟/重组公式的方法依据。[R01] 要求对模型实现可复现、可审查。

### 4.2 必须在最终报告中披露的风险诊断

1. `clipped_fraction`：是否已有大量个体表型被限制在0/1。
2. `mate_failure_rate`：交配空间或雌雄性别比是否造成出生失败。
3. `female_fraction`：小昆虫种群里是否出现雌雄偏斜乃至单性吸收状态。
4. `acceptance_rate` vs `actual feeding_gain`：气味偏好和资源获取是否被混淆。
5. `resource_balance_error`：所有取食资源流量是否守恒。
6. `mutation_events` 和有效 QTL 效应量：突变数据是否和 `G9` 的表型变化匹配。
7. `aborted_run/safety_limit`：运行停止是自然灭绝、到达预定终点还是计算异常；不能混为一种结局。

## 5. 结果何时可写成论文中的“发现”？

只有在 Agent 实现 v0.1.1、所有正确性验收通过、预先定义研究自变量 X 与响应 Y、充分重复和明确参数来源后，才可以填写“结果”章中的具体数量、趋势和统计推断。在模型开发阶段，小规模模拟只用于检查异常和公式实施，不允许为了某个看起来理想的共存或专化趋势持续调参数。

**未来结果报告表（待运行后填写）**：

| 科学问题 | 预定义实验因素（X） | 主要响应（Y） | 当前状态 |
|---|---|---|---|
| 遗传变化能否影响生态互作？ | 抽象等位效应、初始遗传多样性、重组结构 | QTL频率、D/Z/I/A变化、`W_feed` | **待 v0.1.1 运行** |
| 气味选择与实际利用是否一致？ | VOC耦合、感知容忍度（预先声明） | `E_choice`, `W_feed`, `g` | **待 v0.1.1 运行** |
| 交配限制如何参与种群维持？ | 花粉/昆虫交配半径及分布场景 | 无配偶率、雌雄比例、实际出生数、灭绝概率 | **待 v0.1.1 运行** |
| 空间性状差异与网络如何联系？ | 竞争尺度、扩散及初态结构 | 空间共现、网络权重、`N_eff` | **待 v0.1.1 运行** |
| 新物种能否产生？ | **本版未定义** | **无有效成种指标** | **不在 v0.1.1 研究范围** |

## 6. 输出字段与参考依据的对应关系

| 内容 | 主要依据 | 证据性质 |
|---|---|---|
| 个体模型与可复现输出 | Grimm 等 2020 [R01]，Haller & Messer 2023 [R02] | 论文报告/方法学框架 |
| 配子、重组和双亲传递 | Haller & Messer [R02,R03]，Haldane [R04] | 前向遗传范式 + 经典遗传图谱理论 |
| 可遗传数量性状与区分统计效应 | Chu 等 2024 [R05] | 定量遗传理论；不代表现实酶活 |
| 植物交配对种群的影响 | Barrett & Harder 2017 [R06] | 交配系统多样性，有助合理提出理论规则 |
| 寄主选择与 VOC | McFadden 1974 [R07]、Auer 等 2020 [R10] | 概率选择数学形式 + 嗅觉机制实验背景 |
| 化学性状/底物特异性 | Kliebenstein 2001、Köllner 2004、Li 2004 [R08,R09,R11] | 经验现象背景；具体 D/I 函数仍待实证 |
| 竞争、共存和资源反馈 | Chesson 2000、Beverton & Holt 1957 [R12,R13] | 负密度依赖和共存理论背景，非实测系数 |


---

## 参考文献与理论出处（按文内编号）

**[R01]** Grimm V, Railsback SF, Vincenot CE, et al. 2020. The ODD protocol for describing agent-based and other simulation models: A second update to improve clarity, replication, and structural realism. *Journal of Artificial Societies and Social Simulation*, 23(2):7. https://doi.org/10.18564/jasss.4259  
**适用范围**：IBM 的 ODD 模型描述、方法透明度和可复现性；不是具体生态方程的来源。

**[R02]** Haller BC, Messer PW. 2023. SLiM 4: Multispecies eco-evolutionary modeling. *The American Naturalist*, 201(5):E127–E139. https://doi.org/10.1086/723601 （[开放全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10793872/)）  
**适用范围**：多物种显式遗传、QTL、生态互作、nonWF 风格的前向个体模型。

**[R03]** Haller BC, Messer PW. 2019. SLiM 3: Forward genetic simulations beyond the Wright–Fisher model. *Molecular Biology and Evolution*, 36(3):632–637. https://doi.org/10.1093/molbev/msy228  
**适用范围**：个体出生/死亡/适合度的非 Wright–Fisher 正向模拟。evoSim 不是 SLiM 软件的封装。

**[R04]** Haldane JBS. 1919. The combination of linkage values, and the calculation of distances between loci of linked factors. *Journal of Genetics*, 8:299–309. [Haldane 著作目录](https://jbshaldane.org/bib/JBS-Haldane-Bibliography.html)  
**适用范围**：无交换干涉的重组映射 `r=(1−e^(−2d))/2`；适用条件必须说明。

**[R05]** Chu TT, Kristensen PS, Jensen J. 2024. Simulation of functional additive and non-additive genetic effects using statistical estimates from quantitative genetic models. *Heredity*, 133:33–42. https://doi.org/10.1038/s41437-024-00690-5 。同年勘误：https://doi.org/10.1038/s41437-024-00725-x 。  
**适用范围**：个体功能遗传效应与群体统计遗传效应的区分；不等于分子酶功能模拟。

**[R06]** Barrett SCH, Harder LD. 2017. The ecology of mating and its evolutionary consequences in seed plants. *Annual Review of Ecology, Evolution, and Systematics*, 48:135–157. https://doi.org/10.1146/annurev-ecolsys-110316-023021  
**适用范围**：种子植物交配系统多样性和花粉受限的重要性；不证明严格异交适用于所有植物。

**[R07]** McFadden D. 1974. Conditional logit analysis of qualitative choice behavior. In: Zarembka P, ed. *Frontiers in Econometrics*. Academic Press, pp. 105–142. [作者提供的章节](https://eml.berkeley.edu/reprints/mcfadden/zarembka.pdf)  
**适用范围**：多项 Logit / Softmax 离散选择的数学基础；不证明昆虫行为遵从该具体效用核。

**[R08]** Kliebenstein DJ, Kroymann J, Brown P, et al. 2001. Genetic control of natural variation in Arabidopsis glucosinolate accumulation. *Plant Physiology*, 126:811–825. https://doi.org/10.1104/pp.126.2.811  
**适用范围**：植物防御化学谱具有遗传基础；不提供本项目两个 D 轴的校准值。

**[R09]** Köllner TG, Schnee C, Gershenzon J, Degenhardt J. 2004. The variability of sesquiterpenes emitted from two *Zea mays* cultivars is controlled by allelic variation of two terpene synthase genes encoding stereoselective multiple product enzymes. *The Plant Cell*, 16:1115–1131. https://doi.org/10.1105/tpc.019877  
**适用范围**：等位变异改变多组分 VOC 产物谱；不验证本研究的线性混合系数。

**[R10]** Auer TO, Khallaf MA, Silbering AF, et al. 2020. Olfactory receptor and circuit evolution promote host specialization. *Nature*, 579:402–408. https://doi.org/10.1038/s41586-020-2073-7  
**适用范围**：昆虫嗅觉受体及神经环路参与寄主相关行为；不意味着 A 就是某个 OR 的酶学参数。

**[R11]** Li X, Baudry J, Berenbaum MR, Schuler MA. 2004. Structural and functional divergence of insect CYP6B proteins: From specialist to generalist cytochrome P450. *PNAS*, 101:2939–2944. https://doi.org/10.1073/pnas.0308691101  
**适用范围**：昆虫不同解毒酶具有不同的底物广度与效率；不提供本模型的 `q=min(I/D)` 公式。

**[R12]** Chesson P. 2000. Mechanisms of maintenance of species diversity. *Annual Review of Ecology and Systematics*, 31:343–366. https://doi.org/10.1146/annurev.ecolsys.31.1.343  
**适用范围**：稳定化共存机制和种间/种内竞争区别；不证明本项目竞争矩阵。

**[R13]** Beverton RJH, Holt SJ. 1957. *On the Dynamics of Exploited Fish Populations*. Fisheries Investigations, Series II, Vol. 19. [Beverton–Holt 函数与原始文献说明](https://search.r-project.org/CRAN/refmans/LeMaRns/html/rec_BH.html)  
**适用范围**：负密度依赖的有理函数建模传统；本项目局部个体竞争是从该理论形状借鉴的变体。

**证据声明**：引用文献只对相应的“遗传规律、方法形式或经验机制”负责。`μ=0.01`、`σ_mut=0.02`、`λ=2`、`α=0.02/0.005`、化学效用系数、配偶半径、严格异交与 Poisson 出生的特定参数等，凡未直接从独立数据拟合，都必须标记 `theoretical_placeholder` 或 `implementation_choice`。这些假设不能因为列出参考文献而被写成实测事实。


## v0.1.3 的 SLiM 记录变化

`metadata.json.slim_calls` 保存每次成功提交的类群/代际调用身份、后代数、配子和突变 seed、实际执行文件及 binary/script SHA-256。源码指纹包含实际 `.slim` 脚本。无后代时不启动引擎；失败更新不发布该更新的调用或等位记录。

配子表的 `transmitted_homolog_before_mutation` 由同位置中性来源标记判定，`switch_before_locus` 为相邻 QTL 来源切换。压缩图谱没有完整区间交换数，因此 `crossover_count` 在新协议中缺失；旧协议仍按旧定义审计。突变仍按每传递 QTL 拷贝定义，由 Eidos 执行，`mutation_stream` 改为实际 SLiM seed 与代际标识。

引擎缺失、非零退出、超时或非法返回记录为 `backend_error`。该状态保持现有群体和未提交贡献语义，不自动切换至另一遗传器。历史结果、空群体频率、终代和失败更新的缺失值规则保留。
