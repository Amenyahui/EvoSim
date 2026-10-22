# v1.4 遗传机制的假设与边界

## 从基因型到生态

每个个体有形状 `(2, L)` 的整型等位 ID 数组，植物与昆虫有独立位点集合；同一生物类群的各物种可以共享位点坐标，遗传只来自同种本亲代。效应记录为冻结对象，新突变创建新 ID 和上游 ID，旧等位效应不原位改变。

植物效应坐标依次为 `[D_1...D_n, Z_1...Z_K]`，昆虫为 `[I_1...I_n, A_1...A_K]`。每个生物类群有共同背景向量。表型为 `clip(background + 0.5 × 两拷贝所有位点效应之和, 0, 1)`，没有按物种 ID 决定表型的分支。初始化按显式等位频率独立抽取两拷贝和各位点。QTL 模式的旧 `founder_traits` 仅是固定诊断面板，不覆写任何个体。

每个子代先复制本亲代两份等位数组，再在每个拷贝/位点以该位点 `mu` 独立判定突变。新效应为旧效应加零均值正态增量，当前支持声明靶坐标上的独立增量，即对角协方差；未实现任意相关矩阵。非靶坐标效应和标准差必须为零。所有旧直接性状突变函数在 QTL 模式下停用。

[SLiM 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC10793872/)是显式 QTL 与多物种互作的框架先例；[Chebib & Guillaume 2021 的 Model](https://pmc.ncbi.nlm.nih.gov/articles/PMC8664587/)提供连续等位值上的高斯增量先例。具体半和、截断、稀疏靶坐标、参数和克隆方式为项目约定。两个文献均不能验证本项目的分子映射或数值。OBP 的功能标签也不能笼统称为远距嗅觉：[Matsuo et al. 2007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050118)讨论味觉和寄主产卵偏好。其他原设计文献索引和本次核查范围见 evidence_ledger.csv。

## 开发参数与对照

| 项目 | 开发配置 | 解释 |
|---|---|---|
| 遗传表示 | 每个体两拷贝；每生物类群 8 位点 | 每个 D/Z 或 I/A 坐标两个位点；可改维度和位点数 |
| 背景 | 所有坐标为 0 | 物种差异来自等位组合 |
| 创始效应 | 原参考性状的一半；第二等位增加 0.02 | 每个位点各两种等位，频率 0.5/0.5；全部是人为输入 |
| 初始独立性 | 各位点、两拷贝独立抽样 | 仅初态假设；后续不发生随机交配 |
| 突变率 | 0.01 / 拷贝 / 位点 / 代 | 8 位点两拷贝至少一次突变概率约 0.1485；不能解释为 DNA 突变率 |
| 增量 SD | 靶坐标 0.02；其他坐标 0 | 对角协方差；无实测标定 |
| 关闭突变 | `mutation_enabled` 按植物/昆虫分别设定 | 选择和漂变仍运行；不重置等位频率 |
| 遗传随机流 | `[seed,140,phase,taxon,species]` | 与已有生态/人口 130、空间 150 相位独立 |
| 资源、成本、λ、竞争 | 沿用旧配置 | 没有为了延缓灭绝调整 |

分配保护为等位目录上限 100000、累计保存基因型标量数上限 1000000。更新前保守预留所有允许突变拷贝都突变时的空间，故可能早于真正耗尽目录上限停止；这是计算中止，保留全部已抽出生数但不提交下一代。它们不是生态承载量，也不是对进程所有 Python/NumPy 内存的完整上限。当前仍整批保存在内存中；大规模运行前需流式写出与稀疏汇总。

## 输出和独立检查

- `snapshots/g*.npz`：每代个体 ID、物种、格位、基因型和表型，与旧生态审计共用。
- `genotypes.npz`、`phenotypes.npz`：按 `g000000_字段` 导出相同状态，审计逐数组与快照对照。
- `genetic_loci.csv`、`allele_catalog.csv`、`mutations.csv`：位点靶效应、不可变等位谱系和所有真实突变事件；零效应事件也记录。
- `offspring.csv`：同种本亲代 ID、后代 ID、位点拷贝突变次数、空间来源与派生表型。旧 `*_mutated_j` 在 QTL 模式表示影响该坐标的位点发生突变，包括零步长；不意味着表型必然变化。
- `genetic_summary.csv`：每代/物种/位点等位频率，观察杂合度 `Ho`、基因多样性 `1−Σp²`，性状均值、总体方差及原始值超出 `[0,1]` 的截断比例。等位行的位点统计是重复展示，不能再对所有等位行累加突变数。物种灭绝后频率和性状统计为 NA，数量为 0。
- `metadata.json`：配置、源码哈希、遗传架构哈希、独立随机流，以及分开的软件/生态/遗传版本。

审计器不调用运行时 G→P 或遗传函数；从保存的目录、事件、亲子表与基因型重新算表型、重建继承和统计量，再调用原生态审计重算资源、适合度、繁殖和空间。更新输出 SHA 后的人为效应、突变来源和统计量篡改仍会失败。审计是内部一致性与实现检查，不是模型已被自然数据验证的证明。

## 解释限制和后续顺序

二倍体整体克隆意味着完整谱系连锁，无孟德尔分离、重组、交配或有性基因流。`Ho` 与 `1−Σp²` 分开报告，不假定克隆群体处在随机交配平衡。位点效应相加且没有显性/上位性，并不代表实际 OR、P450 或代谢网络机制。

后续依次处理：生态持续性与参数敏感性；明确目标生活史后引入双亲遗传/重组；有经验功能数据后校准 G→P；再扩展多物种网络试验。当前保持固定物种 ID，不能根据性状差异或基因多样性宣称产生新物种。

## 参考文献索引

以下沿用原设计的 R 编号，供 evidence_ledger.csv 对照；当前实现边界以上文为准。



- **[R01] Haller BC & Messer PW (2019).** *SLiM 3: Forward Genetic Simulations Beyond the Wright–Fisher Model*. **Molecular Biology and Evolution**. https://doi.org/10.1093/molbev/msy228 。支持遗传显式前向模拟、非 WF 过程和可定制生活史；不支持具体植物—昆虫生化函数。
- **[R02] Haller BC & Messer PW (2023).** *SLiM 4: Multispecies Eco-Evolutionary Modeling*. **The American Naturalist**. https://doi.org/10.1086/723601 。在共同演化实例（A Coevolutionary Example）直接给出多物种气味匹配、加性 QTL、正态效应及基因型依赖的相互作用；不支持当前 D/Z/A/I 的数字效应。
- **[R03] Haller BC, Ralph PL & Messer PW (2026).** *SLiM 5: Eco-evolutionary Simulations Across Multiple Chromosomes and Full Genomes*. **Molecular Biology and Evolution**. https://doi.org/10.1093/molbev/msaf313 。说明多染色体基因组层建模可扩展，**不要求**v1.4 模拟全基因组。
- **[R04] Chebib J & Guillaume F (2021).** *Pleiotropy or linkage? Their relative contributions to the genetic correlation of quantitative traits and detection by multitrait GWA studies*. **Genetics**. https://doi.org/10.1093/genetics/iyab159 。Methods 明确：双倍体、加性等位值、正态/多元正态突变、连续等位效应增量和多效性—连锁比较；本项目最重要的逐式来源。
- **[R05] Barton NH, Etheridge AM & Véber A (2017).** *The infinitesimal model: Definition, derivation, and implications*. **Theoretical Population Biology**. https://doi.org/10.1016/j.tpb.2017.06.001 。支持多位点定量遗传理论背景；其无穷小极限不等于有限少数大效应基因。
- **[R06] Speed MP, Fenton A, Jones MG, Ruxton GD & Brockhurst MA (2015).** *Coevolution can explain defensive secondary metabolite diversity in plants*. **New Phytologist**. https://doi.org/10.1111/nph.13560 。原模型的毒素—抗性最小瓶颈、成本以及双向选择；不包含我们的显式基因型或VOC映射。
- **[R07] Kliebenstein DJ et al. (2001).** *Genetic Control of Natural Variation in Arabidopsis Glucosinolate Accumulation*. **Plant Physiology**. https://doi.org/10.1104/pp.126.2.811 。植物防御代谢组分受遗传位点影响的原始研究；不规定通用加性形式。
- **[R08] Köllner TG et al. (2004).** *The Variability of Sesquiterpenes Emitted from Two Zea mays Cultivars Is Controlled by Allelic Variation of Two Terpene Synthase Genes Encoding Stereoselective Multiple Product Enzymes*. **The Plant Cell**. https://doi.org/10.1105/tpc.019877 。TPS4/TPS5 等位变异改变挥发性萜烯产物谱；不证明 VOC 与 PSM 线性相关。
- **[R09] Auer TO et al. (2020).** *Olfactory receptor and circuit evolution promote host specialization*. **Nature**. https://doi.org/10.1038/s41586-020-2073-7 。Or22a 等位替换与寄主气味识别、行为因果证据；不提供我们的二维 A 或 Gaussian 选择效用数值。
- **[R10] Matsuo T et al. (2007).** *Odorant-Binding Proteins OBP57d and OBP57e Affect Taste Perception and Host-Plant Preference in Drosophila sechellia*. **PLOS Biology**. https://doi.org/10.1371/journal.pbio.0050118 。OBP 相关基因影响味觉/寄主产卵行为；不等于所有 OBP 只作用于远距嗅觉。
- **[R11] Weng JK, Philippe RN & Noel JP (2012).** *The rise of chemodiversity in plants*. **Science**. https://doi.org/10.1126/science.1217411 。植物特殊代谢、酶多功能性与化学多样化的机制综述；不提供具体 PSM 位点参数。
- **[R12] Li X, Baudry J, Berenbaum MR & Schuler MA (2004).** *Structural and functional divergence of insect CYP6B proteins: From specialist to generalist cytochrome P450*. **PNAS**. https://doi.org/10.1073/pnas.0308691101 。CYP6B1 与 CYP6B8 底物代谢范围有功能表达实验证据；不支持逐毒素一基因的普适映射。
- **[R13] Calla B et al. (2020).** *Substrate-specificity of cytochrome P450-mediated detoxification as an evolutionary strategy for specialization on furanocoumarin-containing hostplants: CYP6AE89 in parsnip webworms*. **Insect Molecular Biology**. https://doi.org/10.1111/imb.12612 。CYP6AE89 的底物特异性与酶活证据；不等于直接测到本项目 I 值。
- **[R14] Moghe GD & Last RL (2015).** *Something Old, Something New: Conserved Enzymes and the Evolution of Novelty in Plant Specialized Metabolism*. **Plant Physiology**. https://doi.org/10.1104/pp.15.00994 。复制与功能分化促进新代谢能力的综述；v1.4 只作为未来扩展理论依据。
- **[R15] Grimm V et al. (2020).** *The ODD Protocol for Describing Agent-Based and Other Simulation Models: A Second Update to Improve Clarity, Replication, and Structural Realism*. **JASSS**. https://doi.org/10.18564/jasss.4259 。模型假设、过程和可复现说明规范；并非生物学公式来源。
