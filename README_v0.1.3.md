# evoSim v0.1.3：在现有 IBM 中接入 SLiM 遗传引擎

**版本：v0.1.3，实际 SLiM 遗传引擎接入。** 更新日期：2026-10-26。起点为 `v0.1.2`（`abb0a63`）。遗传执行协议为 `slim-qtl-1`。本文件替代此前“保留自研遗传、SLiM 仅作对照”的计划；目前采用用户确认的实际接入路线。

## 1. 本版改什么

**在现有 evoSim 中使用真正的 SLiM 5.2。** 每次代际更新，evoSim 把当前亲本基因型和已抽出的双亲名单交给 SLiM，SLiM 用原生染色体和 `addCrossed()` 产生后代；后代返回 evoSim 后继续执行性状计算、扩散及下一代生态互作。

| 职责 | 实际执行方 |
| --- | --- |
| YAML/CSV 输入、创始等位频率抽样 | Python evoSim |
| 空间相遇、寄主接受、取食、资源竞争及存活 | Python evoSim |
| 配偶资格、双亲名单及每母本后代数 | Python evoSim |
| 双亲配子、孟德尔分离、原生染色体重组及后代构建 | **SLiM 原生引擎** |
| 每 QTL 拷贝 Bernoulli 突变、亲本中心正态效应增量 | **在 SLiM 中执行的 Eidos 模型脚本** |
| 后代效应核对、加性性状映射、扩散及账本提交 | Python evoSim |
| 审计、汇总、寄主关系与网络分析 | Python evoSim |

本版删去生产路径中自写的 NumPy 配子/Poisson 重组和突变抽样。缺少 SLiM、引擎退出、超时或返回材料不合法时，记录 `backend_error` 并停止该次更新，不自动改用旧遗传实现。

原有寄主选择、取食、生存、母源出生均值、父本资格、竞争和扩散规则保留，以便辨别遗传执行器更换的影响。v0.1.2 文档里关于父本 `F>0`、母株邻域竞争阶段名称的候选修正仍是未决理论问题，**这次迁移没有同时实施它们**。

## 2. 为什么这样接入

课题仍然是“遗传变异 → 功能性状 → 寄主互作与空间接触 → 生存和双亲贡献 → 遗传组成及网络变化”。成熟遗传引擎负责基础遗传操作，研究工作集中于性状和生态机制。SLiM 支持可配置 QTL、非 WF 过程和多物种生态—演化，适合这条路线。[R1](#r1)、[R2](#r2)、[R3](#r3)

本版保留 Python 调度及输入输出。其代价是每个有后代的类群每代启动一个 SLiM 进程并交换数值材料；当前少量 QTL 的模型可以承受。未来只有实际性能测量表明必要时，才考虑常驻进程或更广泛迁移生态执行。常驻状态需要额外设计失败回滚，当前不增加。

SLiM 不会替我们决定“哪种突变效应合理”“寄主选择如何发生”或“什么时候形成新物种”。这些仍是模型假设，必须有独立理论论证。感知遗传变化参与寄主专化有实证依据，个体过程产生互作网络有模型先例，但它们不直接验证当前全部函数。[R13](#r13)、[R14](#r14)

## 3. 使用与配置

需要 Python ≥3.10 和 **SLiM 5.2 命令行程序**。Python 包不包含 SLiM 二进制。可使用 SLiM 官方安装方式，或在已有 conda 环境中安装：

```powershell
conda install --channel conda-forge slim=5.2
python -m pip install -e ".[dev]"
python example.py
```

公共 YAML 的新增执行设置：

```yaml
parameters:
  model_version: 0.1.3
  genetics:
    schema_version: slim-qtl-1
    recombination_model: slim_haldane_marker_map
    slim:
      executable: slim
      timeout_seconds: 120
```

`executable` 可以填完整路径。也可用 `EVOSIM_SLIM` 环境变量覆盖，方便不同机器指定各自安装位置；实际路径、版本和二进制指纹保存到运行元数据。科学参数继续来自 YAML 和六张 CSV，不内嵌在 Eidos 脚本中。

本机已在被 Git 忽略的 `codex/slim_runtime/` 中安装独立运行环境：`codex/slim_runtime/Library/bin/slim.exe`。默认命令 `slim` 在 PATH 不可用时，开发工作区可找到此安装；其他机器应使用 PATH、YAML 或环境变量配置，不能依赖该本地目录随仓库分发。

现有场景仍可运行：

```powershell
python -m evo_sim --config evo_sim/models/plant_insect/configs/v012_base.yaml --out codex/validation/v013_slim_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/v013_slim_example
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input codex/validation/v013_slim_example
```

`v012_*.yaml` 名称保留其场景来源，展开后的运行版本现为 0.1.3，遗传执行为 SLiM。它们不会重现旧版逐轨迹结果。要重现 v0.1.2，应检出旧 Git 标签并使用当时的环境；历史结果和元数据不修改。

## 4. 遗传规则及理论边界

### 4.1 固定 QTL 与等位身份

仍使用每个体 `(2,L)` 的定相等位 ID。植物和昆虫各自有独立的抽象连锁群及图谱；默认各 8 个 QTL。位点数及功能分配是理论输入，不是实测染色体或已识别的 OR/CYP 等基因。

SLiM 按 `chromosome_id` 建立原生二倍体染色体，导入当前等位状态。共同等位 ID 在多个拷贝出现时共享同一 Mutation 对象。QTL 使用中性直接效应、last-mutation 替换策略、保留固定效应；完整效应向量单独储存，不再作为独立选择系数施压。[R1](#r1)

创始频率独立抽样仍由 Python 完成，作为明确的初始条件导入 SLiM。有限样本可能存在关联；不能把输入频率相同当成相位或实际创始基因型相同。

### 4.2 遗传距离转换与原生重组

对同群相邻 QTL：

\[
d=(x_{j+1}-x_j)/100,\qquad r(d)=\frac{1-\exp(-2d)}2.
\]

x 为 cM，d 为 Morgan。10 cM 对应 r≈0.0906346。不同连锁群对应不同原生染色体，自由组合。每条染色体的 QTL 位置压缩为连续整数，标记边界重组概率设为 r，由 SLiM 原生重组执行。[R1](#r1)、[R6](#r6)

这是标记传递分布的映射。它保留当前功能位点间无干扰模型的遗传关系，**没有保留原图谱区间中的完整交换次数和物理断点位置**。日志中 `switch_before_locus` 表示标记来源是否改变；`crossover_count` 记缺失，不填 0/1 冒充原 Poisson 交换数。每个 QTL 同位置另有中性来源标记，因此纯合亲本也能记录来源同源拷贝。

### 4.3 突变发生与效应

本版保留已声明的 QTL 规则，在 SLiM 的 Eidos 中执行：

\[
B_{\ell h}\sim\operatorname{Bernoulli}(\mu_\ell),\quad
\boldsymbol\delta\sim\mathcal N(\boldsymbol m_\ell,\operatorname{diag}(\boldsymbol\sigma_\ell^2)),\quad
\boldsymbol a'=\boldsymbol a+B_{\ell h}\boldsymbol\delta.
\]

每传递拷贝、位点、代最多一次事件；μ 是抽象 QTL 概率，不能当作每碱基突变率。正态增量以旧效应为中心，而不是重新抽取一个与旧效应无关的绝对值。新等位 ID 储存完整新效应并替换当前拷贝，历史目录保留旧 ID。[R4](#r4)、[R5](#r5)

SLiM 自动核苷酸/基因组突变发生器在本桥接脚本中关闭，以免叠加第二种事件过程。这里使用成熟引擎和自定义模型脚本，不声称所有规则都是 SLiM 默认设置。

同一类群本代先完成全部原生配子抽样，再切换独立 seed 执行突变。给定相同亲本和配对名单时，开关突变不改变本次调用的配子来源。跨代生态状态变化仍可使后续亲本名单和轨迹不同。

### 4.4 性状映射与选择

保留：

\[
\boldsymbol z^{raw}=\boldsymbol b+s_e\sum_{\ell,h}\boldsymbol a_{\ell h},\qquad
\boldsymbol z=\operatorname{clip}(\boldsymbol z^{raw},L,U).
\]

尺度 s_e=0.5 是效应计量约定，不是显性系数。突变只经此映射进入 D/Z/I/A，再经生态生存和双亲后代形成选择；没有为同一效应增加另一项直接遗传适合度。裁剪使完整映射非线性；仍需监测边界比例和潜在效应积累。[R1](#r1)、[R9](#r9)

给定传递基因型，新增突变的原始性状均值和协方差为：

\[
E[\Delta\boldsymbol z^{raw}]=2s_e\sum_\ell\mu_\ell\boldsymbol m_\ell,
\]

\[
M=2s_e^2\sum_\ell[\mu_\ell\Sigma_\ell+
\mu_\ell(1-\mu_\ell)\boldsymbol m_\ell\boldsymbol m_\ell^{\mathsf T}].
\]

这是按当前规则直接推导的突变输入，不能当成种群现存遗传方差或裁剪后表型方差。现有诊断已报告逐坐标输入；完整矩阵及原始遗传值协方差是后续可补充的观测，不是本次已实现功能。

## 5. 实施顺序、出处与完成标准

| 步骤 | 内容 | 理论／方法依据 | 状态与验收 |
| --- | --- | --- | --- |
| 1 | 固定 Python 生态、SLiM 遗传的职责与生命过程顺序 | nonWF 自定义人口学与 ODD：[R2](#r2)、[R10](#r10) | 已接入；本代后代不提前参与本代生态，双亲名单由原生态阶段给出 |
| 2 | 将现有 YAML/CSV 及当前定相基因型导入原生染色体 | 固定 QTL 与连锁图谱：[R1](#r1)、[R6](#r6) | 已实现；不从边际频率重新抽另一份亲本，不改输入效应或单位 |
| 3 | 用原生 `addCrossed` 替换自写配子形成；Eidos 执行明确的突变核 | [R1](#r1)、[R4](#r4)、[R5](#r5) | 已实现；生产路径不调用 NumPy 重组或突变抽样 |
| 4 | 验证来源、效应更新及新 ID，所有材料通过后提交下一代 | 遗传身份及实际双亲贡献：[R1](#r1)、[R9](#r9) | 已实现；失败不产生半代群体或孤儿等位记录 |
| 5 | 升级协议、随机记录、源码身份与独立审计 | 可复现过程描述：[R10](#r10) | 已实现；记录原生引擎版本及 binary/script hash，区分标记切换和交换计数 |
| 6 | 通过实际 SLiM 路径检查分离、连锁、独立染色体、突变及回滚 | 第 4 节解析预期 | 36 组项目检查通过，包含实际 SLiM 路径与失败回滚 |
| 7 | 保存小型实际批次，审计、汇总和安装包验证 | 独立重建与复现：[R10](#r10) | 4 次保存批次及独立审计、汇总通过；wheel 隔离安装后实际运行与审计通过 |
| 8 | 在遗传执行稳定后，另行定稿成本角色和密度作用阶段 | 项目因果推导；防御成本、空间竞争背景：[R11](#r11)、[R12](#r12) | 本次不混入修改；仍按 v0.1.2 理论草案讨论 |

独立与连锁处理应使用相同的实际创始基因型、效应、频率、空间和人口学输入，只改变图谱。偏好与表现的关联有实证依据，也有独立遗传的反例，二者都应作为待研究假说。[R7](#r7)、[R8](#r8)

## 6. 验证什么，不能据此声称什么

当前测试以实际生产入口调用 SLiM，不再用删去的 Python 配子函数证明生产路径正确：

- 同一固定材料检查 Aa 配子 1:1、Aa×Aa 子代 1:2:1。
- 10,000 个三位点后代，共 20,000 份配子，核对 10 cM 的重组频率及相邻区间联合切换；独立染色体检查四类来源组合。
- μ=0/1、一般事件概率、条件正态增量、配置的非零均值、非 target 零效应、旧等位不可变及新效应代数。
- 相同输入与 seeds 可复现，突变开关不改本次调用的遗传来源。
- 实际进程缺失、损坏输出、植物成功而昆虫阶段失败时，目录、后代和贡献不提交。
- 保存结果后独立重建基因型—性状、来源、突变、资源、人口学和贡献账本。

离散随机频率使用预先给定的六标准误或明确绝对容差；这些是实现检查，不是生物学检验。稀有事件不能用一个过宽绝对容差宣称已验证；事件数不足应报告精度不足。不同软件或旧版本同 seed 不要求轨迹相同。

当前局部交配、性别及可变种群数量不满足理想 WF 的全部假设，不能要求整模型漂变方差等于 `p(1-p)/(2N)`。[R2](#r2) 无后代时等位频率和杂合度缺失，不能记零；频率比较需注明非空群体的分母，并另报灭绝。

测试通过意味着已检查场景下实现与规则相容，不证明功能映射、突变率、成本函数已得到实证校准，也不证明已经产生共同演化、长期稳定网络或新物种。[R15](#r15)

## 7. 改动与记录

- `helpers/slim_reproduction.slim`：原生染色体、共享等位对象、来源标记、`addCrossed` 及 Eidos 突变。
- `helpers/slim_backend.py`：执行器定位、版本检查、输入交换、超时及返回材料检查。
- `helpers/genetics.py`：调用实际引擎并衔接现有等位目录、性状和亲子记录。
- `audit_genetics.py`：独立核对 SLiM 调用身份、来源切换、缺失交换次数、突变与效应。
- 配置、运行元数据、源码 hash、包数据及现有测试同步升级。

每次有后代的调用记录 `kind/generation/offspring_count`、两个 SLiM seed、执行文件路径、版本和 binary/script SHA-256。生成的临时数字输入和返回文件在完成后清理；保存的快照、配置、等位目录及亲子记录可重建其科学输入。完整源码身份包含 `.slim` 文件，wheel 包含实际执行脚本。

使用独立代际进程不保存 SLiM tree sequence 或完整物理染色体祖先树。目前保留的是功能位点状态与 evoSim 的双亲/等位谱系账本；不能把它称为全基因组祖先记录。

## 8. 当前交付与后续边界

当前交付是可运行的 SLiM 遗传接入。正式机制实验仍须先完成 IBM 理论定稿。原始遗传值协方差、更广参数域及长程漂变核对可在有明确问题时补充；不为当前接入增加显性、上位性、真实 DNA 或自动成种。

发布前检查命令：

```powershell
python -m ruff check evo_sim tests example.py
python -m ruff format --check evo_sim tests example.py
python -m tests
```

2026-10-26 验证完成：Ruff 检查／格式检查及 **36 组项目测试通过**。测试包含 20 个小型场景各运行两轮的可复现与独立审计检查；没有以是否出现专化或共存作为通过条件。

另外保存了 `results/results_v0.1.3/slim_migration_20261026` 批次：4 次运行、每次 2 次更新，全部完成并通过独立遗传／贡献审计及汇总，共 460 个后代、37 次 QTL 突变，资源超支为零，最大相对平衡误差为 `2.220446049250313e-16`。实际 SLiM 调用数为 16。这些仅是实现与衔接检查。

`evosim-0.1.3-py3-none-any.whl` 已构建，SHA-256 为 `b344ac80526330e022a0c3495e43eca755905006440e2793d5076f39d4c2fa1e`。安装包包含实际 `.slim` 脚本；在仓库外隔离安装后，使用配置的外部 SLiM 成功完成一代并通过独立审计。SLiM 本机二进制 SHA-256 为 `53ce60f903de3443432e90143691ee7ba87b10f9fe859e62c8540227b4cfa7a1`。

### 2026-10-26：200 代探索性运行和图件

原 `sexual_qtl_200gen.yaml` 已完成：突变开启／关闭各 3 个重复，共 6 次 200 代运行，全部审计通过；保存 734,438 个后代、58,611 次 QTL 突变及 1,206 份快照。资源超支为零，最大相对平衡误差 `4.440892098500626e-16`。

A–E 五组图（种群、功能性状、寄主利用与灭绝、空间、遗传多样性）已生成 PNG/PDF/SVG 并目视检查。结果位于 `results/results_v0.1.3/slim_200gen_20261026`，其中 `RUN_REPORT.md` 提供索引，`audits/` 保存审计，`data/summary/` 保存汇总。

H1 在突变开启重复 2、3 中分别于第 78、26 代首次记录为零，在突变关闭重复 1 中于第 121 代首次记录为零；其余重复 H1 保留到第 200 代，两种植物在全部重复中均保留。所有重复均纳入图件，参数没有依据这次结果调整。这是探索性描述，不构成突变导致灭绝或网络长期稳定的证据。

正式机制研究、物种形成和长期网络稳定性验证仍未完成。运行结果、图件和本地 SLiM 安装按仓库约定保存在 Git 忽略目录，不随源码标签分发。

## 9. 文献与具体支持范围

以下文献沿用已核验题录和相应原文／摘要。软件手册说明执行语义；理论和实证文献说明机制依据与边界。文件组织、进程桥接和检验容限属于项目工程／验证选择，不伪称论文结论。

### R1

Haller, B. C., & Messer, P. W. **SLiM Manual, version 5.2**，2026-04-11 修订。[官方 PDF](https://benhaller.com/slim/SLiM_Manual.pdf)；[官方主页](https://messerlab.org/slim/)。手册无 DOI，不编造。

已读取 §13 开头及 §13.1–13.3 相关段落（第 279–286 页）、§14.11 重组概率解释（第 355–356 页）及 `mutationStackPolicy` 说明。支持固定 QTL、效应与适合度分开、等位替换、固定效应保留和重组概率转换；不证明 evoSim 的生物学函数或参数正确。核验 PDF SHA-256：`311fbf067cf6d6974a3dd50bd8af24a58e68f2746d7478c6d73b7121936f0577`。

### R2

Haller, B. C., & Messer, P. W. (2019). **SLiM 3: Forward Genetic Simulations Beyond the Wright–Fisher Model.** *Molecular Biology and Evolution*, 36, 632–637. [DOI: 10.1093/molbev/msy228](https://doi.org/10.1093/molbev/msy228)。

核验原文非 WF 过程说明。支持区分固定规模 WF 与具有显式存活、繁殖、空间和数量变化的模型；不要求 evoSim 改为重叠世代，也不验证当前生活史。

### R3

Haller, B. C., & Messer, P. W. (2023). **SLiM 4: Multispecies Eco-Evolutionary Modeling.** *The American Naturalist*, 201, E127–E139. [DOI: 10.1086/723601](https://doi.org/10.1086/723601)。

核验原论文摘要及方法背景。支持显式遗传与多物种生态互作可以在 SLiM 中结合；不表示必须迁移，也不表示多物种功能等于物种形成。

### R4

Chakraborty, R., & Nei, M. (1982). **Genetic differentiation of quantitative characters between populations or species: I. Mutation and random genetic drift.** *Genetical Research*, 39, 303–314. [DOI: 10.1017/S0016672300020978](https://doi.org/10.1017/S0016672300020978)。

核验出版社全文第 303–305 页。其模型部分明确讨论以原等位效应为中心的连续正态突变传统，并提出自身的离散等位状态模型。支持区分效应状态、突变和漂变；本文并未把该论文的离散模型说成 evoSim 正在采用的模型。该文所述 Kimura (1965) 原始论文 DOI 为 [10.1073/pnas.54.3.731](https://doi.org/10.1073/pnas.54.3.731)，本次核验其题录，未据未读扫描页补充细节。

### R5

Turelli, M. (1984). **Heritable genetic variation via mutation-selection balance: Lerch's zeta meets the abdominal bristle.** *Theoretical Population Biology*, 25, 138–193. [DOI: 10.1016/0040-5809(84)90017-0](https://doi.org/10.1016/0040-5809(84)90017-0)。

本次核验题录，结合既有审查作为突变—选择平衡的理论背景；不据此声称正态效应、μ=0.01 或 σ=0.02 已被校准。第 4.4 节矩阵由本模型独立推导。

### R6

White, I. M. S., & Hill, W. G. (2020；2019 年在线发表). **Effect of heterogeneity in recombination rate on variation in realised relationship.** *Heredity*, 124, 28–36. [DOI: 10.1038/s41437-019-0241-z](https://doi.org/10.1038/s41437-019-0241-z)。

核验原文图谱推导段落及正式题录。支持 Haldane 图谱参数与遗传距离的对应，以及物理长度和图谱长度的区别。本文的奇偶概率直接由 Poisson 分布推导，不将无干扰假设说成普适实证事实。

### R7

Hawthorne, D. J., & Via, S. (2001). **Genetic linkage of ecological specialization and reproductive isolation in pea aphids.** *Nature*, 412, 904–907. [DOI: 10.1038/35091062](https://doi.org/10.1038/35091062)。

核验出版社及原研究摘要。支持特定蚜虫系统中研究生态专化与隔离性状的遗传关联；不能据此固定所有 I/A 位点的连锁距离，不能把关联直接当作成种完成。

### R8

Forister, M. L. (2005). **Independent inheritance of preference and performance in hybrids between host races of Mitoura butterflies (Lepidoptera: Lycaenidae).** *Evolution*, 59, 1149–1155. [DOI: 10.1111/j.0014-3820.2005.tb01051.x](https://doi.org/10.1111/j.0014-3820.2005.tb01051.x)。

核验原研究摘要与既有题录。提供偏好与表现可以独立遗传的实证依据；支持把其关联作为研究问题，而不是先验保证相互匹配。

### R9

Price, G. R. (1970). **Selection and Covariance.** *Nature*, 227, 520–521. [DOI: 10.1038/227520a0](https://doi.org/10.1038/227520a0)。

沿用 v0.1.2 已核验题录和理论背景。支持区分贡献差异与传递变化；本计划的双亲贡献和条件频率期望是针对当前二倍体过程的计数推导，不冒称为该文直接给出的植物—昆虫公式。

### R10

Grimm, V., et al. (2020). **The ODD Protocol for Describing Agent-Based and Other Simulation Models: A Second Update to Improve Clarity, Replication, and Structural Realism.** *Journal of Artificial Societies and Social Simulation*, 23(2), 7. [DOI: 10.18564/jasss.4259](https://doi.org/10.18564/jasss.4259)。

核验原文方法及修订模型描述说明。支持明确实体、状态、过程顺序、假设和模型变化；不证明某个生态函数正确，也不直接规定本计划的统计容限。

### R11

Züst, T., Joseph, B., Shimizu, K. K., Kliebenstein, D. J., & Turnbull, L. A. (2011). **Using knockout mutants to reveal the growth costs of defensive traits.** *Proceedings of the Royal Society B*, 278, 2598–2603. [DOI: 10.1098/rspb.2010.2475](https://doi.org/10.1098/rspb.2010.2475)。

核验原研究摘要。支持防御成本的存在和条件依赖；不支持父源无成本，也不确定本模型成本作用的繁殖角色。第 5.2 节是项目的最小因果取舍。

### R12

Bolker, B. M., & Pacala, S. W. (1999). **Spatial moment equations for plant competition: Understanding spatial strategies and the advantages of short dispersal.** *The American Naturalist*, 153, 575–602. [DOI: 10.1086/303199](https://doi.org/10.1086/303199)。

沿用 v0.1.2 已核验的原研究摘要。支持区分空间竞争与扩散过程；不验证当前 phi 的具体函数，也不要求本版增加落点招募。

### R13

Auer, T. O., Khallaf, M. A., Silbering, A. F., et al. (2020). **Olfactory receptor and circuit evolution promote host specialization.** *Nature*, 579, 402–408. [DOI: 10.1038/s41586-020-2073-7](https://doi.org/10.1038/s41586-020-2073-7)。

沿用 v0.1.2 已核验的原研究摘要及实验说明。支持感知遗传变化参与寄主专化；不支持将抽象 QTL 直接命名为真实 OR，也不验证本模型的偏好函数。

### R14

Maliet, O., Loeuille, N., & Morlon, H. (2020). **An individual-based model for the eco-evolutionary emergence of bipartite interaction networks.** *Ecology Letters*, 23, 1623–1634. [DOI: 10.1111/ele.13592](https://doi.org/10.1111/ele.13592)。

沿用 v0.1.2 对原研究摘要和引言的核验。支持由个体生态—演化过程研究网络形成的总体路线；不直接验证本项目遗传规则、边权或物种定义。

### R15

Dieckmann, U., & Doebeli, M. (1999). **On the origin of species by sympatric speciation.** *Nature*, 400, 354–357. [DOI: 10.1038/22521](https://doi.org/10.1038/22521)。

沿用 v0.1.2 已核验的原文摘要。提供生态分化与交配机制共同参与隔离形成的理论先例；不提供通用成种阈值，不要求植物和昆虫采用相同成种机制。
