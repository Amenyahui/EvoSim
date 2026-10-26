# evoSim v0.1.4 — Model：植物—植食性昆虫显式遗传生态—演化 IBM

## v0.1.4：双向功能观察的实际定义

更新：2026-10-27。本版扩展已保存个体数据的分析，保持遗传与生态规则及其随机抽样顺序。目标是检验双方功能变化怎样影响互作、贡献和寄主网络；不是仅依据性状趋势判定共同演化。

D 为防御功能强度，Z 为独立气味遗传性状，V 为最终气味表征；A 为综合识别/偏好的有效模板，I 为防御对应的利用能力。它们允许讨论 PSM/VOC、嗅觉和解毒的功能假说，但不直接模拟真实化合物浓度、OR/OBP 序列或 P450/GST/UGT 酶动力学。D 同时经防御和 D→V 两条路径起作用，改变 D 不是纯防御干预。[R10,R14]

### 观察性时间交叉试验

从保存的早期 e 和后期 l 状态取真实个体分布，对 (P_e,H_e)、(P_e,H_l)、(P_l,H_e)、(P_l,H_l) 分别评估。每个植物—昆虫个体对都是独立的单寄主试验，统一 R*、d* 及其量纲，使用运行的相同行为和成本参数；不模拟空间、气味打乱、其他消费者、配偶或后代。物种内所有个体对等权平均，保留 A/I、D/Z 的联合分布，不用平均性状代入非线性公式。

设 u 为下文气味效用，外部选项为 u0；允许拒绝时，单寄主接受概率为：

\[
a_{ip}=\frac{\exp(\tau u_{ip})}{\exp(\tau u_{ip})+\exp(u_0)}.
\]

这里外部效用 u0 不再乘 τ，与生产模型一致。τ=0 时 a=1/(1+exp(u0))，并非无条件接受；forced 模式 a=1。该函数是已有离散选择规则的单候选特例，不是新增的生物经验规律。

给定已接受，沿用瓶颈 q 及资源分配函数：

\[
x^+_{ip}=\min(R^*,d^*q_{ip}),\quad
E[x]=a_{ip}x^+_{ip},\quad E[g]=a_{ip}x^+_{ip}/d^*,\quad
E[h]=a_{ip}x^+_{ip}/R^*.
\]

先算接受后的取食，再按接受概率求期望；不能将 a×q 先代入资源截断，这会改变结果。独立试验中 E[S_H]=E[g]，E[S_P]=1−(1−sP0)E[h]。输出 insect_SF=E[S_H]F_H、plant_SF=E[S_P]F_P。F 按各自个体投入计算；先对个体对计算乘积再平均。**SF 是受控条件下的生态代理，不包含密度、性别、配偶、父源成功，不等于实际繁殖适合度。**

固定祖先植物面板的 reference_p 仍是多候选条件相对偏好；它与这里的绝对单寄主接受概率不同。跨时间差异提供功能响应证据，但不足以独立排除漂变或证明共同演化。时间交叉比较的原始实验先例为 Buckling & Rainey (2002), [10.1098/rspb.2001.1945](https://doi.org/10.1098/rspb.2001.1945)，研究细菌—噬菌体；支持检验策略，不直接验证植物—昆虫函数。该试验设计由现有函数推导；不是声称文献验证了本项目的全部参数或函数。

### 性状的双亲贡献与传递残差

对同一物种、已提交的 t→t+1，令 w_i=(母源子代数+父源子代数)/2，z 为已映射表型坐标：

\[
\bar z_w=\frac{\sum_i w_i z_i}{\sum_i w_i},\qquad
\Delta_{contrib}=\bar z_w-\bar z=\frac{\operatorname{Cov}(z,w)}{\bar w},\qquad
\epsilon=\bar z_{t+1}-\bar z_w.
\]

因此 Δz=Δcontrib+ε 是按本项目双亲计数定义的代数恒等式，依据贡献协方差/传递分解的理论思路；不能把 realized_differential 全部归为确定性自然选择。有限种群的抽样也影响 w；ε 汇集遗传传递、突变及截断等差异，不能仅称作突变效应。无后代时加权均值未定义；终代或失败更新未提交时各贡献项为 NA。现有等位频率分解仍用于直接检查遗传传递。

依据：Price (1970), *Selection and Covariance*, [10.1038/227520a0](https://doi.org/10.1038/227520a0)。此文支持分解原则，本项目的双亲权重和截断表型残差由上述定义单独说明。

### 寄主与网络观察

逐代关联机会、实际接受/拒绝、固定面板 p/q、实际资源、母源/父源遗传贡献及双方丰度。资源和贡献两层网络分别按预先指定阈值计边，以相邻代的连续长度定义观察窗口内持续性。所有阈值与窗口来自 `analysis/v014_analysis.yaml` 或批次 `analysis.yaml`，只用于观察。

每只虫的有效寄主数为 exp(−Σr log r)。p 直接作为 r；q 需先除以其跨寄主之和，全零时宽度未定义。逐虫计算后报告均值、方差及正响应个体数；同时保留 q 的绝对水平，避免把低利用误解成专化。该指标是操作统计，不是专化形成的因果判据。

偏好和表现可独立遗传的实证见 Forister (2005), [10.1111/j.0014-3820.2005.tb01051.x](https://doi.org/10.1111/j.0014-3820.2005.tb01051.x)；其产卵/幼虫生活史不能直接验证本模型自主取食。由个体过程形成网络的理论先例见 Maliet et al. (2020), [10.1111/ele.13592](https://doi.org/10.1111/ele.13592)；丰度与机会的解释限制见 Blüthgen (2010), [10.1016/j.baae.2010.01.001](https://doi.org/10.1016/j.baae.2010.01.001)。

完整固定功能分布的四组演化控制、自动判定遗传性寄主转换、扰动后恢复和成种尚未实现。mutation-off 仍可由站立变异发生演化，不能代替冻结对方功能分布。有限窗口的边持续性不能当作永久稳定。

## v0.1.3 SLiM 遗传执行契约

当前在 Python IBM 中实际调用 SLiM 5.2 的原生染色体和 `addCrossed()`。QTL 突变在 Eidos 中按已有离散规则定义，性状映射及生态仍由 Python 执行。协议为 `slim-qtl-1`；引擎身份、来源标记及失败处理见 [v0.1.3 说明](../README_v0.1.3.md)。采用 SLiM 不等于参数获得实证校准。

v0.1.2 的个体联合寄主反应、双亲贡献及网络观察继续使用。父本仍采用 F>0 资格门槛；密度因子仍作用于母源产量，配置阶段标签保留历史名称；当前生态公式未随引擎迁移改变。旧长程模拟只作为历史探索。

## 历史配置说明：2026-10-22 默认实验修正

旧配置的 6 只昆虫在 10×10 网格中只有 1 格交配半径，同种子产生 5 雌 1 雄且无空间合格配偶；每雌 λ=2 在 F<1、雌性比例 0.5 下还低于理想替代水平。此为演示配置不匹配有性生命周期，不能据此宣称昆虫必然快速灭绝。

默认演示场景现改为：植物 P1/P2 各 200 株、昆虫 40 只；昆虫交配半径 2 格；每雌繁殖系数 λ=4；植物完全受食后的基础存活率 sP0=0.5（原 0.1）。单虫需求和单株资源仍均为 1。取食、死亡和交配公式不变，不补种、不强制雌雄比例、不保证存续。sP0 调整是一项显式的植物耐受假设变化，不能当作程序错误修复或实测校准。

依据：近似平衡性别下，E[N(t+1)]/N(t)≈0.5×S×q_mate×λ×F。λ=4 在高存活、配偶充分时允许增长；实际密度仍由资源与空间反馈约束。增加植物初始数量降低初始取食压力；提高 sP0 表示受食后仍有繁殖机会的耐受型理论植物，避免原场景的强烈寄主崩溃反馈。参数选择使用了试运行，因此后续曲线是探索性诊断，不能视为独立验证或正式处理效应检验。

保留全部探索记录：`codex/validation/mating_correction_*pilot.json`。包括只提高昆虫密度、提高植物密度、改变需求尺度仍发生灭绝/暴发的结果，没有删除不利种子。最终场景的三个种子 20261021–20261023 在 50 代时昆虫数为 108/134/127。

200 代配置现为两个处理各 3 个重复。基因型累计存储上限由 1,000,000 提高到 4,000,000 个标量（按 int64 约 32 MB 的纯基因型数据，不含日志与对象）；个体上限 5,000 和 allele 上限 100,000 不变。仍保留安全终止，不截断出生数。完整原始记录可能占用更多内存。

绘图逐重复展示，单个观测使用点标记，时间轴从 0 开始，标注灭绝和实际停止状态；遗传图阴影仍表示种群内 SD，不能当作重复间置信区间。等位频率图每物种/位点最多显示峰值频率最高的 12 个 allele，全部原始数据保留。


> **用途**：论文 *Materials and Methods / Model Description* 的中文初稿；v0.1.4 机制说明；遗传协议 `slim-qtl-1`。
> **模型类型**：显式个体、空间局部互作、双亲有性遗传、非重叠世代的前向随机仿真。  
> **状态**：下列为当前实现的公式契约，不是已获得的生态科学结论；新增观测的字段与缺失语义见 [输出规范](results.md)。
> **证据标记**：`[B]` 经典遗传/数学规律，`[M]` 成熟理论形式，`[H]` 本研究可证伪的机制性假设，`[I]` 实现约定/占位参数；引用的论文只能支撑其相应层次。

## 摘要

本研究构建空间显式的植物—植食性昆虫个体基础模型（individual-based model, IBM），以连接定量遗传变异、化学防御和气味相关功能性状、寄主选择、有限资源取食及生态人口学。模型将植物与昆虫表示为具有二倍体 QTL 基因型的独立个体，允许同种双亲经配子形成、孟德尔分离、连锁重组与新突变产生子代；植物的防御性状及 VOC 性状影响昆虫寄主选择及利用收益，进而影响两类生物的生存与繁殖。时间按非重叠的完整世代推进，不将一代等同于自然年。当前 QTL 与化学维度均为理论性状轴，不是实际 DNA、P450 蛋白活性或单一化合物浓度。该设计沿用 SLiM nonWF 式的出生/死亡分开控制范式 [R02,R03]，并以 ODD 协议界定状态和过程 [R01]。

## 1. 研究对象、空间和时间（ODD：Purpose, Entities, Scales）

- **植物物种**：`P1, P2, ...`；每株植物是独立的二倍体、两性花、严格异交种子母本/花粉父本。
- **植食性昆虫物种**：`H1, H2, ...`；每只昆虫是独立的二倍体雌虫或雄虫，雌雄均经历本代取食与生存检验。当前不模拟真实昆虫染色体性别决定。
- **时间**：离散代数 `t=0,...,T`。从 `t` 到 `t+1` 包括选择、取食、生存、交配、出生及扩散；全部亲代退出种群；终点 `T` 只观察不产生 `T+1`。
- **空间**：大小为 `L_x × L_y` 的二维周期栅格，多株植物可共处同一栅格；格长为抽象斑块，不代表给定米数。昆虫有初始搜寻位置，其本代交配位置按实际取食寄主的位置表示。
- **状态**：植物基因型和性状 `G^P,D,Z,V`；昆虫基因型和性状 `G^H,I,A`；各自的 `species, id, cell`，昆虫稳定 `sex`，生态账本与谱系事件。

**目的**是搭建机制可审计的理论 IBM；是否形成寄主专化、稳定共存、寄主转换或多物种网络，属于模型建成后的研究问题，而不是初始化时预设的结论。当前物种 ID 固定，不包含成种。

## 2. 显式二倍体遗传与表型

### 2.1 基因型及等位基因目录 [B/M]

对类群 `k∈{P,H}`，个体 `i`、同源拷贝 `h∈{0,1}`、位点 `ℓ∈{1,...,L_k}`：

\[
G_i^{(k)}=\{a_{ih\ell}\},\qquad \operatorname{shape}(G^{(k)})=(N_k,2,L_k).\tag{G1}
\]

`a` 是不可变等位基因 ID；目录保存 `allele_id, locus_id, effect_vector, parent_allele_id, birth_generation`。每个位点有 `chromosome_id`（理论连锁群）和 `position_cM`（遗传位置）。目前 `L_P=L_H=8`；每类群的4个性状坐标各由2个抽象 QTL 贡献。**这8个位点数及分配是 [I]，不是基因组测量。**

### 2.2 配子形成、孟德尔分离与连锁重组 [B/M]

在每条连锁群 `c` 上，按遗传位置排序 QTL。形成配子时独立选择起始同源拷贝：

\[
H_{p,c,1}\sim\mathrm{Bernoulli}(1/2),\qquad M_{p,\ell_1}=a_{p,H_{p,c,1},\ell_1}.\tag{G2}
\]

不同连锁群独立选择起始拷贝，代表不连锁位点自由组合 [B]。同一连锁群中，相邻两个 QTL 的遗传距离为：

\[
d_j=(x_{j+1}-x_j)/100\quad \text{Morgan};\qquad C_j\sim\mathrm{Bernoulli}(r(d_j)).\tag{G3}
\]

采用 **Haldane 无交换干扰的标记传递分布**，将 r(d) 输入 SLiM 压缩标记图谱，以原生染色体和 `addCrossed` 生成配子。C 表示标记间来源切换：

\[
H_{p,c,j+1}=H_{p,c,j}\oplus C_j,\qquad
M_{p,\ell_{j+1}}=a_{p,H_{p,c,j+1},\ell_{j+1}}.\tag{G4}
\]

因此标记间重组分数为：

\[
r(d)=\frac{1-\exp(-2d)}2.\tag{G5}
\]

该概率可从区间 Poisson(d) 交换的奇数概率推导；当前执行器不抽取或记录该完整 Poisson 次数。`switch_before_locus` 记录来源切换，`crossover_count` 缺失。相同位置的独立中性来源标记识别纯合亲本的同源拷贝。

例如10 cM 对应 `d=0.10 Morgan`、`r≈0.0906346`，**10 cM 不等于准确 10% 重组配子**；两次交换可能使两端标记的单倍型不变。无经验图谱时，位置与连锁群仅是预定义的理论架构，默认8个 QTL 分属于8个独立的抽象连锁群；在专门的连锁测试或研究情景中再设置多位点连锁群。[R04]

### 2.3 双亲受精与突变 [B/M/I]

同种合格母本 `m` 和父本 `f` 由 SLiM 原生构建后代；随后在本代新生拷贝上执行 Eidos QTL 突变，等价表示为：

\[
G_{o,0,\ell}=M'_{m,\ell},\qquad G_{o,1,\ell}=M'_{f,\ell}.\tag{G6}
\]

SLiM 自动基因组突变率设为零；由 Eidos 对每份传递拷贝在位点 `ℓ` 上独立执行本项目声明的突变核：

\[
U_{o,h,\ell}\sim\mathrm{Bernoulli}(\mu_\ell),\qquad
\boldsymbol e_{a'}=\boldsymbol e_a+\Delta\boldsymbol e_\ell,\quad
\Delta\boldsymbol e_\ell\sim\mathcal N(\boldsymbol m_\ell,\mathrm{diag}(\boldsymbol\sigma_\ell^2)).\tag{G7}
\]

`m_ℓ` 对应每个位点的 `mutation_mean`，由 `inputs/loci.csv` 的 `mutation_mean_j` 列读入，默认全为 0，保留无偏突变假设；`mutation_sd_j` 列指定各坐标标准差。改变均值会引入定向突变，不能解释为选择的结果。新等位 `a'` 使用新的唯一 ID、父等位 ID 和事件 ID；不能原位改变 `a`。`μ_ℓ` 为**每次配子传递一个抽象 QTL 的突变概率**，不是每核苷酸突变率。若各位点相同、两个配子各有 `L` 个位点：

\[
P(\text{子代至少有一个新 QTL 突变})=1-(1-\mu)^{2L}.\tag{G8}
\]

### 2.4 基因型—性状 [M/H]

令个体性状向量 `\boldsymbol T_i` 分别为植物 `(D_1,D_2,Z_1,Z_2)` 或昆虫 `(I_1,I_2,A_1,A_2)`，各等位基因贡献可正可负的功能效应：

\[
T_{iq}^{\rm raw}=b_q+s_e\sum_{\ell=1}^{L_k}\sum_{h=0}^1 e_{a_{ih\ell},q},\qquad
T_{iq}=\operatorname{clip}_{[l_T,u_T]}(T_{iq}^{\rm raw}).\tag{G9}
\]

采用加性 QTL 模型 [M]；`s_e` 读取 `genetics.phenotype_map.effect_scale`，默认 `1/2` 是**本模型对二倍体拷贝效应尺度的约定** [I]，不是显性系数。`[l_T,u_T]` 读取同组 `bounds`，默认 `[0,1]`，且须满足 `0≤l_T<u_T≤1`。显性、上位性与群体方差定标不在本版实现；Chu 等 2024 提供的是个体功能遗传效应与群体统计遗传效应区分的理论，而非 P450 蛋白功能关系 [R05]。必须记录未截断 `T_raw`，否则长期演化的遗传效应隐藏积累无法识别。当前基因突变**仅经 G9**影响性状；性状没有额外随机突变。

## 3. 植物气味、昆虫选择与资源竞争

### 3.1 VOC 合成信号 [H]

植物 `p` 的防御向量为 `\boldsymbol D_p`，独立 VOC 性状为 `\boldsymbol Z_p`，共用 `B` 映射矩阵及耦合系数 `ρ`：

\[
\boldsymbol V_p=\operatorname{clip}_{[l_V,u_V]}[\rho\,B\boldsymbol D_p+(1-\rho)\boldsymbol Z_p].\tag{E1}
\]

`[l_V,u_V]` 从 `voc_bounds` 读取，默认 `[0,1]`，须满足 `0≤l_V<u_V≤1`。这是一种**现象学的信号耦合假设**，不是萜烯合酶反应速率公式。PSM 化学谱与植物 VOC 释放受遗传变异调节的实验背景见 [R08,R09]；具体线性混合形式及 `ρ` 为本研究假设，尚无实验参数估计。

### 3.2 空间候选与寄主效用 [H/M]

设植物位置 `\mathbf x_p`、昆虫搜寻位置 `\mathbf x_i`。周期二维网格上的切比雪夫距离：

\[
d_\infty(\mathbf x,\mathbf y)=\max_{a\in\{x,y\}}\min(|x_a-y_a|, L_a-|x_a-y_a|).\tag{S1}
\]

昆虫只观察 `d_∞≤r_search` 的植物，若可见植物多于 `candidate_size`，则无放回均匀抽取一个候选集合 `\mathcal H_i`；它不是整个理论寄主物种库。[I]

偏好模板 `\boldsymbol A_i` 与候选植物 VOC `\boldsymbol V_p` 的匹配效用：

\[
u_{ip}=b_{\rm odor}-\frac{1}{2\sigma_{\rm odor}^2}
\sum_{k=1}^{n_{VOC}} w_k(V_{pk}-A_{ik})^2.\tag{E2}
\]

若允许拒绝全部候选（outside option），候选 `p` 的被选概率为：

\[
P(H_i=p\mid\mathcal H_i)=
\frac{\exp(\tau u_{ip})}{\exp(u_{\rm out})+\sum_{z\in\mathcal H_i}\exp(\tau u_{iz})},\qquad
P(H_i=\varnothing)=\frac{\exp(u_{\rm out})}{\exp(u_{\rm out})+\sum_z\exp(\tau u_{iz})}.\tag{E3}
\]

`τ≥0` 为模型定义的偏好判别强度；外部选项在此实现中不乘 `τ`，因此 `u_out` 与 `τ` 的尺度关系必须保持固定契约。启用强制选择时省去外部选项并按候选重新归一化。**Softmax 是成熟的离散选择数学形式 [R07]，但不意味着昆虫真实决策必定遵从这个特定效用函数。** 嗅觉受体和环路参与寄主行为演化的经验依据见 [R10]。

寄主分配记为：

\[
H_i\in\mathcal H_i\cup\{\varnothing\},\qquad \text{由 E3 抽取且本代固定。}\tag{E4}
\]

### 3.3 植物利用能力和资源账本 [H/I]

植物防御 `\boldsymbol D_p` 与昆虫利用能力 `\boldsymbol I_i` 共同产生有效需求比例：

\[
q_i=\begin{cases}
\min\left(1,\min_{j=1,\dots,n_D} I_{ij}/D_{H_i,j}\right),&H_i\neq\varnothing;\\
\mathrm{NA},&H_i=\varnothing.
\end{cases}\tag{E5}
\]

约定分母 `D=0` 的该通道不构成限制，按 `+∞` 处理（包括 `I=D=0`），最后由外层最小值限制到 `[0,1]`。`q` 是**相对有效摄食能力**，不是体内酶动力学活性。底物特异性的实证背景见 [R11]；最限制通道及 I/D 摄入预算的直接模型先例见 [R14]。这是有明确条件的理论抽象，不是普适生化定律；本项目的单位推导见 [最小模型规范第 3 节](plan/v0.1.2_minimal_model.md)。

若昆虫种 `s` 的最大单代需求为 `d^{max}_s`、植物 `p` 供给上限 `R_p>0`：

\[
d_i=d^{max}_{s(i)}q_i,\qquad T_p=\sum_{i:H_i=p}d_i.\tag{E6}
\]

\[
\gamma_p=\begin{cases}\min(1,R_p/T_p),&T_p>0,\\1,&T_p=0,\end{cases}
\quad x_i=\gamma_{H_i}d_i,\quad X_p=\sum_{i:H_i=p}x_i\le R_p.\tag{E7}
\]

\[
g_i=x_i/d^{max}_{s(i)},\qquad h_p=X_p/R_p.\tag{E8}
\]

共享资源**同时按比例分配**，不允许昆虫先来先得。`0≤g_i,h_p≤1`；`X_p=min(T_p,R_p)` 在浮点容差内成立。生态适合度根据实际 `x` 与 `g`，不根据未分配的需求 `d`。[I]

### 3.4 防御/气味成本与生存概率 [H]

植物因防御和 VOC 输出支付繁殖成本：

\[
F^P_p=\max\left[0,1-c_D\sum_j D_{pj}-c_V\overline V_p\right],\qquad
\overline V_p=n_{VOC}^{-1}\sum_k V_{pk}.\tag{E9}
\]

昆虫利用能力成本：

\[
F^H_i=\max\left[0,1-c_I\sum_j I_{ij}\right].\tag{E10}
\]

两类生物生存至交配阶段的概率为：

\[
S^P_p=s_{P0}+(1-s_{P0})(1-h_p),\qquad S^H_i=g_i.\tag{E11}
\]

以上线性成本与生存函数均为**理论参数化假设**；它们体现防御、利用、损伤与适合度之间可研究的权衡，但**不是直接摘自某篇实验文章的函数与数值**。`S^H_i=g_i` 的定义意味着没有资源收益则昆虫无法存活至交配阶段；当前不额外叠加虫龄、宿主免疫或生理发育模块。

## 4. 植物密度竞争与空间扩散

令 `\mathcal N_p(r_c)` 为以植物母本 `p` 所在格为中心、半径 `r_c` 的周期邻域，其中同一格只计算一次，且包括焦点植物本身。令 `N^{loc}_{p,s}` 为该邻域植物物种 `s` 的个体数：

\[
\phi_p=\left(1+\sum_s\alpha_{s(p),s}N^{loc}_{p,s}\right)^{-1}.\tag{S2}
\]

这是 Beverton–Holt 型负密度依赖的**局部个体层面变体** [R12,R13]，并非文献直接规定的空间竞争系数。`α` 维度为每个邻域植物个体的效应。`φ_p` 在**种子母本产子数**阶段折减，不是新生种子在目标格的另一次成活抽签。

植物种子与昆虫子代均按独立轴均匀抽取扩散位移：

\[
\mathbf x_o=(\mathbf x_{\rm origin}+\boldsymbol\delta)\bmod(L_x,L_y),\qquad
\delta_a\sim\mathrm{DiscreteUniform}\{-r,\ldots,r\}.\tag{S3}
\]

植物 `origin` 为种子母本格；昆虫 `origin` 为母方本代实际取食寄主格。花粉传播、昆虫交配、昆虫搜寻和子代扩散分别有不同半径，不得混用。扩散范围为理论斑块步数，**无实际距离单位**。

## 5. 有性繁殖、配偶缺失和种群更替

### 5.1 存活和繁殖角色 [B/H]

全体植物和昆虫对 `S_i` 各进行一次独立生存抽签：

\[
R_i\sim\mathrm{Bernoulli}(S_i).\tag{D1}
\]

昆虫出生时性别独立设为：

\[
P(\mathrm{female})=p_f,\qquad P(\mathrm{male})=1-p_f.\tag{D2}
\]

其中 `p_f=0.5` 是默认的理论值，不由当前 QTL 决定。生存抽签 `R_i=1` 仅表示有资格参加交配；无论是否通过抽签，亲代均不会跨世代保留。

### 5.2 空间配偶规则 [H/I]

对种子母本/雌虫 `i`，在当代 `R=1` 的同种个体里寻找满足以下条件的父本：

\[
\mathcal C_i=\{j:\;s(j)=s(i),\ R_j=1,\ F_j>0,\ d_\infty(\mathbf x_i,\mathbf x_j)\le r_{mate},\ \operatorname{compatible}(i,j)\}.\tag{D3}
\]

植物 `compatible` 表示 `j≠i`；昆虫表示母方雌、父方雄。植物花粉距离使用 `r_pollen`，昆虫交配距离使用 `r_insect_mate`；昆虫交配位置为其**本代实际取食寄主格**。父本由子代逐一**均匀有放回**选取：

\[
P(f=j\mid j\in\mathcal C_i)=1/|\mathcal C_i|.\tag{D4}
\]

严格异交是本理论植物系统的设定，文献支持植物交配系统的重要性，但不代表所有植物不能自交 [R06]。无雄性精子数限制、花粉竞争或性选择也属于简化，而不是实验事实。

### 5.3 条件生殖均值、出生及矩 [H/数学推导]

v0.1.2 的母源期望统计覆盖所有植物和所有雌虫，包括未找到寄主的雌虫（期望为零）；雄虫的母源期望缺失。`is_maternal_role` 仅表示角色，不表示已通过生态或交配筛选。按全体昆虫计算母源均值时，雄虫按零贡献计入；按雌虫计算时保留所有雌虫，不能只统计已找到寄主者。此改动修正观测分母，不改变抽样或生殖规则。

植物母本：

\[
\nu^P_i=\lambda_{P,s(i)} F^P_i\phi_i.\tag{D5}
\]

昆虫雌虫母本：

\[
\nu^H_i=\lambda_{H,s(i)} F^H_i.\tag{D6}
\]

这里 `λ_P` 和 `λ_H` 分别为**已经通过生存、已拥有至少一名可交配父本的母本**的基础后代均值；`λ_H` 不是所有雌虫雄虫都产仔的均值。个体子代数：

\[
B_i=\begin{cases}\mathrm{Poisson}(\nu_i),& i\text{ 是合格母本},\ R_i=1,\ |\mathcal C_i|>0;\\
0,&\text{否则}.
\end{cases}\tag{D7}
\]

设 `\mathcal C_i^0` 是生存抽签前、基于固定空间/性别/物种/`F_j>0` 得到的潜在父本集合。假设本代个体 `R` 相互独立，且母本不在自己的潜在父本集合中：

\[
q_i^{mate}=1-\prod_{j\in\mathcal C_i^0}(1-S_j),\qquad
p_i^{birth}=S_i q_i^{mate}.\tag{D8}
\]

在“给定本代取食结果、固定位置与 `S,F,φ`，生存与 Poisson 过程独立”的条件下，母本出生数的期望与方差为：

\[
\mathbb E(B_i)=p_i^{birth}\nu_i,\qquad
\operatorname{Var}(B_i)=p_i^{birth}\nu_i+p_i^{birth}(1-p_i^{birth})\nu_i^2.\tag{D9}
\]

**D9 是 Bernoulli–Poisson 混合分布的概率推导，不是文献中的生物学固定定律。** 在已实现生存和配偶筛选、`R_i=1,|C_i|>0` 的母本条件下，均值和方差均为 `ν_i`。多个母本可能共享父本，其出生数未必无条件独立，不能把所有 D9 的方差直接相加当作种群总方差。雄虫没有“母方子代数”及其矩，应记录父源贡献数量。

### 5.4 下一代组成及世代终止 [I]

\[
N_{s,t+1}=\sum_{i:\,s(i)=s}B_i.\tag{D10}
\]

实际出生子代由两亲本按 G2–G9 形成，按 S3 扩散。模型不引入迁入、种子库、出生保底、环境容量硬上限或成年亲代存续。若某物种 `N_{s,t}=0`，且无迁入/成种，其灭绝是吸收态；终点 `t=T` 记录观察，不计算新的 `B_i`。

## 6. 参数定义、分级及其后续用途

### 6.1 三种不同的参数工作，不能混为一谈

- **C：待校准/不确定参数（后期需要讨论估计、界限与敏感性）**。参考独立文献、实验、可辨识性或预注册敏感性范围；模型未建成前不因结果方向而调节。
- **X：科学情景自变量（未来用于机制比较与图件横轴/分组）**。应在正式假说检验前定义水平，例如化学防御维度、VOC–防御耦合、连锁强度、交配空间结构；这不叫为了结果好看而“调参”。
- **N：数值/配置控制参数（仅工程或实验规模）**。不能用于生物学拟合，如安全阈值、种子、记录步长、重复次数。

另有 **Y：输出指标**（不是参数），见第7节与 `results.md`。某个 C 参数未来也可以被预注册为 X，但同一次分析中必须声明其角色，不能先观察结论再改变分组。

### 6.2 参数总表（未核实数值标明“待定”）

通用参数和 CSV 路径由 [`configs/model_parameters.yaml`](../evo_sim/models/plant_insect/configs/model_parameters.yaml) 声明；物种、竞争、位点、等位效应、初始频率和参考性状分别由 [`inputs/`](../evo_sim/models/plant_insect/inputs/) 中的六张 CSV 提供，读取后展开为现有模型参数结构。Python 不提供模型默认值。所有键的实际位置、CSV 列、继承和覆盖规则见 [参数配置指南](parameters.md)。下表保留展开后参数的科学含义分类；数值均为理论示例。

| 参数 / 配置键 | 数学作用 / 含义 | 单位或约束 | 当前理论示例（如有） | 主分类 | 将来怎么用 / 依据 |
|---|---|---|---|---|---|
| `plant_sizes`, `insect_sizes` | 初始每物种个体数 | 个体，整数≥0 | P1=200, P2=200, H1=40 | X | 初始群落与灭绝风险情景，不作真实密度宣称 |
| `grid_shape` | 空间斑块数量 | 格数，两个正整数 | 10×10 | X | 空间尺度敏感性与空间热图；无物理长度 |
| `search_radius` | 昆虫候选寄主搜索半径 | 格，≥0 | 2 | C | 生态行为参数；理论区间/空间效应检查 |
| `candidate_size` | 每只昆虫最多评估的候选植株数 | 正整数 | 8 | C | 信息限制假设；文献不决定固定值 |
| `sigma_odor` / `sigma` | VOC–偏好距离核宽度 | 性状单位，>0 | 0.2 | C | 选择宽容度/专化潜力；需和 `τ` 联合检查可辨识性 |
| `tau` | Softmax 判别强度 | ≥0 | 1 | C | 选择强度；与 `sigma_odor` 不可任意独立解释 |
| `b_odor` / `b` | 气味效用截距 | 效用单位 | 1 | C | 有外部选项时影响拒绝比例；强制选择时共同截距抵消 |
| `outside_utility` | 不选择任何寄主的效用 | 效用单位 | 0 | C | 拒绝概率、无寄主死亡 [R07] |
| `weights` | VOC 维度距离权重 | 非负；约定和为1 | (0.5,0.5) | X | 化学维度识别/特异性情景；真实值无证据 |
| `B` | D→VOC 映射矩阵 | 2×2，输出限制[0,1] | 单位阵 | X | PSM–VOC 耦合结构理论情景 |
| `rho_voc` | VOC 中 D 映射贡献比例 | [0,1] | 0.5 | X | VOC–防御耦合是否影响寄主选择（绘图因素） |
| `resource_capacity` `R_p` | 每株单代最大可食用资源 | 归一化量，>0 | 各植物1 | C | 资源供给标度，需定义量纲一致性 |
| `max_feeding_demand` `dmax` | 单昆虫潜在摄食上限 | 同 R 的归一化量，>0 | H1=1 | C | 决定资源竞争与虫体收益 |
| `plant_cost` `c_D` | 防御性状繁殖成本 | 每性状单位 | 0.125 | C | 防御—繁殖权衡；可作为后期假设对照 X |
| `voc_cost` `c_V` | VOC 发生成本 | 每性状单位 | 0.05 | C | 气味信号成本，机制依据需限定 |
| `insect_cost` `c_I` | 利用能力繁殖成本 | 每性状单位 | 0.125 | C | 解毒广度与成本权衡，不等于真实代谢成本 |
| `sP0` | 全部资源被摄食时植物残余生存率 | [0,1] | 0.5 | C | 抵抗食害的生存基线 |
| `plant_offspring_lambda` | 可交配植物母本每代基础后代均值 | 子代/母本，>0 | P1/P2=2 | C | 对应 D5，绝非产生的实际种子粒数 |
| `insect_offspring_lambda` | 可交配雌虫基础后代均值 | 子代/雌虫，>0 | H1=4 | C | 对应 D6，不能以旧群体曲线反推 |
| `plant_competition_alpha` | 邻域物种竞争系数矩阵 | 每个邻域个体≥0 | 对角0.02/交叉0.005 | C | 共存机制及密度反馈；可作 X 的理论范围 |
| `competition_radius` | 植物竞争邻域 | 格，≥0 | 1 | X | 空间竞争尺度；绘图比较 |
| `pollen_radius` | 花粉可达范围 | 格，≥0 | 1 | C | 独立于种子扩散的交配空间尺度 |
| `insect_mate_radius` | 昆虫雌雄交配半径 | 格，≥0 | 2 | C | 交配限制与灭绝风险，单独探索 |
| `plant_dispersal_radius` | 种子扩散范围 | 格，≥0 | 1 | X | 空间集聚与迁移情景 |
| `insect_dispersal_radius` | 昆虫子代扩散范围 | 格，≥0 | 1 | X | 虫体空间结构及寄主转换情景 |
| `p_female` | 昆虫子代为雌性概率 | [0,1] | 0.5 | C | 无性染色体时理论性别分配，检验局部小种群效应 |
| `L_P,L_H` | 类群各自抽象 QTL 数 | 正整数 | 各8 | X | 由 `genetics.loci` 列表决定，需覆盖各性状坐标 |
| `chromosome_id`, `position_cM` | 连锁结构、重组几率 | 连锁组、cM | 默认8个独立组 | X | 未来连锁强度比较，**不是经验图谱** |
| `mutation_rate[ℓ]` | QTL/配子传递突变概率 | [0,1] | 每位点0.01 | C | 强理论占位，不能称核苷酸突变率 |
| `mutation_sd[ℓ]` | 突变效应增量标准差 | 性状效应单位 | 目标坐标0.02 | C | 需界限/敏感性；不按进化结果方向微调 |
| `mutation_mean[ℓ]` | 突变效应增量均值 | 性状效应单位；非目标坐标为0 | 全0 | X | 默认无偏；非零为额外的定向突变假设 |
| `genetics.phenotype_map.effect_scale` | 等位效应总和的标度 | 有限非负值 | 0.5 | X | 与等位效应单位共同定义；不是显性参数 |
| `genetics.phenotype_map.bounds`, `voc_bounds` | 遗传表型及 VOC 的截断区间 | 0≤下界<上界≤1 | 均为[0,1] | X | 归一化尺度内的边界假设，不是数值容差 |
| `founder_allele_frequencies` | 初始等位频率 | 每位点和为1 | 例子0.5/0.5 | X | 初始遗传差异对结果的影响；人为适应须明示 |
| `background_traits` | 加性性状基线 `b_q` | 数值，无量纲 | 各0 | C | 初始功能表型均值与截断率 |
| `num_steps_per_episode` | 代数 T | 非负整数 | 200仅示例 | N | 研究观察窗口；不是生理年龄 |
| `seed`, `replicates`, `record_every` | 随机复现/重复数/记录频率 | 整数 | seed按任务设置 | N | 精确复现与统计精度，不是生物性状 |
| `max_individuals_safety` | 内存/计算安全上限 | 个体数 | 5,000仅保护值 | N | 超限**报错不裁断**，不得被解读为承载量 |

### 6.3 参数来源等级及调整纪律

所有运行在 `metadata.json.config` 保存 YAML 与 CSV 展开后的实际参数及来源记录。通过配置文件运行时，`input_sources` 保存输入 CSV 的路径和 SHA-256 指纹；独立审计使用保存的展开数据，不依赖原 CSV。当前 `parameter_provenance` 使用 `kind/source`，遗传及交配来源使用 `source_kind/description`；没有为每个参数自动生成单位、引用或敏感性范围。单位由相关单位字段及本说明定义，科学依据仍需研究者核实并填写。**当前多数数值仍是 theoretical_placeholder。**

- **构建期**：只能因生物/数学逻辑不正确而改变规则；冒烟运行检查守恒、分离率、RNG 和边界；不按趋势方向调数值。
- **模型完成后校准**：用**独立数据**或预定目标统计量估计可辨识参数（C），记录不确定性及参数间耦合；不能把预测图上的“期望趋势”当校准目标。
- **正式机制比较前**：预先定义科学自变量 X 的水平、主要 Y 指标、重复数量与分析计划；随后的敏感性分析必须与定向调参分开。

## 7. 科学输出指标（Y，供后续结果绘图）

| 输出指标 | 定义 | 推荐图件/目的 |
|---|---|---|
| `N_{s,t}` / `extinction_time` | 各物种总数与首次达到0的代数 | 时间序列、灭绝事件分布 |
| 雌雄比、有效母本率、无配偶率 | 生存雌虫占比、母本有父本比例 | 种群维持与繁殖失败分解 |
| 各类 `D,Z,V,I,A` 分布 | 均值、方差、分位数、剪裁率 | 性状演化时序/分布图 |
| QTL 等位频率 `p_{ℓ,a}` | 每个物种、每个位点双拷贝计数/2N | 等位频率轨迹 |
| 杂合度 `H_O`, `H_E` | `H_O=杂合个体比例`, `H_E=1-Σp_a²` | 遗传多样性变化；后者也称 gene diversity |
| 连锁/重组诊断 | 测试配子的重组率、QTL单倍型频率 | 核验 Haldane 实现、连锁结构影响 |
| 实际寄主选择率 | `#(H_i属于植物种s)/#昆虫` | 寄主偏好及接受率 |
| 实际摄食网络权重 | `W_{h,p,t}=Σ_{i:s(i)=h,H_i为p种}x_i` | 加权二部互作网络，基于资源实际流量 |
| 寄主利用集中度 | `π_{h,p}=W_{h,p}/Σ_pW_{h,p}`、`exp[-Σ_pπ lnπ]` | 专化—广食的连续指标，零摄食时 NA |
| 空间占据/重叠 | 每格 P1、P2、昆虫数、同格共现 | 多代空间热图、群落斑块动态 |
| 资源守恒误差 | `max_p(max(0,X_p-R_p))` | 科学核查，不当生态效应 |
| 交配与出生数 | `Σ母方B_i`、`Σ父源贡献`、`no_compatible_mate` | 交配制约/生殖成功的分解 |

*请特别注意：目前固定物种 ID 的模型不能生成真正的新物种，`speciation_rate` 和“物种形成事件”不得作为真实输出。没有植物资源利用时，寄主利用集中度应输出 NA，不能被误判为绝对专食。*

## 8. 检验范围、可识别性与局限

1. **分析级别区分**：代码级验收（分离、连锁、守恒、同种交配）≠科学机制检验（何时产生寄主专化、物种共存）。前者必须在本版完成；后者待模型整体合理性审查完毕后预定义。
2. **概率依赖**：共享同一父本候选集合的母本，会因父本生存状态相关；不可将不同母本 D9 的无条件方差简单相加。
3. **基因效应截断**：`G9` 的 `[0,1]` 会使极端遗传值对应相同表型，需同时报告 raw 值和截断比例。
4. **化学功能缺口**：两个 D/I 与两个 Z/V/A 只是抽象通道；并不等同两种 PSM 或 VOC、真实酶动力学、P450/OR 的基因序列。实验论文支持化学专一性的重要性，不验证本模型的具体效用/取食函数 [R08–R11]。
5. **交配系统边界**：植物为严格异交的两性花理论类群，昆虫雌雄异体；对自交植物、单倍二倍体蜂类、孤雌生殖种类不直接适用。
6. **空间边界**：二维周期格、母源局部扩散、同格多株、局部母本产子折减均为简化；没有把真实面积、发育期、季节变化和种子库接入。

## 9. 论文写作时如何引用“公式来源”

| 公式或机制 | 引用 | 源头性质/不得误写 |
|---|---|---|
| G1/G2/G6：二倍体与双亲配子传递 | 遗传学基本规律；[R02,R03] | SLiM 是实现范式参考，二倍体生殖原则不是 SLiM 首创 |
| G3–G5：Haldane 标记图谱／SLiM 原生重组 | [R04]，SLiM 手册见版本说明 | 无交换干扰的标记分布；压缩边界不记录完整区间交换次数；图谱参数待定 |
| G7–G9：突变效应与加性 QTL | [R02,R05] | 模拟用统计/功能遗传模型；0.01和0.02不是文献实测值 |
| E1：PSM/VOC 遗传功能联系 | [R08,R09] | 实验支持遗传影响化学组成，**不支持线性混合式数值** |
| E2–E4：气味效用与选择 | [R07,R10] | Logit 数学形式+昆虫嗅觉作用的经验背景，函数参数是研究假设 |
| E5–E8：底物瓶颈与资源比例共享 | [R11,R14] + 本研究 H/I | R11 支持底物特异性，R14 提供最限制预算先例；共享资源比例分配仍为本项目假说 |
| E9–E11：生存与生殖成本 | 本研究 H，背景 [R08–R11] | **不得**写成参考文献直接报告的具体适合度公式 |
| S2：负密度依赖 | [R12,R13] | 类 Beverton–Holt 的局部变体，并非实测竞争系数 |
| D3–D4：花粉与异交配偶限制 | [R06] + 本研究 H/I | 植物交配系统有依据；半径/均匀选父是理想化 |
| D5–D10：Poisson、缺配偶及矩 | 基本概率推导 + 本研究 H | 混合分布计算是数学结果；生殖分布选择是理论设定 |


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

**[R14]** Speed MP, Fenton A, Jones MG, Ruxton GD, Brockhurst MA. 2015. *Coevolution can explain defensive secondary metabolite diversity in plants*. New Phytologist, 208:1251–1263. [DOI: 10.1111/nph.13560](https://doi.org/10.1111/nph.13560)。
**适用范围**：最限制抵抗通道、固定取食窗口预算与线性成本的模型先例；不验证本项目父本成本门槛、共享资源分配或具体参数。

**证据声明**：引用文献只对相应的“遗传规律、方法形式或经验机制”负责。`μ=0.01`、`σ_mut=0.02`、`λ=2`、`α=0.02/0.005`、化学效用系数、配偶半径、严格异交与 Poisson 出生的特定参数等，凡未直接从独立数据拟合，都必须标记 `theoretical_placeholder` 或 `implementation_choice`。这些假设不能因为列出参考文献而被写成实测事实。
