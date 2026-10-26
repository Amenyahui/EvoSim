# evoSim 完整模型公式与统计定义（中文）

基准：evoSim v0.1.1（sexual-qtl-1）；当前实验为 run4 的两个处理、各三个重复、200 次世代更新。本文对生成机制、运行时诊断、遗传统计及 analysis 目录中的科学分析指标逐项作数学展开，供 Methods 与 Supplementary Methods 选用。公式编号独立于《Methods_NC_中文初稿》。

“当前实验”表示该机制或设置在本批次中启用；“可选分支”表示代码支持但本批次未启用；“分析指标”表示脚本中的定义，不代表本次写作重新运行了该项分析。文件保存、哈希、绘图尺寸等工程操作不作为生物学公式。

## 1. 符号、状态与时间

令 \(t=0,\ldots,T\)，当前 \(T=200\)。\(\mathcal P_t,\mathcal H_t\) 分别为植物与昆虫个体集合，\(p\) 表示植物，\(i\) 表示昆虫，\(s,r\) 表示物种。用 \(K\in\{P,H\}\) 区分类群。物种身份写为 \(s(u)\)，空间坐标为 \(\boldsymbol x_u\)，唯一身份为 \(\mathrm{id}_u\)。\(\mathbf1(\cdot)\) 为指示函数。

\[
\mathcal P_{s,t}=\{p\in\mathcal P_t:s(p)=s\},\quad
\mathcal H_{s,t}=\{i\in\mathcal H_t:s(i)=s\},\quad
N^P_{s,t}=|\mathcal P_{s,t}|,\quad N^H_{s,t}=|\mathcal H_{s,t}|.
\tag{S1}
\]

植物性状为 \(\boldsymbol T_p^P=(D_{p1},D_{p2},Z_{p1},Z_{p2})\)，昆虫性状为 \(\boldsymbol T_i^H=(I_{i1},I_{i2},A_{i1},A_{i2})\)。\(D\) 为防御，\(Z\) 为独立气味性状，\(I\) 为抗性，\(A\) 为偏好模板。为避免混淆，后文用 \(\nu_u\) 表示条件平均出生数，用 \(\mu_\ell\) 表示突变概率，用 \(q_{ip}\) 表示取食匹配，用 \(\pi_u\) 表示配偶可获得概率。

定义截断函数：
\[
\operatorname{clip}_{[0,1]}(z)=\min\{1,\max(0,z)\},
\tag{S2}
\]
对向量逐坐标应用。除明确说明外，后文所有期望与概率都以当前代生态状态为条件。空集合之和为 0、空集合之积为 1；无有效观测的均值记为 NA，不记为 0。

## 2. 空间、初始化与创始遗传变异

### 2.1 周期距离与随机位置

在 \(L_1\times L_2\) 栅格中，各坐标为 \(0,\ldots,L_a-1\)。周期切比雪夫距离为
\[
d_{\mathrm{per}}(\boldsymbol x,\boldsymbol y)=
\max_{a=1,2}\min\{|x_a-y_a|,L_a-|x_a-y_a|\}.
\tag{S3}
\]

若未指定初始坐标，则
\[
X_{u,a}(0)\sim\operatorname{DiscreteUniform}\{0,\ldots,L_a-1\},
\qquad a=1,2.
\tag{S4}
\]
当前 \(L_1=L_2=10\)，不同个体可占据相同位置。初始 \(N^P_{1,0}=N^P_{2,0}=200,\ N^H_{1,0}=40\)。

### 2.2 创始基因型

每类群具有 \(L_K=8\) 个位点，基因型数组形状为 \(N_K\times2\times L_K\)。令 \(f^0_{s\ell a}\) 为配置中的创始等位基因频率，则
\[
G^K_{u,h,\ell}(0)\sim
\operatorname{Categorical}\left(
\left\{\frac{f^0_{s\ell a}}{\sum_b f^0_{s\ell b}}\right\}_{a}
\right),\qquad h=1,2.
\tag{S5}
\]
同源拷贝与位点分别抽样；当前每物种每位点有两个可用等位基因，各概率 0.5，有限种群实际频率不必恰好等于 0.5。

以 \(\chi_i=1\) 表示雌性、0 表示雄性：
\[
\chi_i(0)\sim\operatorname{Bernoulli}(p_f^0),\qquad p_f^0=0.5.
\tag{S6}
\]
这是本文的性别记号；代码使用 sex=0 表示雌性、sex=1 表示雄性。

### 2.3 基因型—表型映射

令 \(e_{a,q}\) 为等位基因 \(a\) 对性状坐标 \(q\) 的效应，\(b_q^K\) 为类群背景：
\[
T^{K,\mathrm{raw}}_{u,q}
=b_q^K+\frac12\sum_{\ell=1}^{L_K}\sum_{h=1}^{2}
e_{G^K_{u,h,\ell},q},
\qquad
T^K_{u,q}=\operatorname{clip}_{[0,1]}(T^{K,\mathrm{raw}}_{u,q}).
\tag{S7}
\]

当前背景均为零，每个坐标由两个位点贡献。低效应位点的两个等位效应为 0.10、0.12，高效应位点为 0.30、0.32。由式（S5）—（S7）得到创始性状期望：
\[
\begin{aligned}
E[\boldsymbol T^P_{P1}(0)]&=(0.22,0.62,0.22,0.62),\\
E[\boldsymbol T^P_{P2}(0)]&=(0.62,0.22,0.62,0.22),\\
E[\boldsymbol T^H_{H1}(0)]&=(0.22,0.62,0.22,0.62).
\end{aligned}
\tag{S8}
\]
这是根据配置推得的期望，不是每次运行中固定的实际均值。

## 3. 气味生成、候选遭遇与寄主选择

### 3.1 发射气味与感知线索

采用列向量记号：
\[
\boldsymbol V_p^{\mathrm{emit}}
=\operatorname{clip}_{[0,1]}
\{\rho\boldsymbol B\boldsymbol D_p+(1-\rho)\boldsymbol Z_p\}.
\tag{S9}
\]
当前 \(\rho=0.5,\boldsymbol B=\boldsymbol I_2\)。当 \(\rho=1\) 时，映射退化为 \(\operatorname{clip}(\boldsymbol B\boldsymbol D_p)\)。

感知线索为
\[
\boldsymbol V_p^{\mathrm{cue}}=
\begin{cases}
\boldsymbol V_p^{\mathrm{emit}},&\text{no permutation (current experiment)},\\
\boldsymbol V_{\Pi_t(p)}^{\mathrm{emit}},&\text{odor permutation branch}.
\end{cases}
\tag{S10}
\]
其中 \(\Pi_t\) 为当代全部植物的随机排列，可跨物种置换。置换仅改变选择线索，不改变实际防御、植物身份和发射气味成本。

### 3.2 候选集

\[
\mathcal A_i=\{p\in\mathcal P_t:
d_{\mathrm{per}}(\boldsymbol x_i,\boldsymbol x_p)\le r_{\mathrm{search}}\},
\qquad
m_i=\min(M,|\mathcal A_i|).
\tag{S11}
\]
当前 \(r_{\mathrm{search}}=2,\ M=8\)。从 \(\mathcal A_i\) 无放回均匀抽取 \(m_i\) 株：
\[
P(\mathcal C_i=C\mid\mathcal A_i)
={|\mathcal A_i|\choose m_i}^{-1},
\quad C\subseteq\mathcal A_i,\ |C|=m_i.
\tag{S12}
\]
所有个体的候选池基于相同代初状态，不因其他昆虫已选择某植物而缩减。

### 3.3 气味匹配得分

\[
\delta_{ip}^2=\sum_{k=1}^{2}w_k(V^{\mathrm{cue}}_{pk}-A_{ik})^2,
\qquad
u_{ip}=b-\frac{\delta_{ip}^2}{2\sigma^2}.
\tag{S13}
\]
当前 \(w_1=w_2=0.5,\sigma=0.2,b=1\)。中性嗅觉分支将所有候选得分设为 \(u_{ip}=b\)，当前未启用。

### 3.4 外部选项与选择概率

记 \(J_i=0\) 为未选择寄主，\(J_i=p\) 为选择植物 \(p\)：
\[
P(J_i=p\mid\mathcal C_i)=
\frac{\exp(\tau u_{ip})}{\exp(u_0)+\sum_{j\in\mathcal C_i}\exp(\tau u_{ij})},
\quad
P(J_i=0\mid\mathcal C_i)=
\frac{\exp(u_0)}{\exp(u_0)+\sum_{j\in\mathcal C_i}\exp(\tau u_{ij})}.
\tag{S14}
\]
当前 \(\tau=1,u_0=0\)。注意代码中 \(u_0\) 已是外部选项的 logit，不再乘 \(\tau\)。空候选集强制 \(J_i=0\)，不进行选择抽样。

为稳定计算，令
\[
\ell_{ip}=\tau u_{ip},\quad\ell_{i0}=u_0,\quad
a_i=\max_{j\in\mathcal C_i\cup\{0\}}\ell_{ij},\qquad
P(J_i=j)=\frac{\exp(\ell_{ij}-a_i)}
{\sum_k\exp(\ell_{ik}-a_i)}.
\tag{S15}
\]
这是式（S14）的等价数值形式。代码还针对极端参数采用先中心化再缩放，避免中间值溢出。

在强制选择分支中，仅在寄主间归一化：
\[
P^{\mathrm{cond}}_{ip}
=\frac{\exp(\tau u_{ip})}{\sum_{j\in\mathcal C_i}\exp(\tau u_{ij})},
\qquad P(J_i=0)=0.
\tag{S16}
\]
当前实际选择采用式（S14），但式（S16）仍用于条件偏好诊断。当 \(\tau=0\) 或嗅觉中性时，条件偏好为 \(1/m_i\)；在外部选项模型中，\(\tau=0\) 时每候选概率为 \(1/(m_i+\exp u_0)\)。

### 3.5 候选暴露与条件物种偏好

对非空候选集：
\[
a_{ir}=\frac{\sum_{p\in\mathcal C_i}\mathbf1(s(p)=r)}{m_i},
\qquad
c_{ir}=\sum_{p\in\mathcal C_i:s(p)=r}P^{\mathrm{cond}}_{ip}.
\tag{S17}
\]
候选为空时二者为 NA。昆虫物种 \(s\) 对植物物种 \(r\) 的累计候选暴露次数为
\[
E_{sr}=\sum_{i\in\mathcal H_{s,t}}
\sum_{p\in\mathcal C_i}\mathbf1(s(p)=r).
\tag{S18}
\]
这是候选出现次数，同一植株出现在不同昆虫的候选集中可被多次计数。

## 4. 防御匹配、共享资源与竞争诊断

### 4.1 取食能力与需求

对已匹配昆虫：
\[
r_{ipk}=
\begin{cases}
I_{ik}/D_{pk},&D_{pk}>0,\\
+\infty,&D_{pk}=0,
\end{cases}
\qquad
q_{ip}=\min\{1,\min_k r_{ipk}\}.
\tag{S19}
\]
全零防御对应 \(q_{ip}=1\)。未匹配昆虫的 q 字段为 NA，而非零。

\[
d_i=
\begin{cases}
d_i^{\max}q_{iJ_i},&J_i\ne0,\\
0,&J_i=0.
\end{cases}
\tag{S20}
\]

### 4.2 按比例分配

令 \(\mathcal H_p=\{i:J_i=p\}\)，单株资源为 \(R_p>0\)：
\[
D_p^{\mathrm{tot}}=\sum_{i\in\mathcal H_p}d_i,\qquad
\gamma_p=
\begin{cases}
1,&D_p^{\mathrm{tot}}\le R_p,\\
R_p/D_p^{\mathrm{tot}},&D_p^{\mathrm{tot}}>R_p.
\end{cases}
\tag{S21}
\]

\[
x_i=\begin{cases}\gamma_{J_i}d_i,&J_i\ne0,\\0,&J_i=0,\end{cases}
\quad
X_p=\sum_{i\in\mathcal H_p}x_i,\quad
g_i=\frac{x_i}{d_i^{\max}},\quad h_p=\frac{X_p}{R_p}.
\tag{S22}
\]

精确算术下：
\[
X_p=\min(D_p^{\mathrm{tot}},R_p),\quad
0\le g_i,h_p\le1,\quad
\sum_pX_p=\sum_ix_i.
\tag{S23}
\]
当前 \(R_p=d_i^{\max}=1\)。资源不足和近似耗尽诊断分别为
\[
L_p^{\mathrm{limited}}=\mathbf1(D_p^{\mathrm{tot}}>R_p),\qquad
L_p^{\mathrm{exhausted}}=\mathbf1\{X_p\ge(1-10^{-12})R_p\}.
\tag{S24}
\]

若浮点累加产生 \(X_p>R_p\)，代码对该植物所有取食者统一乘
\[
\eta_p=\operatorname{nextafter}(R_p/X_p,0),\quad
x_i\leftarrow\eta_p x_i\ (i\in\mathcal H_p),\quad
\gamma_p\leftarrow\eta_p\gamma_p,
\tag{S25}
\]
随后重算 X，至多修正三次。最终要求 \(X_p\le R_p\)，并以相对容差 \(10^{-12}\)、绝对容差 0 检查式（S23）的第一项。

### 4.3 同种与异种竞争负荷

令
\[
n_{sp}=\sum_{i\in\mathcal H_p}\mathbf1(s(i)=s),\quad
D_{sp}=\sum_{i\in\mathcal H_p:s(i)=s}d_i,\quad
n_p=\sum_sn_{sp}.
\tag{S26}
\]
对选择植物 p 的昆虫 i：
\[
\begin{aligned}
n_i^{\mathrm{same}}&=n_{s(i),p}-1,&
n_i^{\mathrm{other}}&=n_p-n_{s(i),p},\\
d_i^{\mathrm{same}}&=\max(0,D_{s(i),p}-d_i),&
d_i^{\mathrm{other}}&=\max(0,D_p^{\mathrm{tot}}-D_{s(i),p}).
\end{aligned}
\tag{S27}
\]
未匹配个体上述四量均为零。这些是资源函数返回的诊断量，不能假定全部都已写入当前批次的事件 CSV。

## 5. 存活、成本、密度与交配

### 5.1 成本及繁殖因子

\[
C_p^D=c_D\sum_kD_{pk},\qquad
C_p^V=c_V\frac1n\sum_kV^{\mathrm{emit}}_{pk},\qquad
C_i^I=c_I\sum_kI_{ik}.
\tag{S28}
\]
这里 \(n=2\)。VOC 成本使用发射气味，不使用置换后的感知线索。
\[
F_p^P=\max(0,1-C_p^D-C_p^V),\qquad
F_i^H=\max(0,1-C_i^I).
\tag{S29}
\]

\[
S_p^P=s_{P0}+(1-s_{P0})(1-h_p),\quad
S_i^H=g_i,\quad
Y_u^K\sim\operatorname{Bernoulli}(S_u^K).
\tag{S30}
\]
\(Y=1\) 表示存活到交配阶段。当前 \(s_{P0}=0.5\)，无寄主昆虫 \(S_i^H=0\)。

输出中的适合度分量为
\[
W_p^P=S_p^PF_p^P,\qquad W_i^H=S_i^HF_i^H.
\tag{S31}
\]
W 未包括配偶概率、繁殖系数及植物密度因子，不等同于真实的期望子代数。

### 5.2 局部植物竞争

\[
N_r(\mathcal N_p)=\sum_{v\in\mathcal P_t}
\mathbf1(s(v)=r)\mathbf1\{d_{\mathrm{per}}(\boldsymbol x_p,\boldsymbol x_v)\le r_c\},
\quad
\phi_p=\left[1+\sum_r\alpha_{s(p),r}N_r(\mathcal N_p)\right]^{-1}.
\tag{S32}
\]
使用代初数量并包含焦点植物，周期邻域内每个栅格只计一次。当前 \(r_c=1,\alpha_{ss}=0.02,\alpha_{sr}=0.005\ (s\ne r)\)。昆虫 \(\phi_i=1\)。记录的植物物种密度因子为该物种所有个体 \(\phi_p\) 的算术均值。

### 5.3 母方、父方资格与位置

定义交配前的母方角色指示 \(M_u\)、潜在父方资格 \(D_u^{\mathrm{donor}}\) 及交配位置：
\[
\begin{array}{lll}
M_p=1,&D_p^{\mathrm{donor}}=\mathbf1(F_p^P>0),&
\boldsymbol z_p=\boldsymbol x_p,\\
M_i=\mathbf1(J_i\ne0)\chi_i,&
D_i^{\mathrm{donor}}=\mathbf1(J_i\ne0)(1-\chi_i)\mathbf1(F_i^H>0),&
\boldsymbol z_i=
\begin{cases}\boldsymbol x_{J_i},&J_i\ne0,\\\boldsymbol x_i,&J_i=0.\end{cases}
\end{array}
\tag{S33}
\]
未匹配昆虫保留的坐标不使其获得交配资格。

对同一类群的个体 u，潜在父方池为
\[
\mathcal M_u^0=
\{v\ne u:s(v)=s(u),\
D_v^{\mathrm{donor}}=1,\
d_{\mathrm{per}}(\boldsymbol z_u,\boldsymbol z_v)\le r_m^K\},
\tag{S34}
\]
仅在 \(M_u=1\) 时构建，否则为空。当前 \(r_m^P=1,r_m^H=2\)。排除自身实现植物严格异交。

存活后的实际配偶池与母方可繁殖事件为
\[
\mathcal M_u=\{v\in\mathcal M_u^0:Y_v=1\},\qquad
Z_u=M_uY_u\mathbf1(|\mathcal M_u|>0).
\tag{S35}
\]
母方 F=0 时即使 Z=1，其条件出生均值仍为零。父方则必须 F>0。

### 5.4 配偶概率、出生混合分布、期望和方差

由于存活抽样在给定生态状态下独立，且配偶池不含自身：
\[
\pi_u=P(|\mathcal M_u|>0)=
1-\prod_{v\in\mathcal M_u^0}(1-S_v),\qquad
p_u=S_u\pi_u\quad(M_u=1).
\tag{S36}
\]
空配偶池给出 \(\pi_u=0\)。这正是代码 mate_probability 的定义。

\[
\nu_p=\lambda_{s(p)}^PF_p^P\phi_p,\qquad
\nu_i=\lambda_{s(i)}^HF_i^H.
\tag{S37}
\]

\[
B_u=
\begin{cases}
\operatorname{Poisson}(\nu_u),&Z_u=1,\\
0,&Z_u=0.
\end{cases}
\tag{S38}
\]
当前 \(\lambda^P=2,\lambda^H=4\)。昆虫 λ 是每雌的系数。存活、配偶缺失通过 Z 处理，不在 \(\nu\) 中重复乘入。

对具有母方角色的个体：
\[
P(B_u=0)=1-p_u+p_u e^{-\nu_u},\qquad
P(B_u=b)=p_u e^{-\nu_u}\frac{\nu_u^b}{b!},\quad b\ge1.
\tag{S39}
\]

\[
E[B_u]=p_u\nu_u,\qquad
\operatorname{Var}(B_u)=p_u\nu_u+p_u(1-p_u)\nu_u^2.
\tag{S40}
\]
式（S39）是对代码抽样过程的数学展开，式（S40）由 maternal_moments 直接计算。对雄虫等非母方角色，实际母方 B 恒为零，但代码将母方期望与方差字段记为 NA，以区分“不适用”和“适用但不能繁殖”。

不同母方可能共享配偶，导致其繁殖事件相关。因此，不能直接把个体方差相加并宣称为总种群出生方差；代码未实现这一协方差分析。

### 5.5 父方抽样与繁殖贡献

对母方 u 的每个子代 o：
\[
P(F_o=v\mid\mathcal M_u,B_u)=|\mathcal M_u|^{-1},
\quad v\in\mathcal M_u,\qquad M_o=u.
\tag{S41}
\]
按子代有放回抽样。父方 v 的实际贡献数为
\[
B_v^{\mathrm{sired}}=\sum_o\mathbf1(F_o=v).
\tag{S42}
\]
植物既可作为母方也可作为父方，父方贡献不能再次加入总出生数。

## 6. 配子、受精、突变与扩散

### 6.1 连锁群内传递

对于父母任一方 v、某条连锁群 c，按遗传位置排序位点 \(\ell_1,\ldots,\ell_m\)。起始同源拷贝选择：
\[
H_{v,c,1}\sim\operatorname{Bernoulli}(1/2),\qquad
M_{v,\ell_1}=G_{v,H_{v,c,1},\ell_1}.
\tag{S43}
\]
此处 H 采用 0/1 拷贝索引；不同连锁群独立。

\[
d_j=\frac{x_{\ell_{j+1}}^{\mathrm{cM}}-x_{\ell_j}^{\mathrm{cM}}}{100},\quad
K_j\sim\operatorname{Poisson}(d_j),\quad
H_{v,c,j+1}=H_{v,c,j}\oplus(K_j\bmod2),\quad
M_{v,\ell_{j+1}}=G_{v,H_{v,c,j+1},\ell_{j+1}}.
\tag{S44}
\]
\(\oplus\) 表示异或。相邻标记重组分数为
\[
r(d)=P(K\text{ is odd})=\frac{1-e^{-2d}}2.
\tag{S45}
\]
当前每个位点独占一个连锁群，因此只执行独立分离，不产生可观测的群内多位点交换。

### 6.2 子代双亲组合

\[
G^{\mathrm{pre}}_{o,1,\ell}=M_{M_o,\ell},\qquad
G^{\mathrm{pre}}_{o,2,\ell}=M_{F_o,\ell}.
\tag{S46}
\]
母父形成配子的随机流彼此分开。双亲必须为不同的同种个体。

### 6.3 突变与新等位效应

令 \(\epsilon_K\in\{0,1\}\) 指示该类群是否开启新突变：
\[
U_{o,h,\ell}\sim\operatorname{Bernoulli}(\epsilon_K\mu_\ell),
\quad
\Delta\boldsymbol e_{o,h,\ell}\sim
\mathcal N\{\boldsymbol0,\operatorname{diag}(\boldsymbol\sigma_\ell^2)\}.
\tag{S47}
\]

\[
G_{o,h,\ell}=
\begin{cases}
G^{\mathrm{pre}}_{o,h,\ell},&U_{o,h,\ell}=0,\\
a_{\mathrm{new}},&U_{o,h,\ell}=1,
\end{cases}
\quad
\boldsymbol e_{a_{\mathrm{new}}}
=\boldsymbol e_{G^{\mathrm{pre}}_{o,h,\ell}}+\Delta\boldsymbol e_{o,h,\ell}.
\tag{S48}
\]
每次突变都创建新 ID，保留原等位效应。当前指定目标坐标的标准差为 0.02，其余坐标为零，突变率均为 0.01。效应自身不截断，最终表型按式（S7）截断。

每个子代的新突变数及至少一次突变的概率为
\[
M_o^{\mathrm{new}}=\sum_{h=1}^2\sum_{\ell=1}^{L_K}U_{o,h,\ell},
\quad E[M_o^{\mathrm{new}}]=2\epsilon_K\sum_\ell\mu_\ell,
\quad
P(M_o^{\mathrm{new}}>0)=1-\prod_\ell(1-\epsilon_K\mu_\ell)^2.
\tag{S49}
\]
开启处理下，当前均匀突变率给出 \(E[M_o^{\mathrm{new}}]=0.16\)，
\(P(M_o^{\mathrm{new}}>0)=1-0.99^{16}\approx0.14854\)。该概率不等于出现可观测表型变化的概率，因为存在表型截断。

### 6.4 子代性别与空间扩散

\[
\chi_o\sim\operatorname{Bernoulli}(p_f^{\mathrm{offspring}}),
\qquad p_f^{\mathrm{offspring}}=0.5.
\tag{S50}
\]

\[
\boldsymbol x_o^{\mathrm{origin}}=
\begin{cases}
\boldsymbol x_{M_o},&K=P,\\
\boldsymbol x_{J_{M_o}},&K=H,
\end{cases}
\quad
\Delta x_{o,a}\sim\operatorname{DiscreteUniform}\{-r_d^K,\ldots,r_d^K\},
\quad
x_{o,a}=(x_{o,a}^{\mathrm{origin}}+\Delta x_{o,a})\bmod L_a.
\tag{S51}
\]
当前 \(r_d^P=r_d^H=1\)，两坐标独立，包含留在原斑块的可能性。

## 7. 世代替换、灭绝与人口记录

\[
N^K_{s,t+1}=\sum_{u\in\mathcal K_{s,t}}B_u,\qquad
\mathcal K_{t+1}=\{\text{all offspring of taxon }K\text{ produced in this update}\}.
\tag{S52}
\]
亲代不保留，没有迁入与补种。于是
\[
N^K_{s,t}=0\Longrightarrow N^K_{s,t'}=0\quad(t'>t).
\tag{S53}
\]

逐物种记录：
\[
\begin{aligned}
N_s^{\mathrm{surv}}&=\sum_{u:s(u)=s}Y_u,\\
N_s^{\mathrm{mated\ mothers}}&=\sum_{u:s(u)=s}Z_u,\\
N_s^{\mathrm{producing\ mothers}}&=\sum_{u:s(u)=s}\mathbf1(B_u>0),\\
N_s^{\mathrm{producing\ fathers}}&=\sum_{u:s(u)=s}\mathbf1(B_u^{\mathrm{sired}}>0),\\
N_s^{\mathrm{zero\ fitness}}&=\sum_{u:s(u)=s}\mathbf1(W_u=0),\\
N_s^{\mathrm{female}}&=\sum_{i\in\mathcal H_{s,t}}\chi_i,\qquad
N_s^{\mathrm{male}}=N^H_{s,t}-N_s^{\mathrm{female}}.
\end{aligned}
\tag{S54}
\]

最后一代仅观测，不再更新。更新指示为
\[
\mathcal U_t=\mathbf1(t<T)\mathbf1(|\mathcal P_t|+|\mathcal H_t|>0).
\tag{S55}
\]
\(\mathcal U_t=0\) 时不执行繁殖抽样，末代相关出生、存活抽样字段为空，不应解释为这些随机变量已被抽样为零。

## 8. 物种相互作用网络与空间统计

### 8.1 接受数、资源流与每虫资源量

\[
A_{sr}(t)=\sum_{i\in\mathcal H_{s,t}}\mathbf1(J_i\in\mathcal P_{r,t}),\quad
Q_{sr}(t)=\sum_{i\in\mathcal H_{s,t}}x_i\mathbf1(J_i\in\mathcal P_{r,t}),\quad
L_{sr}(t)=\frac{Q_{sr}(t)}{N^H_{s,t}}.
\tag{S56}
\]
当 \(N^H_{s,t}=0\) 时 L=NA。Q=0 不等于从未选择过该植物，可能发生“已选择但取食能力为零”。

### 8.2 资源网络、丰富度与连通度（分析指标）

\[
\mathcal E_t=\{(s,r):Q_{sr}(t)>0\},\quad
S_P(t)=\sum_r\mathbf1(N^P_{r,t}>0),\quad
S_H(t)=\sum_s\mathbf1(N^H_{s,t}>0).
\tag{S57}
\]

\[
C_{\mathrm{initial}}(t)=\frac{|\mathcal E_t|}{S_P(0)S_H(0)},\qquad
C_{\mathrm{extant}}(t)=\frac{|\mathcal E_t|}{S_P(t)S_H(t)}.
\tag{S58}
\]
第一项代码分母为配置中初始物种槽位数之积，当前各物种初始数量均大于零，因此与上式一致。第二项分母为零时记为 NA。

\[
\mathrm{Turnover}(t)=1-
\frac{|\mathcal E_{t-1}\cap\mathcal E_t|}
{|\mathcal E_{t-1}\cup\mathcal E_t|}.
\tag{S59}
\]
t=0 或并集为空时为 NA。相邻均为空的网络不被赋予零更替率。

### 8.3 资源比例、首次使用和占据格数量（分析指标）

\[
Q_{\mathrm{all}}(t)=\sum_{s,r}Q_{sr}(t),\qquad
f_{P2}(t)=\frac{\sum_sQ_{s,P2}(t)}{Q_{\mathrm{all}}(t)},\qquad
t_{P2}^{\mathrm{first}}=\min\{t:\exists s,\ Q_{s,P2}(t)>0\}.
\tag{S60}
\]
总资源量为零时 f=NA；从未使用 P2 时首次使用时间为空。

\[
O^K_s(a,b,t)=
\sum_{u\in\mathcal K_{s,t}}\mathbf1\{\boldsymbol x_u=(a,b)\}.
\tag{S61}
\]
这是空间图每格个体数；初始昆虫位置不等于已发生的取食连边。

## 9. 固定参考面板、标准化偏好与寄主宽度

当前参考性状固定为 P1 的 \((0.2,0.6)\)、P2 的 \((0.6,0.2)\)，D 与 Z 相同，均不是有限创始种群的实际平均值。

\[
\boldsymbol V_r^{\mathrm{ref}}=
\operatorname{clip}\{\rho\boldsymbol B\boldsymbol D_r^{\mathrm{ref}}
+(1-\rho)\boldsymbol Z_r^{\mathrm{ref}}\},\quad
u_{ir}^{\mathrm{ref}}=
b-\frac{\sum_kw_k(V_{rk}^{\mathrm{ref}}-A_{ik})^2}{2\sigma^2},
\quad
p_{ir}^{\mathrm{ref}}=
\frac{\exp(\tau u_{ir}^{\mathrm{ref}})}
{\sum_v\exp(\tau u_{iv}^{\mathrm{ref}})}.
\tag{S62}
\]
各参考物种提供一个等权选项，不含外部选项，也不加现实丰度权重。

\[
\bar p_{sr}^{\mathrm{ref}}=
\frac1{N^H_{s,t}}\sum_{i\in\mathcal H_{s,t}}p_{ir}^{\mathrm{ref}},
\quad
q_{ir}^{\mathrm{ref}}=\min\left\{1,\min_{k:D_{rk}^{\mathrm{ref}}>0}
\frac{I_{ik}}{D_{rk}^{\mathrm{ref}}}\right\},
\quad
\bar q_{sr}^{\mathrm{ref}}=
\frac1{N^H_{s,t}}\sum_{i\in\mathcal H_{s,t}}q_{ir}^{\mathrm{ref}}.
\tag{S63}
\]
全零防御仍规定 q=1，昆虫灭绝则均值为 NA。该诊断不抽样寄主、不分配资源、不改变模拟状态。

物种 s 的有效寄主数采用群体平均偏好的熵指数：
\[
H_s^{\mathrm{ref}}=
-\sum_{r:\bar p_{sr}^{\mathrm{ref}}>0}
\bar p_{sr}^{\mathrm{ref}}\ln\bar p_{sr}^{\mathrm{ref}},
\qquad
B_s^{\mathrm{eff}}=\exp(H_s^{\mathrm{ref}}).
\tag{S64}
\]
它不是“个体有效寄主数的均值”。例如偏好互相相反的专化个体混合后，也可能得到较大的群体有效寄主数。

宏观脚本再对存续昆虫物种等权平均：
\[
\bar B^{\mathrm{eff}}=
\frac1{S_H(t)}\sum_{s:N^H_{s,t}>0}B_s^{\mathrm{eff}},
\qquad
\bar q^{\mathrm{macro}}=
\frac1{S_H(t)}\sum_{s:N^H_{s,t}>0}
\left(\frac1{S_P^{\mathrm{ref}}}\sum_r\bar q_{sr}^{\mathrm{ref}}\right).
\tag{S65}
\]
不是按昆虫丰度加权；当前仅一个昆虫物种。全部昆虫灭绝时为 NA。

## 10. 遗传多样性、表型和图形统计

### 10.1 等位频率与杂合度

\[
f_{s\ell a}(t)=
\frac{\sum_{u\in\mathcal K_{s,t}}\sum_{h=1}^{2}
\mathbf1(G_{u,h,\ell}=a)}{2N^K_{s,t}}.
\tag{S66}
\]
对于非空种群，目录中已有但该种群不存在的等位基因频率为零；整个种群灭绝时频率为 NA。

\[
H^O_{s\ell}(t)=\frac1{N^K_{s,t}}
\sum_{u\in\mathcal K_{s,t}}\mathbf1(G_{u,1,\ell}\ne G_{u,2,\ell}),
\qquad
H^E_{s\ell}(t)=1-\sum_af_{s\ell a}(t)^2.
\tag{S67}
\]
\(H^E\) 没有有限样本无偏校正，不应称为已校正估计。两项在 N=0 时均为 NA。

### 10.2 连锁群内两位点单倍型频率

对同一连锁群上的位点 \(\ell,m\)：
\[
f_{s,\ell m,ab}^{\mathrm{hap}}=
\frac1{2N^K_{s,t}}\sum_{u\in\mathcal K_{s,t}}\sum_{h=1}^{2}
\mathbf1(G_{u,h,\ell}=a,\ G_{u,h,m}=b).
\tag{S68}
\]
代码保存观察到的等位组合；当前每群一个位点，因此该表没有有效的两位点组合。代码未由此计算 \(D'\) 或 \(r^2\)，本文也不把它们写为已有指标。

### 10.3 表型均值、方差与截断比例

\[
\bar T_{sq}=\frac1{N_s}\sum_{u:s(u)=s}T_{uq},\quad
V_{sq}=\frac1{N_s}\sum_{u:s(u)=s}(T_{uq}-\bar T_{sq})^2,\quad
\mathrm{SD}_{sq}=\sqrt{V_{sq}}.
\tag{S69}
\]
方差分母为 N，不是 N−1。

\[
\bar T^{\mathrm{raw}}_{sq}=\frac1{N_s}\sum_{u:s(u)=s}T^{\mathrm{raw}}_{uq},\qquad
f_{sq}^{\mathrm{clip}}=\frac1{N_s}
\sum_{u:s(u)=s}\mathbf1(T^{\mathrm{raw}}_{uq}<0\ \lor\ T^{\mathrm{raw}}_{uq}>1).
\tag{S70}
\]
恰好为 0 或 1 不计入 clipped_fraction。种群为空时以上指标为 NA。

图中阴影为
\[
[\bar T_{sq}-\mathrm{SD}_{sq},\ \bar T_{sq}+\mathrm{SD}_{sq}],
\tag{S71}
\]
这是种群内变异，不是均值置信区间，绘图函数未把阴影重新截断到 [0,1]。

### 10.4 突变与等位图筛选

\[
M_{s\ell}(t)=\#\{e:e.\mathrm{species}=s,\ e.\mathrm{locus}=\ell,\
e.\mathrm{generation}=t\},\quad
M_{\mathrm{run}}=\#\{\text{all mutation events in this run}\},\quad
A_{\mathrm{catalog}}=\#\{\text{founder and newly added allele entries}\}.
\tag{S72}
\]
位点突变数、杂合度等在遗传输出中会随每条等位基因记录重复出现，绘图整理时按“运行—世代—物种—位点”去重，不能逐等位行求和。

每物种每位点的显示排序量为
\[
R_{s\ell a}^{\mathrm{plot}}=\max_{t\ \mathrm{observed}}f_{s\ell a}(t).
\tag{S73}
\]
选取峰值最高的至多 12 个等位基因，且峰值必须大于零；不重新归一化保留曲线，完整频率数据不被删除。

## 11. 灭绝、运行终点与描述性汇总

令 \(\mathcal T_{\mathrm{obs}}\) 为实际记录的世代：
\[
t_s^{\mathrm{ext}}=\min\{t\in\mathcal T_{\mathrm{obs}}:N_{s,t}=0\},\quad
N_s^{\mathrm{end}}=N_{s,\max\mathcal T_{\mathrm{obs}}},\quad
N_s^{\mathrm{peak}}=\max_{t\in\mathcal T_{\mathrm{obs}}}N_{s,t}.
\tag{S74}
\]
无零数量观测时灭绝时间为空。对完成或全部灭绝的有效运行，宏观脚本记录
\[
\delta_s=\mathbf1(t_s^{\mathrm{ext}}\text{ observed}),\qquad
\tilde t_s=
\begin{cases}
t_s^{\mathrm{ext}},&\delta_s=1,\\
T,&\delta_s=0,
\end{cases}
\qquad
\mathrm{right\_censored}_s=1-\delta_s.
\tag{S75}
\]
这是事件/右删失编码，不意味着已拟合生存模型。计算中止或参数无效的运行不进入该脚本的有效生物学汇总。

\[
t^{\mathrm{first\ loss}}=\min_{s:\delta_s=1}t_s^{\mathrm{ext}},\qquad
I^{\mathrm{retained}}=
\mathbf1\{S_P(t_{\mathrm{end}})=S_P(0),\
S_H(t_{\mathrm{end}})=S_H(0)\}.
\tag{S76}
\]
无物种损失时 first loss=NA。总灭绝时间仅在运行状态为 all_species_extinct 时定义。

对于同一处理、同一世代的重复均值，采用
\[
\bar y_{c,t}=
\frac1{|\mathcal R_{c,t}^{\mathrm{valid}}|}
\sum_{r\in\mathcal R_{c,t}^{\mathrm{valid}}}y_{c,r,t}.
\tag{S77}
\]
此式描述空间汇总脚本的 pandas 均值口径：只使用有记录且非缺失的值，若集合为空则为 NA。灭绝种群数量零仍参与种群数量均值；灭绝后的偏好或取食比例 NA 不以零填充。当前遗传图主要按重复分开展示，不报告重复间误差条。

## 12. 初始条件诊断公式

初始诊断脚本在 t=0 汇总：
\[
R_{\mathrm{tot}}=\sum_pR_p,\quad
D_{\max,\mathrm{tot}}=\sum_id_i^{\max},\quad
D_{\mathrm{potential}}=\sum_id_i,\quad
X_{\mathrm{tot}}=\sum_ix_i,\quad
\bar h=\frac1{N_P}\sum_ph_p,\quad
\bar\phi=\frac1{N_P}\sum_p\phi_p.
\tag{S78}
\]

\[
\bar B_P^{\mathrm{expected}}=\frac1{N_P}\sum_pS_p\pi_p\nu_p,\quad
\bar B_H^{\mathrm{maternal}}=
\frac1{|\mathcal H^{M}|}\sum_{i\in\mathcal H^{M}}S_i\pi_i\nu_i,\quad
\bar\nu_P^{\mathrm{no\ damage}}=\frac1{N_P}\sum_p\lambda_{s(p)}^PF_p^P\phi_p,
\tag{S79}
\]
其中 \(\mathcal H^{M}=\{i:M_i=1\}\)，与忽略非母方 NA 的代码均值相符。名称为 expected_without_damage 的最后一项实际上同时不含存活及配偶可获得概率，不能解读为“只移除受食、仍保留其他全部过程”的完整反事实期望。初始诊断脚本还记录全类群个体的 \(\bar S,\bar F\)，均按个体算术平均。

## 13. 代码中的非当前分支与辅助关系

### 13.1 全局密度函数

\[
\phi_s^{\mathrm{global}}=
\left(1+\sum_r\alpha_{sr}N^P_{r,t}\right)^{-1}.
\tag{S80}
\]
该辅助函数存在于代码，空间流程随后使用局部式（S32）覆盖其作用于个体的结果，不应把两者相乘。

### 13.2 理想单物种承载量换算

代码辅助函数给出
\[
\alpha=\frac{\lambda-1}{K}.
\tag{S81}
\]
对应简化确定性关系 \(N_{t+1}=\lambda N_t/(1+\alpha N_t)\) 的正平衡 \(K=(\lambda-1)/\alpha\)，要求 \(\lambda>1,K>0\)。该关系忽略存活损失、配偶限制和空间异质性，不是当前有性模型的实际承载量公式。

### 13.3 继承的昆虫基础存活表达式

资源函数内保留
\[
W_i^{\mathrm{legacy}}=
\begin{cases}
[s_{I0}+(1-s_{I0})g_i]F_i^H,&J_i\ne0,\\
0,&J_i=0.
\end{cases}
\tag{S82}
\]
当前生态子步骤明确覆盖为 \(W_i=g_iF_i^H\)，且当前 \(s_{I0}=0\)。不能据辅助函数声称当前专性寄主模型允许无取食存活。

## 14. 计算边界与随机流

### 14.1 安全阈值

令 \(B_K=\sum_{u\in\mathcal K_t}B_u\)，已保留基因型标量数为 C，等位目录大小为 A。预检要求：
\[
\begin{aligned}
B_P+B_H&\le5000,\\
C+\sum_K2L_KB_K&\le4{,}000{,}000,\\
A+\sum_K2B_K\epsilon_K
\sum_{\ell\in K}\mathbf1(\mu_\ell>0)&\le100{,}000.
\end{aligned}
\tag{S83}
\]
最后一项按所有允许突变拷贝均突变的保守上界预留，而非仅按期望突变数预留。超限时中止更新，既不裁剪出生数，也不提交部分后代。

每次保存代初快照后
\[
C\leftarrow C+2L_PN_P(t)+2L_HN_H(t).
\tag{S84}
\]
个体 ID 另检查 \(\mathrm{nextID}+B_K<2^{63}\)。这些是计算边界，不属于生态密度调节。

### 14.2 随机种子与重复

\[
\mathrm{seed}_r=20261021+(r-1),\quad r=1,2,3,
\tag{S85}
\]
两处理使用相同种子列表。人口/生态和遗传随机流分别通过
\[
\begin{aligned}
\mathcal R_{\mathrm{eco}}&=
\mathrm{default\_rng}(\mathrm{SeedSequence}[
\mathrm{seed}_r,130,\mathrm{phase},K,s]),\\
\mathcal R_{\mathrm{gen}}&=
\mathrm{default\_rng}(\mathrm{SeedSequence}[
\mathrm{seed}_r,140,\mathrm{phase},K,s])
\end{aligned}
\tag{S86}
\]
构建；空间初始化与扩散使用 150 系列专用子流。这里给出算法定义，不将伪随机数序列视为生物学方程。两处理分化后的随机事件不保证逐项配对。

## 15. 当前参数与代码对应

| 参数 | 当前设置 |
|---|---|
| 世代更新次数、重复 | T=200；每处理3次 |
| 初始植物与昆虫 | P1=200，P2=200，H1=40 |
| 空间 | 10×10周期栅格；搜寻半径2，竞争半径1 |
| 配偶半径、扩散半径 | 植物配偶1，昆虫配偶2；两类群扩散均1 |
| 候选数量 | M=8 |
| 气味与偏好 | B=I₂，ρ=0.5，b=1，σ=0.2，τ=1，u₀=0，w=(0.5,0.5) |
| 资源 | R=1，dmax=1 |
| 存活与成本 | sP0=0.5，cD=cI=0.125，cV=0.05 |
| 人口学 | λP=2，λH=4；α同种=0.02，α异种=0.005 |
| 性别 | 创始和子代雌性概率均0.5 |
| QTL | 每类群8位点，每性状2位点；各位点独立连锁群 |
| 突变 | 开启处理 μ=0.01、目标效应SD=0.02；关闭处理 εP=εH=0 |

以下路径相对于 evoSim/evo_sim/models/plant_insect/。

| 内容 | 对应文件/函数 |
|---|---|
| S3–S4，S11，S32，S51，S61 | helpers/spatial.py；analysis/plot_spatial.py |
| S5，S7–S8，S43–S49，S66–S70，S72 | helpers/genetics.py：founders、phenotype、form_gametes、fertilize、summary、haplotype_summary |
| S9–S10，S13–S16，S28–S29 | helpers/voc.py；helpers/traits.py；substeps/ecology.py |
| S11–S12，S17–S27，S56 | helpers/ecology.py；helpers/traits.py；substeps/ecology.py |
| S30–S42，S52–S55 | helpers/demography.py；substeps/demography.py；substeps/reproduction.py |
| S57–S59，S64–S65，S74–S76 | analysis/analyze_macro.py |
| S60–S61，S77 | analysis/plot_spatial.py |
| S62–S63 | helpers/observations.py |
| S69–S73 | helpers/genetics.py；analysis/plot_genetics.py |
| S74及运行汇总 | analysis/summarize_qtl_run.py |
| S78–S79 | analysis/diagnose_initial_conditions.py |
| S80–S82 | helpers/demography.py；helpers/ecology.py；helpers/runtime.py |
| S83–S86 | helpers/genetics.py；helpers/runtime.py；substeps/reproduction.py；batch.py |

## 16. 使用边界

这份公式稿包含 86 组编号表达式，覆盖当前生成过程和上述科学分析代码中的主要定义，同时明确列出了非当前分支。S8、S39、S45、S49 等是代码规则的数学展开或直接推论，不是新增模拟机制。

本稿没有添加代码未实现的 GWAS、QTL 效应估计、选择梯度、连锁不平衡系数、处理效应检验或置信区间。母方后代矩、群体平均偏好和种群内标准差均按实际代码口径定义。论文定稿时可将 S1–S55 的核心机制压缩为主文 Methods，其余作为补充方法，但此处先完整保留。

