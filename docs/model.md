> 本文件说明沿用的科学方程；空间与非空间模式由所选 YAML 配置决定，默认示例为非空间 QTL 模型。初始化空间网络规则见 [host_network.md](host_network.md)。

# 模型公式、来源与 ODD 说明

这里保留当前沿用的生态公式与出处对应；显式 QTL 遗传规则见 assumptions_genetics_v14.md。文献背景不能替代参数校准。

## F01–F16 对应

| 公式 | 实施位置与内容 | 来源与证据级别 | 检查 |
| --- | --- | --- | --- |
| F01 | ecology.settle_feeding→旧consumption，防御正维度min(I/D)，零防御q=1 | [R1]继承 | T01/T02/T15 |
| F02a/b | 旧map_voc、choose_hosts及概率函数；V_emit与V_cue区分、outside保持 | [R2]仅概率选择先例；模板/VOC混合/outside为项目继承假设 | T01/T13/T15 |
| F03 | 旧allocate_resources，d=dmax*q，x=gamma*d | [R6]为竞争背景；比例缩放为继承的项目假设 | T08/T13/T15/T16 |
| F04 | g=x/dmax，h=X/R，拒绝q=NaN/g=0 | [R1/R6]背景，项目账本定义/推导 | T01/T02/T16 |
| F05 | ecology中S_I=g，F_I=max(0,1−c_IΣI)，W=S*F | [R1]β=0函数族；q到实际g再到概率为扩展/未校准近似；[R9]定性支持 | T01–T05 |
| F06 | ecology中S_P=sP0+(1−sP0)(1−h)，F_P沿用真实VOC成本 | [R1]继承，h与VOC项为扩展；[R8]耐受限制 | T01/T15/T16 |
| F07 | density_phi，φ_s=1/(1+Σ_r α_sr N_r)，N为取食前人数 | [R3 Eq1/R4 Eq5]竞争分母形状；应用于个体条件后代为改造 | T06/T07/T16 |
| F08 | alpha_from_K=(λ−1)/K，仅λ>1和K>0 | [R3/R4]理想单物种代数推导，不是硬上限 | T06 |
| F09 | draw_offspring：uniform<S实现Bernoulli | [R1]概率分解/[R4]人口学背景；个体两阶段为项目假设 | T02/T03/T05 |
| F10a/b | μ=λFφ，死者B=0，存活者Poisson(μ)；昆虫φ=1 | [R4 Eq6]种群Poisson背景+[R1]S/F；本个体条件实现为改造/假设 | T03–T05/T08 |
| F11 | expected_offspring=Sμ=λWφ | 全期望定律严格推导，[R1/R4]背景 | T05/T06/T16 |
| F12 | offspring_variance=Sμ+S(1−S)μ² | 全方差定律严格推导；不同于直接Poisson(Sμ) | T05/T16 |
| F13 | inherit按每亲代B重复，N_next=ΣB，同种同亲代遗传/旧突变 | [R1/R2/R4]扩展；无固定规模重采样 | T09/T12/T16 |
| F14 | 没有亲代就没有后代，零为吸收态，保留零网络行 | 无迁入/无成种的项目边界与推导，[R4]人口学背景 | T09–T11 |
| F15 | resource_vector每代按新植物物种赋R，虫赋dmax | [R5]生命周期边界背景+上游项目假设；非同株再生 | T13/T16 |
| F16 | 理想单种S=F=1时E[N_next]=λN/(1+αN) | [R3/R4]种群函数形状，由本个体算法求和推导 | T06 |

没有把两阶段个体抽样、比例共享或气味模板伪称为文献原公式。条件于本代状态，亲代存活与条件产子独立；无环境同步噪声、配偶限制或种子库。这些是建模边界，不是校准结果。

## 文献与本轮核查范围

- [R1 Speed 2015](https://doi.org/10.1111/nph.13560)：使用用户规范与上游明确的q/S/F函数族，本轮未重新取得全文。
- [R2 Ohsaki 2025](https://doi.org/10.1016/j.jtbi.2025.112094)：沿用规范中的概率寄主选择与数量更新背景；项目具体分配/VOC公式不归给原文。
- [R3 Hart 2018](https://www.iee.unibe.ch/unibe/portal/fak_naturwis/d_dbio/b_ioekev/content/e54478/e384613/e639025/e687350/files728857/Hart2_eng.pdf)：本轮可读原始PDF，第2页§2 Eq1核对了竞争分母和系数受体方向；两物种形式推广至矩阵求和已在本项目说明。
- [R4 Shoemaker 2020](https://doi.org/10.1002/ecy.2922)：本轮确认出版信息与人口学/环境随机性主题；作者站PDF初次检索可见，后续403，期刊全文不可访问，Eq5/Eq6的具体编号对应沿用用户规范，未假称本轮完成逐式全文复核。个体Bernoulli→Poisson明确是项目选择。
- [R5 Stouffer 2022](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13965)：本轮可读原文，确认其对离散植物生活史与密度依赖繁殖结构一致性的讨论。该文并非为本项目采用BH函数作普遍生物学背书；保留BH仍是明确近似。
- [R6 Richards 2000](https://doi.org/10.1086/303316)：竞争背景，非当前比例共享原式；本轮未重取全文。
- [R7 Grimm 2020](https://doi.org/10.18564/jasss.4259)：ODD为描述与复现规范，不是生物学验证；本轮使用规范提供的引用。
- [R8 Strauss & Agrawal 1999](https://doi.org/10.1016/S0169-5347(98)01576-6)：耐受/补偿限制的背景，未重取全文。
- [R9 Gols 2022](https://doi.org/10.1111/eea.13192)：数量影响发育/繁殖的定性支持，不标定S=g或λ；未重取全文。

λ与α必须逐物种/逐物种对填写kind/source，固定样地与feeding_window也必须显式提供。开发夹具为theoretical_placeholder且illustrative_only=true；继承生态数值引用上游开发示例，也非实测值。α的单位为1/株，R/dmax为同窗口干质量；不把质量单位缩放用于α。

## ODD：目的、实体、时间顺序与子模型

**目的**：实现可追溯的生态—人口学—遗传流程，检查数学与程序边界；不检验共存、灭绝风险或专食假说。

**实体与状态**：植物D、plant_voc_independent_trait、物种/唯一ID；昆虫I/A、物种/唯一ID。V_emit/V_cue、资源、供需、存活概率和后代均值为当代派生状态。物种槽位固定，人数可为0。空间配置采用周期斑块栅格（局部遭遇、局部竞争与有界扩散，见 host_network.md）；下方公式对应原始非空间样地版本。

**时间顺序**：t代取食前统计N并为新个体赋资源；候选选择；共享取食；S/F；植物φ；Bernoulli；存活者Poisson；安全检查；同亲代遗传与旧突变；提交t+1代。亲代退出，无多年生存续或种子库。最终状态保留生态与条件矩但不额外抽一代。

**随机性**：SeedSequence([seed,130,phase,taxon,species])，phase0初始化、1生态、2存活、3生育、4突变。旧选择行为计算、q和突变函数复用。只要求本版本同配置同种子位级复现；不要求旧人口学版本整体轨迹相同。

**观察与初始化**：初始物种人数可为非负整数；空物种保留身份。初始性状用旧initialize，支持测试传入D/I/A/植物气味性状。每代变长NPZ、个体亲代表、子代追踪、生态账本与所有物种人口行均保存。未实现适应性行为调整、营养同化或碳质量到后代的守恒；空间落位已有有界扩散（周期栅格，见 host_network.md），不是连续移动轨迹。

**计算限制**：max_individuals_safety作用于下一代植物+昆虫合计人数，在抽样后、分配后代前检查；超过停止并保存全体抽样计数，N_next空且update_applied=false，不把它伪装为生态承载量。异常前已完成的原始状态和账本保留，失败阶段记录在status。为完整审计，所有世代保留内存，长运行需未来流式写出。
