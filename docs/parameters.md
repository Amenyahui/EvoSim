# v0.1.12 参数、输入与证据

本版参数为理论假设，未按自然系统标定，也未按模拟共存、专化或成种结果拟合。
方程、适用范围见[模型](model.md)，科学修改理由见[设计](ecology_v0.1.12_design.md)。
原0.1.11逐字段表及文档完整保存于Git提交fd9ba73和冻结源码包；当前[参数清单](parameter_inventory.csv)
从当前默认配置展开，包含新参数、单位、来源、证据地位和完整输入CSV的哈希；
其他配置的差异以各YAML及其继承关系为准。

## 输入与遗传

默认 `configs/current.yaml`：十植物类型各20株、200只H1、五轮、有虫与无虫各4重复，
双方突变关闭；`current_base.yaml`保留启用突变与关闭突变的对照；
`review_diagnostics.yaml`是300轮独立设计，配置存在不表示已运行。

输入CSV位于 `inputs/current/`。植物顺序PSM1/PSM2/Z/VOC1/VOC2（5维），
昆虫顺序EI1/EI2/muI1/muI2/EA1/EA2/muA1/muA2/G1raw/G2raw（10维）。
CSV effect/target/mutation_mean/mutation_sd均有10列，植物只填前5列，其余留空。
原位点ID `plant_x_N_*`、`plant_x_S_*`、`plant_x_V_*`分别作用PSM1、PSM2、VOC2，
沿用ID不代表旧三轴映射；新增Z、VOC1及G位点有明确名称。
默认每轴2QTL、同轴0/10cM、轴间独立；targets可包含多个轴，不强制遗传相关。

| 参数 | 当前值 | 单位/意义 | 证据地位 |
| --- | --- | --- | --- |
| phenotype_map.effect_scale | 0.5 | 等位效应到归一化坐标的缩放 | 理论假设 |
| phenotype_map.bounds | [0,1] | 原发坐标边界 | 理论假设 |
| background | plant五个0、insect十个0 | 加性背景 | 理论假设 |
| 各位点mutation_rate | 0.01 | 每传递QTL拷贝、每出生后代 | 理论假设，非核苷酸率 |
| mutation_mean / mutation_sd | 0 / 0.02 | 目标坐标等位效应增量 | 理论假设 |
| 新Z、VOC1、G原始中心 | 0.5 | 归一化原发坐标 | 理论假设 |
| 创始频率、效应偏差 | 各0.5，中心/2±0.025 | 每轴两QTL、独立抽拷贝 | 理论假设 |
| G功能权重 | 2*raw-1 | [-1,1]综合接触响应 | 理论假设，非GR表达 |

保留旧PSM1/PSM2、VOC2和八个昆虫坐标的输入效应中心，没有根据新模拟结果反向调整。
SLiM 5.2双亲遗传及标记重组算法不变；新增位点导致抽样和轨迹不与旧协议相同。

## 摄入与化学处理

| 配置项 | 当前值 | 单位与含义 | 来源/限制 |
| --- | --- | --- | --- |
| chemical_processing.dose_scale | [1,1] | 抽象化学负荷/可食资源量；与归一化PSM相乘 | theoretical_placeholder |
| processing_scale | [1,1] | 每处理表型、每窗口的抽象负荷容量 | theoretical_placeholder |
| toxicity_weights | [1,1] | 每轴残余负荷的逆单位危害系数 | theoretical_placeholder |
| antinutrition_weights | [0.1,0.1] | 每轴摄入负荷的逆单位危害系数 | theoretical_placeholder |
| nutrition_exponent | 1 | 正摄入的g指数；0为移除营养不足对照 | theoretical_placeholder |
| contact_response.intercept | 0 | 接触logit截距，无量纲 | theoretical_placeholder |
| contact_response.slope | 1 | 签名接触权重与归一化PSM的logit缩放 | theoretical_placeholder |

剂量为x*PSM*dose_scale，已处理负荷为min(dose,I*processing_scale)。
抗营养作用按全部摄入负荷计算，不默认被处理模块消除。不同轴可独立指定权重。
所有上述数值均要求非负有限值（接触截距允许负值），provenance必须声明理论来源。
相关生理实验只支持区分终点，不标定此处数值，参考文献及解释见模型文档。

## 植物与成本

| 配置项 | 当前值 | 单位/解释 | 证据地位 |
| --- | --- | --- | --- |
| background_mortality_rate | 0.01005033585350145 | 每昆虫代背景危害 | 理论假设 |
| defoliation_mortality_rate | 0.13862943611198905 | 每昆虫代、单位实际损伤的附加危害 | 数字保留，定义改为X/B_before，未标定 |
| resource_recovery_rate | 0.6931471805599453 | Z=0时每代恢复率 | 理论假设 |
| tolerance_recovery_rate | 0.6931471805599453 | 每代、单位Z的恢复率增量 | 新增理论假设 |
| newborn_resource_fraction | 0.25 | 新生株初始资源/K | 理论假设 |
| plant_cost | 0.125 | 每单位PSM合成投入的母源折减 | 理论假设 |
| voc_cost | 0.05 | 每单位VOC合成投入的母源折减，两轴均收费 | 理论假设 |
| insect_cost / sensory_cost | 0.125 / 0.05 | 处理表达/嗅觉表达母源成本 | 理论假设 |
| 耐受固定成本、额外处理成本 | 无 | 不自动收费，不重复收费 | 本版范围 |

取食前资源、取食后资源和恢复后资源分开记录；恢复后资源进入当轮植物繁殖。
无食害的低资源新株只受背景死亡，不被误判为食害；Z不直接改变存活或PSM。

## 保留的人口、空间与实现参数

默认十类型K=10、dmax=1，单位为归一化可食干生物量/幼虫窗口；
基准CSV的双植物小示例K=1。植物lambda=2提案/合格幸存母株/轮，
昆虫lambda=4潜在卵/有配偶雌虫/轮。10×10周期网格、每格3个建立位点、
寄主搜索/交配半径2、植物/昆虫扩散半径1、候选上限8、性别概率0.5，均为理论设置。
嗅觉宽度及处理宽度0.25、idle=.1、sensory_half=.2、valence=[1,-1]，仍为固定功能近似。
资源容量、配偶资格、建立抽签、无移民/种子库及世代重叠规则保留。
max_individuals=5000、max_alleles=100000、累计基因型值4000000是计算保护，不是生态密度调节。
历史结果不重算、不覆盖；参数改变须更新配置、本文及对应测试。
