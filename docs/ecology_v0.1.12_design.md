# v0.1.12 生态机制修改方案

本次为用户明确授权的 Scientific Modification，以提交 fd9ba73 的 v0.1.11 为历史基线。
完整源码、输入、依赖与执行说明及逐文件 SHA256 保存在
`codex/history/before_v0.1.12_fd9ba73/`；Git 基线提交不改写。

## 当前问题与替代机制

v0.1.11 将三维植物投入固定读为 D=(xN,xS)、V=(eta*xS,xV)，昆虫
min(I/D) 同时限制需求与食物相关存活。低资源新株被计为受害，繁殖使用恢复前资源。
新版本取消这些关系，保留空间、配偶、SLiM 原生减数分裂、突变及原子提交。

植物原发坐标依次为 (PSM1,PSM2,Z,VOC1,VOC2)，均在 [0,1]。
D 为前两轴，V 为后两轴，Z 为第三轴。默认每轴两个独立组内连锁的 QTL；
已有 targets/effect 向量允许多效性，默认不引入相关方向。
昆虫前八轴不变，后两轴为接触坐标；响应权重 G=2*coordinate-1。
这些是功能表型，不是具体 GR/P450/GST/UGT 基因或表达量。

## 明确的理论方程

每卵在空间候选内按现有 VOC 嗅觉效用抽取未检查植物，接触接受概率为
sigmoid(contact_intercept + contact_slope * sum(G_k*PSM_k))。
接触拒绝移除当前植物并继续；最多检查候选数次。outside 模式仍允许嗅觉外部选项
终止搜索，forced 只移除该选项，并不跳过接触检查。每次检查记录植物、概率与结果。
候选与逐卵随机流的种子派生规则不变；新行为消耗额外抽样，新旧版本轨迹不要求相同。

幼虫实际需求 d=dmax，与 I/PSM 无关。沿用同时按需求比例分享资源，
x=gamma*d，g=x/dmax。每轴摄入剂量 dose=x*PSM*dose_scale；
有效处理容量 capacity=I*processing_scale，processed=min(dose,capacity)，
residual=dose-processed。抗营养负荷 antinutrient=dose，作为独立的食物质量作用，
本版不假定处理模块消除抗营养作用。营养、毒性与抗营养存活分量分别为
g^nutrition_exponent、exp(-sum(toxicity_weights*residual))、
exp(-sum(antinutrition_weights*antinutrient))，总存活为三者乘积。
零摄入的营养存活定义为0，即使 nutrition_exponent=0；正摄入且 exponent=0
是移除摄入不足效应的理论对照。没有额外重复的解毒成本。

植物实际取食损伤 h=X/B_before（B_before=0 时 h=0），
S=exp(-mu0-muD*h)。恢复率 r=resource_recovery_rate+tolerance_recovery_rate*Z，
B_recovered=K-(K-B_post)*exp(-r)。母源乘子为
max(0,1-cPSM*sum(PSM)-cVOC*sum(VOC))*B_recovered/K。
恢复只对幸存株实现；潜在恢复量用于幸存后的条件繁殖。不增加耐受固定成本。
新株仍以配置比例的 K 建立，背景死亡、空位抽签与双亲遗传不变。

## 参数、证据与预期后果

新增剂量缩放、处理缩放、毒性权重默认各轴1，抗营养权重各轴0.1，
营养指数1，接触截距0、斜率1，耐受恢复率系数 ln(2)。均为显式的
theoretical_placeholder，未用模拟结果拟合。剂量缩放单位为抽象负荷/资源量，
处理缩放为抽象负荷/处理表型/窗口，毒性及抗营养权重为负荷的倒数，
恢复率为每昆虫代。浓度与功能坐标为归一化量，不宣称对应实测生化单位。
新 Z、VOC1、G 原始中心均为0.5；保留旧 PSM、VOC2 及八个昆虫坐标中心，
每轴位点效应偏差沿用 ±0.025，属于新版本理论初始条件。

[Kaplan 等 2014](https://doi.org/10.1371/journal.pone.0093714) 的摄入、增长与转化实验
支持区分这些终点，但不确定本版乘积存活方程。
[Rayapuram & Baldwin 2006](https://doi.org/10.1111/j.1365-3040.2006.01534.x)
提供营养指标与植物抗性案例；[Strauss & Agrawal 1999](https://doi.org/10.1016/S0169-5347(98)01576-6)
为耐受概念综述。接触 logistic、线性容量、独立指数危害及 Z 线性恢复率是
本项目理论闭合，不是上述文献的实测定律。

预期：处理能力改变残余负荷但不改变当轮消费；PSM/VOC 可独立进化；
接触可接受或拒绝与嗅觉不同的植物；高 Z 更快恢复，同损伤死亡率不受 Z 直接改变。
无虫新株的低资源不再被计为食害。恢复后繁殖会改变人口反馈，但不保证共存或专化。

验收分别覆盖八项用户要求、参数边界、固定种子重复、独立重建全部新记录、
重新散列后的篡改拒绝以及失败时不提交半轮数据。旧读者拒绝新协议，旧输出按旧基线读取。

## 实施文件与原因

下表路径均相对 `evo_sim/models/plant_insect/`，测试与文档除外。

| 文件组 | 修改及原因 |
| --- | --- |
| `biology/plants/chemistry.py`、`biology/genetics/qtl.py`、`biology/agents.py`、`config/input_loader.py`、`inputs/current/*.csv` | 植物5维、昆虫10维；独立PSM/VOC与Z/G位点、效应、频率和参考表型保持一致，保留可配置多效性。 |
| `biology/insects/host_selection.py`、`simulation/stages/oviposition.py` | 在嗅觉接近后执行有符号接触响应；拒绝后搜索剩余候选，避免把嗅觉选择当作最终接受。 |
| `biology/insects/feeding.py`、`biology/traits.py` | 食物需求与处理能力分开；摄入后逐轴记账并组合三种存活分量，保留原功能模块表达成本。 |
| `biology/plants/demography.py`、`biology/plants/defense.py`、`simulation/stages/plant_dynamics.py` | 实际损伤与资源亏损分开；Z控制恢复，恢复量参与繁殖；两维VOC成本显式相加，不新增Z固定成本。 |
| `config/schema.py`、`configs/current_base.yaml` | 所有新系数显式配置、校验并标注理论假设；升级协议以拒绝混用历史含义。 |
| `records/events.py`、`records/observations.py`、`simulation/runtime.py`、`batch.py` | 保存接触轨迹、化学负荷与恢复阶段；记录接收计算结果，元数据明确旧列的新含义。 |
| `audit/lifecycle.py`、`audit/lifecycle_checks/` | 独立重放接触抽样、化学剂量和资源恢复，核对新记录及原有遗传/事务约束。 |
| `analysis/summaries/lifecycle.py`、`analysis/plots/diagnostics.py`及分析协议入口 | 用实际总存活计算预期招募，区分损伤和资源亏损，更新5/10维图轴；生殖隔离计算规则不变。 |
| `tests/scientific/test_ecology_v012.py`及相关单元/集成/验收测试 | 验证八项要求和边界；替换与旧机制绑定的断言，保留历史测试副本。 |
| `README.md`、`docs/model.md`、`docs/parameters.md`、`docs/results.md`及版本声明 | 同步方程、假设、单位、协议和实际验证记录。 |

最终116组测试、8次各3轮的小规模运行与隔离安装包校验通过；
执行证据与尚未完成的长期/经验验证见 [results.md](results.md#本版验证)。
