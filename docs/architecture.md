# evoSim v0.1.5 — 当前实现架构

更新：2026-10-28。化学协议shared-chemical-1，遗传协议slim-qtl-2，SLiM仍为5.2。当前因果链为“基因型→原发化学模块/昆虫性状→个体互作→存活繁殖→双亲遗传→网络观察”。

## 1. 数据与模块职责

| 位置（相对evo_sim/models/plant_insect） | 当前职责 |
| --- | --- |
| configs/*.yaml、inputs/*.csv | 参数、物种、遗传结构、创始频率和参考面板；无隐藏科学默认值 |
| helpers/config.py、inputs.py | 新旧配置读取、字段/域/协议检查；新输入第四植物坐标必须空 |
| helpers/genetics.py | 类群分别3/4维加性QTL表型、原生遗传桥接、来源和统计 |
| helpers/chemistry.py | 从同一三模块x确定D/V，验证域、维度与一致性；历史映射仅显式读取旧协议 |
| helpers/voc.py、traits.py、ecology.py | 原气味选择、成本与瓶颈、同时资源分配；植物成本读x而非V均值 |
| helpers/agents.py | PlantPopulation保存chemical_traits及派生defense；旧voc_trait仅历史读取，不在新状态生成Z |
| substeps/ecology.py | 每代先保存状态、计算D/V、选择、资源、存活概率和参考观察 |
| substeps/demography.py | 生存抽样、不含母源Fm门槛的父本池、原母源Poisson产量 |
| substeps/reproduction.py | 按母源计数生成后代并均匀选父本，交给SLiM，验证后原子提交 |
| helpers/slim_backend.py、slim_reproduction.slim | 保留SLiM5.2、原生addCrossed、Haldane标记图谱及既有Eidos突变核 |
| helpers/contributions.py | 从真实双亲计数提交贡献；不增加SFm加权 |
| audit.py、audit_genetics.py | 独立重建化学、遗传与资源/生殖账本；按原协议审计历史记录 |
| analysis/ | 适配x的只读汇总和原有观察；plot_diagnostics为唯一绘图入口 |

## 2. 状态和维度

植物原发性状为N_P×3，D/V为N_P×2；植物基因型默认N_P×2×6。昆虫I/A各N_H×2，基因型N_H×2×8。数量与位点来自CSV，当前化学构型固定三模块/两条防御和气味轴。快照保存plant_chemical_traits及D、V_emit/V_cue，不生成plant_voc_independent_trait。

SLiM效应矩阵按调用类群建立，不把全局异宽等位目录强行转成矩形。另一类群的未使用等位ID行置零，仅为保持全局ID索引；被传递基因型的位点/类群归属必须验证，不能利用占位行制造性状。

## 3. 一代与衔接

初始化按CSV频率产生二倍体，由加性效应计算x和I/A。每代按x生成D/V，再按原空间候选和效用选择；资源同时结算，形成S。植物母源成本仅按x_N+x_S及x_V计算，共享模块不重复收费。母本条件产量为λFmφ，存活且有父本才抽样；父本资格不含Fm，植物同一株可承担两种角色，昆虫母雌父雄。

出生数和双亲索引由Python确定，SLiM不再施加成本选择。遗传、性别与母源扩散完成并验证后，整体替换亲代并提交Bm/Bf及(Bm+Bf)/2。失败或终代没有提交贡献，不能记零冒充实际繁殖。

## 4. 版本与观察隔离

当前新运行要求软件/模型0.1.5及slim-qtl-2。旧配置保留原含义，历史协议只读分析与审计仍可使用，重建旧模拟使用匹配代码/配置。不能只凭相同0.1.5标识认定旧四轴草稿与新三轴模型兼容。

元数据保存chemical_schema、坐标顺序、观察窗口、完整配置/架构哈希和SLiM身份。旧B、rho_voc、voc_bounds在新配置明确为null；eta_shared为独立参数。历史结果及源码哈希不改写。

观察计算只读取结果，不消耗生态/遗传随机数或反向控制网络。结构性迁移使祖先气味关系改变，旧场景名只表示输入来源，不构成新旧生态处理等价保证。完整公式、文献及迁移清单见 [model.md](model.md) 和 [实施记录](plan/v0.1.5_chemical_shared_function.md)。
