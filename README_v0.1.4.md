# evoSim v0.1.4：双向功能与寄主网络观察

更新：2026-10-27。实现依据为 [版本计划](docs/plan/v0.1.4_plan.md)，遗传协议继续使用 `slim-qtl-1`。

## 本版完成的内容

围绕“遗传变异—功能性状—个体互作—多物种网络”，增加五类数据分析：

1. 用保存的早期/后期真实植物和昆虫个体，组成四种时间交叉组合，比较单寄主接受、利用、资源和损伤及双方 SF 生态代理。
2. 将植物 D/Z、昆虫 I/A 与实际双亲遗传贡献按 ID 连接，分解亲代贡献差异和下一代表型残差。
3. 联合报告候选暴露、实际接受/拒绝、固定面板 p/q、双方丰度、实际资源及寄主相关繁殖贡献。
4. 先逐个体计算参考反应的有效寄主宽度，再汇总均值、方差和有响应个体数。
5. 分别跟踪资源边和昆虫繁殖贡献边，记录首次正边、连续代数及预定窗口内的持续性。

所有新分析只读取已有数据，不影响模拟状态或随机数。生态及遗传规则沿用 v0.1.3，版本和默认输出目录升级为 0.1.4。两组历史 mutation-on/off 仍只比较新突变输入，不能当作演化与不演化对照。

## 使用方法

SLiM 5.2 的安装与引擎配置沿用 [v0.1.3 说明](README_v0.1.3.md)。在项目根目录运行：

```powershell
# 新批次示例；目标目录必须不存在
python -m evo_sim --config evo_sim/models/plant_insect/configs/v012_base.yaml --out codex/validation/v014_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/v014_example
python -m evo_sim.models.plant_insect.analysis.summarize_coevolution --input codex/validation/v014_example
```

五张新 CSV 及 `coevolution_metadata.json` 写入批次 `data/summary/`。常规汇总和 A–E 绘图入口继续可用；本版没有另建绘图入口。默认新批次目录为 `results/results_v0.1.4/runN/`。

默认观察参数位于 [analysis/v014_analysis.yaml](evo_sim/models/plant_insect/analysis/v014_analysis.yaml)。需要自定义时复制到批次根目录并命名 `analysis.yaml`，提供全部键。可指定两时期、统一试验资源/需求、资源和贡献边阈值及持续窗口；定义见 [参数说明](docs/parameters.md)。

分析逐项检查所读取原始文件的保存哈希，并记录当前分析代码、设计、输入和原审计的指纹。历史结果保留原版本与源码身份；允许用新分析读取，但不能把重分析写成重新完成的演化实验。

## 理论依据与解释边界

| 项目 | 依据 | 支持范围 |
| --- | --- | --- |
| 时间交叉评估 | Buckling & Rainey (2002), [10.1098/rspb.2001.1945](https://doi.org/10.1098/rspb.2001.1945) | 细菌—噬菌体原始实验支持跨时期比较策略；不验证植物—昆虫函数 |
| 双亲贡献分解 | Price (1970), [10.1038/227520a0](https://doi.org/10.1038/227520a0) | 贡献协方差与传递项的区分；本模型具体权重和表型残差另行定义 |
| 偏好与利用分别观察 | Forister (2005), [10.1111/j.0014-3820.2005.tb01051.x](https://doi.org/10.1111/j.0014-3820.2005.tb01051.x) | 两类性状可独立遗传的实验支持；生活史不同，不直接验证自主取食 |
| 个体互作产生网络 | Maliet et al. (2020), [10.1111/ele.13592](https://doi.org/10.1111/ele.13592) | IBM 网络涌现的理论先例，不预设本模型网络形状 |
| 丰度与机会控制解释 | Blüthgen (2010), [10.1016/j.baae.2010.01.001](https://doi.org/10.1016/j.baae.2010.01.001) | 限定网络指标的生态解释，不支持任意阈值为生物常数 |

单寄主试验与野外多候选选择不同；固定面板 reference_p 是条件相对偏好，不能等同绝对接受率。SF 不含配偶、性别、密度与父源成功，不能等同真实繁殖适合度。实际贡献差异含随机性，传递残差不只是突变。全部具体公式见 [model.md](docs/model.md)。

固定功能分布的四组演化控制、自动因果归类寄主转换、扰动后恢复、成种和真实分子基因家族映射尚未实现。物种形成继续是解释新增网络节点与多样性产生的后续里程碑，不用物种标签阈值强制产生隔离。

## 同步文档

`docs/` 根目录的 [architecture](docs/architecture.md)、[assumptions](docs/assumptions.md)、[methods](docs/methods.md)、[model](docs/model.md)、[parameters](docs/parameters.md)、[results](docs/results.md) 均按当前规则、公式、配置和输出更新；历史设计保存在 `docs/plan/`，论文材料在 `docs/paper/`。

## 验证记录（2026-10-27）

- `python -m ruff check evo_sim tests example.py` 和 `python -m ruff format --check evo_sim tests example.py` 通过。
- `python -m tests`：40 组全部通过，包含真实 SLiM 遗传、资源账本、终点/失败原子性，以及新增的解析、双向功能、贡献恒等式、持续窗口和保存结果被修改后的拒绝检查。
- 原 `results/results_v0.1.3/slim_200gen_20261026` 六次 200 代运行已经用新入口完成重分析，未重跑演化过程。新表分别为 cross_assay 48 行、contribution_differentials 12,960 行、host_relations 2,412 行、host_edges 4,824 行、host_breadth 2,412 行。1,242 份实际输入的保存哈希核对通过。
- 其中 12 行交叉组合因昆虫灭绝而不可用；60 行终代性状贡献记录未提交，贡献差异为 NA。它们保留在输出中，没有删去不利重复或补成零。

这些检查证明实现按定义执行，不构成适应、共同演化、稳定性或成种的科学结论。
