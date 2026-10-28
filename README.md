# evoSim 0.1.6

植物—植食性昆虫共同演化个体模型。以雌成虫择主、幼虫植食的生活史为抽象原型，将双方 QTL 遗传、植物防御/VOC、母体产卵、子代取食、资源竞争和双亲繁殖连成代际过程。

本版保留 SLiM 5.2、`slim-qtl-2` 和植物三模块化学映射，改变昆虫生活史：母亲的 A 决定逐卵寄主选择，SLiM 生成的子代以自己的 I 取食，幼虫存活后才进入下一代成虫。参数是未校准的理论示例；OR/OBP/P450 等真实分子功能和自动成种尚未实现。

## 运行与核查

要求 Python ≥ 3.10 和 SLiM 5.2。可通过 PATH、`genetics.slim.executable` 或 `EVOSIM_SLIM` 指定引擎；引擎缺失或失败会明确停止，不回退到替代遗传实现。

```powershell
python -m pip install -e ".[dev]"
python example.py
python -m evo_sim --config evo_sim/models/plant_insect/configs/v012_base.yaml --out results/v016_example
python -m evo_sim.models.plant_insect.audit --input results/v016_example
python -m evo_sim.models.plant_insect.analysis.summarize_lifecycle --input results/v016_example
```

省略 `--out` 时创建 `results/results_v0.1.6/run<N>`，已有目录不覆盖。`v012_` 配置文件名沿用历史名称，当前继承本版配置，用于少量机制检查，不代表执行旧生活史。

新汇总只整理出生、成年补充、寄主来源和事件网络账本，不新增宏观指标。旧 `summarize_qtl_run`、`summarize_coevolution`、`plot_diagnostics` 明确拒绝 v0.1.6，避免把旧阶段含义套到新结果；历史数据仍按原协议读取、审计。当前版本只开展小规模实现核查，不按专化、物种数或网络形态调参。

```powershell
python -m tests
python -m ruff check evo_sim tests example.py
python -m ruff format --check evo_sim tests example.py
```

实际检查记录见 [v0.1.6 实施记录](docs/plan/v0.1.6_implementation.md)。

## 一轮发生什么

1. 已完成前轮幼虫存活的成虫，在当前坐标寻找同种异性配偶。母源潜在卵数服从 Poisson，母亲自身抗性成本只折减产卵，不限制雄虫父源资格。
2. 每枚潜在卵依据母亲的气味偏好，在可达植物中选择寄主或拒绝。成功产卵后，从合格父本池逐子代均匀抽父本，批量交给 SLiM 生成基因型。
3. 同批幼虫以自己的利用性状与出生寄主防御匹配，共享有限植物资源；真实消费同时决定植物损伤与幼虫一次存活概率。
4. 幸存幼虫从出生寄主格扩散，成为下一代成虫。植物受食后经历存活、密度制约、严格异交和后代扩散；双方后代整体替换亲代。

终点仅保存成虫/植物快照及明确标注的参考诊断，不生成额外幼虫事件。`lambda` 对昆虫现表示潜在卵基准，不能与旧版有效后代数直接比较。

## 输出与科学边界

- `offspring.csv`、`gamete_transmissions.csv` 保留全部出生者，包括死亡幼虫；成年快照只保留成功补充者。
- `parents.csv` 分别记录出生和成年补充的双亲贡献；两套总和分别等于出生数和下一代成虫数。
- `oviposition.csv` 通过事件 ID 连接真实子代；拒绝没有子代或植物损伤。
- `host_flows.csv` 记录亲本出生背景到子代出生背景的定向来源贡献，包括双亲均来自 P1、子代在 P2 成活的路径；创始背景为 unknown，不推测补填。
- `networks.csv` 记录实际消费和成功补充；`lagged_networks.csv` 按亲本自己的出生背景记录真实繁殖贡献，两者不混用。
- 固定物种、可消失重现的寄主背景群、双亲系谱图分别保存。背景变化不制造生殖隔离，也不生成新物种或候选独立谱系。

审计从已保存数据独立重建方程、身份和遗传账本，只证明实现一致性，不等于生物学验证。

## 文件与修改入口

| 入口 | 内容 |
| --- | --- |
| `evo_sim/models/plant_insect/substeps/lifecycle.py` | 母体产卵—幼虫取食—成年补充及植物更新 |
| `evo_sim/models/plant_insect/helpers/` | 资源、化学、遗传、空间和配置检查 |
| `evo_sim/models/plant_insect/configs/`、`inputs/` | YAML 参数及物种/QTL/等位效应 CSV |
| `evo_sim/models/plant_insect/audit_v016.py` | 新生活史独立审计；旧审计路径继续保留 |
| [正式模型](docs/model.md)、[假设](docs/assumptions.md)、[方法](docs/methods.md) | 方程、生命阶段和证据边界 |
| [参数](docs/parameters.md)、[输出](docs/results.md)、[架构](docs/architecture.md) | 数据定义、配置和执行流程 |
| [修订后计划](docs/plan/v0.1.6_plan.md)、[Astra 审核](docs/plan/v0.1.6_astra_review_2026-10-28.md) | 理论实现合同和审查依据 |
| [Sol 功能文献核实](docs/plan/v0.1.6_literature_functions_sol_2026-10-28.md)、[Sol 谱系文献核实](docs/plan/v0.1.6_literature_lineages_sol_2026-10-28.md) | 各引用支持范围与核实深度 |

历史说明保留于 [v0.1.4](README_v0.1.4.md)、[v0.1.3](README_v0.1.3.md) 和 [v0.1.5 化学迁移](docs/plan/v0.1.5_chemical_shared_function.md)。开发约定见 [CONTRIBUTING.md](CONTRIBUTING.md)。
