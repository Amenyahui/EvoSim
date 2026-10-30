# evoSim 0.1.8（开发实现）

植物—植食性昆虫共同演化个体模型。以雌成虫择主、幼虫植食的生活史为抽象原型，将双方 QTL 遗传、植物防御/VOC、母体产卵、子代取食、资源竞争和双亲繁殖连成代际过程。

软件 v0.1.8 保留 v0.1.7 科学模型协议，当前使用 `slim-qtl-4` / `expression-tuning-2`，保留 SLiM 5.2 和原成年空间交配：植物三个化学投入轴、昆虫感觉/处理模块的八个可遗传表达与调谐坐标，以及只读候选谱系和动态事件网络。昆虫从出生寄主扩散决定成年位置，在局部兼容异性池中交配，再由母体逐卵择主；未来产卵寄主不能决定此前交配。植物沿用存活、局部同祖先兼容池异交、均匀花粉供源与母株后代扩散，不增加花期。植物化学性状仍参与双向选择和遗传演化；本版未提供额外的植物适应性成种机制。参数是理论示例，抽象模块不是真实OR/OBP/P450家族，候选谱系不等于确证物种。

此前含植物花期的 `slim-qtl-3` 草稿、含交配前定居的 `expression-tuning-1` 仅保留历史回读；不允许作为当前运行协议。交配前定居及其开关、专属输出已删除。已有草稿CSV、原始运行数据及source/metadata哈希保持原身份，不把移除机制后的文档套回旧结果。

新版本需显式选择 `v017_modules.yaml` 或 `v017_single_ancestor.yaml`。默认旧YAML及 `example.py` 保留v0.1.6模型，软件版本与模型版本分别记录，避免把旧输入静默换成新语义。完整路线见 [v0.1.7计划](docs/plan/v0.1.7_plan.md)，当前交付与验收见 [基线修订记录](docs/plan/v0.1.7_baseline_revision_2026-10-29.md)；[旧实施记录](docs/plan/v0.1.7_implementation.md)保留历史事实。

## 运行与核查

要求 Python ≥ 3.10 和 SLiM 5.2。可通过 PATH、`genetics.slim.executable` 或 `EVOSIM_SLIM` 指定引擎；引擎缺失或失败会明确停止，不回退到替代遗传实现。

```powershell
python -m pip install -e ".[dev]"
python -m evo_sim --config evo_sim/models/plant_insect/configs/v017_single_ancestor.yaml --out codex/validation/v017_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/v017_example
python -m evo_sim.models.plant_insect.analysis.summarize_lifecycle --input codex/validation/v017_example
python -m evo_sim.models.plant_insect.analysis.summarize_lineages --input codex/validation/v017_example --design evo_sim/models/plant_insect/analysis/v017_lineages.yaml
# A–E 图比较要求两处理均完整且至少运行4代；正式200代批次例如：
# python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input results/results_v0.1.7/200gen_H400_corrected_20261029
```

省略 `--out` 时按软件版本创建 `codex/validation/v0.1.8/run<N>`，已有目录不覆盖。开发测试、短程验收及其汇总统一放在 `codex/validation/`；`results/` 仅用于显式指定的正式科研实验。`v012_` 配置沿用历史名称且当前继承v0.1.6配置，不代表执行v0.1.2生活史。

生命周期汇总整理出生、成年补充、寄主来源和事件网络账本；新谱系汇总另输出观察分群、分裂/合并/终止、交配机会和两种时标明确的网络。旧 `summarize_qtl_run`、`summarize_coevolution` 命令继续拒绝新生活史。唯一绘图入口 `plot_diagnostics` 支持v0.1.6和当前v0.1.7：A含终点成年数，B直接展示植物三化学轴及昆虫八模块轴，C展示幼虫寄主资源份额和灭绝前或末轮繁殖阶段，D为空间快照，E为QTL遗传多样性；不将派生I/A误画成旧偏好模板。按用户安排完成的400只昆虫、200代正式批次，完整性和结果见[基线修订记录](docs/plan/v0.1.7_baseline_revision_2026-10-29.md)；不按专化、物种数或网络形态调参。

```powershell
python -m tests
python -m ruff check evo_sim tests example.py
python -m ruff format --check evo_sim tests example.py
```

当前检查见 [v0.1.7基线修订记录](docs/plan/v0.1.7_baseline_revision_2026-10-29.md)；原版本检查保留于 [v0.1.6实施记录](docs/plan/v0.1.6_implementation.md)。

生活史与人口学基线的70次小规模机制对照见 [基线检查报告](docs/plan/v0.1.6_demographic_baseline_2026-10-28.md)。可用 `configs/v016_demographic_controls.yaml` 复用“初始识别匹配×配偶范围”的20代检查，输出显式放在 `codex/validation/`；此配置保留原默认参数作为对照，不替换默认初始条件。

研究已有寄主关系的变化时，显式选择 `configs/v016_ancestral_matched.yaml`：初始识别和利用均偏向P1，匹配写在可遗传的A等位效应中。研究识别不匹配下的建立/失败时，选择 `configs/v016_recognition_mismatch.yaml`，保留原输入。两者均是理论情景，不能把前者称为稳定共存或把后者称为完全未适应祖先。新增50代核定中，适宜寄主单独存在时4/4重复存续，双寄主时3/4存续；完整结果和存活函数证据见 [核定结论](docs/plan/v0.1.6_reference_qualification_2026-10-28.md)。

```powershell
# 50代、关闭新突变、两种寄主配置各4重复的开发核定；目录须不存在
python -m evo_sim --config evo_sim/models/plant_insect/configs/v016_reference_validation.yaml --out codex/validation/v016_reference_new
python -m evo_sim.models.plant_insect.audit --input codex/validation/v016_reference_new
```

## v0.1.6保留生活史与v0.1.7扩展

下列四步保留母体择主—子代利用基础。成虫位置由此前从出生寄主的扩散决定，直接用于空间交配；本轮逐卵选择不移动成虫或重建配偶池。母体和幼虫使用各自基因型的功能模块；新协议I/A是派生诊断，不能当作旧直接遗传坐标。`shuffle_voc`只打乱逐卵产卵线索。

1. 已完成前轮幼虫存活的成虫，在当前坐标寻找同祖先兼容池异性配偶。母源潜在卵数服从 Poisson；v0.1.7母亲自身处理/感觉表达成本只折减产卵，不限制雄虫父源资格（v0.1.6仍按母亲I收费）。
2. 每枚潜在卵依据母亲的气味偏好，在可达植物中选择寄主或拒绝。成功产卵后，从合格父本池逐子代均匀抽父本，批量交给 SLiM 生成基因型。
3. 同批幼虫以自己的利用性状与出生寄主防御匹配，共享有限植物资源；真实消费同时决定植物损伤与幼虫一次存活概率。
4. 幸存幼虫从出生寄主格扩散，成为下一代成虫。植物受食后抽存活，在局部同池非自身的存活个体中逐后代均匀抽父本；合格母株按Poisson(lambda_P Fm phi)产出后代，经SLiM遗传并从母株格扩散。双方后代整体替换亲代。

终点仅保存成虫/植物快照及明确标注的参考诊断，不生成额外幼虫事件。`lambda` 对昆虫现表示潜在卵基准，不能与旧版有效后代数直接比较。

## 输出与科学边界

- `offspring.csv`、`gamete_transmissions.csv` 保留全部出生者，包括死亡幼虫；成年快照只保留成功补充者。
- `parents.csv` 分别记录出生和成年补充的双亲贡献；两套总和分别等于出生数和下一代成虫数。
- `oviposition.csv` 通过事件 ID 连接真实子代；拒绝没有子代或植物损伤。
- `host_flows.csv` 记录亲本出生背景到子代出生背景的定向来源贡献，包括双亲均来自 P1、子代在 P2 成活的路径；创始背景为 unknown，不推测补填。
- `networks.csv` 记录实际消费和成功补充；`lagged_networks.csv` 按亲本自己的出生背景记录真实繁殖贡献，两者不混用。
- 原始账本分别保存祖先兼容池、可消失重现的寄主背景群和双亲系谱图。v0.1.7后处理另追踪操作谱系；背景或谱系标签不制造生殖隔离，不自动宣称新生物物种。

审计从已保存数据独立重建方程、身份和遗传账本，只证明实现一致性，不等于生物学验证。

## 文件与修改入口

| 入口 | 内容 |
| --- | --- |
| `evo_sim/models/plant_insect/substeps/lifecycle.py` | 母体产卵—幼虫取食—成年补充及植物更新 |
| `evo_sim/models/plant_insect/helpers/` | 资源、化学、遗传、空间和配置检查 |
| `evo_sim/models/plant_insect/configs/`、`inputs/` | YAML 参数及物种/QTL/等位效应 CSV |
| `evo_sim/models/plant_insect/audit_v016.py` | 新生活史独立审计；旧审计路径继续保留 |
| [现行模型](docs/model.md)、[假设](docs/assumptions.md) | v0.1.7 方程、生命阶段和证据边界 |
| [历史方法](docs/methods.md)、[历史模型](docs/history/model_v0.1.6.md) | v0.1.6 方法与原始理论合同 |
| [参数](docs/parameters.md)、[输出](docs/results.md)、[架构](docs/architecture.md) | 数据定义、配置和执行流程 |
| [修订后计划](docs/plan/v0.1.6_plan.md)、[Astra 审核](docs/plan/v0.1.6_astra_review_2026-10-28.md) | 理论实现合同和审查依据 |
| [Sol 功能文献核实](docs/plan/v0.1.6_literature_functions_sol_2026-10-28.md)、[Sol 谱系文献核实](docs/plan/v0.1.6_literature_lineages_sol_2026-10-28.md) | 各引用支持范围与核实深度 |

历史说明保留于 [v0.1.4](README_v0.1.4.md)、[v0.1.3](README_v0.1.3.md) 和 [v0.1.5 化学迁移](docs/plan/v0.1.5_chemical_shared_function.md)。开发约定见 [CONTRIBUTING.md](CONTRIBUTING.md)。
