# evoSim

植物—植食性昆虫共同演化个体模型：母体择主、幼虫资源竞争、成年空间交配与双亲 QTL 遗传共同决定代际变化。软件0.1.8，当前工作模型协议0.1.7，遗传/功能协议 `slim-qtl-4` / `expression-tuning-2`；版本号分别记录，不强行统一。

## 运行与阅读

要求 Python ≥3.10、SLiM 5.2。SLiM通过 PATH、`EVOSIM_SLIM` 或配置指定，缺失或失败时停止。默认Runner、example.py与命令行示例共用当前多虫理论配置。

SLiM 5.2可按[官方安装说明](https://github.com/MesserLab/SLiM/blob/v5.2/README.md)从conda-forge安装。项目根目录运行：

```powershell
python -m pip install -e ".[dev]"
conda create --prefix ./codex/slim_runtime --override-channels -c conda-forge slim=5.2 -y
& ./codex/slim_runtime/Library/bin/slim.exe -v
$env:EVOSIM_SLIM = (Resolve-Path ./codex/slim_runtime/Library/bin/slim.exe).Path
python -m tests
```

运行时不随Python包发布；测试会先核实SLiM存在。`audit`默认检查账本完整性，科研运行验收须加示例中的成功状态和最小转移数要求。

```powershell
python -m evo_sim --config evo_sim/models/plant_insect/configs/current.yaml --out codex/validation/multilarva_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/multilarva_example --require-success --min-transitions 1
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/multilarva_example --sensitivity
```

该入口比较R=10/30/100，各4个种子、5代；200代延长配置为 `configs/current_200gen.yaml`。指定目录须不存在；省略 `--out` 默认放入 `codex/validation/v<version>/run<N>`。正式科研实验须由用户明确安排并显式指定 `--out results/...`。

| 文档 | 内容 |
| --- | --- |
| [模型与理论](docs/model.md) | 方程、生活史、假设边界、寄主隔离与谱系网络 |
| [参数与文献](docs/parameters.md) | 输入、取值、单位、分系统证据及未标定项；附[逐字段清单](docs/parameter_inventory.csv) |
| [输出与绘图](docs/results.md) | 文件字段、审计、汇总及A–H/敏感性图 |
| [开发规则](AGENTS.md) | 代码入口、检查和结果保护 |

<a id="current-work"></a>

## 当前状态与下一步

当前P40/H40，单株归一化资源R=10、单幼虫需求dmax=1；关闭新突变，保留站立变异与重组。创始处理表达按固定核反解至平均防御阈值附近，保留原站立方差，作为机制参考；未做实测标定，初始化元数据记录全部初始植物—昆虫组合的q与防御响应范围。R并非实测容量，具体依据集中见[参数R1](docs/parameters.md#r1)。人口学基线尚未证明可长期存续，200代是运行目标，实际停止代和状态须从输出确认。

匹配空间机会、固定参考代有限QTL拷贝祖源、候选谱系及网络观察已有最小实现。受控F1/F2/回交执行器、特定跨群队列的持续贡献和独立物种验证仍待完成；候选群不自动成为物种。定义和反例集中见[隔离理论](docs/model.md#isolation)。

当前主要人口学风险：满食害时植物存活sP0=0.5、有效后代基准λP=2，即使必有父本，每株期望后代也至多 `0.5×2×Fm×phi≤1`，成本/竞争使其低于替代；局部严格异交还可能限制父本。这是满食害条件下的诊断，不是所有轨迹必然灭绝的证明；只提高资源不足以解决更新问题。

1. 分开核对无虫植物更新、充分供食幼虫招募、局部配偶与食害反馈，记录条件期望及实际数量，区分生态灭绝、数值错误和计算停止。
2. 按共同阶段、损伤定义、资源与空间尺度审查植物存活、繁殖、竞争及花粉范围。真实系统证据分开列；理论选值明确标记，不能把为存续选择的数字称为实测标定。
3. 预定多种子、对照和停止规则，先做单侧及短程检查，再做200代验证；报告全部重复、删失和不确定性，绘制人口、阶段损失及资源轨迹，不增加人口保底或筛掉灭绝重复。
4. 人口学条件合格后检验共同环境表现、分群稳定性和持续遗传流动；独立性证据不足保留uncertain，不用标签禁配或额外杂交死亡制造物种。

以上人口学基线修改和正式实验尚待实施。早期R=1、20代敏感性及多虫5代的结果属于历史报告记载，本次未找到所指v0.1.8原始验证批次，未重新复核；旧测试数不作为当前检查结论。

本轮实施验收（非正式实验）：89组测试与Ruff检查通过；R=10/30/100各4个种子，共12次运行均完成5代，严格执行验收通过，累计60次已提交转移，资源超用和相对守恒误差均为0。第5代植物76–143、昆虫135–291，3张敏感性诊断图已生成核对。该短程检查不代表长期存续或200代已验证。

文档只维护本页、AGENTS及docs中的三份说明。合并前原文及审计记录保存于 [codex/history](codex/history/context_merge_2026-11-04/manifest.json)，历史材料按需追溯；历史源码与输入已归档，原始结果身份不改，当前创始表达已按参数说明修订。

本轮核实确认旧创始组合的q全饱和，并修订当前理论初态；寄主死亡与同轮幼虫独立、资源上限及拒绝边界仍是明确的模型假设。执行失败与最小转移数现纳入严格验收，累计账本容量按真实保存量预算，当前入口不再注册退役生活史。

当前配置链为`current_base.yaml`→`current_single_ancestor.yaml`→`current.yaml`；两池当前协议与单池干预设计分别保留，所有CSV位于`inputs/current/`。`B`、`rho_voc`、`sigma`等保留字段不是额外的当前化学机制。旧substeps、审计、分析和配置输入已冻结至[codex/history/claude_review](codex/history/claude_review/before_cleanup_20261104T080008Z/RUN_HISTORY.md)，按原路径/SHA256追溯；旧结果须用冻结源码读取，原始结果身份不改。
