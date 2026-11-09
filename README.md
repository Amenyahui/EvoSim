# evoSim

植物—植食性昆虫共同演化个体模型。当前科学版本0.1.12：遗传slim-qtl-6、功能expression-tuning-contact-3、生活史dose-tolerance-lifecycle-2。植物五维PSM1/PSM2/Z/VOC1/VOC2，昆虫十维含两个接触响应坐标。实际摄入后计算化学处理、毒性与抗营养存活，Z促进植物叶片恢复，恢复后资源参与繁殖。新系数均为显式理论假设，未按共存或专化结果拟合。

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
python -m evo_sim --config evo_sim/models/plant_insect/configs/current.yaml --out codex/validation/resource_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/resource_example --require-success --min-transitions 1
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/resource_example --sensitivity
```

默认入口使用10种初始化学类型，每型20株，共200株；10×10格、每格3个建立位点，H1初始200只。指定输出目录须不存在；省略 `--out` 默认放入 `codex/validation/v<version>/run<N>`。正式科研实验须由用户明确安排并显式指定 `--out results/...`。

| 文档 | 内容 |
| --- | --- |
| [模型与理论](docs/model.md) | 方程、生活史、假设边界、寄主隔离与谱系网络 |
| [参数与文献](docs/parameters.md) | 输入、取值、单位、分系统证据及未标定项；附[逐字段清单](docs/parameter_inventory.csv) |
| [输出与绘图](docs/results.md) | 文件字段、审计、汇总及A–H/敏感性图 |
| [开发规则](AGENTS.md) | 代码入口、检查和结果保护 |
| [代码架构](docs/architecture.md) | 模块职责、依赖边界、迁移前后对应关系 |
| [开发与验证](docs/development.md) | 完整测试、固定种子比较、SLiM 与打包检查 |

代码按 `config/`、`biology/`、`simulation/`、`records/`、`audit/` 和
`analysis/` 分层，生命周期在 `simulation/stages/` 中依次执行。测试按用途分为
`unit/`、`integration/`、`scientific/` 和 `acceptance/`，仍使用 `python -m tests`
运行全部检查。本次是科学机制修改；固定种子在0.1.12内可复现，新旧版本轨迹不要求相同。SLiM核心脚本、随机流种子派生、空间与事务提交机制保持。原分析命令继续可用，历史数据用冻结源码读取。

<a id="current-work"></a>

## 当前状态与下一步

本版取消PSM与VOC固定关联，母虫按空间候选→VOC嗅觉→接触PSM→接受/拒绝的顺序产卵，
接触拒绝后可继续剩余候选。I改变摄入后的残余毒性，不限制实际食物需求。
实际取食损伤为X/B_before，Z仅改变恢复速率，不自动承担固定成本。
类型仍为母系初始来源标签；植物全局异交、昆虫同兼容池局部交配和观察谱系规则保持。
不增加物种形成或隔离机制。默认十类型对照仍关闭双方新突变。

详细方程见[模型](docs/model.md)，参数值、单位及证据见[参数](docs/parameters.md)，
修改论证见[设计](docs/ecology_v0.1.12_design.md)，实际测试和小规模运行见[输出与验证](docs/results.md)。
历史0.1.11的验收不构成0.1.12验证，固定期限存续也不证明长期共存或共同演化。

0.1.10及以前的结果不迁移、不改哈希；简化前源码、配置、输入、测试、路径映射和SHA256见[冻结包](codex/history/before_simplification_0.1.11_20261106T034842Z/README.md)。历史0.1.9/0.1.10验证及读取边界见[输出](docs/results.md#feedback-validation)。
