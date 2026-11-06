# evoSim

植物—植食性昆虫共同演化个体模型。当前开发协议为0.1.11：遗传 `slim-qtl-5`、功能 `expression-tuning-2`、生活史 `resource-coupled-plants-larval-recruitment-1`。植物每轮取食后存活、繁殖和恢复，植物全局随机异交并通过SLiM双亲遗传；昆虫保持SLiM有性遗传及非重叠世代。数值为理论探索，不是自然系统标定。

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

<a id="current-work"></a>

## 当前状态与下一步

植物每轮按食后资源亏损抽存活，幸存株按表达成本和剩余资源产生扩散提案；空位抽签决定建立，资源逐轮恢复。类型是母系初始来源标签，不表示植物间生殖隔离或固定后代化学类型；每格3个位点只限制建立，不表示共享光、水或养分竞争。每株幸存植物均可作为其他幸存植物的父本，不按距离或初始类型限制，排除自交；后代继承双亲基因并可重组和突变。默认十类型配置关闭双方新突变作为机制对照，没有新等位变异，选择和重组仍可改变性状；进化速率从结果估计，不设置双方速度参数。昆虫的空间交配和双亲遗传保持不变。公式与边界见[模型](docs/model.md)，默认值和证据见[参数](docs/parameters.md)。

植物繁殖基准λP仍为2，但现在以每轮为单位，不能与旧版每5轮的数值等同；未通过拟合参数确保共存。**0.1.11已通过109组测试及8次五轮开发验收，40次转移独立审计通过，资源超用及收支误差为0。** 记录、绘图和安装包校验见[当前验收](docs/results.md#current-validation)，历史结果不构成本版验证。 固定期限存续不能证明长期共存或共同演化。

0.1.10及以前的结果不迁移、不改哈希；简化前源码、配置、输入、测试、路径映射和SHA256见[冻结包](codex/history/before_simplification_0.1.11_20261106T034842Z/README.md)。历史0.1.9/0.1.10验证及读取边界见[输出](docs/results.md#feedback-validation)。
