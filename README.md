# evoSim

植物—植食性昆虫共同演化个体模型：母体择主、幼虫资源竞争、成年空间交配与双亲 QTL 遗传共同决定代际变化。软件与模型0.1.9，生活史协议`perennial-plants-larval-recruitment-1`，遗传/功能协议`slim-qtl-4` / `expression-tuning-2`。植物跨昆虫代留存、恢复叶量并按独立脉冲繁殖；昆虫保持非重叠代。

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
python -m evo_sim --config evo_sim/models/plant_insect/configs/current.yaml --out codex/validation/perennial_example
python -m evo_sim.models.plant_insect.audit --input codex/validation/perennial_example --require-success --min-transitions 1
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input codex/validation/perennial_example --sensitivity
```

该入口仅检查最大叶量K=10的基线，4个种子、5个昆虫世代；当前资源存量不再逐代重置。200代延长配置为 `configs/current_200gen.yaml`，同样仅运行K=10、4个重复。K=30/100退出默认验证，既有记录保留。指定目录须不存在；省略 `--out` 默认放入 `codex/validation/v<version>/run<N>`。正式科研实验须由用户明确安排并显式指定 `--out results/...`。

| 文档 | 内容 |
| --- | --- |
| [模型与理论](docs/model.md) | 方程、生活史、假设边界、寄主隔离与谱系网络 |
| [参数与文献](docs/parameters.md) | 输入、取值、单位、分系统证据及未标定项；附[逐字段清单](docs/parameter_inventory.csv) |
| [输出与绘图](docs/results.md) | 文件字段、审计、汇总及A–H/敏感性图 |
| [开发规则](AGENTS.md) | 代码入口、检查和结果保护 |

<a id="current-work"></a>

## 当前状态与下一步

当前单池P40/H40，10×10格、每格1个共享建立位点；公共两池基础P1/P2各40株。植物最大叶量K=10、单幼虫需求dmax=1；创始满叶，新株以.25K建立。植物幸存者保留ID、位置和基因型，食后叶量按有界恢复式更新；每5个昆虫世代发生一次植物繁殖脉冲，提案只在真正空位中建立。新株初始少叶也会承担额外建立风险；该叶量亏损不全来自虫害，尚未区分成熟体型和幼株储备。原sP0、alpha和母源密度φ退役，空位不强填，植物和昆虫均允许灭绝。数值均为理论参考，机制依据与局限见[模型](docs/model.md)和[参数P1/P2](docs/parameters.md#p1-多年生存活恢复和独立繁殖时间)。

创始处理表达仍使用固定核下反解的阈值参考，关闭新突变、保留站立变异与重组；并未重新用种群存续结果拟合。叶量恢复、背景及损伤危害率、5轮繁殖钟均无实测标定。匹配空间机会、有限QTL祖源及候选谱系网络是只读观察；本轮先解决人口学时间尺度，不加入额外寄主禁配、杂交死亡或新的双寄主实验。

**软件实现验收通过，长期共存人口学基线未合格。** 101组测试及71个Python文件的Ruff检查通过；最终安装包60项文件、43个模块和6个配置通过隔离校验。6次10代短程运行通过严格审计，共60次转移。无虫植物对照2次均完成200代，期末96/97株，通过每次至少200转移的严格审计；这支持无虫植物基线可执行，不代表双方长期共存。

此前200代耦合验证共12次，结果如下；各格中的顺序对应4个重复。该批次保留原配置与结果身份，后续默认仅验证K=10：

| 最大叶量 | 实际结局 |
| --- | --- |
| K=10 | 昆虫分别在64/96/166/87代灭绝；前三次继续至200代，期末植物1/96/1株；第四次在87代双方全灭。 |
| K=30 | 分别在112/119/157/83代双方全灭。 |
| K=100 | 分别在37/14/13/71代计算截断：37/71代达到400万累计基因型值限制，14/13代潜在卵超过5000；不是已知灭绝。 |

因此没有一次已验证双方共同存续至200代。主批次12次运行的完整账本审计已通过，共1293次已提交转移，资源超用及相对守恒误差均为0；但3次完成200代、5次全灭终止、4次计算截断，未满足全部运行至少200次转移的要求，`execution_success=false`、`acceptance_passed=false`。账本完整与200代执行验收不等价，不能据审计通过宣称长期共存。

当前读取器为重放生产端已有的资源舍入预算纠偏而修复了`audit_lifecycle.py`，因此`source_matches_current=false`；生态、遗传及空间执行核未变。20次运行均与[冻结执行包](codex/history/perennial_run_0.1.9_20261105/READ_FROZEN_SOURCE.txt)核对，原metadata、执行源码和输出哈希保持原样。终点、完整审计和只读人口图见[输出结果](docs/results.md#perennial-validation)。

人口学条件例提示时标紧张：若持续5轮满叶量亏损、F固定.86且配偶充分，旧株留存概率约.4755；以建立比例a计，每周期期望替代因子为`.4755×(1+2×.86×a)`，维持期望替代至少需a约.6413。实际K=10/30在昆虫存在期间建立比例约.41–.54、亏损均值约.74–.91，且配偶机会不足；这些均值不能直接代入条件式证明因果，也不是通用存续阈值。它只为下一步预先敏感性设计提供线索，推导及限制见[模型](docs/model.md#52-植物存活繁殖脉冲及空位建立)。

1. 预先声明时标、背景/亏损死亡、恢复、繁殖脉冲及空间配偶机会的敏感性范围、重复与停止规则，检验限制来自哪一阶段；不按存活结果挑参数，不增加人口保底或把计算截断当灭绝。
2. 为真实寄主差异设计分别明确化学、可食存量、取食窗口和人口学证据，理论值与实测值分开；不将不同自然系统拼成一套标定。
3. 当前不进入成种实验证据阶段。人口学基线合格后，再检验共同环境表现、分群稳定性和持续遗传流动；不通过额外禁配或杂交死亡制造独立群体。

文档只维护本页、AGENTS及docs中的三份说明。合并前原文及审计记录保存于 [codex/history](codex/history/context_merge_2026-11-04/manifest.json)，历史材料按需追溯；历史源码与输入已归档，原始结果身份不改，当前创始表达已按参数说明修订。

本轮核实确认旧创始组合的q全饱和，并修订当前理论初态；寄主死亡与已完成摄入幼虫独立、资源期望上限及默认拒绝边界仍为明确假设；植物自身改为跨轮存活与有界叶量恢复。执行失败与最小转移数现纳入严格验收，累计账本容量按真实保存量预算，当前入口不再注册退役生活史。

当前配置链为`current_base.yaml`→`current_single_ancestor.yaml`→`current.yaml`；两池当前协议与单池干预设计分别保留，所有CSV位于`inputs/current/`。`B`、`rho_voc`、`sigma`等保留字段不是额外的当前化学机制。旧substeps、审计、分析和配置输入已冻结至[codex/history/claude_review](codex/history/claude_review/before_cleanup_20261104T080008Z/RUN_HISTORY.md)，按原路径/SHA256追溯；旧结果须用冻结源码读取，原始结果身份不改。
