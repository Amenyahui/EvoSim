# 项目工作规则与开发约定

## 现行文档

项目介绍、运行与当前计划见 README.md；公式、假设和隔离理论见 docs/model.md；参数与文献见 docs/parameters.md；输出和绘图见 docs/results.md。每份原位更新，不另建版本化 README、计划、Methods 或审计说明。被替代文件归档到 `codex/history/`，保留原文、校验值与路径映射并修复现行引用；归档仅按需追溯，不自动批量加载。配置和入口只维护当前协议；退役源码、配置、输入和测试须先冻结为可独立读取历史结果的完整副本，保留原路径、运行说明与SHA256。仍参与当前实现的文件不可仅按旧版本名删除，实验原始结果不按版本号清理。

## 输出与结果保护

- 开发测试、冒烟、短程机制检查、实施验收、审计及其日志和汇总统一写入 `codex/validation/`，自动测试也可用系统临时目录。
- `results/` 仅供用户明确安排的正式科研实验，显式传入 `--out results/...`；不能因运行模拟就视为正式实验。批量默认路径为 `codex/validation/v<version>/run<N>`，目录不得覆盖。
- 示例和复现命令遵守同样路径约定；迁移产物同步现行引用，保留原始数据、源码指纹和历史结果身份。不得修改历史元数据哈希掩盖差异，不自动搬迁已有实验数据。
- `results/`、`tests/reports/`、`codex/` 被Git忽略；必要使用说明保留在 README/docs，不能只存本地日志。以上约定适用后续版本，无需用户重复提醒。

## 开发与检查

Python ≥3.10、UTF-8、4空格、100字符行宽；Ruff配置以 pyproject.toml 为准。变量/函数 snake_case、类 PascalCase、常量 UPPER_SNAKE_CASE；复用既有辅助函数，注释解释科学假设、单位及限制，公开入口用简短docstring。空 __init__.py 可保留。

改科学规则或默认参数时同步模型、参数说明及相关验收测试。纯整理不改变公式、随机抽样顺序或原始结果；汇总/绘图调整用已有数据验证。测试使用普通assert，不以 `python -O` 执行。

```powershell
python -m pip install -e ".[dev]"
python -m ruff check evo_sim tests example.py
python -m ruff format --check evo_sim tests example.py
python -m tests
```

需要整理时使用 `ruff check --fix` 和 `ruff format`，随后审阅差异。文档合并检查链接、参数等价和原始文件保护即可，不为整理重复运行科研模拟。

## 代码入口

`evo_sim/core/` 负责通用调度。下表路径相对于 `evo_sim/models/plant_insect/`：

| 位置（相对 evo_sim/models/plant_insect） | 职责 |
| --- | --- |
| configs/*.yaml、inputs/current/*.csv | 生活史、版本、参数、物种/QTL及参考面板；数值为理论示例 |
| helpers/config.py、inputs.py | 当前配置展开、输入/域/协议检查；历史合同由冻结副本读取 |
| helpers/agents.py | 植物三模块、资源存量和建立位点；昆虫模块及派生 I/A、性别、出生寄主、队列及基因型 |
| helpers/chemistry.py | 当前 x→D/V 的确定性映射；历史映射由冻结副本读取 |
| helpers/voc.py、ecology.py、demography.py、spatial.py | 选择/拒绝、摄入与成本、植物死亡/资源恢复、空间候选和空位建立 |
| substeps/lifecycle.py | AdvanceGeneration 完成成年交配→逐卵择主→SLiM幼虫→取食/存活→植物保留及繁殖脉冲→成年补充 |
| helpers/genetics.py、slim_backend.py、slim_reproduction.slim | 当前3/8维加性 QTL、原生 addCrossed、标记重组、项目突变及位点传递 |
| helpers/contributions.py | 阶段双亲计数、等位贡献分解；本版生命周期分别调用出生/补充两次 |
| batch.py、simulator.py | 调度、批次、状态、文件及指纹；历史结果不覆盖 |
| audit.py、audit_lifecycle.py | 核对当前遗传/资源/阶段来源/事件映射，区分账本完整性与执行验收 |
| analysis/summarize_lifecycle.py | 新阶段原始账本只读汇总；不增加宏观网络指标 |
| analysis/common.py | 只接收当前协议的常规审计批次，历史结果交由冻结副本读取 |

`helpers/modules.py` 实现功能核和表达成本；`analysis/summarize_lineages.py` 只读计算候选群、来源、事件和网络，不回写交配或遗传。植物跨昆虫世代保留身份、基因型及位置，资源逐轮恢复，繁殖只在预定脉冲发生。`offspring` 和遗传贡献只计新生后代，存活老株另外计数；`plant_proposals` 保存所有传播提议及空位建立结果。整轮资源/遗传/空间检查成功后才提交账本，失败不提交半轮；计算保护停止不截断人口。

`analysis/summarize_*.py` 输出汇总；`plot_diagnostics.py` 是统一绘图入口，提供A–E诊断、F–H谱系及敏感性图。各汇总共享Batch；当前命令及适用协议见 docs/results.md，`summarize_qtl_run`仅提供阶段无关的性状/遗传汇总，阶段生态量由`summarize_lifecycle`读取幼虫账本。汇总写批次 `data/summary/`，绘图可指定 `--out`。

默认Runner和example.py使用`configs/current.yaml`；`current_base.yaml`明示当前协议的完整参数，不继承退役生活史配置。当前输入集中在`inputs/current/`，打包时须包含全部YAML、CSV及SLiM脚本；外部SLiM运行时按README安装。
