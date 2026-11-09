# EvoSim 架构与渐进迁移计划

## 范围与基线

2026-11-09：当前任务为 **Scientific Modification**，以main提交fd9ba73为基线，
科学版本升级为0.1.12。五维植物、十维昆虫继续使用既有QTL/SLiM与分层包。
行为、取食、资源方程分别属于biology/insects和biology/plants；simulation和records接收结果，
audit/lifecycle_checks独立重建，analysis使用新版本语义。植物阶段内部改为
实际食害存活→Z恢复→按恢复资源繁殖，保持原子提交及阶段顺序。
冻结源码和输入清单见[科学设计](ecology_v0.1.12_design.md)。

下文记录2026-11-06已完成的 **Refactoring**。当时基线为 v0.1.11，提交
`20d2496f6e85efe37c02cbf495ea50613f0ec0f8`。2026-11-06 检查时远端
`main` 指向同一提交；工作区已有用户更新的 `AGENTS.md`。

科学方程、参数、配置协议、生命周期、随机流及抽样顺序、SLiM 5.2 脚本、
输出列名和标识符均须保持。结构变更不会提升科学模型版本。

## 当前架构

2026-11-06，用户已批准本文计划；源码迁移与测试组织已完成，测试已归入四个用途目录。
`helpers/` 和 `substeps/` 下的实现已迁出；`core/helpers.py` 继续保存原有通用操作。

```text
evo_sim/core/                     通用调度、文件写入与计算保护异常
evo_sim/models/plant_insect/
    config/                      schema、YAML 入口、CSV 加载
    biology/                     群落数组、共享性状、空间
        genetics/                QTL 状态与原生 SLiM 脚本/传输
        plants/                  化学读出、表达成本、资源与存活
        insects/                 感觉、逐卵择主、取食、繁殖机制
    simulation/                  runtime、固定随机流、世代控制器
        stages/                  七个有明确输入/输出的执行阶段
    records/                     观察、事件、遗传贡献；只记录已有结果
    audit/                       文件/批次验收与独立生命周期重建
        lifecycle_checks/        独立遗传、人口、取食、资源、事件检查
    analysis/
        summaries/               生命周期、谱系、QTL 汇总
        plots/                   诊断绘图
tests/{unit,integration,scientific,acceptance}/
```

当前执行依赖为 `simulator → simulation → biology`；执行层将已有阶段结果交给
`records` 格式化。记录中的阶段类型仅在 `TYPE_CHECKING` 下导入，不构成运行时环。
`audit` 只共享配置和文件设施，独立重建方程与抽样；不导入生产 `biology` 计算。
`ComputationalLimit` 放入通用 `core/runner.py`，避免遗传层依赖运行层。

七阶段依次为 `prepare_adults → oviposit → inherit_offspring → develop_larvae →
update_plants → reproduce_plants → recruit`。阶段用简单 dataclass 显式传递结果，
检查点和原子提交仍由 `AdvanceGeneration.forward()` 管理。原约 450 行的世代函数
现在负责顺序调用和未提交结果的返回。

旧 `analysis/summarize_*.py`、`analysis/plot_diagnostics.py` 和 `audit_lifecycle.py`
保留薄兼容入口，不保留两份机制。旧审计模块现在是 `audit` 包，原模块命令和导出仍可用。
内部导入及测试补丁位置已更新；未承诺继续支持内部 `helpers.*` 导入路径。

暂保留 QTL 状态类、较长的审计逐代重建循环和谱系汇总：这些代码需要共享明确的
跨事件检查状态，继续拆分应另选有独立职责的块，不为缩短行数制造动态状态容器。
遗传后端已有原生突变/重组实现，不另设空的 Python 模块。

2026-11-09 的 review 诊断属于 **Analysis Extension**，继续使用现有结构：
`summaries/lifecycle.py` 读取已审计账本并汇总有限期限结局、资源耗尽与繁殖失败原因；
`summaries/qtl.py` 读取遗传汇总和 D/I 快照，以个体数加权分解类型内与类型间方差。
两者不改变模拟状态、不重新计算生物输出，保留原始记录与既有汇总列。
`configs/review_diagnostics.yaml` 是独立实验配置，默认入口和科学协议仍为 v0.1.11。

## 原始模块地图（迁移前）

下表路径相对于 `evo_sim/models/plant_insect/`，除非另有说明。

| 模块 | 实际职责 | 目标与处理方式 |
| --- | --- | --- |
| `evo_sim/core/` | YAML 合并、注册、调度、状态更新、原子写入、源码摘要 | 保持模型无关；已有 `helpers.py` 暂不因名字单独迁移 |
| `simulator.py`、`batch.py`、`main.py` | 注册模型、增量执行、批次与持久化、入口 | 保持入口与运行语义 |
| `helpers/config.py`（407 行） | 参数 dataclass、域验证、序列化，并延迟调用生物配置验证 | `config/schema.py`，保留延迟导入避免环 |
| `helpers/inputs.py`（221 行） | YAML 入口、CSV 解析、指纹、替换及合并规则 | 分为 `config/loader.py` 与 `config/input_loader.py` |
| `configs/*.yaml`、`config.yaml`、`inputs/current/*.csv` | 当前协议配置及输入 | 保留原内容与相对路径；`config.yaml` 仍是便捷 API 的入口 |
| `helpers/agents.py` | 植物、昆虫和群落 dataclass、数组转换与快照 | `biology/agents.py` |
| `helpers/genetics.py`（554 行） | QTL 输入验证、等位目录、创始抽样、表型、SLiM 传递、遗传汇总 | `biology/genetics/qtl.py` 起步；仅在职责可独立时继续拆分 |
| `helpers/slim_backend.py`、`slim_reproduction.slim` | 引擎发现、输入传输、原生交配及输出验证 | 一起迁入 `biology/genetics/`；修正源码目录层级对应的本地引擎查找路径 |
| `helpers/modules.py`、`traits.py` | 表达/调谐核、功能读出、成本；偏好与限制通道摄入函数 | 功能读出集中到 `biology/traits.py`；偏好/取食函数按实际过程拆出 |
| `helpers/chemistry.py` | 共享化学模块到防御/VOC、参考面板 | `biology/plants/chemistry.py` |
| `helpers/voc.py` | 偏好概率与植物表达成本，实际跨两个职责 | 概率进入昆虫 `sensory.py`，植物成本进入植物机制模块 |
| `helpers/ecology.py`（256 行） | 创始取食诊断、同时择主、比例资源分配、幼虫存活及账本量 | `biology/insects/host_selection.py` 与 `feeding.py`；诊断不新增随机抽样 |
| `helpers/demography.py` | 植物资源恢复/存活、性繁殖配置验证、候选配偶、母本矩 | 植物与昆虫人口过程分别归属；计算保护异常保持共同身份 |
| `helpers/spatial.py` | 周期空间、传播、候选寄主、空位建立 | `biology/spatial.py` |
| `helpers/runtime.py`（287 行） | 随机流、资源向量、初始化、状态转换、失败回滚、最终输出 | `simulation/runtime.py`；随机流可原样提取到 `simulation/random.py` |
| `helpers/observations.py`、`contributions.py` | 参考观察、双亲贡献和等位贡献记录 | `records/`；保留既有输出含义 |
| `substeps/lifecycle.py`（705 行） | 亲本抽样、快照、一次世代全部过程、记录与提交 | `simulation/lifecycle.py` 与真实阶段、记录模块 |
| `audit.py`（166 行） | 文件完整性、批次审计、执行验收、命令行 | `audit/run.py`，包导出既有 API 与模块命令 |
| `audit_lifecycle.py`（1411 行） | 独立重建遗传、取食、人口、空间、事件与摘要 | `audit/lifecycle.py` 与独立检查模块；不导入被检查的生产生物方程 |
| `analysis/common.py`、`lineage_observations.py` | 只读批次接口、机会基线和祖源分析 | 保持职责；祖源函数较长但没有证据表明须立即拆分 |
| `analysis/summarize_*.py` | 生命周期、谱系与 QTL 汇总 | `analysis/summaries/`；旧命令使用薄兼容入口 |
| `analysis/plot_diagnostics.py`（1174 行） | A–H 与敏感性绘图 | `analysis/plots/diagnostics.py`；暂不按图逐一碎片化 |
| `analysis/diagnose_initial_conditions.py`、`lineages.yaml` | 初始条件批次诊断、谱系分析配置 | 保留入口与配置路径 |
| `tests/test_*.py`（26 个模块） | 单函数、跨模块、科学边界、完整运行与绘图混合 | 按测试目的渐进归类；保留所有仍检查当前机制的测试 |
| `tests/support.py`、`fixtures/`、`__main__.py` | 共用配置、断言、连锁图、完整测试入口 | 保持 `python -m tests`；归类时改为递归发现，避免漏跑 |

## 依赖与实际问题

迁移前入口链为：`batch/simulator → helpers/runtime + substeps/lifecycle → helpers 生物机制`；
审计由 `audit → audit_lifecycle` 独立重建；分析读取批次和审计文件。
`core/` 当前没有导入植物昆虫模型。生物机制当前未依赖执行、分析或审计层。
配置验证对遗传/空间/繁殖验证使用函数内导入，迁移不能直接改成全部顶层导入。

必须解决的问题：

1. `helpers/` 的命名掩盖了机制、记录和执行三种不同依赖层。
2. `transition()`（225–674 行）同时执行所有过程并构造全部事件记录，难以核对阶段边界。
3. `audit_run_lifecycle()`（378–1411 行）跨多种检查；嵌套亲本/后代重建还捕获大量局部状态，不能机械切割。
4. 文件移动会影响 SLiM 本地开发回退路径、包数据、源码指纹和测试中的补丁位置，不能只替换导入字符串。
5. 当前测试发现仅扫描顶层；直接移动测试会静默漏测。

可选改善：逐图拆分绘图模块、进一步缩短谱系汇总、重命名 `core/helpers.py`。
这些没有改变科学行为的必要性，不与第一批迁移混做。

重复计算检查：偏好概率在 `traits` 与 `voc` 中分别实现 outside/forced 两种合同，
不能因相似而合并改变边界行为；`modules.effective()` 是昆虫功能读出的共同实现，
应继续复用。审计中的距离、化学读出、随机过程及资源方程重建是刻意独立，必须保留。
本轮检查未确认需要修正的生产方程重复或科学 bug。

## 必须保持的执行细节

实际代码不是七段可以任意调换的模板。其顺序为：成年母本成本/产卵计数 →
逐卵择主 → 为被接受卵抽父本 → 昆虫 SLiM 遗传 → 同时取食 →
幼虫存活、性别、成年传播 → 植物资源/存活与传播子计数 →
传播及空位抽签 → 只对建立成功的植物抽父本并执行 SLiM →
拼接保留植株和新生个体 → 检查/构造全部记录 → 原子提交。

`sample_parents()` 同时用于植物和昆虫，植物分支还抽存活随机数；拆阶段不能重复调用。
观察在世代转移前执行，终末世代仅观察。`AdvanceGeneration.forward()` 保留随机流检查点，
失败走现有回滚路径，禁止提交半轮账本。

## 分批迁移计划

本计划先提交供审阅，已获用户批准。依据 `AGENTS.md` 第 10 节，大批文件迁移在批准后执行。
每批完成相关测试和固定种子比较后再进入下一批，不一次改写整个仓库。

| 批次 | 具体变更 | 验收 |
| --- | --- | --- |
| 0：检查与基线 | 本文、开发说明、模块/依赖清单；记录原提交、源码 SHA256、环境和基线输出 | 原有 lint、format、完整测试；真实 SLiM 固定种子多场景基线 |
| 1：遗传 | 迁入 `biology/genetics/`，保留 `GeneticState`/`Allele` 身份与接口；同步打包、引擎路径和调用方 | 遗传单元/集成、突变、连锁与 SLiM 传输测试；完整轨迹比较 |
| 2：配置与性状 | YAML/CSV 与 schema 分离；移动 agents、化学、共享性状；保留配置协议和输入文件 | 配置/CSV/共享化学/模块测试；基线比较 |
| 3：生态与空间 | 依次拆择主、取食、植物/昆虫人口与空间；不制造空模块 | 边界、资源守恒、配偶、空间、科学合同测试；基线比较 |
| 4：执行与记录 | runtime、随机流、观察/贡献归位；按真实边界逐段提取阶段，最后简化控制器 | 增量步进、失败回滚、存储保护、ID 与全部中间记录逐项比较 |
| 5：独立审计与分析 | 独立检查先提取可独立子块，再组织 audit 包；汇总/绘图归类；保留旧命令 | 重新散列的篡改测试、独立审计测试、历史输出列兼容、绘图 CLI |
| 6：测试组织与交付 | 递归发现后渐进移动；更新 README/results/development/architecture；清理已确认冗余文件 | 三项规定的完整检查、构建包检查和固定种子最终比较；审阅 diff 并提交 |

兼容策略：保留公开 `Config/get_runner/simulate`、CLI、文件格式；旧模块若被现有调用方
或测试补丁引用，使用显式薄导出或更新内部引用，不保留两份生物实现。
SLiM 已实现突变和重组，不为目标目录树另造 Python `mutation.py/recombination.py`。

## 历史与清理边界

用户已授权删除不需要的历史文件，并说明 Git 有备份，因此不强制额外复制全部退役源文件。
Git 中的基线提交和 tag 不删除、不重写。名称带旧版本号的测试仍有当前机制覆盖价值，
不能按文件名删除。`codex/slim_runtime/` 是可运行依赖，必须保留。

`build/` 是构建副本，可重新生成；Python/Ruff 缓存可清理。
`codex/history/` 含 Git 忽略的本地快照和分析，远端提交一致并不证明其中全部内容都已备份；
只有逐项核实被 Git 包含或内容冗余后才删除。`results/`、原始 CSV 和当前配置不作为废旧文件处理。

本轮批量删除与单独删除 `build/` 均被自动审批拒绝（`blocked by policy`），因此
未清理本地历史目录或生成目录。已生成源码清单供核查：部分本地历史文件只存在于
本地 Git/tag 或忽略目录，不能根据 `origin/main` 一致就断言已远端备份。
退役的受版本控制实现通过迁移删除，仍可从原 Git 提交取回。
