# evoSim 0.1.2

植物—植食性昆虫的个体基（agent-based）共进化模拟器。每一代中，昆虫在周期空间栅格里局部遇见植物、按气味线索选择寄主、共享取食资源，植物与昆虫再经历存活、密度依赖繁殖、QTL 遗传、突变与扩散，提交下一代；可连续运行数百代，观察种群动态与性状的演化响应。

技术要点：

- **纯 NumPy + YAML / CSV 输入驱动**：YAML 声明流程和通用参数，CSV 保存物种、遗传结构与初始组成，代码实现方程和输入检查。
- **位级可复现**：随机流按 `阶段 × 类群 × 物种` 分层播种，同配置同种子完全复现。
- **不是深度学习框架**：代码组织参考 AgentTorch 风格（Runner / Controller / Registry，组件实现 `forward()`），但没有 PyTorch 后端，不提供 GPU、自动求导或参数校准。
- **参数为理论示例**：所有数值配置均标注 illustrative，不是实测标定；模拟输出不能作为校准数据。

本版接通个体联合寄主反应、双亲遗传贡献、寄主贡献网络和等位频率分解。保留既有生态与遗传生成规则，明确当代昆虫自主选择、取食和繁殖的生活史；物种标签仍固定。依据与边界见 [v0.1.2 计划](docs/v0.1.2_plan.md)、[最小模型规范](docs/v0.1.2_minimal_model.md) 及 [实施检查报告](docs/v0.1.2_causal_check_report.md)。

## 快速上手

在本目录打开终端：

```powershell
# 1. 安装（Python ≥ 3.10，依赖：numpy / pandas / matplotlib / PyYAML）
python -m pip install -e .

# 2. 最小示例：默认配置（sexual_qtl_development，5 代）跑一遍并打印状态
python example.py

# 3. 小规模批次：2 代、独立遗传/无新突变，各两个重复
python -m evo_sim --config evo_sim/models/plant_insect/configs/v012_base.yaml --out results/v012_example
```

省略 `--out` 时输出到 `results/results_v0.1.2/run<N>` 的下一个空闲编号；显式输出目录必须不存在，已有结果不会覆盖。

实验跑完后，建议依次做独立审计和汇总作图：

```powershell
# 4. 独立重建取食/人口/遗传账本，核对数值一致性（不依赖运行时函数）
python -m evo_sim.models.plant_insect.audit --input results/v012_example

# 5. 汇总逐代人口/性状/等位基因
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input results/v012_example

# 6. 生成已有诊断图（PNG + PDF）
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input results/v012_example

# 7. 本版小规模因果对照：20 次运行，逐次独立审计
python -m evo_sim.models.plant_insect.analysis.check_causal_chain --out results/v012_causal_example
```

`results/results_v0.1.0/` 和 `results/results_v0.1.1/run4` 属于历史版本，不属于本版验证。历史结果应使用生成它们时的代码审计；不改旧元数据以通过新审计。小规模因果检查只验证规则与对照设计，不证明专化、共存或共同演化。

其他常用命令：

```powershell
python -m tests                                        # 全部自动检查（无需 pytest）
```

`python -m evo_sim` 等价于安装后提供的 `evosim` 命令。本版只接受 `sexual-qtl-1`；旧协议会报错，不自动迁移。

## 模拟内容：一代发生了什么

每代按三个子步骤执行（顺序与状态路径声明在 `evo_sim/models/plant_insect/config.yaml`）：

1. **选寄主与取食**：植物按防御性状 $D$ 和独立气味性状 $Z$ 释放 VOC 线索；昆虫在候选池（空间模式下为局部邻域随机抽取）中按化学得分 softmax 选择寄主，允许"不选择"；被接受的个体按 $d = d_{\max} q$ 产生需求，共享植株资源并按比例缩放实际转移量。
2. **存活与繁殖**：植物存活受取食负荷 $h$ 影响，昆虫存活取决于实际取食满足率 $g$；成本乘子 F 连续作用于母源后代均值。存活且有配偶的植物母本/雌虫按条件 Poisson($\mu = \lambda F \varphi$) 抽后代，植物繁殖另受局部密度限制。父本经 `F>0` 等资格筛选后均匀抽样，该成本门槛仍是未获直接实证支持的假设。
3. **遗传与扩散**：真实双亲各形成一份配子，同连锁群按 Haldane–Poisson 重组、不同组独立分离；逐 QTL 突变后由新基因型重算性状，再从母方位置扩散。植物严格异交，昆虫仅雌虫产子；最终代不再繁殖。

完整参数说明与理论公式（含文献出处）见 [docs/model.md](docs/model.md)；假设与边界见 [docs/assumptions.md](docs/assumptions.md)。`S×F` 仅为表现代理量；实际一代遗传贡献为母源、父源后代数之和的一半。

## 内置实验配置

通用参数和输入路径在 [model_parameters.yaml](evo_sim/models/plant_insect/configs/model_parameters.yaml)，物种数量、资源、繁殖率、竞争、位点、等位效应及初始频率在 [inputs/](evo_sim/models/plant_insect/inputs/) 的六张 CSV 表中。主 YAML 从 908 行缩减到 139 行，原数值不变。Python 不提供模型参数默认值；缺少必填项会报错。实验 YAML 继承公共配置，只写差异，用 `--config` 选择。CSV 列含义、单位、约束和覆盖示例见 [参数配置指南](docs/parameters.md)。

| 配置 | 代数 | 重复 | 处理 | 用途 |
| --- | --- | --- | --- | --- |
| `sexual_qtl_smoke.yaml` | 1 | 1 | 突变开 / 关 | 最小批次与审计检查 |
| `sexual_qtl_development.yaml` | 5 | 1 | 突变开 / 关 | 最小实现示例：2 种植物（各 200 株）+ 1 种昆虫（40 个），QTL 遗传 |
| `sexual_qtl_200gen.yaml` | 200 | 3 | 突变开 / 关 | 上例的长程延伸；不用于假说检验，安全上限可提前停止 |
| `initial_conditions.yaml` | 0 | 4 | 7 组初始种群数量 | 初始资源与密度诊断；原 Python 中的实验设定已移入 YAML |
| `v012_base.yaml` / `v012_linked.yaml` | 2 | 2 | 独立/无新突变/连锁 | 同一创始输入的架构与变异来源检查 |
| `v012_mirror.yaml` / `v012_formation.yaml` / `v012_multi.yaml` | 2 | 2 | 重标记/形成/2P×2H | 比较边界见参数指南 |
| `v012_*intervention*.yaml` | 1 | 2 | 对方性状组成、搜寻尺度 | 一代路径检查，保持焦点初始基因型一致 |
| `v012_causal.yaml` | 1–2 | 2 | 汇集上述 10 个场景 | 使用 `analysis.check_causal_chain` 入口 |

批量运行器（`batch.py`）按"处理 × 重复"逐一模拟，输出逐代 CSV、NPZ 快照、等位基因目录和运行状态；**已有结果目录不会被覆盖**，换输出目录即可。

## 目录结构

```text
evoSim/
├── evo_sim/
│   ├── core/                  通用调度：Runner / Controller / Registry、配置与状态路径
│   └── models/plant_insect/   植物—昆虫模型
│       ├── config.yaml        子步骤顺序、状态与输入输出路径声明
│       ├── configs/           实验参数 YAML（每个实验一个文件）
│       ├── inputs/            物种、竞争、位点、等位效应、初始频率和参考性状 CSV
│       ├── simulator.py       组件注册与运行生命周期
│       ├── substeps/          各阶段 forward()：生态 / 人口 / 遗传与扩散
│       ├── helpers/           方程、个体数组、QTL 规则、参数校验
│       ├── analysis/          汇总、绘图、初始条件诊断（统一 --input/--out）
│       ├── batch.py           批量实验与数据保存
│       └── audit*.py          独立账本重建与完整性检查
├── docs/                      方程、遗传假设、空间规则、数据字典
├── tests/                     断言式自动检查（无 pytest 依赖）与检查报告
├── results/results_v<版本>/run<N>/     模拟输出：审计表 CSV、snapshots/、settings.json、metadata.json及 audits/ 审计报告
├── example.py                 Runner API 最小示例
└── pyproject.toml             安装信息与命令入口
```

## 结果如何验证

- **`audit.py` / `audit_genetics.py`**：从保存的快照、亲子表、等位基因目录独立重算表型、资源、适合度、繁殖与空间账本，不调用运行时方程函数；篡改任一输出都会失败。审计是内部一致性与实现检查，不是生物学验证。
- **`python -m tests`**：布局、资源守恒、空间交配、重组、双亲谱系、回滚与可复现性断言检查；报告在 `tests/reports/`。

## 修改入口

- **改公共模型参数**：编辑 `evo_sim/models/plant_insect/configs/model_parameters.yaml`，换新输出目录重跑。
- **改物种与遗传输入**：编辑 `inputs/species.csv`、`competition.csv`、`loci.csv`、`alleles.csv`、`founder_frequencies.csv`；观测参考面板在 `reference_traits.csv`。CSV 可用 Excel 编辑并保存为 UTF-8。
- **改实验差异**：编辑 `configs/` 下对应实验 YAML；嵌套字典逐项覆盖，列表整体替换。
- **改流程顺序 / 状态路径**：`evo_sim/models/plant_insect/config.yaml`。
- **改个体行为**：`substeps/` 中对应类的 `forward()`；通用数值计算在 `helpers/`。
- **新增阶段**：实现一个 Substep 类 → 在 `simulator.get_registry()` 注册 → 在 YAML 中引用 generator 名称。

## 文档索引

| 文档 | 内容 |
| --- | --- |
| [docs/model.md](docs/model.md) | 论文式模型说明：ODD、全部参数、理论公式（F01–F16 + QTL）与参考文献 |
| [docs/assumptions.md](docs/assumptions.md) | 全部建模假设与边界（生态/取食/嗅觉/空间/遗传/标定验证），逐条标注类别 |
| [docs/parameters.md](docs/parameters.md) | YAML 通用参数、CSV 列说明、输入路径与实验覆盖规则 |
| [docs/v0.1.2_plan.md](docs/v0.1.2_plan.md) | 六步计划、13 篇参考文献及证据边界 |
| [docs/v0.1.2_minimal_model.md](docs/v0.1.2_minimal_model.md) | 生活史、预算、双亲贡献与寄主关系规范 |
| [docs/v0.1.2_causal_check_report.md](docs/v0.1.2_causal_check_report.md) | 实现范围、小规模检查及未解决问题 |
| [docs/results.md](docs/results.md) | 输出内容、指标定义与 `results/results_v<版本>/run<N>` 目录文件契约 |
| [docs/architecture.md](docs/architecture.md) | AgentTorch 风格框架的组件对应关系与每代调度顺序 |

默认生态数值沿用 v0.1.1 的探索性配置，历史修改见 [修正说明](README_v0.1.1.md)。本次未按预期演化方向调参，也未重跑 200 代实验。

代码格式、目录职责与检查命令见 [开发约定](CONTRIBUTING.md)。
