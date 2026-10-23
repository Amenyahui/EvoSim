# evoSim 0.1.1

植物—植食性昆虫的个体基（agent-based）共进化模拟器。每一代中，昆虫在周期空间栅格里局部遇见植物、按气味线索选择寄主、共享取食资源，植物与昆虫再经历存活、密度依赖繁殖、QTL 遗传、突变与扩散，提交下一代；可连续运行数百代，观察种群动态与性状的演化响应。

技术要点：

- **纯 NumPy + YAML 配置驱动**：流程顺序、状态路径和实验参数全部由 YAML 声明，代码只实现方程。
- **位级可复现**：随机流按 `阶段 × 类群 × 物种` 分层播种，同配置同种子完全复现。
- **不是深度学习框架**：代码组织参考 AgentTorch 风格（Runner / Controller / Registry，组件实现 `forward()`），但没有 PyTorch 后端，不提供 GPU、自动求导或参数校准。
- **参数为理论示例**：所有数值配置均标注 illustrative，不是实测标定；模拟输出不能作为校准数据。

## 快速上手

在本目录打开终端：

```powershell
# 1. 安装（Python ≥ 3.10，依赖：numpy / pandas / matplotlib / PyYAML）
python -m pip install -e .

# 2. 最小示例：默认配置（sexual_qtl_development，5 代）跑一遍并打印状态
python example.py

# 3. 完整实验：200 代、两个处理（QTL 突变开/关），输出到新的 runN（当前最终批次为 run4）
python -m evo_sim --config evo_sim/models/plant_insect/configs/sexual_qtl_200gen.yaml
```

输出目录统一命名为 `results/results_v<版本>/run<N>`：一个版本一个文件夹，省略 --out 时自动取下一个 run 编号；--out 可指定新目录。当前只保留最终批次 `run4`（200 代，两处理 × 3 重复）；此前的 run1–run3 已按用户要求删除，其状态与删除原因记录在本地修正报告（`codex/mating_correction_report.md`，不随 Git 发布）。下次默认创建 run5。

实验跑完后，建议依次做独立审计和汇总作图：

```powershell
# 4. 独立重建取食/人口/遗传账本，核对数值一致性（不依赖运行时函数）
#    run4 的绘图源码在运行期间更新过，仅对两个绘图文件显式放行并记录
python -m evo_sim.models.plant_insect.audit --input results/results_v0.1.1/run4 --allow-plot-source-drift

# 5. 汇总逐代人口/性状/等位基因，输出到 results/results_v0.1.1/run4/data/summary/
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input results/results_v0.1.1/run4

# 6. 生成结果图（PNG + PDF）到 results/results_v0.1.1/run4/figures/
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input results/results_v0.1.1/run4
```

历史克隆模型结果仍保留在 `results/results_v0.1.0/run1/`（不同版本，不适用新审计，不可与新曲线混用）。本版的最终验证输出为 `results/results_v0.1.1/run4`（200 代，六个重复全部通过审计）。run1–run3 已删除，状态记录见本地修正报告（`codex/mating_correction_report.md`，不随 Git 发布），原始日志仍在 `codex/validation/`。

其他常用命令：

```powershell
python -m tests                                        # 全部自动检查（无需 pytest）
```

`python -m evo_sim` 等价于安装后提供的 `evosim` 命令。本版只接受 `sexual-qtl-1`；旧协议会报错，不自动迁移。

## 模拟内容：一代发生了什么

每代按三个子步骤执行（顺序与状态路径声明在 `evo_sim/models/plant_insect/config.yaml`）：

1. **选寄主与取食**：植物按防御性状 $D$ 和独立气味性状 $Z$ 释放 VOC 线索；昆虫在候选池（空间模式下为局部邻域随机抽取）中按化学得分 softmax 选择寄主，允许"不选择"；被接受的个体按 $d = d_{\max} q$ 产生需求，共享植株资源并按比例缩放实际转移量。
2. **存活与繁殖**：植物存活受取食负荷 $h$ 与防御/VOC 成本影响，昆虫存活取决于实际取食利用率 $g$；存活且有合格同种配偶的植物母本/雌虫按条件 Poisson($\mu = \lambda F \varphi$) 抽后代，植物繁殖受密度因子 $\varphi = (1 + \sum_r \alpha_{sr} N_r)^{-1}$ 抑制。
3. **遗传与扩散**：真实双亲各形成一份配子，同连锁群按 Haldane–Poisson 重组、不同组独立分离；逐 QTL 突变后由新基因型重算性状，再从母方位置扩散。植物严格异交，昆虫仅雌虫产子；最终代不再繁殖。

新版本设计见 [README_v0.1.1.md](README_v0.1.1.md)。完整参数说明与理论公式（含文献出处）见 [docs/model.md](docs/model.md)；全部建模假设与边界见 [docs/assumptions.md](docs/assumptions.md)。

## 内置实验配置

实验参数在 `evo_sim/models/plant_insect/configs/`，一个 YAML 一个实验，用 `--config` 选择：

| 配置 | 代数 | 重复 | 处理 | 用途 |
| --- | --- | --- | --- | --- |
| `sexual_qtl_smoke.yaml` | 1 | 1 | 突变开 / 关 | 最小批次与审计检查 |
| `sexual_qtl_development.yaml` | 5 | 1 | 突变开 / 关 | 最小实现示例：2 种植物（各 200 株）+ 1 种昆虫（40 个），QTL 遗传 |
| `sexual_qtl_200gen.yaml` | 200 | 3 | 突变开 / 关 | 上例的长程延伸；不用于假说检验，安全上限可提前停止 |

批量运行器（`batch.py`）按"处理 × 重复"逐一模拟，输出逐代 CSV、NPZ 快照、等位基因目录和运行状态；**已有结果目录不会被覆盖**，换输出目录即可。

## 目录结构

```text
evoSim/
├── evo_sim/
│   ├── core/                  通用调度：Runner / Controller / Registry、配置与状态路径
│   └── models/plant_insect/   植物—昆虫模型
│       ├── config.yaml        子步骤顺序、状态与输入输出路径声明
│       ├── configs/           实验参数 YAML（每个实验一个文件）
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

- **改实验参数**：编辑 `evo_sim/models/plant_insect/configs/*.yaml`，换新输出目录重跑。
- **改流程顺序 / 状态路径**：`evo_sim/models/plant_insect/config.yaml`。
- **改个体行为**：`substeps/` 中对应类的 `forward()`；通用数值计算在 `helpers/`。
- **新增阶段**：实现一个 Substep 类 → 在 `simulator.get_registry()` 注册 → 在 YAML 中引用 generator 名称。

## 文档索引

| 文档 | 内容 |
| --- | --- |
| [docs/model.md](docs/model.md) | 论文式模型说明：ODD、全部参数、理论公式（F01–F16 + QTL）与参考文献 |
| [docs/assumptions.md](docs/assumptions.md) | 全部建模假设与边界（生态/取食/嗅觉/空间/遗传/标定验证），逐条标注类别 |
| [docs/results.md](docs/results.md) | 输出内容、指标定义与 `results/results_v<版本>/run<N>` 目录文件契约 |
| [docs/architecture.md](docs/architecture.md) | AgentTorch 风格框架的组件对应关系与每代调度顺序 |

默认实验已修正有性生殖的密度、每雌繁殖率与植物耐受假设，详见 [修正说明](README_v0.1.1.md)。当前配置以 run4 为准；run1–run3 已删除，历史 v0.1.0 结果不代表当前配置。

修正后的长程批次：`results/results_v0.1.1/run4`，两个处理各 3 个重复，均完整运行 200 代。昆虫有晚期灭绝，也有存续至终点的重复；完整结果见 本地修正报告（`codex/mating_correction_report.md`，不随 Git 发布）。

代码格式、目录职责与检查命令见 [开发约定](CONTRIBUTING.md)。
