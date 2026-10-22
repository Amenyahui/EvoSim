# evoSim

植物—植食性昆虫个体模拟。按照本地 AgentTorch 的代码组织和调用方式重构：Registry 注册组件，Runner 管理运行，Controller 根据 YAML 执行 Observation → Policy → Transition，各组件实现 `forward()`。

数值计算采用 NumPy，保持原模型方程、随机流和数据语义。这是独立实现的运行框架，参考 `E:/myh/03PhD/ABM/AgentTorch/agent_torch`；未接入 AgentTorch/PyTorch 后端，不提供 GPU、自动求导或参数校准。

## 文件结构

```text
evoSim/
├── evo_sim/
│   ├── core/                       Runner、Controller、Registry、配置和状态路径
│   └── models/plant_insect/
│       ├── config.yaml             子步骤、状态、输入输出路径
│       ├── configs/                实验参数和处理组合（YAML）
│       ├── simulator.py            模型注册和运行生命周期
│       ├── substeps/               forward()：生态、繁殖、遗传与扩散
│       ├── helpers/                方程、个体数组、遗传规则和参数检查
│       ├── analysis/               汇总、绘图和初始条件分析
│       ├── batch.py                批量实验与数据保存
│       └── audit*.py               独立重建账本、检查完整性
├── docs/                           方程、数据字段和结构说明
├── tests/                          自动检查与小型测试输入
│   └── reports/                    检查报告
├── results/<实验名>/
│   ├── settings.json               实际参数和源文件哈希
│   ├── data/                       模拟数据、运行状态和衍生汇总
│   └── figures/                    PNG / PDF 结果图
├── example.py                      Runner 调用示例
├── pyproject.toml                  安装信息和命令入口
└── requirements.txt               依赖
```

`results` 不存放代码、源码副本、验证报告或历史备份。已有结果目录不会被覆盖。

## 运行

在本目录打开终端：

```powershell
python -m pip install -e .
python example.py
python -m evo_sim --config evo_sim/models/plant_insect/configs/qtl_200gen.yaml --out results/qtl_200gen
python -m evo_sim.models.plant_insect.audit --input results/qtl_200gen
python -m evo_sim.models.plant_insect.analysis.plot_genetics --input results/qtl_200gen
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input results/qtl_200gen
```

完整检查：`python -m tests`。可选的新旧代码对照：

```powershell
python -m tests.compare_reference --reference ../simu_v1.4
```

新项目运行不依赖原工程；只有可选迁移对照会读取原工程。

## 修改入口

- 实验参数：`evo_sim/models/plant_insect/configs/*.yaml`。
- 阶段顺序和状态路径：`evo_sim/models/plant_insect/config.yaml`。
- 个体行为：`substeps` 中对应类的 `forward()`；通用数值计算在 `helpers`。
- 新增行为：实现 Substep 类，在 `simulator.get_registry()` 注册，在 YAML 引用 generator 名称。

详见 [结构说明](docs/architecture.md)、[科学方程](docs/model.md) 和 [数据字段](docs/data_dictionary.md)。
