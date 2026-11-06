# 开发与验证

项目规范见根目录 `AGENTS.md`，当前模块地图和分批迁移计划见
[architecture.md](architecture.md)。结构重构保持 v0.1.11 科学协议，不能同时修改方程、
默认参数、随机数消耗、生活史顺序或输出含义。

## 环境与命令

Python >= 3.10；Python 依赖和开发依赖由 `pyproject.toml` 定义：

```powershell
python -m pip install -e ".[dev]"
python -m ruff check evo_sim tests example.py
python -m ruff format --check evo_sim tests example.py
python -m tests
```

完整测试使用现有无框架断言接口，必须启用 Python 断言，禁止使用 `python -O`。
测试入口首先检查真实 SLiM 5.2。引擎使用 `EVOSIM_SLIM`、配置中的 executable 或 PATH；
开发工作区也可使用 `codex/slim_runtime/Library/bin/slim.exe`，不要作为历史文件删除。
安装说明见 [README](../README.md)。

## 重构验证与记录

1. 先阅读实现及全部调用点、相关测试，记录固定种子原始结果和随机流状态。
2. 一次只迁移一个职责；不改算式、数组顺序、抽样函数或检查点/提交边界。
3. 比较前后快照、事件、遗传传递、资源记录、阶段摘要及最终状态，不能只比较最终人口。
4. 独立审计继续根据配置和原始账本重建，不能调用生产生物方程核验自身。
5. 所有相关检查通过后审阅完整 diff，再提交；报告实际执行的测试和限制。

验证目录使用 `codex/validation/v0.1.11/run<N>/`，每次新建目录，不能覆盖旧运行。
正式研究输出使用 `results/`，开发验证不写入该目录。源码 SHA256 因迁移改变属于预期，
但必须单独记录，不能通过改写旧 metadata 或旧审计结果掩盖来源差异。

## 历史清理

仅清理已确认无运行依赖、无独有科研信息的历史副本和可再生成缓存。
保留 Git 提交/tag、当前输入/配置和研究结果。忽略目录不等于远端已备份；
删除前必须检查具体内容。不要按测试文件中的旧版本号判断其是否过时。
