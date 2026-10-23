# 开发约定

Python ≥ 3.10。先在项目根目录安装开发依赖：

```powershell
python -m pip install -e ".[dev]"
```

## 代码与目录

- 使用 UTF-8、4 空格缩进、100 字符行宽。函数和变量采用 snake_case，类采用 PascalCase，常量采用 UPPER_SNAKE_CASE。
- 格式及导入顺序由 Ruff 统一管理，以 pyproject.toml 为准。
- `evo_sim/core/` 管理通用调度；`models/plant_insect/helpers/` 实现模型计算；`substeps/` 组织每代流程。
- `analysis/summarize_*.py` 和 `analyze_macro.py` 只输出数据汇总；`plot_diagnostics.py` 是唯一绘图入口，负责 A–E 图。
- 使用已有辅助函数，避免重复实现。注释说明科学假设、单位和限制；公开入口用简短 docstring 说明输入或用途。
- 改动科学规则或默认参数时，同步更新 docs/model.md、docs/assumptions.md 和相关验收测试。
- 空的 __init__.py 是包标记，可以保留。

## 修改后检查

```powershell
python -m ruff check evo_sim tests example.py
python -m ruff format --check evo_sim tests example.py
python -m tests
```

需要自动整理时运行 `python -m ruff check evo_sim tests example.py --fix` 和 `python -m ruff format evo_sim tests example.py`，随后检查差异。测试采用普通 assert 函数，不使用 `python -O`。

数据汇总或绘图调整应使用既有结果验证；纯整理不得改变模拟公式、随机抽样顺序或已有原始结果。不要修改历史元数据中的源码哈希来掩盖差异。

## 本地产物

`results/`、`tests/reports/` 和 `codex/` 已被 Git 忽略。codex 中的开发建议、验证日志和临时报告仅保存在本地；对使用者必要的说明应写入 README 或 docs，不能只链接到被忽略的文件。

汇总命令仅接受 `--input`，统一写入批次的 `data/summary/`。`summarize_qtl_run` 一次分块读取同时生成性状和遗传多样性表；原独立 `summarize_genetics` 入口已合并。只有 `plot_diagnostics` 接受 `--out` 并创建图目录。内部汇总函数共用 `Batch` 对象。
