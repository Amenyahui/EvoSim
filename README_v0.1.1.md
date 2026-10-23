# evoSim 0.1.1：有性QTL运行说明

## 2026-10-22：默认实验修正

旧配置的 6 只昆虫在 10×10 网格中只有 1 格交配半径，同种子产生 5 雌 1 雄且无空间合格配偶；每雌 λ=2 在 F<1、雌性比例 0.5 下还低于理想替代水平。此为演示配置不匹配有性生命周期，不能据此宣称昆虫必然快速灭绝。

默认演示场景现改为：植物 P1/P2 各 200 株、昆虫 40 只；昆虫交配半径 2 格；每雌繁殖系数 λ=4；植物完全受食后的基础存活率 sP0=0.5（原 0.1）。单虫需求和单株资源仍均为 1。取食、死亡和交配公式不变，不补种、不强制雌雄比例、不保证存续。sP0 调整是一项显式的植物耐受假设变化，不能当作程序错误修复或实测校准。

依据：近似平衡性别下，E[N(t+1)]/N(t)≈0.5×S×q_mate×λ×F。λ=4 在高存活、配偶充分时允许增长；实际密度仍由资源与空间反馈约束。增加植物初始数量降低初始取食压力；提高 sP0 表示受食后仍有繁殖机会的耐受型理论植物，避免原场景的强烈寄主崩溃反馈。参数选择使用了试运行，因此后续曲线是探索性诊断，不能视为独立验证或正式处理效应检验。

保留全部探索记录：`codex/validation/mating_correction_*pilot.json`。包括只提高昆虫密度、提高植物密度、改变需求尺度仍发生灭绝/暴发的结果，没有删除不利种子。最终场景的三个种子 20261021–20261023 在 50 代时昆虫数为 108/134/127。

200 代配置现为两个处理各 3 个重复。基因型累计存储上限由 1,000,000 提高到 4,000,000 个标量（按 int64 约 32 MB 的纯基因型数据，不含日志与对象）；个体上限 5,000 和 allele 上限 100,000 不变。仍保留安全终止，不截断出生数。完整原始记录可能占用更多内存。

绘图逐重复展示，单个观测使用点标记，时间轴从 0 开始，标注灭绝和实际停止状态；遗传图阴影仍表示种群内 SD，不能当作重复间置信区间。等位频率图每物种/位点最多显示峰值频率最高的 12 个 allele，全部原始数据保留。


唯一运行协议是 `sexual-qtl-1`。旧配置、单亲克隆、直接性状突变及旧版本分支已移除；历史 results 只读保留。参数仍为理论占位值，不用于假说检验。

## 运行

```powershell
python example.py
python -m tests
python -m evo_sim --config evo_sim/models/plant_insect/configs/sexual_qtl_development.yaml
python -m evo_sim.models.plant_insect.audit --input results/results_v0.1.1/run4 --allow-plot-source-drift
python -m evo_sim.models.plant_insect.analysis.summarize_qtl_run --input results/results_v0.1.1/run4
```

输出目录统一命名为 `results/results_v<版本>/run<N>`，省略 --out 时自动取下一个 run 编号；--out 可指定新目录。当前只保留最终批次 `run4`（200 代，两处理 × 3 重复）；run1–run3 已删除，其状态记录见 本地修正报告（`codex/mating_correction_report.md`，不随 Git 发布），下一次默认创建 run5。已有输出目录不可覆盖。smoke/development/200gen分别请求1/5/200次更新；提前灭绝和计算保护须报告真实状态。

## 本版机制

1. 保留原VOC、寄主选择、比例共享资源、S/F、植物局部竞争与λ/α。
2. 植物为严格异交两性花；昆虫0=雌、1=雄，性别只在出生时赋值，默认雌性概率0.5。所有个体先抽一次存活；同种、距离合格、F>0的存活供体形成候选池；仅有配偶的种子母本/雌虫抽Poisson出生数。
3. 每个实际子代均匀有放回选择父本。父母分别按连锁群形成配子，组内相邻间隔K~Poisson(ΔcM/100)，奇数切换同源拷贝，组间独立分离。默认每类群8个独立抽象连锁群。
4. 配子逐QTL突变产生新不可变allele，copy0母源/copy1父源。D/Z/I/A由子代基因型半和加性计算并剪裁；不直接修改性状。
5. 植物从种子母本位置、昆虫从雌虫当代实际寄主位置扩散；整代替换，亲本全部退出。终代只观测，不抽繁殖随机数。

## 审计和边界

`offspring.csv`记录双亲；`gamete_transmissions.csv`逐位点记录来源、交换数和突变；`haplotype_summary.csv`仅在图谱含同组位点对时记录频率。安全检查后临时生成所有子代，失败不提交新allele、子代或ID，并恢复繁殖随机流；当代生态和尝试计数可保留作失败证据。所有结果累积在内存，不承诺大规模全基因组吞吐。

Haldane只是在无交换干涉条件下的理论模型。SLiM源码用于核对双亲边界和前向模拟设计，本项目不依赖、不复制SLiM引擎，也不将cM当碱基。单性、无配偶、零食物可自然灭绝，无补种或补雄。

详见 本地开发规格（`codex/evoSim_v0.1.1_Agent_README.md`）、[模型公式](docs/model.md)、[输出规范](docs/results.md)、本地实施报告（`codex/mating_correction_report.md`）。

修正后的长程批次：`results/results_v0.1.1/run4`，两个处理各 3 个重复，均完整运行 200 代。昆虫有晚期灭绝，也有存续至终点的重复；完整结果见 本地修正报告（`codex/mating_correction_report.md`，不随 Git 发布）。

## A–E 综合诊断图

现有汇总表生成后，使用统一入口：

```powershell
python -m evo_sim.models.plant_insect.analysis.plot_diagnostics --input results/results_v0.1.1/run4 --out results/results_v0.1.1/run4/figures
```

借鉴用户提供的 revised A–E 脚本，修正 I/A 生物学标签、空间行列方向与写死的范围/重复数。要求完整且等长的两个处理，拒绝将计算停止的数据混入均值；D 默认展示 mutation_off 的 rep01/rep02，可用 --spatial-reps 显式指定。旧绘图已删除；新入口自动更新汇总。
