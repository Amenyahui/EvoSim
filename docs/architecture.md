# AgentTorch 风格对应关系

参考本地 AgentTorch 的 core/runner.py、controller.py、substep.py、registry.py 和 models/predator_prey，采用相同职责划分和 forward 调用约定。

| AgentTorch 风格 | evoSim 实现 |
|---|---|
| Runner.init / step / reset | 通用 Runner；PlantInsectRunner 处理终止代与模型状态 |
| Controller.observe / act / progress | 根据 YAML 的 active_agents、generator 和路径调度 |
| Registry.register | 每个运行器独立注册初始化、观察、策略和转换 |
| Observation.forward(state) | ObserveCues：本代快照、植物 VOC 线索 |
| Action.forward(state, observation) | ChooseHosts：选择寄主；SampleOffspring：存活与后代抽样 |
| Transition.forward(state, action) | SettleFeeding、CollectOffspring、InheritAndDisperse |
| state / agents / environment | 对齐的植物和昆虫个体数组，environment/runtime 保存运行上下文 |
| 模型局部配置、helpers、substeps | models/plant_insect 下集中组织 |

YAML 的 input_variables 映射状态路径，output_variables 声明输出；Controller 检查并提交转换返回的新状态。parameters 存模型参数，simulation_metadata.num_steps_per_episode 存代数，batch 定义处理与重复。实验通过 extends 继承公共流程。

## 每代顺序

1. ObserveCues 保存当前个体状态、计算 V_emit / V_cue；ChooseHosts 按线索和空间候选池选寄主。
2. SettleFeeding 结算共享资源、适合度和植物密度，记录个体数据与网络。
3. SampleOffspring 对植物、昆虫分别先抽存活，再抽条件 Poisson 后代；CollectOffspring 汇集结果。
4. InheritAndDisperse 检查分配上限，完成遗传、突变与扩散后返回下一代 agents。

初态为 t=0，完成 G 次更新后记录 t=G 的生态状态和条件矩，不再抽后代。step(n) 可以分段运行；reset() 使用相同种子重新初始化。随机流继续按阶段、类群、物种划分。快照复制数组，避免受后续状态修改影响。

记录在阶段内累积；下一代 agents 仅在全部生成成功后替换。灭绝、内存上限与算术异常分别记录终止原因。

保留原 NumPy 数值内核与独立审计方程，重构调度和模块边界。此实现不是 Torch 后端；学习参数与 GPU 请求会明确报错。迁移到 Torch 需要另行验证随机流和科学等价性。
