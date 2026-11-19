# genetic-relatedness-in-PRS

> ⚠️ 本仓库代码位于 `reference_code/`，是 2020 年前后读博期间的历史实现，仅作参考与溯源，不代表当前可用状态；新的开发请另建目录（如 `src/`），不要直接改这里的文件。

## 这个项目研究什么

这个项目研究当 PRS 的目标样本与发现样本存在亲缘关系或样本重叠时，预测 R² 为什么会被高估。工作重点是区分样本重叠与遗传亲缘这两种来源。

## 方法/数据

合作者提供的 R 脚本使用亲缘关系矩阵、Wishart 矩计算和 UKBB 模拟，比较不同相关结构下的 PRS 过拟合程度。

## 当前状态

已有理论推导和模拟脚本；运行仍依赖合作者的启动环境与 UKB 基因型数据，本仓库不包含这些数据，尚未形成可直接复跑的流程。

## 代码地图

- `reference_code/tim_mak_scripts/code/2019-01/`：合作者提供的早期 UKBB 亲缘关系与模拟脚本。
- `reference_code/tim_mak_scripts/code/relatedness-simulations/`：样本相关与重叠的模拟及 PBS 批处理脚本。
- `reference_code/tim_mak_scripts/moments/`：矩与 Wishart 相关计算脚本。
- `reference_code/tim_mak_scripts/paper/20201218/`：对应版本的 moments 与 simulations 计算脚本。
