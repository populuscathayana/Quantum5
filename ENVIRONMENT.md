# 运行环境

电子结构：CP2K 2026.2原生GFN1-xTB，中性单重态，XYZ周期；晶胞和全部坐标见每项 `calc.inp`。裸表面固定底32Cu；吸附、路径及频率固定64Cu。路径为5副本CI-NEB，频率中心差分步长0.01 bohr，仅O的3个活动坐标。

离线复核建议Python3.11或3.12，依赖见 `requirements.txt`。绘图使用Matplotlib Agg、DejaVu Sans；报告使用XeLaTeX、ctex、Fandol。结果图像可因库版本不同出现少量渲染差别，数值核验有显式容差。

`python reproduce.py --tex` 使用PATH中的 `xelatex`；编译不执行任何量子计算。无需集群账号、SSH或原作者本机路径。

若独立运行CP2K，须自行配置CP2K及计算资源，在各输入目录执行。仓库没有自动提交、重试或远程连接程序；特别是 `failed-initial` 与 `symmetry-baseline` 仅为历史对照，不是推荐表面。
