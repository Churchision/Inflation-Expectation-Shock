# Inflation-Expectation-Shock
## Fasani2023.mod
原来的代码添加中文注释，添加了不确定性变量
## expectationshock.mod
添加并只设置了expectation shock
## revise_expectationshock.mod
校准了参数
## only_revise_expectationshock.mod
删除所有冲击，只保留通胀预期冲击和技术冲击（因为发现别的zeps也会动）
在IRF_mat增加了对不确定性的稳态偏离定义
## China_only_revise_expectationshock.mod
校准了中国的参数。
* xie仍为6.51（校准后解不出稳态，其有可能能被稳态反解）
* ec，xc未校准，缺乏相应数据
* Taylor规则中phiy变化（0.01变为0.5）可能会导致剧烈变化
* 其他参数没发现如此重要的作用
结果不收敛
