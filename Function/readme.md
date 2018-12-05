```matlab
function [x,O_min] = min_convex_newton(O_fun, dd_f, d_f, A_eq,b_eq, O_neq,dd_neq,d_neq, x0)
```
主要针对牛顿下降法，其中各个变量含义如下：
- O_fun: 目标函数
- dd_f: 目标函数的二阶导数
- d_f: 目标函数的一阶导数
- A_eq: 等式约束的A
- b_eq: 等式约束的b
- O_neq: 不等式约束写成罚函数的形式
- dd_neq: 罚函数对应的二阶导数
- d_neq: 罚函数对应的一阶导数
- x0: 初值(满足不等式约束和等式约束)
- x, O_min: 最优结果
