> [← 返回总目录](MCM_A_备赛学习指南.md)　·　本册：第四章 P1 重点专项

## 第四章 P1 重点专项

P1 的要求不是"能推完所有定理"，而是：**知道它解决什么问题、能写出标准形式、能在两个候选模型之间说出选择理由**。每节给最小例子与真题落点。

### 4.1 最优控制（PMP 与变分）

**定位**：P1｜难度 ★★★★｜6 h｜前置：3.2、3.9、变分法｜2022 的核心工具。

#### 1. 定义

一句话：**系统怎么演化已经由 ODE 定死，剩下要选的是"每一刻推多大劲"这条曲线。**

$$
\min_{u(\cdot)\in\mathcal U}\ J=\phi(\mathbf{x}(T))+\int_0^T L(\mathbf{x},\mathbf{u},t)\,dt
\quad\text{s.t.}\quad \dot{\mathbf{x}}=\mathbf{f}(\mathbf{x},\mathbf{u},t),\ \ \mathbf{x}(0)=\mathbf{x}_0,\ \ \mathbf{x}(T)\in\mathcal X_T .
$$

已知：向量场 $\mathbf{f}$、运行代价 $L$、终端代价 $\phi$、初值 $\mathbf{x}_0$、时域 $T$（可固定也可自由）、控制集 $\mathcal U$。未知：整条控制曲线 $\mathbf{u}(\cdot)$、由它决定的轨线 $\mathbf{x}(\cdot)$，以及作为副产品的协态 $\boldsymbol\lambda(\cdot)$。注意未知量是**函数**而不是有限维向量，所以这是无穷维优化（变分问题），不是 3.9 里的 $\min f(\mathbf{x})$。

#### 2. 核心思想与大白话

笨办法是"把 $u(t)$ 打成 $N$ 段，丢给 `fmincon`"。它能算，但有三个毛病：$N$ 大时变量多、梯度贵；解出来是一串数，看不出结构（该不该开满、何时切换）；更没法说明最优性。

PMP（Pontryagin 最小值原理，Pontryagin's minimum principle）换了个问法：不问"每时刻该给多少"，而问"**多一单位状态值多少钱**"。这个问题交给一个随时间倒着演化的价格 $\boldsymbol\lambda(t)$，再让每一刻都买"最便宜的即时代价加按价折算的状态变化"。于是一个无穷维问题变成**两点边值问题**（boundary value problem, BVP）：$2n$ 个一阶 ODE，$n$ 个初值条件、$n$ 个终端条件。

> **大白话。** 你开车从 A 到 B，油钱按油门算，迟到要罚款。你不需要事先规划整段油门曲线，只需要一张"时间价格表"：晚到一分钟罚多少、早到一分钟省多少。这张表随时间倒着推（越靠近终点，剩余时间的价值越清楚），每一刻你按当前价格挑最划算的油门。$\boldsymbol\lambda(t)$ 就是这张表——**此刻状态多 1 单位，未来总代价能降多少**。价格表的单位是"代价/状态"，所以 $\boldsymbol\lambda^\top\mathbf{f}$ 才和 $L$ 一样是"代价/时间"。

这个类比与数学一致：$\boldsymbol\lambda(t)=\partial V/\partial\mathbf{x}$（值函数的梯度），见第 6 小节第十二步。

#### 3. 数学基础

- **泛函与一阶变分（first variation）**：$J$ 的定义域是函数空间，导数取 Gateaux 方向导数 $\delta J=\frac{d}{d\varepsilon}J[\mathbf{u}+\varepsilon\boldsymbol\eta]\big|_{\varepsilon=0}$，最优的必要条件是对所有容许 $\boldsymbol\eta$ 有 $\delta J=0$。用来把"找最优曲线"化成逐点方程。
- **分部积分 + 变分法基本引理**：$\int(\cdots)\boldsymbol\eta\,dt=0$ 对任意 $\boldsymbol\eta$ 成立，推出被积函数恒为零。这是从积分条件得到微分方程的唯一通道。
- **Lagrange 乘子（函数版）**：ODE 在每一时刻都是一个等式约束，所以乘子不是向量而是函数 $\boldsymbol\lambda(t)$。用来把带约束变分变成无约束变分。
- **Hamiltonian 与 Legendre 变换**：$H=L+\boldsymbol\lambda^\top\mathbf{f}$ 把二阶的 Euler–Lagrange 方程降成一阶 Hamilton 系统，既便于数值求解，也暴露结构（线性于 $u$ 与否）。
- **两点边值问题**：$n$ 个条件在 $t=0$、$n$ 个条件在 $t=T$，中间用打靶（shooting）或配点（collocation）连接。
- **二阶条件**：$H_{uu}\succ0$（Legendre–Clebsch 条件）是局部最小的必要条件；$H$ 对 $(\mathbf{x},\mathbf{u})$ 联合凸且 $\mathbf{f}$ 线性时，一阶条件升格为充分条件。
- **值函数与 HJB 方程**：$V(\mathbf{x},t)$ 满足一阶 PDE，$\boldsymbol\lambda=\nabla_{\mathbf{x}}V$。用来给 PMP 一个全局解释，也用来判断数值解是否合理。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| $\mathbf{f},L,\phi$ 对 $(\mathbf{x},\mathbf{u})$ 连续可微 | 一阶变分存在 | 机理模型光滑 | 开关、干摩擦、饱和处 $\boldsymbol\lambda$ 跳变，需分段拼接；漏掉跳变使终端误差达 10%—30% | 沿最优轨线检查 $\mathbf{f}$ 的 Lipschitz 常数与分段点 |
| 控制集 $\mathcal U$ 紧（闭且有界） | 逐点最小值存在 | 执行器有物理上限 | 无界时"最优"可能是脉冲（Dirac），数值解发散 | 画可行域；看解是否贴边 |
| 无状态约束，或约束不起作用 | 只需 ODE 乘子 | 状态自然满足安全边界 | 状态越界却未加约束，$J$ 被低估且方案不可行 | 解后逐时刻检查 $\mathbf{x}(t)\le\mathbf{x}_{\max}$ |
| 时域与终端条件已给定 | 横截条件封闭 BVP | 题目给了截止时间 | 自由 $T$ 却漏掉 $H(T)=-\phi_T$，时间不最优，成绩差 5%—20% | 检查 $T$ 是题目给定还是由最优性决定 |
| 控制无延迟、无磨损 | $\mathbf{u}(t)$ 即时生效 | 理想执行器 | 实际有 1—2 s 延迟时策略不可实现，峰值超限 | 对比指令与实际响应的互相关 |
| 参数确定 | 优化用名义值 | 参数已估计 | 参数偏 10% 可使最优策略的优势从 5% 掉到 1% | 做参数敏感性（3.10）与情景重解 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $t$ | 时间 | s | $0$—$10^4$ | 自变量；$T$ 越大，终端代价相对越轻 |
| $\mathbf{x}$ | 状态向量 | 各分量不同（m、J、K） | 2—10 维 | 系统"存量"；越大表示可用资源或惯性越多 |
| $\mathbf{u}$ | 控制向量 | W、N、K/s、— | 1—3 维 | 每时刻的"用力"；受执行器上限约束 |
| $\boldsymbol\lambda$ | 协态（costate） | 代价/状态单位 | 与 $\partial\phi/\partial\mathbf{x}$ 同量级 | 状态的边际价格；绝对值大表示该状态很值钱 |
| $H$ | Hamiltonian | 代价/时间 | 与 $L$ 同量级 | 即时代价加状态的"贴现"价值 |
| $L$ | 运行代价率 | 代价/时间 | $u^2$、$1$、$c(u)$ | 每一秒的花钱速度 |
| $\phi$ | 终端代价 | 代价 | $(x(T)-x_d)^2$ | 只惩罚终点偏差 |
| $T$ | 时域 | s | 题目给定 | 决策视野长度 |
| $\mathcal U$ | 控制可行集 | — | $[u_{\min},u_{\max}]$ | 执行器的物理边界 |

代价若用能量，1 Wh = 3600 J；论文里必须写明 $L$ 的单位，否则 $\boldsymbol\lambda$ 的数值无法解释。

#### 6. 从基本原理推导核心公式

**第一步，泛函与扰动。** 固定 $\mathbf{x}_0,T$，每个容许控制 $\mathbf{u}(\cdot)$ 经 ODE 决定一条轨线，从而给出实数 $J[\mathbf{u}]$。设 $\mathbf{u}^*$ 最优，取任意扰动 $\boldsymbol\eta$，令 $\mathbf{u}_\varepsilon=\mathbf{u}^*+\varepsilon\boldsymbol\eta$，则 $\frac{d}{d\varepsilon}J[\mathbf{u}_\varepsilon]\big|_{\varepsilon=0}=0$。

**第二步，变分的动力学。** 令 $\delta\mathbf{x}=\frac{d}{d\varepsilon}\mathbf{x}[\mathbf{u}_\varepsilon]\big|_{\varepsilon=0}$，对 ODE 求导得线性化方程
$$\delta\dot{\mathbf{x}}=\mathbf{f}_{\mathbf{x}}\,\delta\mathbf{x}+\mathbf{f}_{\mathbf{u}}\,\boldsymbol\eta,\qquad \delta\mathbf{x}(0)=\mathbf{0}.$$

**第三步，乘子与分部积分。** 把约束乘进目标，乘子 $\boldsymbol\lambda$ 取 $\boldsymbol\lambda^\top(\mathbf{f}-\dot{\mathbf{x}})$ 这一符号（这样 $\boldsymbol\lambda$ 与 $H=L+\boldsymbol\lambda^\top\mathbf{f}$ 同号）：
$$\delta J=\phi_{\mathbf{x}}\delta\mathbf{x}(T)+\int_0^T\Big[L_{\mathbf{x}}\delta\mathbf{x}+L_{\mathbf{u}}\boldsymbol\eta+\boldsymbol\lambda^\top\big(\mathbf{f}_{\mathbf{x}}\delta\mathbf{x}+\mathbf{f}_{\mathbf{u}}\boldsymbol\eta-\delta\dot{\mathbf{x}}\big)\Big]dt .$$
对 $-\int\boldsymbol\lambda^\top\delta\dot{\mathbf{x}}\,dt$ 分部积分得 $-\big[\boldsymbol\lambda^\top\delta\mathbf{x}\big]_0^T+\int\dot{\boldsymbol\lambda}^\top\delta\mathbf{x}\,dt$（用 $\delta\mathbf{x}(0)=\mathbf{0}$）。合并同类项：
$$\delta J=\big[\phi_{\mathbf{x}}-\boldsymbol\lambda^\top\big]\delta\mathbf{x}(T)+\int_0^T\Big[\big(L_{\mathbf{x}}+\dot{\boldsymbol\lambda}^\top+\boldsymbol\lambda^\top\mathbf{f}_{\mathbf{x}}\big)\delta\mathbf{x}+\big(L_{\mathbf{u}}+\boldsymbol\lambda^\top\mathbf{f}_{\mathbf{u}}\big)\boldsymbol\eta\Big]dt .$$

**第四步，选 $\boldsymbol\lambda$ 消去状态变分。** 令
$$\dot{\boldsymbol\lambda}=-\big(L_{\mathbf{x}}+\mathbf{f}_{\mathbf{x}}^\top\boldsymbol\lambda\big)=:-H_{\mathbf{x}},\qquad \boldsymbol\lambda(T)=\phi_{\mathbf{x}}^\top .$$
这两条正是协态方程与横截条件（transversality condition）。含 $\delta\mathbf{x}$ 的项全消，只剩
$$\int_0^T\big(L_{\mathbf{u}}+\boldsymbol\lambda^\top\mathbf{f}_{\mathbf{u}}\big)\boldsymbol\eta\,dt=0\quad\text{对任意 }\boldsymbol\eta,$$
由变分法基本引理，被积函数恒零：$L_{\mathbf{u}}+\boldsymbol\lambda^\top\mathbf{f}_{\mathbf{u}}=H_{\mathbf{u}}=\mathbf{0}$。

**第五步，Euler–Lagrange 是它的特例。** 若 $\mathbf{f}=\mathbf{u}$（即 $\dot{\mathbf{x}}=\mathbf{u}$）、$L=L(\mathbf{x},\dot{\mathbf{x}},t)$，则 $H=L+\boldsymbol\lambda^\top\dot{\mathbf{x}}$，$H_{\mathbf{u}}=L_{\dot{\mathbf{x}}}+\boldsymbol\lambda=\mathbf{0}$ 给出 $\boldsymbol\lambda=-L_{\dot{\mathbf{x}}}$；代入协态方程得 $\frac{d}{dt}L_{\dot{\mathbf{x}}}-L_{\mathbf{x}}=0$，就是经典 Euler–Lagrange 方程。所以 **PMP 是带微分约束的 Euler–Lagrange，$\boldsymbol\lambda$ 是广义动量**。

**第六步，有界控制转成变分不等式。** 当 $\mathbf{u}^*$ 贴在边界上，$\varepsilon\boldsymbol\eta$ 只有单侧可行，$\delta J\ge0$ 而非 $=0$，于是
$$\big(L_{\mathbf{u}}+\boldsymbol\lambda^\top\mathbf{f}_{\mathbf{u}}\big)(\mathbf{u}-\mathbf{u}^*)\ge0\ \ \forall\mathbf{u}\in\mathcal U\ \Longleftrightarrow\ \mathbf{u}^*=\arg\min_{\mathbf{u}\in\mathcal U}H .$$
这就是"最小值原理"这个名字的来历。英文文献常写成 $\max$，因为有人把 $\boldsymbol\lambda$ 定义成收益的符号——**符号约定必须在论文里声明**，否则 $\lambda=-1$ 会被读反。

**第七步，为什么 $\partial H/\partial\mathbf{u}=0$ 给出 $\mathbf{u}^*$。** $H$ 中 $\mathbf{u}$ 只通过 $L$ 与 $\mathbf{f}$ 出现，而 $\boldsymbol\lambda$ 是乘子、不是 $\mathbf{u}$ 的函数，对 $\mathbf{u}$ 求偏导时当常数。内点解满足 $H_{\mathbf{u}}=\mathbf{0}$；若 $H_{\mathbf{uu}}\succ0$（例如 LQR 中 $H=\mathbf{u}^\top R\mathbf{u}+\cdots$，则 $H_{\mathbf{uu}}=2R\succ0$），这是严格全局最小，解出 $\mathbf{u}^*=-\frac{1}{2}R^{-1}B^\top\boldsymbol\lambda$；若把运行代价写成 $\frac12\mathbf{u}^\top R\mathbf{u}$，系数 $2$ 被吸收，就得到常见的 $\mathbf{u}^*=-R^{-1}B^\top\boldsymbol\lambda$。两式只是 $R$ 的定义差一个因子，混用会差 2 倍的增益。

**第八步，终端约束的乘子。** 若终点受 $\boldsymbol\psi(\mathbf{x}(T))=\mathbf{0}$ 约束，$\delta\mathbf{x}(T)$ 不再任意，需再加乘子 $\boldsymbol\nu$：
$$\boldsymbol\lambda(T)=\phi_{\mathbf{x}}^\top+\boldsymbol\psi_{\mathbf{x}}^\top\boldsymbol\nu .$$
$\boldsymbol\nu$ 与 $\boldsymbol\lambda$ 同量纲，是"终端约束的影子价格"。自由终端时 $\boldsymbol\nu=\mathbf{0}$，回到 $\boldsymbol\lambda(T)=\phi_{\mathbf{x}}^\top$；固定终端 $\mathbf{x}(T)=\mathbf{x}_T$ 时 $\delta\mathbf{x}(T)=\mathbf{0}$，$\boldsymbol\lambda(T)$ 不被约束。

**第九步，自由终端时间。** 若 $T$ 也可选，对 $T$ 变分多出一项，条件为 $H(T)=-\phi_T$（$\phi$ 不显含 $T$ 时 $H(T)=0$）。自治系统（$\mathbf{f},L$ 不显含 $t$）沿最优轨线 $H\equiv$ 常数，这是可用来降阶的首次积分。

**第十步，bang-bang 与奇异弧。** 若 $H$ 对 $\mathbf{u}$ 线性，$H=H_0(\mathbf{x},\boldsymbol\lambda,t)+\boldsymbol\sigma^\top(t)\mathbf{u}$，其中 $\boldsymbol\sigma=H_{\mathbf{u}}$ 称**切换函数**（switching function）。逐点最小化给出
$$u_i^*=u_{i,\min}\ (\sigma_i>0),\qquad u_i^*=u_{i,\max}\ (\sigma_i<0),$$
即 **bang-bang 控制**；$\sigma_i$ 的零点就是切换时刻。若某段区间上 $\boldsymbol\sigma\equiv\mathbf{0}$，$H_{\mathbf{u}}=0$ 给不出 $\mathbf{u}$，该段称**奇异弧**（singular arc）：必须继续对 $t$ 求导 $\frac{d}{dt}H_{\mathbf{u}}=0,\ \frac{d^2}{dt^2}H_{\mathbf{u}}=0,\dots$ 直到 $\mathbf{u}$ 显式出现，再用 Kelley 条件（$\frac{d^{2k}}{dt^{2k}}H_{uu}$ 的符号）判断它是最小还是最大。

**第十一步，量纲与极限检查。** $H=L+\boldsymbol\lambda^\top\mathbf{f}$ 两项必须同量纲，故 $[\boldsymbol\lambda]=$ 代价/状态；算例中 $J,x,u$ 全无量纲，$\lambda=-1$ 无量纲，自洽。极限：控制权重 $R\to0$（用力不要钱）时解趋向 bang-bang 脉冲；$R\to\infty$ 时 $\mathbf{u}^*\to\mathbf{0}$、$J\to\phi(\mathbf{x}_0)$；$T\to0$ 时只剩终端代价，退化为静态投影。

**第十二步，与 HJB 的关系。** 值函数 $V(\mathbf{x},t)$ 是最优剩余代价，满足 $-V_t=\min_{\mathbf{u}}H(\mathbf{x},\mathbf{u},\nabla_{\mathbf{x}}V)$。沿最优轨线 $\boldsymbol\lambda(t)=\nabla_{\mathbf{x}}V(\mathbf{x}(t),t)$，于是 $\dot{\boldsymbol\lambda}=-H_{\mathbf{x}}$ 就是 HJB 对 $\mathbf{x}$ 求全导数（envelope theorem）。算例可验证：$V=a(t)(x-1)^2$、$a(t)=1/(2-t)$，则 $\lambda(0)=2a(0)(0-1)=-1$、$V(0,0)=0.5$，与第 12 小节一致。

#### 7. 求解方法

- **解析（LQR）**：线性系统加二次目标时，值函数是二次型，Riccati 方程 $\dot P=-(A^\top P+PA+Q-PBR^{-1}B^\top P)$、$P(T)=Q_f$，$u^*=-R^{-1}B^\top Px$。自治情形可令 $\dot\lambda=0$ 化成代数方程（ARE）。
- **间接法**：由 PMP 得 BVP，用 `bvp4c` 或 `bvp5c`；初值猜测难，标准做法是**延拓**（homotopy）：先解 $R$ 很大（近似无控制）的问题，再逐步减小到目标值，每一步用上一步的解当初值。
- **直接法**：单射击（只把 $\mathbf{u}$ 的节点当变量，$\mathbf{x}$ 由积分得到，无约束但目标非光滑）或配点（collocation，把 $\mathbf{x}$ 也当变量，动力学写成等式约束，规模大但稳定）。
- **多重射击**：把 $[0,T]$ 分段，段间连续性作约束，条件数从 $e^{LT}$ 降到 $e^{L\Delta t}$，是长时域的标准解法。
- **状态约束**：罚函数或松弛变量（简单但可能不可行），或间接法加跳变条件（精确但难实现）。
- **陷阱**：网格太粗会把 bang-bang 抹成斜坡，切换时刻偏 1—2 个网格；`fmincon` 的有限差分梯度在 $N>200$ 时极慢，应给解析梯度或稀疏结构；终端条件写错一个符号，解会漂到完全不同的轨线。

#### 8. 适用场景

**该用（判据）。**

- 控制是连续可调的量：功率、力、温度设定、施药速率；题目出现"每一时刻""连续调整"。
- 有可写的 ODE 机理（3.2），状态大致 2—10 维。
- 目标能写成"运行代价积分 + 终端惩罚"，且能说清代价单位。
- 需要结构性结论：该不该开满、何时切换、阈值在哪、$\boldsymbol\lambda$ 的经济含义。

**不该用。**

- 决策是离散档位/开关且状态少 → 用 MDP（4.2）。
- 状态维数超过 10 且强非线性不可微 → 直接法也吃力，考虑 MPC 或 RL。
- 只有输入输出数据、没有机理 → 先做系统辨识（3.7），别硬套 PMP。
- 要在线实时重规划 → 滚动时域 MPC 比一次性开环最优更合适。

#### 9. 局限

- **只是必要条件**：非凸时可能有 2—3 个局部极小，单次打靶只给一个；不加多初值就宣称全局最优是常见错误。
- **打靶条件数爆炸**：误差随 $T$ 指数放大（约 $e^{LT}$），$LT>10$ 时初值误差放大约 $10^4$ 倍，必须多重射击或配点。
- **状态约束要额外条件**：实现错会让解"穿墙"，且穿墙后目标反而更小，看起来"更优"。
- **开环最优对参数敏感**：参数偏 $\pm10\%$ 时，实际代价可比名义最优高 5%—15%。
- **奇异弧上二阶信息退化**：$H_{uu}=0$，数值解震荡，需正则化或解析延拓。
- **目标不可微**（$\max$、分段常数）时 $H_{\mathbf{u}}$ 不存在，一阶条件失效，只能走直接法加非光滑优化。

#### 10. 改进

1. **无量纲化 + 换时间尺度**（代价几乎为零）：$t\mapsto t/T$、状态除以典型量级，直接改善 BVP 条件数。
2. **多重射击 + 延拓**（代码约 2 倍，计算 5—20 倍）：解决长时域与初值敏感。
3. **直接配点 + 解析梯度/稀疏 Jacobian**（实现复杂）：工业界主流，绕开 BVP 不稳。
4. **iLQR / DDP**（需迭代，只在邻域最优）：高维非线性下收敛快，适合做精细解。
5. **鲁棒或机会约束优化**（规模乘情景数）：把参数不确定写进约束，解决"最优但脆弱"。

#### 11. 与其他模型的组合

- **+ 参数估计（3.7）**：先用数据拟合 $\mathbf{f}$ 中的 $\theta$，再在名义 $\theta$ 上求 $u^*$。接口是 $\theta$ 进入 $\mathbf{f}$，控制问题把它当常数；必须补一步"$\theta$ 误差 → $J^*$ 误差"的敏感性。
- **+ Sobol / 蒙特卡洛（3.10、3.11）**：对 $\theta$ 采样重解，得到 $J^*$ 的分布而不是一个数，这样"最优策略"才带误差条。
- **+ MDP（4.2）**：把时间离散成 $N$ 段、状态离散成网格，DP 给全局最优，可用来验证 PMP 解没有掉进局部极小。

#### 12. 手算算例

取 $J=\int_0^1u^2\,dt+(x(1)-1)^2$、$\dot x=u$、$x(0)=0$、$u$ 无界。

**第一步，写 $H$。** $H=u^2+\lambda u$。

**第二步，驻点条件。** $\frac{\partial H}{\partial u}=2u+\lambda=0\Rightarrow u^*=-\frac{\lambda}{2}$；$H_{uu}=2>0$，确是最小。

**第三步，协态方程。** $\dot\lambda=-\frac{\partial H}{\partial x}=0$，故 $\lambda$ 为常数，从而 $u^*$ 也是常数。

**第四步，终端条件。** $\lambda(1)=\frac{\partial\phi}{\partial x(1)}=2(x(1)-1)$，而 $x(1)=\int_0^1u\,dt=-\frac{\lambda}{2}$，代入：
$$\lambda=2\left(-\frac{\lambda}{2}-1\right)=-\lambda-2\ \Longrightarrow\ \lambda=-1,\quad u^*=0.5,\quad x^*(t)=0.5t .$$

**第五步，代回目标。** $J=\int_0^1(0.5)^2dt+(0.5-1)^2=0.25+0.25=0.5$。

**第六步，扰动检验（确认是极小）。** 令 $u=0.5+\varepsilon$，则 $x(1)=0.5+\varepsilon$，
$$J(\varepsilon)=(0.25+\varepsilon+\varepsilon^2)+(0.25-\varepsilon+\varepsilon^2)=0.5+2\varepsilon^2 ,$$
在 $\varepsilon=0$ 取最小，二阶系数 $2=H_{uu}\cdot T$ 与理论一致。

**第七步，值函数检验。** $V(x,t)=a(t)(x-1)^2$、$a(t)=1/(2-t)$ 满足 HJB，$V(0,0)=0.5=J^*$，$\lambda(0)=V_x(0,0)=-1$。这个例子适合当代码标尺：数值解应给出 $J\approx0.5$、$u\approx0.5$、$x(1)\approx0.5$。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% 最优控制直接法：min ∫_0^1 u^2 dt + (x(1)-1)^2,  x'=u,  x(0)=0
N  = 50; dt = 1/N;
u0 = 0.5*ones(N,1);
obj = @(u) sum(u.^2)*dt + (sum(u)*dt - 1)^2;   % x(1)=Σu·dt（矩形积分）
opts = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',1e4);
lb = -2*ones(N,1); ub = 2*ones(N,1);
[uopt, J] = fmincon(obj, u0, [], [], [], [], lb, ub, [], opts);
x1 = sum(uopt)*dt;
fprintf('J=%.4f (解析 0.5), mean u=%.4f (解析 0.5), x(1)=%.4f (解析 0.5)\n', ...
        J, mean(uopt), x1);
```

输入：段数 $N$、初值 $u_0$、控制界。输出：最优节点 `uopt`、代价 `J`、终端状态 `x1`。判据：三者同时接近解析值才算对（只对上 $J$ 可能是"少花能量换终端误差"的巧合）；把 $N$ 从 50 加到 500，$J$ 应单调下降收敛到 0.5，误差量级 $O(dt)$，这是矩形积分的正常表现。

#### 14. 真实美赛 A 题应用

- **2022（核心）**：赛道、风、疲劳已知，控制是功率 $P(t)$，状态含速度、能量、疲劳；最小化完赛时间等价于"$L=1$、终端时间自由"的最优控制。第 9.4 节给了 $\dot v$ 与 $\dot F$ 的耦合方程；PMP 在这里的用处是判断何时该超阈值、何时该恢复，并给出 $\lambda_F$（疲劳的边际时间成本）。
- **2026**：连续时间 SOC 模型加使用建议；若要给"最优使用策略"，可把负载/亮度当控制，最小化"给定时间内的耗电"或最大化续航。
- **2020**：把捕捞强度 $E(t)$ 当控制、生物量当状态、利润为目标，就是经典生物经济最优控制；$\lambda$ 是"多一吨鱼值多少钱"。

#### 15. 常见错误

1. **符号约定不声明**：把 $\lambda$ 当收益的影子价格却用 $\min H$，结论正好相反；论文里先写清"$\lambda$ 为成本的边际"。
2. **必要条件当充分**：非凸问题必须多初值或网格验证，否则报的是局部最优。
3. **忘记控制上下界**：解出 $u=3$ 而电机只能到 1.5，方案不可执行。
4. **横截条件写错**：终端自由却固定 $\lambda(T)$，或终端固定却加 $\lambda(T)=\phi_x$。
5. **单位混用**：$L$ 用 Wh、$t$ 用 s，$\lambda$ 的数值差 3600 倍，解释全错。
6. **网格太粗**：bang-bang 被抹成斜坡，切换时刻偏 1—2 个网格，成绩差几个百分点。
7. **把开环最优当反馈策略**：参数一变就失效，需 MPC 或鲁棒化。
8. **不查状态约束**：解穿越物理边界还在报"最优"，评委一眼看出不可行。

#### 16. 学习要求（验收）

- 能对给定的 $\mathbf{f},L,\phi$ 写出 $H$ 与 PMP 四式，并逐条说清它来自"变分、分部积分、逐点最小、终端变分"中的哪一步。
- 能独立推出算例的 $\lambda=-1$、$u^*=0.5$、$J=0.5$，并用 $H_{uu}>0$ 与扰动检验确认是最小而非最大。
- 能判断一道题该用 PMP、直接法还是 MDP，并说出判据（连续/离散、维数、是否有转移概率）。
- 能用直接法复现解析解，并说明误差随 $N$ 的收敛阶。

### 4.2 动态规划与 MDP

**定位**：P1｜难度 ★★★｜5 h｜前置：3.4、概率｜分阶段/随机决策。

#### 1. 定义

一句话：**把"从头算到尾"的一次大决策，拆成"每一步只看当前状态"的一串小决策。**

$$
V^*(s)=\max_{\pi}\ \mathbb E_{\pi}\Big[\sum_{t=0}^{\infty}\gamma^{t}r(s_t,a_t)\ \Big|\ s_0=s\Big],
$$

马尔可夫决策过程（Markov decision process, MDP）是一个五元组 $(\mathcal S,\mathcal A,P,r,\gamma)$：状态集、动作集、转移核 $P(s'\mid s,a)$、奖励 $r(s,a)$、折扣因子 $\gamma\in(0,1)$。已知 $P,r,\gamma$（或可从数据估计），未知的是最优值函数 $V^*$ 与最优策略 $\pi^*$。有限时域版本把 $\sum_{t=0}^{T-1}$ 求和、不加折扣，从 $t=T$ 倒推到 $0$。

注意未知量是一个**策略函数** $\pi:\mathcal S\to\mathcal A$，不是一条轨迹；这就是它比"搜轨迹"便宜的原因。

#### 2. 核心思想与大白话

笨办法有两个：枚举全部策略（$|\mathcal A|^{|\mathcal S|}$ 个，稍大就爆），或对每条轨迹做搜索（指数爆炸）。动态规划（dynamic programming, DP）利用一个结构性事实——**最优策略的后半段也必须是最优的**（Bellman 最优性原理，principle of optimality）——把指数复杂度降成多项式：每个状态只算一次，算完就复用。

代价是必须知道转移概率，而且要对所有状态都算一遍，包括永远到不了的状态。

> **大白话。** 从家去机场，你不必枚举所有路线；只要知道"从每个路口到机场最少还要多久"，站在当前路口就只需比较"下一条路的路况 + 下个路口的值"。折扣 $\gamma$ 是"未来的钱打折"：$\gamma=0.9$ 时，10 步之后的 1 块钱只当 $0.35$ 块（$0.9^{10}\approx0.349$）。倒着算，越靠近终点越清楚，一路推回起点，整张"值表"就出来了。

一致性检查：$|r|\le R_{\max}$ 时 $|V|\le R_{\max}\sum_t\gamma^t=R_{\max}/(1-\gamma)$ 有限——这正是要求 $\gamma<1$（或存在吸收态）的技术理由，不是建模偏好。

#### 3. 数学基础

- **马尔可夫性（Markov property）**：给定当前状态，未来与历史条件独立，轨迹概率是转移核的连乘。没有它，值函数不能只写成 $V(s)$。
- **全期望公式（law of total expectation）**：把无穷和拆成"第一步 + 剩余"，这是 Bellman 方程唯一的推导工具。
- **压缩映射与 Banach 不动点定理**：$\mathcal T$ 在 sup 范数下是 $\gamma$-压缩，故有唯一不动点且迭代几何收敛。值迭代能收敛全靠这一条。
- **策略评估的线性方程**：$V^\pi=(I-\gamma P^\pi)^{-1}r^\pi$；$|\mathcal S|$ 大时用迭代解而不用求逆（求逆是 $O(|\mathcal S|^3)$）。
- **策略改进定理**：对 $V^\pi$ 贪心得到的 $\pi'$ 满足 $V^{\pi'}\ge V^\pi$，保证策略迭代单调不降。
- **线性规划视角**：$V^*$ 是 $\min\sum_s V(s)$ s.t. $V(s)\ge r(s,a)+\gamma\sum_{s'}P(s'\mid s,a)V(s')$ 的解，说明"DP 与 LP 是同一个问题的两种算法"。
- **部分可观测（POMDP）**：状态看不见时，把信念 $b(s)$ 当状态，问题变成连续状态的 MDP——这是它的自然延伸，也是维数灾难的入口。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| 马尔可夫性 | 当前状态含全部历史信息 | 状态变量选得够 | 漏掉关键状态（年龄、累积损伤）使策略次优 10%—30% | 把一阶历史并进状态，看 $V$ 是否显著变化 |
| 状态/动作有限可枚举 | 能建表 | 题目是档位、年龄段、库存档 | 连续状态必须离散化，粗网格让策略在阈值附近抖动 | 网格减半重解，看最优动作是否稳定 |
| 转移概率已知或可估 | 能算期望 | 有历史数据或机理 | $P$ 的误差 $\varepsilon_P$ 使 $V$ 误差约 $\varepsilon_P R_{\max}/(1-\gamma)^2$ | 用 $P$ 的置信区间扰动后重解 |
| 奖励只依赖 $(s,a)$ | 收益当步结算 | 题目按年/按次结算 | 路径依赖（学习效应、习惯）被忽略，长期值系统性偏低 | 检查奖励里是否含"累计量" |
| $\gamma<1$ 或存在吸收态 | 级数收敛 | 有折扣或终止状态 | $\gamma=1$ 且无吸收时值发散，迭代永不收敛 | 检查 $\gamma^kR_{\max}/(1-\gamma)$ 是否可控 |
| 转移平稳（不随 $t$ 变） | 一套 $P$ 用到底 | 环境稳定 | 政策、气候、技术变化时策略迅速过时 | 分时段估 $P$，比较前后差异 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $s$ | 状态 | 各自（尾、kg、元） | 10—$10^4$ 个 | 决策时所处的"格位"；越细策略越准、表越大 |
| $a$ | 动作 | — | 2—10 个 | 可选的操作；越多单步计算越贵 |
| $P(s'\mid s,a)$ | 转移概率 | — | 每行和为 1 | 动作的"可靠性"；越接近 1 越确定 |
| $r(s,a)$ | 即时奖励 | 元/步、尾/年 | 与 $V$ 同单位 | 当步收益；尺度决定停止阈值 |
| $\gamma$ | 折扣因子 | — | 0.9—0.99 | 越接近 1 越看重长期，收敛越慢 |
| $V(s)$ | 值函数 | 与 $r$ 同单位 | $\le R_{\max}/(1-\gamma)$ | 从 $s$ 出发的最优总收益 |
| $Q(s,a)$ | 动作值 | 与 $r$ 同单位 | 同量级 | 先做 $a$ 再最优的总收益；$V=\max_aQ$ |
| $\pi$ | 策略 | — | 表或函数 | 每个状态选哪个动作 |
| $T$ | 时域 | 步 | 10—100 | 有限期问题的倒推长度 |

#### 6. 从基本原理推导核心公式

**第一步，定义回报。** 对策略 $\pi$，$V^\pi(s)=\mathbb E_\pi\big[\sum_{t\ge0}\gamma^tr(s_t,a_t)\mid s_0=s\big]$。由 $|r|\le R_{\max}$ 得 $|V^\pi|\le R_{\max}/(1-\gamma)$，绝对收敛，可以随意交换求和与期望。

**第二步，拆出第一步（全期望公式）。** 对 $t=0$ 的动作和转移取条件期望：
$$V^\pi(s)=\sum_{a}\pi(a\mid s)\Big[r(s,a)+\gamma\sum_{s'}P(s'\mid s,a)V^\pi(s')\Big].$$
括号里出现的是同一个 $V^\pi$——这一步只用了马尔可夫性，没有用别的。

**第三步，取 max 得 Bellman 最优方程。** 定义 $\mathcal T$ 为右端的 max 版本。两个方向合起来给出
$$V^*(s)=\max_{a\in\mathcal A}\Big\{r(s,a)+\gamma\sum_{s'}P(s'\mid s,a)V^*(s')\Big\}.$$
论证要点：任意 $\pi$ 的 $V^\pi$ 满足 $V^\pi\le\mathcal TV^\pi$，故 $V^\pi\le V^*$（$V^*$ 是 $\mathcal T$ 的不动点且不等式迭代保号）；而由 $V^*$ 取贪心得到的策略恰好达到 $V^*$。

**第四步，$\mathcal T$ 是 $\gamma$-压缩。** 对任意 $V_1,V_2$，
$$\big|\mathcal TV_1(s)-\mathcal TV_2(s)\big|\le\gamma\max_a\sum_{s'}P(s'\mid s,a)\big|V_1(s')-V_2(s')\big|\le\gamma\Vert V_1-V_2\Vert_\infty .$$
第一个不等号用 $|\max_af(a)-\max_ag(a)|\le\max_a|f(a)-g(a)|$，第二个用行和为 1。

**第五步，Banach 不动点 → 几何收敛。** 压缩映射有唯一不动点，且
$$\Vert V_k-V^*\Vert_\infty\le\gamma^k\Vert V_0-V^*\Vert_\infty\le\frac{\gamma^kR_{\max}}{1-\gamma}.$$
实用停止准则：若 $\Vert V_{k+1}-V_k\Vert_\infty\le\frac{\varepsilon(1-\gamma)}{2\gamma}$，则 $\Vert V_{k+1}-V^*\Vert_\infty\le\varepsilon/2$。达到 $\varepsilon=0.01$、$R_{\max}=5$、$\gamma=0.9$ 约需 $90$ 轮（由 $\gamma^kR_{\max}/(1-\gamma)\le\varepsilon/2$ 解出 $k\approx88$）；$\gamma=0.99$ 约需 $1150$ 轮——这就是"$\gamma$ 越接近 1 收敛越慢"的定量版本。

**第六步，策略迭代为什么通常更快。** 先解线性方程做策略评估，再对 $V^\pi$ 贪心改进；策略改进定理保证单调，策略数有限故有限步终止。它是 Newton 型方法：一轮的信息量远大于值迭代的一轮，通常 5—20 轮收敛，代价是每轮要解 $|\mathcal S|$ 维线性系统。

**第七步，与最优控制的关系。** 令 $\gamma=e^{-\rho\Delta t}$，把 Bellman 方程在 $\Delta t\to0$ 展开：
$$V(s)=\max_a\Big\{r(s,a)\Delta t+e^{-\rho\Delta t}\,\mathbb E[V(s')]\Big\}\ \Longrightarrow\ \rho V=\max_a\big\{r+\nabla V^\top\mathbf f\big\},$$
右端正是把 $\boldsymbol\lambda=\nabla V$ 代入后的 Hamiltonian。所以：**HJB（4.1）是连续时间的 Bellman 方程，PMP 是 HJB 的特征线**；离散线性二次情形下，Riccati 递推就是"值函数为二次型"这一假设代入 Bellman 方程的结果。

**第八步，量纲与极限。** $V$ 与 $r$ 同单位（$r$ 是每步奖励时，$V$ 是总奖励），$\gamma$ 无量纲，$\rho$ 的单位是 1/时间。$\gamma\to0$：$V\to\max_ar(s,a)$，策略变短视；$\gamma\to1$ 且有吸收态：折扣和趋于平均奖励准则，必须换目标；$R_{\max}\to\infty$ 且无折扣：值发散，模型无意义。

#### 7. 求解方法

- **值迭代**：$V_{k+1}\leftarrow\max_a\{r+\gamma PV_k\}$，每轮 $O(|\mathcal S|^2|\mathcal A|)$（稀疏 $P$ 可降一个量级）；同步更新稳定但慢，Gauss–Seidel（就地更新）常快 2—5 倍。
- **策略迭代**：评估 + 改进交替，轮数少但每轮贵；二者折中是**修正策略迭代**（只评估 $m$ 步再改进）。
- **优先扫描（prioritized sweeping）**：只更新"值变化大"的状态，适合稀疏转移。
- **LP 形式**：把 Bellman 不等式当约束求解，$|\mathcal S|$ 小的时候可用，好处是直接给对偶变量（状态的"约束价值"）。
- **连续状态**：网格离散（光滑 $V$ 时误差 $O(h^2)$，但策略在切换面附近误差 $O(h)$）；或函数逼近（线性基、神经网络）做拟合值迭代。
- **陷阱**：$\gamma$ 接近 1 时迭代次数按 $1/(1-\gamma)$ 增长；奖励尺度大时停止阈值要跟着缩放；不可达状态要从状态集里剔掉，否则它们把 sup 范数撑大、拖慢收敛。

#### 8. 适用场景

**该用（判据）。**

- 阶段天然清晰：年、月、批次、工序、充放电循环。
- 状态可枚举且维数低（1—3 维）：年龄段、库存档、电量档、健康档。
- 转移概率能从历史数据估计，或有可运行的随机模拟器。
- 需要**全局**最优，而不是局部改进；策略要能写成一张可执行的表。

**不该用。**

- 连续高维：每维 10 格，4 维就是 $10^4$ 个状态，再乘动作和迭代轮数就不可接受。
- 状态不可观测又没有信念建模：先用滤波（3.11）或直接改 POMDP。
- 转移概率完全未知且没有模拟器：应走强化学习，但要准备大量样本。
- 需要解释"为什么这个动作好"：DP 只给数值，PMP 的 $\lambda$ 更能讲故事（可用 $Q$ 的差补一部分）。

#### 9. 局限

- **维数灾难**：状态每加一维，规模乘 10—100；4 维以上基本不可行。
- **离散化误差**：网格粗会使策略在真实阈值附近抖动，收益差 2%—10%；网格细则规模爆炸。
- **转移概率估计误差被放大**：$V$ 的误差约 $\varepsilon_PR_{\max}/(1-\gamma)^2$，$\gamma=0.95$ 时放大约 400 倍。
- **折扣改变目标**：$\gamma$ 等价于每步有 $1-\gamma$ 的终止概率，报"长期最优"时必须说明这一点，否则结论无法与别的队比较。
- **平稳性假设**：政策、气候、技术一变，$P$ 与 $r$ 都变，旧策略失效。
- **不提供灵敏度**：不像 PMP 有 $\lambda$ 这种边际量，需额外算 $Q$ 的差或做扰动重解。

#### 10. 改进

1. **状态聚合/特征压缩**（代价小）：把相近状态合并，用聚合误差界说明损失。
2. **拟合值迭代**（代价：需要采样分布）：用回归近似 $V$，把状态从离散表变成特征。
3. **修正策略迭代 / 优先扫描**（代价：实现复杂）：在同样精度下省 2—10 倍计算。
4. **近似动态规划（ADP）/ 强化学习**（代价：样本量与调参）：模型未知或状态连续时的唯一出路。
5. **鲁棒 MDP / 分布鲁棒**（代价：规模乘不确定集大小）：解决 $P$ 估不准的问题。
6. **半马尔可夫/分层 MDP**（代价：建模复杂）：解决"每一步时间长度不同"和"多时间尺度"。

#### 11. 与其他模型的组合

- **+ 种群/资源模型（3.2）**：状态是资源量或年龄结构，动作是捕捞配额，转移来自随机增长模型；DP 给出"按资源量分档的最优配额表"。
- **+ 贝叶斯（4.3）**：参数 $\theta$ 不确定时，把后验当信念状态做信念 MDP，接口是"后验 → 转移核的混合"。
- **+ 最优控制（4.1）**：连续极限下 DP 收敛到 HJB；反过来可用 DP 的网格解给 PMP 的 BVP 提供全局初猜。

#### 12. 手算算例

两状态 $\{L,H\}$、两动作 $\{S,P\}$、$\gamma=0.9$：

| $s,a$ | $r$ | 转移到 |
|---|---|---|
| $L,S$ | 1 | $H$:0.8, $L$:0.2 |
| $L,P$ | 3 | $L$:0.7, $H$:0.3 |
| $H,S$ | 2 | $H$:0.9, $L$:0.1 |
| $H,P$ | 5 | $H$:0.6, $L$:0.4 |

**第一轮**（以即时奖励为 $V_1$）：$V_1(L)=\max(1,3)=3$，$V_1(H)=\max(2,5)=5$。

**第二轮**：
$$Q(L,S)=1+0.9(0.8\times5+0.2\times3)=5.14,\qquad Q(L,P)=3+0.9(0.7\times3+0.3\times5)=6.24,$$
$$Q(H,S)=2+0.9(0.9\times5+0.1\times3)=6.32,\qquad Q(H,P)=5+0.9(0.6\times5+0.4\times3)=8.78.$$
故 $V_2(L)=6.24$、$V_2(H)=8.78$，策略为"两状态都选 $P$"。

**第三轮**（补充，用来看单调性）：$Q(L,S)=8.4448$、$Q(L,P)=9.3018$、$Q(H,S)=9.6734$、$Q(H,P)=11.9876$，故 $V_3=(9.3018,\,11.9876)$，仍在上升且策略不变。

**不动点**（补充）：若两状态都选 $P$，解二元一次方程组 $V(L)=3+0.9(0.7V(L)+0.3V(H))$、$V(H)=5+0.9(0.6V(H)+0.4V(L))$，得 $V^*\approx(37.40,\,40.14)$。检查贪心：在 $V^*$ 处 $Q(L,P)=37.40>Q(L,S)=36.63$、$Q(H,P)=40.14>Q(H,S)=37.88$，策略自洽。注意 $V^*$ 远大于 $V_2$：$\gamma=0.9$ 时值迭代约 45 轮误差才降到 1 以内、约 90 轮才达到 $\varepsilon=0.01$（第 6 小节第五步），手算两轮只够看方向。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% 值迭代：两状态 {L,H}、两动作 {S,P}、gamma=0.9
Ps = [0.2 0.8; 0.1 0.9];      % 动作 S：行=当前 L,H，列=下一 L,H
Pp = [0.7 0.3; 0.4 0.6];      % 动作 P
r  = [1 3; 2 5];              % r(s,a)，行 s，列 a
gamma = 0.9; V = zeros(2,1);
for k = 1:500
    Q = [r(:,1) + gamma*Ps*V, r(:,2) + gamma*Pp*V];   % Bellman 右端
    Vn = max(Q, [], 2);                               % 逐状态取最优
    if norm(Vn - V, inf) < 1e-10, V = Vn; break; end
    V = Vn;
end
[~, pol] = max([r(:,1)+gamma*Ps*V, r(:,2)+gamma*Pp*V], [], 2);
fprintf('V = [%.4f; %.4f], 策略 = [%d; %d], 轮数 = %d\n', V(1), V(2), pol(1), pol(2), k);
```

输入：两张转移矩阵、奖励矩阵、$\gamma$。输出：$V^*$ 与策略。判据：跑两轮后手工核对 $V=[6.24;8.78]$（把循环上限临时改成 2）；收敛后应得 $V\approx[37.40;40.14]$、策略 $[2;2]$、轮数约 250（`tol` 取 $10^{-10}$，比第 6 小节的 $\varepsilon=0.01$ 严得多，故轮数按 $\ln 10^{-10}/\ln 0.9$ 量级增加）。

#### 14. 真实美赛 A 题应用

- **2020（渔业）**：分阶段投资与配额决策——状态是船队规模/资源量，动作是"买船/卖船/配额档"，转移来自资源随机增长，是标准的多期 MDP；折扣对应"50 年经营"的贴现。
- **2026（电池）**：若把使用策略写成离散模式（档位、充电时段），状态是 SOC 档与温度档，动作是模式选择，可用 MDP 求"按档位的最优策略表"；连续时间部分仍应交给 ODE/最优控制（4.1）。
- **对照**：2022 的功率分配是连续控制，不要硬套 MDP；离散化会抹掉 bang-bang 结构。

#### 15. 常见错误

1. **折扣语义不清**：把 $\gamma=0.9$ 当成"打九折的收益"，实际它等价于每步 10% 的终止概率，长期结论完全不同。
2. **状态漏变量**：只用当前资源量、不含年龄结构，马尔可夫性被破坏，策略次优还找不到原因。
3. **网格太粗**：阈值附近的动作来回跳，报出的策略表无法执行。
4. **忘记检查策略自洽**：算出 $V$ 后不回代贪心验证，报出的策略与 $V$ 不匹配（异步更新或提前停止时尤其），自己却以为已收敛。
5. **$\gamma=1$ 直接迭代**：值和发散，程序跑满迭代次数却"看起来收敛"。
6. **把 DP 用于连续高维**：状态 5 维还坚持建表，内存先爆。
7. **奖励单位不统一**：生态量（尾）与经济量（元）直接相加，$V$ 没有意义。
8. **只报策略不报值**：没有 $V$ 就无法比较"这个策略比基线好多少"。

#### 16. 学习要求（验收）

- 能写出 Bellman 最优方程，并说明它由"全期望公式 + 逐状态取 max"两步得到。
- 能独立手算本例两轮值迭代，得到 $V_2=(6.24,\,8.78)$，并解释为什么 $V^*$ 会大一个量级。
- 能证明 $\mathcal T$ 是 $\gamma$-压缩，并由此给出值迭代的误差界与停止准则。
- 能判断一道题该用 DP 还是最优控制，判据说到"时间连续/离散、状态能否枚举、转移是否已知"三条。

### 4.3 参数反演与贝叶斯估计

**定位**：P1｜难度 ★★★★｜7 h｜前置：3.7、概率｜2025 的核心工具。

#### 1. 定义

一句话：**正问题是"给定参数算结果"，反问题是"看着结果倒推参数"，贝叶斯则把"倒推"变成一个分布而不是一个数。**

正问题：$y=G(\boldsymbol\theta)+\boldsymbol\varepsilon$，给定 $\boldsymbol\theta\in\mathbb R^d$ 算 $y\in\mathbb R^n$。反问题：给定观测 $y$，求 $\boldsymbol\theta$。两种标准提法：

$$
\hat{\boldsymbol\theta}_{\text{LS}}=\arg\min_{\boldsymbol\theta}\ \Vert y-G(\boldsymbol\theta)\Vert^2,\qquad
p(\boldsymbol\theta\mid y)=\frac{p(y\mid\boldsymbol\theta)\,p(\boldsymbol\theta)}{\int p(y\mid\boldsymbol\theta')p(\boldsymbol\theta')\,d\boldsymbol\theta'} .
$$

已知：前向映射 $G$（可能是一次 ODE/PDE 数值解）、观测 $y$、噪声模型 $\boldsymbol\varepsilon$、先验 $p(\boldsymbol\theta)$。未知：$\boldsymbol\theta$ 的点估计、区间、以及"哪些方向根本估不出来"。注意 $G$ 通常没有解析逆，每次求值都要解一次正问题，所以反演的代价主要由 $G$ 的次数决定。

反问题的三种典型病症（英文术语首现即给出）：**不适定**（ill-posed，Hadamard 意义下解不存在、不唯一，或解不连续依赖数据）；**病态**（ill-conditioned，条件数大，小噪声被放大成大误差）；**不可辨识**（non-identifiable，多组参数给出同一输出，Fisher 信息矩阵奇异）。三者层次不同：前两者是"算得准不准"，可以靠正则化、更好的算法或更多数据缓解；后者是"根本算不出来"，任何算法都无济于事，只能改模型结构、补独立观测或引入先验（第 6 小节第五步给出判据）。

#### 2. 核心思想与大白话

笨办法是"试参数、对曲线、手调"：改一个参数、跑一次模型、看拟合好不好。它的问题不是慢，而是**不告诉你唯一性和误差**——参数改 5% 曲线几乎不变时，你无法判断是数据不够还是模型结构就不允许。

贝叶斯（Bayesian）把这件事做成概率运算：先验是"看数据前的常识"，似然是"这个参数下数据出现的可能性"，后验是"看完数据后的全部知识"。三者由贝叶斯公式连接，输出的不是一个点，而是一张"可信度地图"——地图上哪里窄（估得准）、哪里是一条长谷（估不准），一眼可见。

> **大白话。** 侦探破案。先验是案发前的嫌疑人名单和初始怀疑度；似然是每个嫌疑人与现场证据的吻合程度；后验是结案时每个人"还剩多少嫌疑"。如果两个嫌疑人永远同进同出（比如只知道两人身高之和），证据再强也分不开他们——这就是**不可辨识**（non-identifiability）。此时正确的结论不是"两人各多高"，而是"只能定出他们的和"，并说明还需要什么独立证据。

这个类比与数学一致：后验在不可辨识方向上是一条平台（平坦方向），沿该方向移动参数似然不变，数据不提供任何信息。

#### 3. 数学基础

- **贝叶斯公式与乘积法则**：$p(\boldsymbol\theta\mid y)\propto p(y\mid\boldsymbol\theta)p(\boldsymbol\theta)$，分母只是归一化常数（与 $\boldsymbol\theta$ 无关），实际计算只需右端。
- **条件独立**：观测独立时似然分解为 $\prod_ip(y_i\mid\boldsymbol\theta)$，负对数似然变成求和——这是所有最小二乘/极大似然的来源。
- **高斯噪声 → 加权最小二乘**：$y_i\sim N(G_i,\sigma_i^2)$ 时 $-\log p=\frac12\sum(y_i-G_i)^2/\sigma_i^2+\text{const}$，权重 $1/\sigma_i^2$ 就是"这条数据多可信"。
- **Fisher 信息矩阵**：$I(\boldsymbol\theta)=\mathbb E[(\nabla\log p)(\nabla\log p)^\top]$，衡量数据对参数的局部信息量；奇异（有零特征值）意味着某个方向完全没信息。
- **共轭先验（conjugate prior）**：先验与似然同族，后验仍属该族，可解析写出（正态—正态、Beta—二项、Gamma—泊松）。
- **Laplace 近似**：后验在众数附近用二次型近似，协方差取负 Hessian 的逆，把"算积分"降级为"算 Hessian"。
- **MCMC（Markov chain Monte Carlo）**：MH（Metropolis–Hastings）、Gibbs、HMC，用一条马尔可夫链的样本近似后验；收敛诊断看 $\hat R$、有效样本量 ESS、迹图。
- **剖面似然（profile likelihood）**：固定一个参数、对其余参数重新优化，得到"实际可辨识"的置信区间，比 Wald 区间可靠（非线性、边界时差别很大）。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| 前向模型 $G$ 结构正确 | 误差只来自噪声 | 机理已核对（3.12） | 模型误差被吸进参数，参数偏 10%—20%，区间还不含真值 | 残差图看结构；留出数据做预测检验 |
| 噪声模型正确（独立、同方差/已知 $\sigma$） | 似然可写 | 测量过程清楚 | 权重错、$\sigma$ 低估使 95% 区间实际覆盖只有 60% | 残差的正态性/异方差检验；重复测量 |
| 先验不主导后验 | 数据信息强于先验 | 先验来自独立知识 | 后验≈先验，形成"确认偏误"，结论是自证 | 换宽先验重算，看后验是否移动 |
| 参数可辨识 | 似然有唯一峰 | 观测设计合理 | 后验长谷，点估计随初值/先验漂移，区间无意义 | 看 Fisher 矩阵最小特征值、后验相关矩阵 |
| 观测独立或相关结构已知 | 似然可分解 | 采样设计独立 | 有效样本量被高估，区间缩窄 $\sqrt{n/n_{\text{eff}}}$ 倍 | 残差自相关；按时间分块 |
| 参数在观测期内不变 | 静态参数 | 过程平稳 | 时变参数被平均掉，反演结果只代表"某段平均" | 分段反演，比较前后 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $\boldsymbol\theta$ | 待反演参数 | 各自（1/s、Pa、—） | $d=1$—$10$ | 决定机理的常数；$d$ 越大越难辨识 |
| $y$ | 观测向量 | 各自（m、K、个） | $n=10$—$10^3$ | 数据；$n/d$ 越大信息越多，但不解决结构不可辨识 |
| $G(\boldsymbol\theta)$ | 前向映射 | 与 $y$ 同 | 一次 ODE/PDE 解 | "给定参数会看到什么" |
| $\sigma$ | 噪声标准差 | 与 $y$ 同 | 测量精度的量级 | 越大后验越宽，$\propto1/\sigma^2$ 的权重越小 |
| $p(\boldsymbol\theta)$ | 先验 | — | $N(\theta_0,\tau_0^2)$ | 看数据前的知识；$\tau_0$ 越小越"固执" |
| $p(\boldsymbol\theta\mid y)$ | 后验 | — | 近似高斯或长谷 | 全部结论的来源 |
| $I(\boldsymbol\theta)$ | Fisher 信息 | 1/参数单位$^2$ | 矩阵 $d\times d$ | 最小特征值 $\approx0$ 即不可辨识 |
| $n_{\text{eff}}$ | 有效样本量 | — | $\le n$ | 相关数据下真实信息量 |

#### 6. 从基本原理推导核心公式

**第一步，贝叶斯公式。** 由条件概率定义 $p(\boldsymbol\theta\mid y)=p(y,\boldsymbol\theta)/p(y)$ 与 $p(y,\boldsymbol\theta)=p(y\mid\boldsymbol\theta)p(\boldsymbol\theta)$，得
$$p(\boldsymbol\theta\mid y)=\frac{p(y\mid\boldsymbol\theta)p(\boldsymbol\theta)}{p(y)},\qquad p(y)=\int p(y\mid\boldsymbol\theta)p(\boldsymbol\theta)\,d\boldsymbol\theta .$$
分母与 $\boldsymbol\theta$ 无关，所以"正比于先验乘似然"不是近似，而是精确到常数。

**第二步，高斯情形化成最小二乘。** 取 $y_i\sim N(G_i(\boldsymbol\theta),\sigma^2)$ 且独立，则
$$-2\log p(y\mid\boldsymbol\theta)=\frac{1}{\sigma^2}\sum_{i=1}^n\big(y_i-G_i(\boldsymbol\theta)\big)^2+n\log(2\pi\sigma^2).$$
再取先验 $\boldsymbol\theta\sim N(\boldsymbol\theta_0,\Sigma_0)$，$-2\log p(\boldsymbol\theta)=(\boldsymbol\theta-\boldsymbol\theta_0)^\top\Sigma_0^{-1}(\boldsymbol\theta-\boldsymbol\theta_0)+\text{const}$，于是
$$-2\log p(\boldsymbol\theta\mid y)=\frac{\Vert y-G(\boldsymbol\theta)\Vert^2}{\sigma^2}+(\boldsymbol\theta-\boldsymbol\theta_0)^\top\Sigma_0^{-1}(\boldsymbol\theta-\boldsymbol\theta_0)+\text{const}.$$
**MAP（最大后验）就是"数据项 + 正则项"的最小化**：$\sigma\to0$ 退化为最小二乘（MLE），$\Sigma_0\to\infty$ 先验消失，$\Sigma_0$ 有限则等价于 Tikhonov 正则化。三个方法不是三套理论，是同一式子的三个极限。

**第三步，共轭先验给出解析后验（本小节的算例就用它）。** 设先验 $\theta\sim N(\mu_0,\tau_0^2)$、观测 $y_i\sim N(\theta,\sigma^2)$、$n$ 个独立样本、$\sigma$ 已知。配方：
$$-\frac{(\theta-\mu_0)^2}{2\tau_0^2}-\frac{\sum_i(y_i-\theta)^2}{2\sigma^2}=-\frac{(\theta-\mu)^2}{2\tau^2}+\text{const},$$
比较 $\theta^2$ 系数得**精度相加** $\frac{1}{\tau^2}=\frac{1}{\tau_0^2}+\frac{n}{\sigma^2}$，比较 $\theta$ 系数得 $\mu=\tau^2\big(\frac{\mu_0}{\tau_0^2}+\frac{n\bar y}{\sigma^2}\big)$。取 $\mu_0=0,\tau_0=1,\sigma=1,n=10$：精度 $=1+10=11$，故 $\mu=10\bar y/11$、标准差 $=1/\sqrt{11}=0.3015$——第 12 小节的 $1.7727$ 与 $0.3015$ 就是这两个式子的数值。

**第四步，Laplace 近似。** 在众数 $\hat{\boldsymbol\theta}$ 处展开 $\log p(\boldsymbol\theta\mid y)$：
$$\log p\approx\log p(\hat{\boldsymbol\theta}\mid y)-\tfrac12(\boldsymbol\theta-\hat{\boldsymbol\theta})^\top H(\boldsymbol\theta-\hat{\boldsymbol\theta}),\qquad H=-\nabla^2\log p(\boldsymbol\theta\mid y)\big|_{\hat{\boldsymbol\theta}},$$
故 $p(\boldsymbol\theta\mid y)\approx N(\hat{\boldsymbol\theta},H^{-1})$。误差是 $O(\Vert\boldsymbol\theta-\hat{\boldsymbol\theta}\Vert^3)$：后验单峰、近似对称时很好；偏斜、多峰、或参数贴边界时不能用。$H$ 的最小特征值直接对应"最不确定的方向"。

**第五步，结构性不可辨识：Fisher 信息奇异。** 取 $y_i=(a+b)x_i$、$\sigma=1$，参数 $\boldsymbol\theta=(a,b)^\top$。灵敏度矩阵
$$S=\frac{\partial G}{\partial\boldsymbol\theta}=\begin{pmatrix}x_1 & x_1\\ \vdots & \vdots\\ x_n & x_n\end{pmatrix},\qquad I=S^\top S=\Big(\sum_i x_i^2\Big)\begin{pmatrix}1&1\\1&1\end{pmatrix}.$$
$I$ 的秩为 1（两列相同），零空间方向 $\mathbf v=(1,-1)^\top/\sqrt2$：$S\mathbf v=\mathbf 0$。沿 $\mathbf v$ 移动参数（即 $a+\delta,b-\delta$），**所有预测值一字不变**，似然完全平坦——数据里根本没有 $a-b$ 的信息。这不是"数据少"，加多少观测都一样：只要 $G$ 只通过 $a+b$ 依赖参数，$I$ 就永远奇异。这个例子里 $a+b=2.997$（真值 3）可估，$a-b$ 不可估；任意分裂 $a=b=1.499$ 给出同样的拟合。

**第六步，长谷的几何图像。** 后验等高线是椭球，主轴方向由 $H$ 的特征向量给出，轴长 $\propto1/\sqrt{\lambda_i}$。这里 $\lambda_1=2\sum x_i^2/\sigma^2$（横向，很窄）、$\lambda_2=0$（纵向，无限长），轴比 $\sqrt{\lambda_1/\lambda_2}\to\infty$：后验是一条**长谷**（ridge）。后果有三个：(i) 点估计沿谷底任取，随初值/先验变；(ii) MCMC 在谷里做随机游走，横向步长受窄方向限制，纵向几乎不动，混合极慢，有效样本量 $\propto\lambda_2/\lambda_1$；(iii) 任何"$a$ 的 95% 区间"都会宽到无意义。**出路是重参数化**：改用 $(\phi_1,\phi_2)=(a+b,\,a-b)$，则 $\phi_1$ 有信息、$\phi_2$ 完全由先验决定，问题立刻变成"一个可估参数 + 一个纯先验参数"。

**第七步，量纲与极限检查。** $a$ 与 $b$ 必须同量纲，否则 $a+b$ 本身量纲不成立——**量纲不齐就是结构不可辨识的物理信号**（详见 2.9、2.14）。极限：$\sigma\to0$ 时后验在可辨识方向上收缩成点质量；$n\to\infty$ 时后验方差趋于 $I^{-1}$（Bernstein–von Mises）；$\tau_0\to\infty$ 时后验趋于似然；$\tau_0\to0$ 时后验锁死在先验上，与数据无关。

#### 7. 求解方法

- **最小二乘反演**：`lsqnonlin`（信赖域/Levenberg–Marquardt），给初值与边界，多初值排除局部解；每次迭代解一次正问题。
- **MAP**：把负对数后验写成一个函数，用 `fminunc`/`fmincon`；先验项是二次的，梯度好写。
- **Laplace 区间**：在 MAP 处数值求 Hessian（`hessian` 或有限差分），协方差取逆，再取对角线的平方根。
- **剖面似然**：对每个 $\theta_i$ 扫描网格，其余参数重新优化，得到 $\Delta(-2\log L)=3.84$ 的区间；非线性/边界时这是唯一可信的区间。
- **MCMC**：随机游走 MH（一维目标接受率约 0.44，高维约 0.234）；共轭块用 Gibbs；维数高用 HMC/NUTS（需梯度）。
- **诊断**：多链 $\hat R<1.01$、ESS $>400$、迹图无趋势、后验预测检验（PPC）能复现观测的分布特征。
- **陷阱**：$\sigma$ 与 $\boldsymbol\theta$ 同时估会出现尺度不可辨识（$\sigma$ 与参数整体缩放互补）；硬边界会让后验在边界堆积，应换成物理上合理的先验；不重参数化时长谷里 MCMC 几乎不动。

#### 8. 适用场景

**该用（判据）。**

- 有可信的正模型（$G$ 的机理已被独立验证），参数有物理意义。
- 数据少、但先验知识可用（文献值、同类实验、物理约束）。
- 题目要"区间/可信度/哪些结论不可辨识"，而不只是"给个数"。
- 需要把多源数据（不同精度）融合：似然里按 $1/\sigma_i^2$ 加权即可。

**不该用。**

- 正模型本身不可信：此时反演只是把模型误差塞进参数，结论不可辩护。
- 纯预测、不关心参数：直接用回归/机器学习，别绕道反演。
- 参数维数极高又无结构（像素级反演）：需要正则化/降维/伴随法，超出 5 天赛程。

#### 9. 局限

- **计算贵**：每次似然一次正问题求解；MCMC $10^4$—$10^6$ 次采样 → 小时到天，PDE 情形必须用代理模型。
- **先验主观**：先验强度用 $\tau_0$ 度量，若 $\tau_0$ 小于后验宽度，后验被先验拉走，偏差可达 1—2 个 $\sigma$。
- **模型误差被吸收**：结构错 5% 可让参数偏 10%—20%，且 95% 区间覆盖真值的概率接近 0。
- **不可辨识时输出无意义**：长谷里"$a$ 的区间"可以任意宽，报出来等于没报。
- **高维混合慢**：$d>20$ 时随机游走 MH 的 ESS 随维数指数下降，必须用梯度方法。
- **局部极值**：非线性最小二乘多解，单初值的结果可能只是局部最优。

#### 10. 改进

1. **重参数化**（代价几乎为零）：换成可辨识组合（如 $a+b$、$a-b$），长谷立刻变成"一窄一平"，诊断也清楚。
2. **剖面似然替代 Wald 区间**（代价：每个参数 $O(20)$ 次重优化）：解决非线性与边界附近的区间失真。
3. **分层贝叶斯（partial pooling）**（代价：模型复杂度）：多个个体共享超参数，小样本个体向总体收缩，解决"每个个体都估不准"。
4. **加模型偏差项（Kennedy–O'Hagan）**（代价：需额外先验）：用高斯过程吸收结构误差，避免它被塞进物理参数。
5. **代理模型 + MCMC**（代价：训练一次，采样便宜 $10^2$—$10^3$ 倍）：PDE 反演的标准做法。
6. **伴随法算梯度**（代价：实现一次伴随方程）：把梯度代价从 $O(d)$ 次正问题降到 $O(1)$ 次。

#### 11. 与其他模型的组合

- **+ ODE（3.2）**：正模型是 ODE，参数是速率常数；外层 `lsqnonlin`、内层 `ode45`，接口是"参数 → 轨迹 → 残差"。
- **+ 敏感性分析（3.10）**：先用 Sobol 看哪些参数主导输出，再决定优先补测哪个——比盲目加数据省一个量级。
- **+ 最优控制 / MDP（4.1、4.2）**：先反演参数，再在名义值上优化（certainty equivalence）；严谨做法是从后验采样多组参数、重做优化，得到"策略性能的分布"。

#### 12. 手算算例

**算例 A（不可辨识，已核对）。** 模型 $y=(a+b)x$，真值 $a+b=3$。用数据拟合得 $a+b=2.997$，且 $a=b=1.499$。

验证"任意分裂同样好"：取 $a=1.499+\delta$、$b=1.499-\delta$，则 $a+b=2.998$ 不变（$\delta$ 任意），所有预测 $y_i=(a+b)x_i$ 完全相同，残差平方和一模一样。故 $\delta$ 无法由数据确定。

再看 Fisher：$I=(\sum x_i^2)\begin{pmatrix}1&1\\1&1\end{pmatrix}$，特征值 $2\sum x_i^2$ 与 $0$，零空间方向 $(1,-1)^\top$。结论：**只有和可辨识**；若物理上 $a,b$ 都重要，必须补充独立观测（例如另一组不同 $x$ 的实验，或直接测其中一个）或加先验。

**算例 B（MH 与解析后验对照，已核对）。** 10 个观测 $y_i\sim N(\theta,1)$，先验 $\theta\sim N(0,1)$，样本均值 $\bar y=1.95$。

后验（由第三步的共轭公式）：精度 $1+10=11$，故
$$\theta\mid y\sim N\!\left(\frac{10\times1.95}{11},\ \frac{1}{11}\right)=N(1.7727,\ 0.3015^2),$$
即后验均值 $1.7727$、标准差 $1/\sqrt{11}=0.3015$。

用随机游走 MH（提议 $\sigma=0.5$）采样 2 万次（丢前 1000 步预热），得均值 $1.777$、标准差 $0.303$、接受率约 $0.67$——与解析解在蒙特卡洛误差内一致，说明采样器写对了。两点旁注：(i) 一维随机游走的最优接受率约 $0.44$，$0.67$ 偏高说明步长偏小、样本相关性强，ESS 会明显小于 2 万，但均值仍无偏；(ii) 均值差 $1.777-1.7727=0.0043$，若 ESS 约 2000，则标准误约 $0.3015/\sqrt{2000}\approx0.0067$，差值在 1 个标准误内，正常。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% Metropolis-Hastings：y_i ~ N(theta,1)，先验 theta ~ N(0,1)
rng(7);
y = 1.95 + randn(10,1);  y = y - mean(y) + 1.95;    % 强制样本均值 = 1.95
logpost = @(t) -0.5*t^2 - 0.5*sum((y - t).^2);       % 先验 + 似然（常数略）
N = 20000; sig = 0.5; th = zeros(N,1); acc = 0;
for k = 2:N
    prop = th(k-1) + sig*randn;                      % 对称随机游走提议
    if log(rand) < logpost(prop) - logpost(th(k-1))
        th(k) = prop; acc = acc + 1;                 % 接受
    else
        th(k) = th(k-1);                             % 拒绝则复制上一样本
    end
end
th = th(1001:end);                                   % 丢掉预热段
fprintf('mean=%.4f (解析 1.7727), sd=%.4f (解析 0.3015), acc=%.2f\n', ...
        mean(th), std(th), acc/N);
```

输入：观测 $y$、先验、提议步长 `sig`、链长 `N`。输出：后验均值/标准差/接受率。判据：均值应接近 $1.7727$（误差 $\le2$ 个标准误），标准差接近 $0.3015$，接受率落在 $0.2$—$0.7$；把 `sig` 改成 $0.1$ 或 $2$，接受率会明显偏离、ESS 变差，但均值仍应稳定——这是"链收敛但效率不同"的标准现象。

#### 14. 真实美赛 A 题应用

- **2025（核心）**：由楼梯磨损深度场反演使用强度与方向偏好，必须讨论可辨识性——由 2.14 的 $\pi$ 群可知只有 $nK$ 可辨识，所以"每天多少人"必须附上"$K$ 取什么值"的说明；这正是贝叶斯/剖面似然要回答的问题。
- **2020（渔业）**：由渔获量与捕捞努力量反演资源量与可捕系数 $q$，是典型的"观测少、参数多"反演；必须给区间并说明 $q$ 与资源量的权衡。
- **2026（电池）**：由充放电曲线反演内阻 $R$、热容 $C_{th}$、散热系数 $hA$；OCV 平台段 SOC 不可辨识，需在方法里明说并只用有信息的区段。

#### 15. 常见错误

1. **不检查可辨识性就报"唯一历史"**：$a+b$ 型模型给出漂亮拟合，却把不可辨识的分裂当成结论。
2. **先验过强**：为"稳定"把先验设得很窄，后验等于先验，形成确认偏误。
3. **只报点估计**：没有区间，评委无法判断结论的可靠度。
4. **不做后验预测检验**：后验能拟合历史，但预测新情景系统性偏离，说明模型结构有问题。
5. **把模型误差当噪声**：残差有结构还硬说"噪声大"，导致 $\sigma$ 被高估、参数区间虚宽。
6. **Wald 区间用于非线性/边界**：区间可能跨到物理上不可能的区域，应改用剖面似然。
7. **MCMC 不查收敛**：单链、无 $\hat R$、迹图有趋势，却把前 100 步当后验。
8. **单位不统一**：$y$ 用 mm、$G$ 用 m，反演出的参数差 1000 倍，且拟合残差看起来"很大"。

#### 16. 学习要求（验收）

- 能写出正问题、似然、先验、后验四件事，并说明 MLE、MAP、贝叶斯三者的关系（同一式子的三个极限）。
- 能对正态—正态共轭情形独立推出"精度相加"公式，并算出 $1.7727$ 与 $0.3015$。
- 能用 Fisher 矩阵证明 $y=(a+b)x$ 只有和可辨识，并画出后验长谷的几何图像、说明 MCMC 为何在谷里混合慢。
- 能给出一个"点估计 + 区间 + 明确写出哪些结论不可辨识"的完整反演结论。

### 4.4 随机过程与 Monte Carlo

**定位**：P1｜难度 ★★★｜6 h｜前置：概率、3.11｜随机输入、风险与区间。

#### 1. 定义

随机过程（stochastic process）是概率空间 $(\Omega,\mathcal F,\mathbb P)$ 上的一族随机变量 $\{X_t\}_{t\in T}$，$T$ 为时间指标集。建模时我们只钉死它的**转移规律**，不钉死单条轨迹：马尔可夫过程满足 $\mathbb P(X_{t+\Delta}\in A\mid \mathcal F_t)=\mathbb P(X_{t+\Delta}\in A\mid X_t)$；泊松过程是计数型马尔可夫过程；扩散过程写成随机微分方程（SDE）
$$dX_t=a(X_t,t)\,dt+b(X_t,t)\,dW_t,$$
其中 $a$ 是漂移（drift）、$b$ 是扩散（diffusion）、$W_t$ 是标准布朗运动。Monte Carlo（MC）则是用样本均值估计期望：
$$\theta=\mathbb E[g(X)]\ \longrightarrow\ \hat\theta_M=\frac1M\sum_{m=1}^{M}g(X^{(m)}),$$
已知量是 $a,b$ 或转移核与样本量 $M$，未知量是 $\theta$、分位数与路径分布。

#### 2. 核心思想与大白话

确定性模型只给一条轨迹；可现实里"明年是旱是涝"本身就有概率。随机过程把这种"可能性"写进状态转移，Monte Carlo 则把"算不出来的积分/期望"换成"大量独立抽样后的平均数"——它替代的笨办法是解析积分或高维网格求积，后者在维数 $d$ 大时误差按 $M^{-1/d}$ 衰减，几乎不可用。

为什么可行：大数定律保证 $\hat\theta_M\to\theta$，中心极限定理给出误差尺度 $\sigma_g/\sqrt M$。所以 MC 的精度**只取决于样本量和方差**，与维数无关——这是它最重要的结构性优势。

**大白话：** 你算不出"所有可能的明天"平均是什么样，那就掷骰子造出一万个平行的明天，各跑一遍，再数一数。想知道"50 年内崩溃的概率"，就是数一万条平行时间线里有几条崩了。掷的骰子必须真的按你假设的分布掷（种子要固定、样本要够多），否则数出来的只是"你编的那几个明天"。

#### 3. 数学基础

1. **条件期望与马尔可夫性**：把"已知现在，未来与过去无关"写成 $\mathbb E[f(X_{t+s})\mid\mathcal F_t]=\mathbb E[f(X_{t+s})\mid X_t]$，这是所有转移核推导的起点。
2. **独立增量与稀有事件极限**：$P(\text{一步内发生两次})=o(\Delta t)$，由此得到泊松分布与指数间隔（第 6 小节）。
3. **Itô 微积分**（Itô calculus）：$(dW_t)^2=dt$ 导致普通链式法则失效，多出 $\frac12b^2f''$ 项；GBM 的显式解全靠它。
4. **大数定律与中心极限定理**：给 MC 的一致性收敛速度，$\sqrt M(\hat\theta-\theta)\Rightarrow N(0,\sigma_g^2)$。
5. **全方差公式**：$\mathrm{Var}(Y)=\mathbb E[\mathrm{Var}(Y\mid X)]+\mathrm{Var}(\mathbb E[Y\mid X])$，用来把输出方差拆成"参数"与"结构"两部分（接 3.11）。
6. **测度变换**：重要性抽样用 Radon–Nikodym 导数 $w=p/q$ 换分布，是稀有事件估计的核心工具。
7. **混合与遍历性**：马尔可夫链的收敛速度由谱隙（spectral gap）控制，决定"跑多少步才忘掉初值"。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| 增量独立、平稳 | 未来增量不依赖历史，率 $\lambda$ 不随时间变 | 泊松/指数模型的公理 | 计数方差大于均值（过散布），区间偏窄且随过散布程度加剧 | 用方差/均值比与自相关函数检验 |
| 马尔可夫性 | 下一步只依赖当前状态 | 状态变量选得足够全 | 遗漏状态导致转移概率有偏、长期概率漂移 | 加滞后状态做似然比检验 |
| 参数时不变 | $\mu,\sigma,\lambda$ 在模拟期内不变 | 短时段近似 | 趋势被当成噪声，长期预测系统性偏低或偏高 | 分段估计参数，看是否越出置信带 |
| 噪声为高斯/白噪声 | $dW$ 独立、$N(0,dt)$ | 中心极限、聚合效应 | 厚尾事件被低估，VaR/CVaR 被系统性低估 | 残差 QQ 图、超额峰度 |
| 样本独立同分布 | MC 的 $X^{(m)}$ 互不影响 | 独立抽样 | 标准误被低估，区间过窄 | 批次均值法、自相关检查 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $t$ | 时间 | s（生态用 yr） | $0\sim T$ | 过程推进的自变量 |
| $X_t$ | 状态量（数量、价格、水位） | 随问题 | — | 我们关心的随机输出 |
| $W_t$ | 标准布朗运动 | —（$\sqrt{\mathrm{s}}$ 量纲的载体） | 增量 $N(0,\Delta t)$ | 累积的微小随机扰动 |
| $\lambda$ | 事件到达率 | $\mathrm{s^{-1}}$ | $0.01\sim10$ | 越大事件越密，间隔均值 $1/\lambda$ |
| $\mu$ | 漂移率 | $\mathrm{s^{-1}}$ | $10^{-3}\sim10^{-1}$ | 越大趋势越强 |
| $\sigma$ | 波动率 | $\mathrm{s^{-1/2}}$ | $0.01\sim1$ | 越大路径越毛糙、尾部越厚 |
| $\Delta t$ | 时间步长 | s | $T/10^3\sim T/10^4$ | 越小离散偏差越小、计算越贵 |
| $M$ | MC 样本数/路径数 | —（条） | $10^3\sim10^6$ | 越大标准误按 $1/\sqrt M$ 降 |
| $\alpha$ | VaR 置信水平 | — | $0.95,0.99$ | 越大越关注尾部 |

#### 6. 从基本原理推导核心公式

**（1）泊松过程与指数间隔互推。** 设 $N(t)$ 有独立增量，且 $P(N(t+h)-N(t)=1)=\lambda h+o(h)$、$P(\ge2)=o(h)$。先看"零事件"概率：
$$p_0(t+h)=p_0(t)\big(1-\lambda h+o(h)\big)\ \Longrightarrow\ \dot p_0=-\lambda p_0\ \Longrightarrow\ p_0(t)=e^{-\lambda t}.$$
首达时刻 $T_1$ 满足 $\{T_1>t\}=\{N(t)=0\}$，故 $P(T_1>t)=e^{-\lambda t}$，即 $T_1\sim\mathrm{Exp}(\lambda)$，密度 $\lambda e^{-\lambda t}$；并且
$$P(T_1>s+t\mid T_1>s)=\frac{e^{-\lambda(s+t)}}{e^{-\lambda s}}=e^{-\lambda t},$$
这就是无记忆性——"已经等了 $s$ 年"不改变还要等多久的分布。

反方向：若间隔 $T_k$ 独立同分布于 $\mathrm{Exp}(\lambda)$，则 $S_n=\sum_{k\le n}T_k\sim\mathrm{Gamma}(n,\lambda)$，且
$$P(N(t)=n)=P(S_n\le t<S_{n+1})=\int_0^t\frac{\lambda^n s^{n-1}e^{-\lambda s}}{(n-1)!}\,e^{-\lambda(t-s)}\,ds=\frac{(\lambda t)^n}{n!}e^{-\lambda t},$$
最后一步用 $\int_0^t s^{n-1}ds=t^n/n$。量纲检查：$\lambda$ 的单位是 $\mathrm{s^{-1}}$，$\lambda t$ 无量纲，这是泊松分布唯一可能的自变量组合。极限检查：$\lambda t\ll1$ 时 $P(N=1)\approx\lambda t$、$P(N\ge2)\approx(\lambda t)^2/2=o(\lambda t)$，与公理自洽。

**（2）生灭过程的平衡方程与稳态。** 状态 $n$，出生率 $\lambda_n$、死亡率 $\mu_n$（单位 $\mathrm{s^{-1}}$）。概率流：$n\to n+1$ 为 $\pi_n\lambda_n$，$n+1\to n$ 为 $\pi_{n+1}\mu_{n+1}$。稳态时每个状态的净流入为零：
$$\pi_n(\lambda_n+\mu_n)=\pi_{n-1}\lambda_{n-1}+\pi_{n+1}\mu_{n+1}.$$
一维生灭过程可逆，细节平衡（detailed balance）逐条成立：$\pi_n\lambda_n=\pi_{n+1}\mu_{n+1}$，迭代得
$$\pi_n=\pi_0\prod_{k=0}^{n-1}\frac{\lambda_k}{\mu_{k+1}},\qquad \pi_0=\Big(1+\sum_{n\ge1}\prod_{k=0}^{n-1}\frac{\lambda_k}{\mu_{k+1}}\Big)^{-1}.$$
M/M/1：$\lambda_n=\lambda$、$\mu_n=\mu\ (n\ge1)$，故 $\pi_n=(1-\rho)\rho^n$，$\rho=\lambda/\mu$；归一化要求 $\rho<1$，否则级数发散、不存在稳态（队列无限增长）。平均队长 $L=\sum_n n\pi_n=\rho/(1-\rho)$，Little 律给 $W=L/\lambda=1/(\mu-\lambda)$。量纲：$\rho$ 无量纲，$W$ 的单位是时间 ✓。再看两个容易写错的变体：常数迁入 + 人均死亡（$\lambda_n=\lambda$、$\mu_n=\mu n$）给泊松分布 $\pi_n=e^{-a}a^n/n!$、$a=\lambda/\mu$；而线性生灭（分支过程，$\lambda_n=\lambda n$、$\mu_n=\mu n$）的 $0$ 是**吸收态**、根本没有稳态分布，其关键结论是灭绝概率 $q=\min(1,\mu/\lambda)$，$\lambda>\mu$ 且条件于未灭绝时准平稳分布才是几何分布（比值 $\mu/\lambda$），不是 $(\lambda/\mu)^n$——这一点常被写错。

**（3）GBM 的 Itô 解。** 对 $f(S)=\ln S$ 用 Itô 公式（关键是 $(dW)^2=dt$）：
$$d\ln S=\Big(\mu S\cdot\frac1S+\frac12\sigma^2S^2\cdot\big(-\frac1{S^2}\big)\Big)dt+\sigma S\cdot\frac1S\,dW=\big(\mu-\tfrac12\sigma^2\big)dt+\sigma\,dW.$$
积分得
$$S_t=S_0\exp\Big[\big(\mu-\tfrac12\sigma^2\big)t+\sigma W_t\Big].$$
检查：$\mathbb E[S_t]=S_0e^{\mu t}$（用 $\mathbb E e^{\sigma W_t}=e^{\sigma^2t/2}$），中位数为 $S_0e^{(\mu-\sigma^2/2)t}$，二者相差因子 $e^{\sigma^2t/2}$——这就是"均值被少数暴涨路径拉高"的数学来源。极限：$\sigma\to0$ 退回 $S_0e^{\mu t}$；$\mu=0$ 时 $e^{-\sigma^2t/2+\sigma W_t}$ 是鞅。量纲：$\sigma$ 的单位是 $\mathrm{s^{-1/2}}$，$\sigma W_t$ 无量纲，指数整体无量纲 ✓。

**（4）Euler–Maruyama 的强/弱收敛阶。** 格式为 $X_{k+1}=X_k+a\Delta t+b\sqrt{\Delta t}Z_k$。单步局部误差的主项来自 Itô 修正 $\frac12bb'(X_k)\big((\Delta W_k)^2-\Delta t\big)$，其均方根为 $O(\Delta t)$；$T/\Delta t$ 步近似独立累积，方差相加：
$$\mathrm{Var}\big(X_T-X_T^{\Delta t}\big)\sim\frac{T}{\Delta t}\cdot C\Delta t^2=CT\Delta t\ \Longrightarrow\ \big(\mathbb E|X_T-X_T^{\Delta t}|^2\big)^{1/2}\le C\Delta t^{1/2}.$$
弱误差看的是 $\mathbb E g(X_T)$：同一项的期望为零（$\mathbb E(\Delta W)^2=\Delta t$），一阶偏差抵消，剩 $O(\Delta t^2)$，累积 $T/\Delta t$ 步得 $O(\Delta t)$。所以 EM 的**强收敛阶为 $1/2$、弱收敛阶为 $1$**（原文作"强收敛一阶"，此处按标准结果给出，待复核）。实践含义：估期望用弱阶，$\Delta t$ 减半偏差减半；估路径、首达时间、极值必须用强阶，此时改用 Milstein 格式（加 $\frac12bb'\big((\Delta W)^2-\Delta t\big)$，强阶升到 1）。检验办法：$\Delta t$ 逐次减半，强误差应按 $\sqrt2$ 下降、弱误差应按 2 下降；不按这个斜率走，说明步长还没进入渐近区。

**（5）MC 的方差与 $O(M^{-1/2})$。** 若 $X^{(m)}$ 独立同分布、$\mathrm{Var}(g(X))=\sigma_g^2<\infty$，则
$$\mathrm{Var}(\hat\theta_M)=\frac{\sigma_g^2}{M},\qquad \mathrm{SE}(\hat\theta_M)=\frac{\sigma_g}{\sqrt M}.$$
中心极限定理给 $\sqrt M(\hat\theta_M-\theta)\Rightarrow N(0,\sigma_g^2)$，故 95% 区间半宽约 $1.96\sigma_g/\sqrt M$。两个直接推论：误差与维数 $d$ 无关（对比网格求积的 $O(M^{-1/d})$）；每多要一位有效数字，样本量要乘 100。这就是"跑 100 次不够、跑 1 万次才有两位"的定量理由。

**（6）方差缩减为什么有效。** 目标都是"不增加 $M$，只降低 $\sigma_g^2$"。
- **对偶变量（antithetic variates）**：用同一组 $Z$ 与 $-Z$，
$$\hat\theta=\frac12\big[g(Z)+g(-Z)\big],\qquad \mathrm{Var}(\hat\theta)=\frac{\sigma_g^2}{2}\big(1+\rho\big),\quad \rho=\mathrm{Corr}\big(g(Z),g(-Z)\big).$$
$g$ 单调时 $g(Z)$ 与 $g(-Z)$ 反向，$\rho<0$，方差严格小于 $\sigma_g^2/2$；$g$ 接近线性时 $\rho\to-1$，几乎零成本地砍掉一半方差。若 $g$ 非单调，$\rho$ 可能为正，反而更差——必须先小样本测 $\rho$。
- **控制变量（control variates）**：找一个期望已知的伴随量 $C$，
$$\hat\theta_\beta=\hat\theta-\beta(\bar C-\mathbb E C),\qquad \beta^*=\frac{\mathrm{Cov}(g,C)}{\mathrm{Var}(C)},\qquad \mathrm{Var}_{\min}=\sigma_g^2\big(1-\rho_{gC}^2\big).$$
即"用已知答案的相似问题当脚手架"，残差相关性越高省得越多；$\rho_{gC}=0.9$ 时方差降到 19%，等效样本量放大约 5 倍。
- **重要性抽样（importance sampling）**：换分布采样，
$$\mathbb E_p[g]=\mathbb E_q\Big[g\frac{p}{q}\Big],\qquad \mathrm{Var}_q=\mathbb E_p\Big[g^2\frac pq\Big]-\theta^2.$$
取 $q^*(y)\propto|g(y)|p(y)$ 时方差为零（$g\ge0$），实际中用一个把采样集中在"重要区域"的近似 $q$（如把均值平移到稀有事件阈值 $a$）。代价是 $q$ 的尾部不能比 $p$ 轻，否则权重 $p/q$ 无界、方差爆炸——这是 IS 最常见的翻车方式。
- **分层/拉丁超立方（LHS）**：把 $[0,1]$ 切成 $M$ 层各抽一个，去掉"层间"那部分方差，同样 $M$ 下通常更稳。

实践中这几种手段在样本量不变时通常可降方差 2—10 倍（原文数字）；但收益依赖问题结构，必须报告"缩减前后各跑一次"的实测方差，不能只宣称"用了方差缩减"。

**量纲与极限总检查：** 漂移项 $a\,dt$ 与状态同量纲；扩散项 $b\,dW$ 中 $b$ 带"状态/$\sqrt{\text{时间}}$"的量纲，乘 $\sqrt{\Delta t}$ 后与状态同量纲，故 EM 两项可相加。$\sigma\to0$ 时 SDE 退化为 ODE；$\Delta t\to0$ 时 EM 收敛到真解；$M\to\infty$ 时 $\hat\theta\to\theta$。三条极限任一不成立，模型或代码就有问题。

#### 7. 求解方法

1. **解析优先**：泊松计数、指数间隔、M/M/1 稳态、GBM 的期望/分位数都有闭式，先写出来当基准。
2. **路径模拟**：定 $\Delta t$ 使 $T/\Delta t\ge10^3$，用 EM；对价格/浓度等必须为正的量改用 log-Euler（GBM 下它在格点上无离散偏差）。
3. **收敛检验**：$\Delta t$ 减半看弱误差是否减半；$M$ 加倍看标准误是否降 $\sqrt2$、分位数变化是否 $<5\%$。
4. **MATLAB 工具**：`randn`（正态）、`rand`（均匀）、`exprnd`（指数）、`prctile`（分位数）、`histcounts`（分布）；SDE 无通用内置求解器，自己写时间循环最稳（Econometrics Toolbox 的 `simulate` 面向 ARMA/GARCH，不是通用 SDE）。
5. **常见数值陷阱**：负值越界（截断会改变分布，优先换对数变量）；用同一随机数流比较不同参数却不声明（这会引入相关性，除非你**故意**用公共随机数做配对比较）；只报均值不报区间；把 $\Delta t$ 造成的偏差误当成模型误差。

#### 8. 适用场景

**该用：** 输入本身是随机变量或随机过程（降雨、到达、故障）；结论需要概率、区间或风险度量（VaR/CVaR、超限概率、崩溃概率）；系统维数高、解析积分不可行；需要把参数不确定性一起传播（接 3.11）。

**不该用：** 系统本质确定、只是参数未知——那是认知不确定性，应做区间/贝叶斯（4.3），不是给它硬套随机过程；只要均值、且 ODE 已够用；稀有事件概率小于 $10^{-6}$ 而又没有解析尾部或重要性抽样方案（朴素 MC 要 $10^8$ 量级样本）；数据不足以标定转移概率或波动率。

#### 9. 局限

1. **分布假设主观**：把"干旱年"设成两状态，本身丢掉了连续干旱强度；假设错了，区间再窄也是假的。可用拟合检验（KS、QQ 图）兜底，但无法证明假设正确。
2. **稀有事件方差大**：$p=10^{-4}$ 时相对标准误约 $1/\sqrt{Mp}$，要 10% 就得 $M\sim10^6$。
3. **SDE 参数难估**：连续时间路径不可全观测，$\mu,\sigma$ 常不可辨识（漂移与扩散在离散采样下混在一起）。
4. **离散化偏差**：EM 弱阶只有 1，$\Delta t=0.01$ 时期望的相对偏差约 $10^{-2}$ 量级，做精细比较时必须减半验证。
5. **混合慢**：接近临界（$\rho\to1$、$\lambda\approx\mu$）时链的混合时间发散，有限步模拟没有代表稳态。

#### 10. 改进

1. **换格式**（代价最低）：Milstein 提强阶到 1；log-Euler 保正；隐式格式处理刚性。
2. **换采样**（中等）：对偶、控制变量、分层/LHS、重要性抽样；稀有事件用交叉熵自适应找 $q$。
3. **拟蒙特卡洛（QMC）**：Sobol/Halton 序列把误差从 $M^{-1/2}$ 改善到接近 $M^{-1}$（低有效维时），代价是实现复杂、高维收益衰减。
4. **多层蒙特卡洛（MLMC）**：把不同 $\Delta t$ 的层级配对，把总代价从 $\varepsilon^{-3}$ 降到 $\varepsilon^{-2}\log\varepsilon$ 量级。
5. **贝叶斯校准**：给 $\mu,\sigma,\lambda$ 先验，用观测做后验抽样（4.3），把参数不确定性与过程噪声分开传播。

#### 11. 与其他模型的组合

- **随机过程 + UQ（3.11）**：过程给出"偶然不确定性"，参数后验给出"认知不确定性"，两者叠加才是完整预测区间。
- **随机过程 + 参数反演（4.3）**：2025 磨损、2020 渔业都可用"确定性正模型 + 随机观测噪声"写似然，再用 MCMC 抽样。
- **随机过程 + 动态规划（4.2）**：状态转移概率就是 MDP 的 $P(s'\mid s,a)$，MC 用来评估策略。
- **随机过程 + PDE（3.5）**：把 PDE 的系数/初值当随机场，得到随机 PDE；MC 是唯一的通用传播手段。

#### 12. 手算算例

**（a）SDE 单步（原文数字）。** 取 $S_0=1$、$\mu=0.1\ \mathrm{yr^{-1}}$、$\sigma=0.3\ \mathrm{yr^{-1/2}}$、$\Delta t=0.01\ \mathrm{yr}$、$Z_1=1$：
$$\underbrace{0.1\times1\times0.01}_{\text{漂移}=0.001}+\underbrace{0.3\times1\times\sqrt{0.01}\times1}_{\text{扩散}=0.03}=0.031,\qquad S_1=1+0.001+0.03=1.031.$$
第二步取 $Z_2=-1$：漂移 $0.1\times1.031\times0.01=0.001031$，扩散 $0.3\times1.031\times0.1\times(-1)=-0.03093$，故 $S_2=1.031+0.001031-0.03093=1.001101$。量纲核对：$0.03$ 与 $1$ 同量纲，可直接相加 ✓。

**（b）泊松计数。** $\lambda=0.2\ \mathrm{yr^{-1}}$，5 年内无事件概率 $e^{-0.2\times5}=e^{-1}=0.368$；至少一次为 $0.632$；等待时间均值 $1/\lambda=5$ yr。

**（c）M/M/1 稳态。** $\lambda=3\ \mathrm{d^{-1}}$、$\mu=4\ \mathrm{d^{-1}}$，$\rho=0.75$，$L=0.75/0.25=3$ 人，$W=1/(4-3)=1$ d。校验：$L=\lambda W=3\times1=3$ ✓。

**（d）2023 干旱的两状态链（示意数值）。** 设年类型为"正常/干旱"，转移矩阵（行：今年，列：明年）
$$P=\begin{pmatrix}0.7&0.3\\0.6&0.4\end{pmatrix},$$
即正常年有 0.3 转干旱，干旱年有 0.6 转正常。平稳分布解 $\pi_D=0.3\pi_N/0.6=0.5\pi_N$，归一化得 $\pi_N=2/3$、$\pi_D=1/3$；平均干旱持续 $1/0.6=1.67$ yr、正常持续 $1/0.3=3.33$ yr。50 年内期望干旱年数 $50/3\approx16.7$；若每个干旱年独立地以 $q=0.05$ 触发崩溃，则 $P(\text{50 年内崩溃})\approx1-(1-0.05)^{16.7}=0.575$。用 $M=10^4$ 条路径估计它，标准误 $\sqrt{0.575\times0.425/10^4}=0.0049$——这就是报告里必须附上的"概率的误差棒"。敏感性：若 $p_{10}$ 由 0.6 估成 0.5，$\pi_D=0.3/0.8=0.375$，比 1/3 高 12.5%；所以原文强调"转移概率的不确定性要一起传播（Beta 后验抽样）"，否则风险被系统性低估。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% GBM：EM 与精确解对比 + MC 区间（4.4）
mu = 0.10; sig = 0.30; S0 = 1; T = 1; dt = 0.01;
M = 2000; rng(7);                       % 固定种子，结果可复现
Nt = round(T/dt); S = zeros(M,1); Sx = zeros(M,1);
for m = 1:M
    z = randn(Nt,1);
    S(m)  = S0*prod(1 + mu*dt + sig*sqrt(dt)*z);            % EM 路径
    Sx(m) = S0*exp((mu-sig^2/2)*T + sig*sqrt(dt)*sum(z));   % 精确解（格点上）
end
fprintf('E[S_T]: EM=%.4f exact=%.4f theory=%.4f\n', ...
        mean(S), mean(Sx), S0*exp(mu*T));   % 理论值 e^{mu T}
ci = prctile(S,[2.5 97.5]);
fprintf('95%% interval=[%.3f,%.3f], median=%.3f\n', ci, median(S));
% 收敛检验：dt 减半重跑，|mean(S)-theory| 应按 dt 缩小（弱阶 1）
```

输入是参数与步长，输出是 $\mathbb E[S_T]$、区间与中位数。**怎么判断对**：`exact` 与 `theory` 应吻合到 $10^{-2}$ 量级（差的是 MC 标准误，本例约 $0.008$）；EM 的均值应接近 `theory` 但略有偏差，且 $\Delta t$ 减半后偏差减半；若 EM 出现负值，说明 $\sigma\sqrt{\Delta t}$ 太大，应换 log-Euler。

#### 14. 真实美赛 A 题应用

- **2023 干旱**：降雨/干旱年是典型的随机过程，题目要求"群落存续概率"，必须用 Markov 链或泊松降雨驱动，并给概率区间。
- **2020 渔业**：recruitment（补充量）受温度年际波动影响，把"平均 $r,K$"换成"带噪声的 $r_t,K_t$"才能解释丰年与荒年。
- **2026 电池**：使用情景（亮屏时长、温度、应用负载）本身随机，续航结论要给分布而非单值。

#### 15. 常见错误

1. 把"随机（aleatory）"与"无知（epistemic）"混为一谈——前者不能靠更多数据消除，后者可以；混了就会把参数区间错报成预测区间。
2. 只跑一条路径就下结论：单条轨迹的"看起来像趋势"可能只是噪声。
3. 不报告样本量与标准误：读者无法判断 0.57 是真概率还是 $M=20$ 的抖动。
4. 用 EM 模拟必须为正的量且不截断，出现负种群/负价格还继续算。
5. 把 $P(T_1>t)=e^{-\lambda t}$ 与"5 年内至少一次"写成 $1-e^{-5\lambda}$ 之外的式子（漏乘 $t$）。
6. 重要性抽样用了比目标更轻的尾部，权重爆炸还当成"方差缩减"。
7. 拿同一随机种子跑所有情景，却把"差异很小"解读为"模型不敏感"——那可能只是公共噪声。
8. 用 $M=100$ 估 $10^{-3}$ 概率，然后宣称区间"很窄"。

#### 16. 学习要求（验收）

1. 能不看笔记推出泊松分布与指数间隔的**双向**关系，并说明无记忆性。
2. 能写出 M/M/1 的平衡方程、解出 $\pi_n$ 并解释 $\rho<1$ 的必要性。
3. 能用 Itô 公式推出 GBM 显式解，并解释均值与中位数为何不同。
4. 能说明 EM 的强阶 $1/2$、弱阶 $1$，并设计一个数值实验验证斜率。
5. 能对给定稀有事件写出重要性抽样的估计量并说明其失效条件。

### 4.5 元胞自动机（CA）

**定位**：P1｜难度 ★★★｜4 h｜前置：3.4、编程概念｜空间蔓延、局部交互。

#### 1. 定义

元胞自动机（cellular automaton, CA）是四元组 $(\mathcal L,\mathcal S,\mathcal N,f)$：格点集 $\mathcal L$（有限或无限）、有限状态集 $\mathcal S$、邻域模板 $\mathcal N(i)\subset\mathcal L$、局部更新规则 $f:\mathcal S^{|\mathcal N|+1}\to\mathcal S$。所有格点**同步**更新：
$$s_i^{t+1}=f\big(s_i^{t},\{s_j^{t}\}_{j\in\mathcal N(i)}\big),\qquad i\in\mathcal L.$$
若 $f$ 带随机性（以概率 $p$ 取某状态），称为概率元胞自动机（probabilistic CA）。已知量是网格、邻域、规则与初值，未知量是 $t$ 步后的格局 $\{s_i^t\}$ 及其统计量（占据率、斑块大小分布、前沿位置）。CA 不是"一个方程"，而是一台**把局部规则迭代成全局格局**的机器。

#### 2. 核心思想与大白话

PDE 假设空间连续、密度可微，用微分算子描述"邻居之间的净流入"；可当机制本身是**离散的**——一棵树要么活着要么烧了、一个格点要么被占据要么空着、传播需要"邻居有火"这个逻辑条件——微分算子就使不上劲。CA 的做法是把"邻居影响"写成一个查表规则，让时间一层层迭代。

它替代的笨办法是：为了表达阈值/离散状态，硬把 PDE 的分段项写成示性函数，结果方程不可微、数值解在界面处抖动。CA 则天然处理"有/无"和"达到阈值才发生"。

**大白话：** PDE 像给每个地点都记一本账，还要记邻居之间转了多少账；CA 像给每个格子发一张"下一步该变成什么"的卡片，卡片只看自己和隔壁。每家只看邻家，整片小区最后却自己长出了斑块、边界和蔓延前沿——这就是"局部规则产生全局格局"。

关键取舍：CA 赢在机制表达力，输在可解析性与网格伪影；它能逼近反应—扩散方程，但**只在特定离散化条件下**（第 6 小节给出条件），不能随口说"CA 就是 PDE"。

#### 3. 数学基础

1. **格图与邻域模板**：$\mathcal L$ 是 $\mathbb Z^2$ 的子集，von Neumann 邻域 $|N|=4$、Moore 邻域 $|N|=8$；模板决定离散拉普拉斯的各向异性。
2. **同步映射**：$F:\mathcal S^{|\mathcal L|}\to\mathcal S^{|\mathcal L|}$，$s^{t+1}=F(s^t)$；CA 的一切长期行为都是这个高维映射的轨道。
3. **离散拉普拉斯与 Taylor 粗粒化**：$\frac1{2d}\sum_{j\in N(i)}(s_j-s_i)$ 一阶逼近 $\frac{a^2}{2d}\nabla^2s$，是把 CA 与 PDE 接起来的桥。
4. **渗流（percolation）与临界指数**：占据率超过阈值 $p_c$ 才出现跨越全场的连通簇；二维方格点渗流 $p_c\approx0.5927$，关联长度 $\xi\propto|p-p_c|^{-\nu}$，$\nu=4/3$。
5. **平均场近似**：把邻居状态替换成全局均值，得到 ODE；它给出定性趋势，但会抹掉聚集与有限尺寸效应。
6. **有限尺寸标度**：$L\times L$ 网格把临界点附近的转变抹宽约 $L^{-1/\nu}$，这是"换分辨率结果就变"的定量解释。
7. **遍历性与吸收态**：火灾、感染类 CA 有吸收态（全空/全烧），长期行为要么收敛到吸收态，要么进入统计稳态。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| 空间可离散为规则格 | 现象尺度远大于格距 $a$ | 分辨率足够、机制尺度明确 | 前沿速度与斑块形状依赖格子朝向，误差可达 20% | 把 $a$ 减半，比较前沿速度与占据率 |
| 状态集有限且局部足够 | 每个格点只需自身+邻域状态 | 机制是局部的 | 遗漏长程作用（风、水系）使蔓延被低估 | 加长程耦合项看结果是否显著变化 |
| 规则时不变 | $f$ 在模拟期内不随季节/政策变 | 短时段、稳定环境 | 参数漂移被当成噪声，长期预测偏移 | 分段标定，检查参数漂移 |
| 同步更新 | 同一时刻全体更新 | 世代离散（植物、年际） | 出现棋盘振荡等非物理伪影 | 改成随机异步更新，比较稳态 |
| 邻域各向同性 | 4/8 邻域不偏向某方向 | 方格近似 | 对角方向传播快 $\sqrt2$ 倍，锋面变方 | 旋转初始条件，比较各方向速度 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $a$ | 格距（一个格点的物理边长） | m | $0.1\sim10$ | 越小越接近连续极限，计算量按 $a^{-2}$ 涨 |
| $\tau$ | 时间步长（一代） | s（生态常用 yr） | $1$ d $\sim1$ yr | 与 $a$ 必须满足扩散缩放 $\tau\propto a^2$ |
| $L$ | 网格边长（格点数） | —（格） | $50\sim500$ | 决定有限尺寸抹宽 $L^{-1/\nu}$ |
| $p_g$ | 生长/恢复概率 | 每格每步 | $10^{-3}\sim10^{-1}$ | 越大恢复越快、占据率越高 |
| $f$ | 雷击/自发点燃概率 | 每格每步 | $10^{-5}\sim10^{-3}$ | 越大扰动越频繁、稳态占据率越低 |
| $p_s$ | 蔓延概率（概率 CA） | — | $0.2\sim1$ | 小于 1 时前沿速度按 $p_s$ 折减 |
| $\phi$ | 树/占据格点占比 | — | $0\sim1$ | 宏观序参量，$0.751$ 即一例 |
| $D$ | 有效扩散系数 | $\mathrm{m^2/s}$ | $a^2/(2d\tau)$ | 由格距与步长诱导，不是独立参数 |
| $n_{\text{seed}}$ | 随机种子数 | —（个） | $\ge30$ | 决定统计量的误差棒 |

#### 6. 从基本原理推导核心公式

**（1）更新规则与粗粒化：CA 何时收敛到反应—扩散方程。** 设格点 $i$ 位于 $x$，格距 $a$、步长 $\tau$，状态 $s$ 可视为光滑密度。取一类"平均+反应"的规则：
$$s_i^{t+1}=s_i^{t}+\tau R(s_i^{t})+\frac{2d\,\tau D}{a^{2}}\Big[\frac1{2d}\sum_{j\in N(i)}s_j^{t}-s_i^{t}\Big],$$
其中 $2d$ 是邻居数（正方格 4 邻域 $d=2$），这个因子必须留着，否则粗粒化后的扩散系数会差 $2d$ 倍。对邻居做 Taylor 展开：
$$\frac1{2d}\sum_{j\in N(i)}s_j=s_i+\frac{a^2}{2d}\nabla^2s_i+O(a^4),$$
代回并两边除以 $\tau$：
$$\frac{s_i^{t+1}-s_i^{t}}{\tau}=D\nabla^2s_i+R(s_i)+O(Da^2)+O(\tau).$$
$O(Da^2)$ 来自刚才丢掉的 $O(a^4)$ 邻居项：它被系数 $2d\tau D/a^2$ 放大、再除以 $\tau$ 得到。在扩散缩放 $D=\kappa a^2/\tau$（即 $\tau\propto a^2$）下，$Da^2=\kappa a^4/\tau=O(a^2)$，故截断误差为 $O(a^2)+O(\tau)$。
要让 $a,\tau\to0$ 时极限存在且扩散项不被吃掉，必须保持
$$\frac{D\tau}{a^2}=\kappa=\text{常数},\qquad D=\frac{\kappa a^2}{\tau},$$
即**扩散缩放**（diffusive scaling）$\tau\propto a^2$。此时极限就是
$$\frac{\partial s}{\partial t}=D\nabla^2s+R(s).$$
三个必须同时成立的条件：① 规则的平均漂移与 $R(s)$ 一致（一致性）；② $\tau/a^2$ 固定（否则 $\tau\to0$ 快于 $a^2$ 会让 $D\to0$，只剩 ODE，空间耦合消失）；③ 格点对称性足够（4 邻域的截断误差含 $\cos4\theta$ 项，见下文各向异性）。显式格式的稳定性还要求
$$\frac{D\tau}{a^2}\le\frac1{2d}=\frac14\quad(d=2),$$
正好对应"一步之内信息最多传到一格"，超过就会出现非物理振荡。

**（2）前沿速度：把格子单位换算成物理单位。** 取 $R(s)=rs(1-s/K)$（Fisher–KPP），连续极限的行波速度
$$c=2\sqrt{Dr}.$$
在格子上数出来的是"每步走几格"：$c_{\text{cell}}=c\,\tau/a=2\tau\sqrt{Dr}/a$。代入 $D=\kappa a^2/\tau$：
$$c_{\text{cell}}=2\sqrt{\kappa\,r\tau}.$$
量纲检查：$r\tau$ 是"每步的增长概率"，无量纲，故 $c_{\text{cell}}$ 无量纲 ✓。极限检查：$r\to0$ 时前沿不动；$\kappa=1/4$（稳定性上限）时 $c_{\text{cell}}=\sqrt{r\tau}$，即每步约走 $\sqrt{r\tau}$ 格。这条式子给出一个强检验：改 $a$ 而保持物理 $D,r$ 不变时，$c_{\text{cell}}$ 应按 $\sqrt{\tau}$ 变，乘 $a/\tau$ 后回到同一个 $c$；若对不上，说明离散化没进入渐近区。

**（3）渗流阈值：为什么占据率 0.751 意味着大火。** 方格子点渗流的临界占据率 $p_c\approx0.5927$：$\phi<p_c$ 时树簇是指数小的孤岛，火最多烧掉一个小簇；$\phi>p_c$ 时存在跨越全场的连通簇，一次雷击就可能烧掉大片。平均场（Bethe 格，配位数 $z=4$）给 $p_c=1/(z-1)=1/3$，明显低于真实二维值——因为平均场忽略了"绕路"和"闭环"，这也解释了为什么用 ODE 平均场估计火灾规模会系统性偏大。$0.751>0.5927$，所以该稳态下树是连通的，火灾规模分布会出现重尾；这正是"不能只看平均树占比"的原因。注意连通只是**必要**条件：蔓延概率 $p_s$ 若低于键渗流阈值（方格 4 邻域约 $0.5$），火仍会被掐断在小簇里，见第 12 小节（b）。

**（4）平均场平衡与稳态占据率。** 忽略空间关联，设树占比 $\phi$，每步生长 $p_g(1-\phi)$、雷击损失 $f\phi$、蔓延损失 $-\beta\phi^2$（火与树相遇的概率正比于两者密度），则
$$\dot\phi=p_g(1-\phi)-f\phi-\beta\phi^2.$$
若蔓延项可忽略（$\beta\to0$），稳态 $\phi^*=p_g/(p_g+f)$：$f\ll p_g$ 时 $\phi^*\to1$，$f\gg p_g$ 时 $\phi^*\to0$，两端极限都合理。反过来，若实测 $\phi^*=0.751$ 且蔓延可忽略，则 $f/p_g\approx0.33$；实际有蔓延时 $\beta>0$，同样的 $\phi^*$ 需要更小的 $f/p_g$——这就是"平均场反推参数会偏高"的方向。$\beta$ 不能由 $\phi^*$ 单独定出，必须再匹配火灾规模分布，这是参数不可辨识的典型例子。

**（5）有限尺寸标度：换分辨率为什么结果会变。** 临界点附近关联长度 $\xi\propto|\phi-p_c|^{-\nu}$（二维 $\nu=4/3$）。当 $\xi$ 超过网格尺寸 $L$ 时，系统"看不见"临界性，转变被抹宽：
$$\Delta\phi\sim L^{-1/\nu}=L^{-3/4}.$$
$L=120$ 时 $\Delta\phi\approx120^{-0.75}=0.0276$；$L=240$ 时约 $0.0164$，抹宽缩到 $0.59$ 倍。所以"$L=120$ 的稳态占比"与"$L=240$"的差若在 $0.03$ 量级，是有限尺寸效应而不是模型错了；若差出 $0.1$，才说明规则或参数有问题。这条式子把"格点各向异性检验"从口号变成了可算的容差。

**（6）同步与异步：一个能手算的反例。** 两个相邻格点，规则都是"下一时刻取邻居的值"。同步更新：$(0,1)\to(1,0)\to(0,1)\to\cdots$，永远周期 2，永远不收敛；随机异步更新（每步只随机挑一个格点更新）：第一步就变成 $(1,1)$ 或 $(0,0)$，之后不动。同样的局部规则，只因调度不同，长期行为从"振荡"变成"固定点"。这就是第 4 小节把"同步更新"列为假设、第 15 小节把"同步/异步混用"列为错误的原因。

**（7）各向异性：4 邻域到底偏了多少。** 离散拉普拉斯的截断误差含 $\cos4\theta$ 项，比例约 $(a/\ell)^2$，$\ell=\sqrt{D/r}$ 是前沿宽度。若 $\ell\approx10a$，方向间速度差约 $1\%$；若 $\ell\approx2a$（陡前沿、粗网格），差可达 $20\%$ 以上，锋面会呈现明显的方形/菱形。检验方法：同一初始条件旋转 $45^\circ$ 再跑，比较两个方向的前沿位置；差异超过 $5\%$ 就应加密网格或改用 9 点模板。

#### 7. 求解方法

1. **双缓冲同步**：先算新状态数组，再整体替换；绝不能在遍历时原地改（否则变成异步，且顺序依赖）。
2. **向量化**：用 `conv2` 做邻域计数（Moore 用 `ones(3)`，von Neumann 用十字核），比双重循环快 1—2 个数量级；边界用 `'same'` 加零填充，等价于"外界无树"。
3. **重复与统计**：每个参数点跑 $n_{\text{seed}}\ge30$ 个种子，报告均值与标准差；只跑一次得到的 $0.751$ 没有误差棒。
4. **分辨率检验**：$L=60,120,240$ 各跑一遍，看稳态占比与前沿速度是否按 $O(a^2)$ 或 $L^{-3/4}$ 收敛；不收敛就先怀疑规则依赖格距。
5. **时间步检验**：概率 CA 的每步概率 $p$ 应与 $\tau$ 成正比（$p=r\tau$）；若改 $\tau$ 而固定 $p$，物理时间就被偷偷改了。
6. **MATLAB 工具**：`conv2`（邻域求和）、`rand`（概率判定）、`imshow`/`imagesc`（看格局）、`regionprops`（斑块统计）、`rng`（种子）。

#### 8. 适用场景

**该用：** 机制含离散状态与阈值（活/死、占/空、达标/未达标）；传播依赖"邻居是否有某状态"（火、病、菌丝）；个体位置与空间格局本身是结论（斑块、边界、通道）；需要加异质格点（地形、肥力）而 PDE 参数化困难。

**不该用：** 关注平滑梯度与守恒通量（用 PDE，3.5）；尺度远大于格距且场光滑（CA 只是更贵的差分格式）；需要解析解或严格的参数置信区间；格点朝向会污染结论且无法加密（此时换不规则网格或 ABM）。

#### 9. 局限

1. **格点各向异性**：4 邻域在陡前沿下方向速度差可达 20%，锋面形状被格子"带偏"。
2. **参数多且难标定**：$p_g,f,p_s$ 加上邻域、更新顺序、初值，至少 5—6 个自由度，而可用宏观统计量往往只有 2—3 个，参数不可辨识。
3. **有限尺寸效应**：$L=120$ 的转变抹宽约 $0.028$，临界附近结论对 $L$ 敏感。
4. **同步伪影**：同步更新可产生周期 2 振荡或棋盘格局，被误读成"真实的波动"。
5. **无解析解**：拿不到闭式解与严格的敏感性导数，只能靠重复模拟与统计。
6. **计算随分辨率增长很快**：$a$ 减半，格点数 $\times4$、步数 $\times4$（扩散缩放），总代价 $\times16$，即 $\propto a^{-(d+2)}$（二维 $a^{-4}$）。

#### 10. 改进

1. **换邻域/格子**（最便宜）：9 点模板或六边形格降低各向异性，代价是实现复杂一点。
2. **改调度**：随机异步或 Gillespie 型连续时间更新，消除同步伪影，代价是每步只能改一个格点、更慢。
3. **连续状态 CA**：格点存连续密度而非离散状态，介于 CA 与 PDE 之间，可用差分格式加速。
4. **自适应/多尺度网格**：只在界面附近加密，代价是数据结构复杂。
5. **CA + 参数校准**：用近似贝叶斯（ABC）或模拟矩估计，把 $p_g,f,p_s$ 的后验与预测区间一起给出（接 4.3）。
6. **代理模型**：先用 CA 生成训练集，再拟合高斯过程/神经网络代理，把敏感性分析与优化搬到代理上（接 4.7）。

#### 11. 与其他模型的组合

- **CA + PDE（3.5）**：CA 提供离散机制与阈值，PDE 提供粗粒化极限与解析前沿速度；两者互检——CA 数出的前沿速度应接近 $2\sqrt{Dr}$，对不上说明格距太大。
- **CA + ABM（4.6）**：CA 是"格点状态"的特例，ABM 把格点换成有属性的个体；当个体差异重要时从 CA 升级到 ABM。
- **CA + 渗流理论**：用 $p_c$ 与 $\xi$ 预测"何时从局部小火变成全场大火"，比盲跑参数更省。
- **CA + 参数反演（4.3）**：以斑块统计量（占据率、斑块数、边界长度）为摘要统计量做校准。

#### 12. 手算算例

**（a）森林火灾基线（原文数字）。** $120\times120$ 格、60 步后树占比约 $0.751$。格点总数 $120^2=14400$，对应树格约 $0.751\times14400=10814$ 个。若 $a=1$ m、$\tau=1$ d，则模拟域 $120$ m $\times120$ m、物理时间 60 d；这个换算必须在报告里写清，否则"60 步"没有物理含义。复核提示：用第 13 小节的规则与 $p_g=0.02$、$p_f=10^{-4}$、$p_s=0.3$ 跑，$t=60$ 一般得 $0.65$—$0.68$，$t>300$ 后稳定在 $0.70$ 上下；$0.751$ 与之同量级（单次实现加有限尺寸涨落），但精确复现需要原文未列出的完整规则，故该数字标为待复核、此处不改动。

**（b）与渗流阈值对照。** 树的**静态**连通性由占据率决定：$0.751>p_c=0.5927$，几何上存在跨越全场的树簇。但"连通"不等于"一定烧成大火"：火能不能沿这条通路走过去，还取决于蔓延概率 $p_s$ 是否越过键渗流的临界值（方格 4 邻域 $p_c^{\text{bond}}=0.5$，Moore 8 邻域会把它压低）。$p_s$ 低于阈值时火被掐断在小簇里，$0.751$ 这种高占比才能维持；$p_s$ 接近 1 时连通就意味着一次雷击能烧掉全场，稳态树占比会显著下降。原文只写了占据率阈值、没写这道"蔓延闸门"，这一点待补。若某组参数给出 $\phi^*=0.45<p_c$，火只能烧小簇，火灾规模分布会从重尾变成指数尾——这是判断"模型处于哪个 regime"的快速检查。

**（c）平均场反推。** 忽略蔓延，$\phi^*=p_g/(p_g+f)=0.751$ 给 $f/p_g=0.331$；若 $p_g=0.02\ \mathrm{d^{-1}}$，则 $f\approx0.0066\ \mathrm{d^{-1}}$。但 60 步内树从初值长到 0.751，时间尺度 $\sim1/p_g=50$ d，与 60 步同量级，自洽。注意：这是"若忽略蔓延"的反推，真实 $f$ 应更小（蔓延也消耗树），不能把它当成参数估计。

**（d）分辨率容差。** $\Delta\phi\sim L^{-3/4}$：$L=120$ 给 $0.0276$、$L=240$ 给 $0.0164$。所以两次模拟的稳态占比差在 $0.03$ 以内属正常；若 $L=120$ 与 $L=240$ 差 $0.10$，先查规则是否隐含依赖 $a$（例如"每步固定烧一层"在粗网格上物理上更快）。

**（e）前沿速度量纲检查。** 取 $r=0.02\ \mathrm{d^{-1}}$、$D=0.1\ \mathrm{m^2/d}$（此时 $\kappa=D\tau/a^2=0.1\le1/4$，满足第（1）条的稳定性），则 $c=2\sqrt{0.1\times0.02}=0.089$ m/d；60 d 走约 $5.4$ m，即 $5.4$ 格。若代码数出的前沿 60 d 走了 60 格，说明每步无条件烧一层（$c_{\text{cell}}=1$），物理上等价于无限扩散，必须改规则或缩步长。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% 森林火灾 CA：生长—雷击—Moore 蔓延（4.5）
L = 120; pg = 0.02; pf = 1e-4; ps = 0.3; T = 60; rng(3);
tree = 1; fire = 2; empty = 0;
S = tree*ones(L); S(rand(L) < 0.2) = empty;   % 初值：八成树
for t = 1:T
    S(S==empty & rand(L) < pg) = tree;         % 生长（先判定后写，保持同步）
    S(S==tree  & rand(L) < pf) = fire;         % 雷击
    oldFire = (S==fire);
    S(oldFire) = empty;                        % 上一代火熄灭
    nb = conv2(double(oldFire), ones(3), 'same') > 0;   % Moore 邻域
    S((S==tree) & nb & (rand(L) < ps)) = fire; % 蔓延受 p_s 限制
end
fprintf('tree fraction = %.3f\n', mean(S(:)==tree));
imagesc(S); axis image; colormap([1 1 1; 0.2 0.7 0.2; 1 0.3 0]);
```

输入是 $L,p_g,p_f,p_s,T$ 与种子；输出是稳态树占比与格局图。**怎么判断对**：$p_g=0$ 时树只减不增、占比单调降到 0；$p_f=0$ 时树单调涨到 1；两个极限都对，才说明生长/熄灭的符号没写反。再把 $p_f$ 与 $p_s$ 两个极端各跑一遍：$p_s=1$、$p_f=5\times10^{-4}$ 是"火势活跃"regime（占比掉到 0.3 以下，格局碎成小斑块），$p_s=0.3$ 是本例的"火被掐断"regime（占比维持在 0.6 以上，能看见连片大树簇）。最后把 $L$ 改成 60 与 240，看占比是否落在 $\pm0.03$ 的有限尺寸容差内；换 5 个种子，看标准差是否远小于均值。

#### 14. 真实美赛 A 题应用

- **2021 真菌**：菌丝在枯落物上的扩展是"局部接触 + 资源阈值"的过程，CA 的格点状态（定殖/未定殖）与 Moore 邻域正好对应，参数可用分解速率与扩展速度标定。
- **2023 植物群落**：干旱年的空间斑块与局部竞争适合用 CA 表达"邻居促进/竞争 + 随机扰动"，输出斑块大小分布而不只是平均盖度。

#### 15. 常见错误

1. 同步与异步混用：代码里一边遍历一边改状态，结果依赖遍历顺序，却宣称"同步 CA"。
2. 不报告种子与分辨率：$0.751$ 没有误差棒，读者无法判断它是否稳定。
3. 把 CA 当精确预测：CA 给的是机制下的格局分布，不是"明年这块地一定烧"。
4. 改 $a$ 却不同步改 $\tau$（或反之），破坏扩散缩放，物理极限直接跑偏。
5. 用平均场公式反推参数，却忽略蔓延项，$f$ 被系统性高估。
6. 概率参数不随步长缩放：$p=0.02$ 每步，改 $\tau$ 后物理速率变了三倍还不知道。
7. 边界条件不说明：零填充、周期、反射给出的前沿速度不同。
8. 用单次运行的"漂亮图案"当验证，而不做极限检查（$p_g=0$、$p_f=0$）。

#### 16. 学习要求（验收）

1. 能写出 CA 四元组并用粗粒化推出 $\partial_t s=D\nabla^2s+R(s)$，说清 $\tau\propto a^2$ 与稳定性上限 $\kappa\le1/4$。
2. 能用 $c=2\sqrt{Dr}$ 把格子单位的前沿速度换算成物理速度，并设计分辨率检验。
3. 能解释 $p_c\approx0.5927$ 与 $\phi=0.751$ 的关系，说明为什么平均场会高估火灾规模。
4. 能给出同步/异步的反例，并说明调度选择如何写进模型假设。
5. 能报告"种子数、网格尺寸、边界条件"三件套并给出有限尺寸容差。

### 4.6 主体建模（ABM）

**定位**：P1｜难度 ★★★★｜5 h｜前置：概率、CA｜异质个体 + 局部交互 + 政策仿真。

#### 1. 定义

基于主体的模型（agent-based model, ABM）是五元组 $(\mathcal A,\mathcal E,\{\pi_i\},\sigma,\mathcal I)$：主体集 $\mathcal A=\{1,\dots,N\}$，每个主体有状态 $x_i\in\mathcal X_i$（位置、年龄、库存、信念）；环境 $\mathcal E$ 存资源与空间；决策规则 $\pi_i$ 把"局部观测"映成动作；调度器 $\sigma$ 规定谁在何时更新；交互结构 $\mathcal I$ 规定谁能看到谁（邻域、网络、市场）。一次运行是从初值到终点的**一条随机轨迹**，$R$ 次独立运行给出输出的**分布**：
$$\{x_i(t+1)\}_{i=1}^{N}=\Phi\big(\{x_i(t)\},\mathcal E(t),\xi(t)\big),\qquad \xi(t)\sim\text{随机数流}.$$
已知量是主体数 $N$、规则 $\pi_i$、调度 $\sigma$、交互图 $\mathcal I$ 与种子；未知量是宏观输出（总量、分布、峰值、空间格局）及其不确定性。

#### 2. 核心思想与大白话

ODE 把一屋子人压缩成一个"平均人"：所有个体同参数、充分混合、只关心总量。可很多 A 题结论恰恰藏在**差异**里——谁先跑、谁走得慢、谁结伴、谁靠右——这些差异经过局部交互与反馈，会放大成宏观格局（拥堵、斑块、两极分化）。ABM 就是保留每个个体的状态与邻居，让宏观现象**涌现**（emergence）出来，而不是事先写进方程。

它替代的笨办法是：为了解释"为什么磨损中间偏右更重"，在 ODE 里加一个手工的对流项；ABM 则让"靠右行走"这条个体规则自己生成侧向分布。

**大白话：** ODE 像统计一个班的平均身高，然后假设全班都是这个身高；ABM 给每个同学量身高、排座位、规定谁跟谁说话，最后看班级队形怎么变。平均身高解释不了"为什么前排矮后排高"，ABM 可以。

关键分寸：ABM 的产出是**机制下的分布**，不是精确预测；规则越多越像"用代码编故事"，所以第 3 小节那条建模纪律（最小规则、校准、验证）比代码本身重要。

#### 3. 数学基础

1. **乘积空间上的马尔可夫链**：全体状态 $(x_1,\dots,x_N)$ 构成高维状态，一次更新是转移核 $K$；ABM 的长期行为由这个链的平稳分布决定。
2. **主方程（master equation）**：$p_{t+1}(y)=\sum_x p_t(x)K(x,y)$，是"分布如何演化"的精确描述；ODE 只跟踪其一阶矩。
3. **平均场极限与混沌传播**：$N\to\infty$、交互全局混合时，经验分布收敛到确定性 ODE；有限 $N$ 的涨落为 $O(N^{-1/2})$。
4. **网络与图**：交互结构用图 $G=(V,E)$ 表示，度分布、聚类系数、路径长度决定传播速度与公平性。
5. **方差分解**：$\mathrm{Var}(Y)=\mathbb E[\mathrm{Var}(Y\mid \text{seed})]+\mathrm{Var}(\mathbb E[Y\mid\text{seed}])$，区分"运行内噪声"与"运行间差异"。
6. **伪随机数生成**：Mersenne Twister 周期 $2^{19937}-1$，种子决定整条流；公共随机数（common random numbers, CRN）用于配对比较。
7. **近似贝叶斯计算（ABC）与模拟矩估计**：ABM 似然写不出，用摘要统计量的距离做校准。
8. **模式导向建模（pattern-oriented modeling）**：用多个尺度上的可观察模式同时约束模型，而不是只拟合一条曲线。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| 规则抓到了主导机制 | 每条规则都有行为/物理依据 | 文献、观测、访谈 | 结论是"代码的结论"而非"世界的结论"，外推全错 | 删掉一条规则看输出是否变；无依据的规则不得保留 |
| 主体异质性可参数化 | 个体差异用分布（速度、阈值）描述 | 有样本或合理假设 | 分布选错使宏观分布的峰位与宽度都偏，结论从"单峰"变"双峰"或反之 | 用观测分布做拟合检验；对分布做敏感性 |
| 交互结构正确 | 邻域/网络/市场范围写对 | 空间数据、社交网络 | 局部聚类被抹平，传播被高估 | 比较全局混合与局部交互两版结果 |
| 调度不影响定性结论 | 同步/异步、更新顺序不改变结论 | 规则无强次序依赖 | 出现周期振荡或顺序偏差，结论不可复现 | 换 3 种调度重跑，比较分布 |
| 随机数独立且可复现 | 种子固定、流不重叠 | PRNG 质量 | 结果不可复现；并行流相关导致误差棒虚假变窄 | 固定种子重跑；换种子看分布；检查流偏移 |
| 时间步足够小 | 一步内至多一个事件 | 事件率 $\times\Delta t\ll1$ | 并发事件被吞，动力学被人为减慢 | $\Delta t$ 减半，比较轨迹与总量 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $N$ | 主体数 | —（个） | $10^2\sim10^6$ | 越大越接近平均场，涨落按 $1/\sqrt N$ 缩小 |
| $L$ | 空间尺度（楼梯宽、地块边长） | m | $1\sim10^3$ | 与交互半径 $r$ 的比决定"局部"还是"全局" |
| $v_i$ | 个体速度 | m/s | $0.5\sim1.5$ | 越大越容易超越，拥堵模式改变 |
| $r$ | 交互/感知半径 | m | $0.5\sim10$ | 越大越接近全局混合，空间格局越弱 |
| $\beta$ | 单次接触的传播/采纳概率 | — | $0.01\sim0.5$ | 越大级联越快，阈值行为越强 |
| $\theta_i$ | 个体采纳阈值 | —（比例） | $0\sim1$ | 分布越宽，宏观结果越可能双峰 |
| $\Delta t$ | 时间步长 | s | $0.1\sim10$ | 需满足"事件率 $\times\Delta t\ll1$" |
| $R$ | 独立运行（种子）数 | —（次） | $\ge30$ | 决定输出分布与误差棒的精度 |
| $T$ | 模拟时长 | s（或 d、yr） | 与过程同量级 | 太短只看到瞬态，太长白烧算力 |

#### 6. 从基本原理推导核心公式

**（1）从个体规则到平均场 ODE。** 设 $N$ 个主体充分混合，每个感染者单位时间随机接触 $c$ 个主体，每次接触传播概率 $\beta$，恢复率 $\gamma$。在 $[t,t+\Delta t)$ 内，期望新增感染为
$$\mathbb E[\Delta I\mid S,I]=S\cdot(c\Delta t)\cdot\beta\cdot\frac{I}{N}+o(\Delta t).$$
两边除以 $N\Delta t$，令 $\Delta t\to0$，记 $s=S/N$、$i=I/N$：
$$\frac{ds}{dt}=-\beta c\,s\,i,\qquad \frac{di}{dt}=\beta c\,s\,i-\gamma i,\qquad R_0=\frac{\beta c}{\gamma}.$$
量纲：$\beta$ 无量纲、$c$ 与 $\gamma$ 的单位是 $\mathrm{s^{-1}}$，故 $R_0$ 无量纲 ✓。极限检查：$i=0$ 或 $s=0$ 时无新增；$R_0<1$ 时 $i$ 单调衰减；$N\to\infty$ 且初始比例固定时，大数定律使 ODE 成为精确极限，有限 $N$ 的偏差为 $O(N^{-1/2})$。注意这一步偷偷用了"每次接触的对象在人群中均匀随机"——这正是 ABM 要打破的假设。

**（2）局部交互为什么让传播变慢。** 在格子/网络上，每个主体只有 $n_h$ 个邻居。设 $\chi=\mathbb P(\text{邻居被感染}\mid\text{主体被感染})$ 为配对感染率。由于聚集，$\chi$ 大于全局感染率 $i$；写成一阶形式 $\chi=g\,i$，$g\ge1$ 是配对相关因子（早期同一簇内 $g$ 可接近 $n_h$）。则有效传播率为 $\beta(1-\chi)$ 而非 $\beta(1-i)$，增长率之比
$$\frac{\beta(1-g i)}{\beta(1-i)}\approx1-(g-1)i.$$
所以局部交互让 ODE 的感染力 $\beta(1-i)$ 偏高：ABM 里真正还能被感染的邻居比例是 $1-\chi=1-gi$，比"全局易感比例" $1-i$ 更小，差距随 $i$ 与 $g$ 一起增大。结论：**同样的 $\beta,\gamma$，ODE 会高估增长率、最终规模与峰值**，且高估幅度随聚集度单调增大；当 $gi$ 接近 1 时有效易感比例已压到 0 附近，传播实际停滞，而 ODE 仍按 $1-i$ 继续推——这是高估的极端情形。这就是"ABM 与 ODE 给出不同结论"的第一个可解释来源，不是随机噪声。

**（3）异质阈值为什么产生双峰分布（ODE 做不到）。** 设每个主体有采纳阈值 $\theta_i\sim F$，当"已采纳邻居比例"$\ge\theta_i$ 时采纳。平均场下，若全局采纳比例为 $\phi$，则下一步
$$\phi_{t+1}=F(\phi_t)=:T(\phi_t).$$
不动点满足 $\phi^*=F(\phi^*)$；稳定性由 $|F'(\phi^*)|<1$ 判定。取光滑阶跃 $F(\phi)=3\phi^2-2\phi^3$（把 $[0,1]$ 映到 $[0,1]$，$F(0)=0,F(1)=1$），则
$$F(\phi)-\phi=\phi(1-\phi)(2\phi-1),$$
三个不动点 $0,\ \tfrac12,\ 1$；$F'(\phi)=6\phi(1-\phi)$，在 $\phi=\tfrac12$ 处 $F'=1.5>1$——中间点**不稳定**，两端稳定。这正是"临界质量"：初始采纳率略高于 1/2 就滚到 1，略低就退回 0。

ODE 从 $\phi_0=1/2$ 出发会永远停在 1/2（不稳定平衡，数值上被舍入误差推走，但推去哪一侧是任意的）。ABM 不同：有限 $N$ 带来人口噪声，一步的涨落约 $\sqrt{\phi(1-\phi)/N}$，它把系统推离 1/2，之后非线性放大，最终落在 0 或 1。于是**跨种子看，输出是双峰分布**：一部分运行全员采纳，一部分运行无人采纳，比例由初值和噪声强度决定。这就是"ABM 能产生 ODE 不能产生的分布"的准确含义——不是 ABM 更准，而是 ODE 只给一个点，丢掉了"结局有两种"这一结构性信息。同理，若输出是随机输入的凸函数，Jensen 不等式给 $\mathbb E[f(X)]\ge f(\mathbb E X)$，用均值代入会系统性偏低。

**（4）随机种子、样本量与公共随机数。** 固定种子只是保证可复现，不保证代表性。$R$ 次独立运行的样本均值标准误为 $s/\sqrt R$；要以 $1-\alpha$ 置信度分辨差异 $\delta$，需
$$R\ge\Big(\frac{z_{\alpha/2}\,s}{\delta}\Big)^2.$$
例：$s=0.1$、$\delta=0.05$、$z=1.96$，得 $R\ge16$；若 $s=0.3$，则 $R\ge139$——"跑 10 次"通常连误差棒都画不出来。

比较两个情景时，用**同一组种子**（CRN）可消掉共同噪声：
$$\mathrm{Var}(Y_1-Y_2)=2\sigma^2(1-\rho)\ \Longrightarrow\ \frac{\mathrm{SD}_{\text{CRN}}}{\mathrm{SD}_{\text{indep}}}=\sqrt{1-\rho}.$$
$\rho=0.9$ 时标准差降到 $0.316$ 倍，达到同样精度所需运行数减少到 $1/10$。但 CRN 只对"同种子下两情景轨迹正相关"有效；若两情景的随机消耗顺序不同，相关性消失，收益归零。

**（5）调度顺序的偏差。** 同步更新可能出现周期振荡：两个相邻主体都"复制对方状态"，同步下 $(0,1)\to(1,0)\to(0,1)\cdots$；随机异步下第一次更新就固定。随机顺序更新等价于每步随机选一个主体，单步顺序偏差 $O(1/N)$，$N$ 步累积 $O(1)$——所以小 $N$ 时调度真的会改变结论，大 $N$ 时影响变小但不会消失。事件驱动的 Gillespie 更新在"事件率恒定、指数等待"时是精确的，固定步长格式则有 $O(\Delta t)$ 偏差。

**（6）量纲与极限总检查。** $N\to\infty$ 且 $r\to L$（全局混合）时 ABM 应退化为 ODE，可作为代码正确性检验；$r\to0$ 时退化为局部 CA；$\Delta t\to0$（事件驱动极限）时轨迹应不随步长变；$R\to\infty$ 时输出分布应收敛。四条极限任一不成立，先查代码而不是查模型。

#### 7. 求解方法

1. **时间驱动**：固定 $\Delta t$，每步按调度更新主体；简单、易向量化，但要求事件率 $\times\Delta t\ll1$。
2. **事件驱动**：给每个主体抽指数等待时间，用优先队列推进；精确、步长无关，实现更复杂（Gillespie 型）。
3. **邻域搜索**：朴素两两距离 $O(N^2)$；用网格分箱（cell list）或 `rangesearch`/kd-tree 降到近 $O(N)$。
4. **集成运行**：$R\ge30$ 个种子，报告均值、标准差、分位数；用 `parfor` 并行时给每个 worker 独立子流（`rng(seed+k)`），避免流重叠。
5. **校准**：选 2—3 个宏观统计量（总量、分布的分位数、峰值位置），用 ABC 或模拟矩估计；不要用全部数据拟合，否则过拟合。
6. **验证**：留出独立年份/情景，做模式检验（空间格局、分布形状、时间相关性），而不只是总量吻合。
7. **MATLAB 工具**：`rng`（种子）、`rand`/`randn`/`randsample`（抽样）、`rangesearch`（邻域）、`histcounts`（分布）、`parfor`（并行）、结构体数组或 `handle` 类（主体）。

#### 8. 适用场景

**该用：** 个体差异与局部决策共同决定宏观结果（拥堵、踩踏、抢购、扩散）；政策/规则仿真（靠右行走、隔离、捕捞配额）需要"如果……会怎样"；空间使用模式本身就是结论（磨损分布、斑块）；交互网络不均匀且网络结构是机制的一部分。

**不该用：** 只关心总量且个体近似同质（ODE 更可验证）；数据只能支撑 2—3 个参数而 ABM 要 10 个；无法做多种子与独立验证（结论不可信）；需要解析灵敏度或严格置信区间；实时性要求高（ABM 通常太慢）。

#### 9. 局限

1. **参数不可辨识**：规则多、数据少，多组参数给同一宏观输出；报告单点估计是虚假的精确。
2. **计算贵**：$N$ 大、$R$ 大、$T$ 长时，一次校准可能上千 CPU 小时；朴素邻域搜索还是 $O(N^2)$。
3. **验证困难**：没有似然、没有闭式解，只能靠模式检验；通过检验不代表模型正确。
4. **过拟合风险**：规则越加越多，总能拟合历史，但外推崩溃。
5. **代码即模型**：实现 bug（更新顺序、边界、抽样）会静默改变结论，且难以 review。
6. **随机噪声掩盖信号**：$R$ 不足时，情景差异被淹没在种子噪声里，被误读为"无影响"。
7. **"讲故事"风险**：每条规则都"看起来合理"，但合起来没有可检验的预测。

#### 10. 改进

1. **加纪律**（零成本）：规则表 + 依据 + 删规则检验；固定种子与版本；报告 $R$ 与误差棒。
2. **换调度/步长**（低成本）：随机异步、事件驱动，做 $\Delta t$ 减半检验。
3. **加速**（中等）：cell list、向量化、并行种子、代理模型（高斯过程/神经网络）替代内层仿真。
4. **校准升级**（中等）：ABC-SMC、模拟矩估计、历史匹配（history matching），给出参数后验而非点估计。
5. **混合建模**（高）：把"平均场部分"写成 ODE、只把异质/局部部分留给 ABM，兼顾可解释与速度。
6. **敏感性**（高）：Morris 筛选 + Sobol 分解，识别哪个规则参数真正驱动输出分布。

#### 11. 与其他模型的组合

- **ABM + ODE（3.2）**：ODE 当"平均场骨架"，ABM 当"异质性修正"；两者在同一参数下对比，差异就是异质性的贡献。
- **ABM + CA（4.5）**：CA 是格点上的同质主体，ABM 给格点加属性与决策；空间蔓延用 CA 的邻域，个体差异用 ABM 的属性。
- **ABM + 贝叶斯反演（4.3）**：2025 楼梯题的正模型可用 ABM 生成 $n(x,y)$，再用观测磨损场做反演；ABM 提供机制，反演提供参数。
- **ABM + 网络（5.3）**：接触网络决定传播阈值，网络结构本身就是政策变量（如社交距离）。

#### 12. 手算算例

**（a）2025 楼梯：ABM 输出与反演对账（原文数字）。** 主体 = 行人（速度、方向、是否结伴）；环境 = 楼梯宽度与踏步；规则 = 靠右行走、超越条件、避让；输出 = 每个位置的踩踏次数 $n(x,y)$，磨损 $h=n\,\delta$，$\delta=0.21$ nm/次（2.14 的 Archard 算例）。若反演给出 $1700$ 人次/天、使用 $30$ 年，则总踩踏
$$1700\times365\times30=1.8615\times10^{7}\ \text{次},\qquad h=1.8615\times10^{7}\times2.1\times10^{-10}\ \mathrm{m}=3.91\times10^{-3}\ \mathrm{m}\approx4\ \mathrm{mm},$$
与 2.14 的实测 4 mm 自洽 ✓。ABM 的任务不是重新算这个总量，而是解释**形状**：靠右规则使上行/下行流各偏一侧，叠加后中央偏右出现峰值；两人并行时横向分布变宽。把 ABM 的 $n(x,y)$ 归一化后与反演得到的相对分布对比，形状吻合才说明"靠右 + 并行"是合理机制。

**（b）种子与误差棒。** 单次运行 $N=200$ 个行人，右侧行走比例真值 $0.8$，则单次估计的标准误 $\sqrt{0.8\times0.2/200}=0.028$；跑 $R=30$ 个种子，均值的标准误 $0.028/\sqrt{30}=0.0051$。若两情景的真实差异只有 $0.01$，30 个种子给出的 95% 区间半宽约 $1.96\times0.0051=0.010$——刚好卡在能否分辨的边缘，必须加种子到 $\ge100$。

**（c）双峰的临界质量。** 用（3）的 $F(\phi)=3\phi^2-2\phi^3$：初值 $\phi_0=0.55$，迭代 $\phi_1=0.575\to0.611\to0.664\to0.737\to0.829\to0.923\to0.983\to0.999$，加速奔向 $1$；初值 $0.45$ 则 $0.425\to0.389\to0.336\to0.263\to0.171\to0.077\to0.018\to0.001$，滑向 $0$。验证一步：$F(0.55)=3(0.3025)-2(0.166375)=0.9075-0.33275=0.57475$。同一个模型、同一组参数，只因初值差 $0.1$，结局相反——这就是"临界质量"与双峰的来源，也是政策仿真里"小幅干预可能翻转结局"的数学依据。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% 楼梯人流 ABM：靠右规则 → 侧向踩踏分布（4.6）
N = 200; T = 300; rng(11);
up = rand(N,1) < 0.5;                    % 上行/下行
x  = 0.2*(2*up-1) + 0.05*randn(N,1);     % 初始侧向位置（m）
edges = linspace(-1,1,42); H = zeros(1,41);
for t = 1:T
    b = 5e-4*(2*up-1);                   % 靠右偏置（上行向右、下行向左）
    x = x + b + 0.02*randn(N,1);         % 偏置 + 避让噪声
    x = max(-1,min(1,x));                % 楼梯边界
    H = H + histcounts(x,edges);         % 累计每个位置的踩踏次数
end
bar(edges(1:end-1)+0.025,H); xlabel('lateral x (m)'); ylabel('visits');
fprintf('total visits=%d (should be %d)\n',sum(H),N*T);
```

输入是 $N,T$、靠右偏置与噪声；输出是侧向踩踏直方图 $n(x)$。**怎么判断对**：`sum(H)` 必须等于 $N\times T$（每次每人记一次，守恒）；偏置取 0 时直方图应近似对称，偏置加大时峰值向两侧分开；把 $N$ 加倍，直方图应更平滑但峰位不变（$1/\sqrt N$ 噪声缩小）。这三条都过，才说明"靠右规则 → 双峰"的因果链成立。

#### 14. 真实美赛 A 题应用

- **2025 楼梯磨损**：ABM 用"行人 + 靠右 + 避让"生成踩踏分布，解释磨损为何中间偏右，并与 2.14 的反演结果对照；注意原文的定位——ABM 在这里生成**机制假设**，不是替代反演。
- **2023 植物群落**：个体竞争（光照、水分）与干旱耐受差异决定斑块与存续，ABM 能给出"部分个体死亡、部分存活"的分布，而 ODE 只给平均盖度。
- **2021 真菌**：菌丝个体在枯落物上的扩展可用 ABM 表达"局部接触 + 资源阈值"，与 CA 版对照可区分"格点效应"与"个体差异"。
- **2020 渔业**：船队主体的捕捞决策（去哪里、捕多少）与资源反馈构成 ABM，可模拟配额政策下的分布效应。

#### 15. 常见错误

1. 规则无依据：每条规则都"看起来合理"，但拿不出文献/观测，最后是用代码编故事。
2. 单次运行下结论：把一次轨迹的波动当成机制，不报种子数与误差棒。
3. 调度不说明：同步还是异步、更新顺序如何，直接决定结果，必须写进假设。
4. 用 ABM 输出当精确预测：ABM 给分布，不给"明天几点几分"。
5. 校准用了全部数据：没有留出验证集，过拟合无法发现。
6. 并行时种子流重叠：不同 worker 用了相关随机数，误差棒虚假变窄。
7. 邻域半径 $r$ 与域长 $L$ 不分：$r\approx L$ 时已经是平均场，却声称"局部交互"。
8. 时间步太大：一步内多个事件被吞，动力学被人为减慢，却归因于参数。
9. 只比较均值：两情景均值相同但分布不同（一个窄一个双峰），结论可能完全相反。

#### 16. 学习要求（验收）

1. 能写出 ABM 五元组规格表（主体/规则/环境/调度/交互），并说明每条规则的依据。
2. 能推出均匀混合下的平均场 ODE，并说明局部交互在什么条件下使它高估传播。
3. 能用阈值模型 $F(\phi)=3\phi^2-2\phi^3$ 手算三个不动点与稳定性，解释双峰分布从何而来。
4. 能设计 $R\ge30$ 的种子实验，给出误差棒与"分辨 $\delta$ 所需 $R$"的公式。
5. 能用守恒量（总踩踏数、总人数）与极限情形（$r\to L$ 退化为 ODE）验证代码。

### 4.7 空间插值与 Kriging

**定位**：P1｜难度 ★★★｜4 h｜前置：3.8 概率统计、线性代数｜测点稀疏时的空间场重建。

#### 1. 定义

给定区域 $\mathcal D\subset\mathbb R^d$ 上的测点 $\{(\mathbf x_i,z_i)\}_{i=1}^n$，把待重建的场看成**随机场（random field）** $Z(\mathbf x)$：均值 $\mathbb E Z(\mathbf x)=m$（普通 Kriging 假定 $m$ 未知但为常数），空间相关结构由变异函数 $\gamma(\mathbf h)$ 描述。目标是在未测点 $\mathbf x_0$ 处给出**线性**预测 $\hat z_0=\sum_iw_iz_i$，并同时要求**无偏** $\mathbb E(\hat Z_0-Z_0)=0$ 与**最小方差** $\mathrm{Var}(\hat Z_0-Z_0)$ 最小；未知量是 $n$ 个权重 $w_i$。

按对均值的假定分三档：均值已知为常数是**简单 Kriging（simple Kriging, SK）**；均值未知、靠 $\sum w_i=1$ 消掉是**普通 Kriging（ordinary Kriging, OK）**；均值含趋势 $\sum_k\beta_kf_k(\mathbf x)$ 是**泛 Kriging（universal Kriging, UK）**。美赛用 OK 足够，趋势明显时升 UK。

#### 2. 核心思想与大白话

**为什么会有它。** 反距离加权（inverse distance weighting, IDW）规定权重是 $d_i^{-p}$：形式人为指定，与"这片场在一公里内到底还相不相关"无关；测点一密一疏，IDW 会在每个测点周围鼓起一个"牛眼（bull's eye）"。Kriging 换了个问法：不问"该用哪个公式加权"，而问"在无偏前提下哪组权重使误差方差最小"——权重由数据自己的空间相关性解出来，不用调参，还顺带给出预测方差。

**大白话。** 把每个测点想成一个消息源，它对自己周围的可信度随距离衰减，衰减曲线就是变异函数；预测未知点时按可信度加权平均，但强制"权重和为 1"，否则平均会整体偏高或偏低。IDW 只按距离给权重，Kriging 还看"这个距离上数据实际差多少"：两点观测本来就接近，权重可以给大些；差得离谱，说明这段距离上变异剧烈，权重该收小。

**一个不能省的附加条件。** 衰减曲线必须"合法"，即对应某个真实存在的协方差函数。随手用最小二乘拟合一条曲线当变异函数，算出的预测方差可能是负数——相当于自称可信度超过 100%，说明这组权重在数学上根本不可能是一个方差最小解。

#### 3. 数学基础

- **二阶平稳（second-order stationarity）**：$\mathbb E Z(\mathbf x)=m$、$\mathrm{Cov}(Z(\mathbf x),Z(\mathbf x+\mathbf h))=C(\mathbf h)$ 与位置无关。作用：把 $n\times n$ 协方差矩阵写成只依赖点对距离的形式。
- **内蕴假设（intrinsic hypothesis）**：只要求增量 $Z(\mathbf x+\mathbf h)-Z(\mathbf x)$ 二阶平稳。比二阶平稳弱，是变异函数存在的条件，也是"用差值不用绝对值"的理由。
- **$\gamma(\mathbf h)=C(\mathbf 0)-C(\mathbf h)$**：把协方差翻译成变异函数。作用：方程组可只用 $\gamma$ 写出，而 $\gamma$ 由差值估计，不必先知道均值。
- **条件负定性（conditional negative definiteness）**：$\sum_i\sum_ja_ia_j\gamma(\mathbf x_i-\mathbf x_j)\le0$ 对所有 $\sum_ia_i=0$ 成立。作用：这是"存在合法协方差"的充要条件，也是 Kriging 方差非负的保证。
- **拉格朗日乘子**：等式约束 $\sum w_i=1$ 下的二次型最小化。乘子 $\mu$ 有单位（场值²），是"无偏约束的价格"，不是可随手调的惩罚系数。
- **LU 分解**：解对称系统，$O(n^3)$，是 Kriging 的计算内核；但系数矩阵不是正定矩阵（见第 6 小节），不能直接 Cholesky。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
|---|---|---|---|---|
| 内蕴/二阶平稳 | 均值常数，协方差只依赖 $\mathbf h$ | 同一气候/地质单元内 | 有趋势时 OK 系统性削峰，山脊处低估 10—30% | 残差对坐标回归，看趋势项是否显著 |
| 各向同性 | $\gamma$ 只依赖距离 $\Vert\mathbf h\Vert$ | 距离度量与方向无关 | 变程被方向平均，长轴方向偏差可达一个变程量级 | 分 0°、45°、90° 画实验变异函数 |
| 变异函数模型合法 | 拟合模型属条件负定族 | 只用球状/指数/高斯/Matérn | 出现负预测方差，权重正负振荡（±5 以上） | 检查所有 $\sigma_0^2\ge0$、$\Gamma$ 特征值 |
| 增量近似对称、无离群 | 无单点主导变异函数 | 测量误差小、场连续 | 一个离群点把基台抬高近 1 倍，整片预测被拉偏 | 变异函数云图、箱线图、重复测量 |
| 观测误差可忽略 | $z_i$ 近似真值（块金小） | 仪器精度高于空间变异 | 预测退化为平滑曲面而非过点插值，极值被低估 | 用重复测量估计 $c_0$ |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
|---|---|---|---|---|
| $\mathbf x$ | 空间位置 | m（或经纬度） | 研究区尺度 | 决定点对距离，进而决定 $\gamma$ 的取值 |
| $Z(\mathbf x)$ | 场值（随机变量） | 与 $z$ 同（°C、mg/L、mm） | — | 待重建的量；其空间变异被 $\gamma$ 概括 |
| $z_i$ | 第 $i$ 个测点观测 | 同上 | 20—500 个点 | 数据唯一来源；点数决定变异函数能否估出 |
| $\gamma(\mathbf h)$ | 变异函数（半方差） | 场值² | $0\sim(c_0+c)$ | 越大表示相隔 $\mathbf h$ 的两点差异越剧烈 |
| $c_0$ | 块金（nugget） | 场值² | 基台的 0—30% | 测量误差加微尺度变异；变大则预测不再过点 |
| $c$ | 偏基台（partial sill） | 场值² | 与 $c_0$ 同量级 | 结构性方差；$c_0+c$ 即基台（sill） |
| $a$ | 变程（range） | m | 几十米—几十千米 | 超过此距离近似不相关；变大则场更平滑、更"全局" |
| $w_i$ | 权重 | — | $[0,1]$，$\sum w_i=1$ | 该测点对本次预测的贡献份额 |
| $\mu$ | 拉格朗日乘子 | 场值² | 可正可负 | 无偏约束的价格；为正时抬高方差，表示"被约束住"的代价 |
| $\sigma_0^2$ | 预测方差 | 场值² | $\ge0$，远点趋于基台量级 | 该位置的局部不可信程度 |

#### 6. 从基本原理推导核心公式

**第一步：把"无偏"变成约束。** 记 $C_{ij}=C(\mathbf x_i-\mathbf x_j)$、$C_{i0}=C(\mathbf x_i-\mathbf x_0)$、$C_{00}=C(\mathbf 0)$。误差 $e=\hat Z_0-Z_0=\sum_iw_iZ_i-Z_0$，故 $\mathbb E e=m\sum_iw_i-m=m(\sum_iw_i-1)$。因为 $m$ 未知而结果要对任意 $m$ 无偏，只能要求 $\sum_iw_i=1$——这一步就是 OK 与 SK 的分界：若硬把 $m$ 当 0，约束可以省掉，但均值估错会直接变成偏差。

**第二步：写出方差。** 由线性组合方差公式 $\sigma^2(\mathbf w)=\sum_i\sum_jw_iw_jC_{ij}-2\sum_iw_iC_{i0}+C_{00}$。量纲：$C$ 是场值²，$w$ 无量纲，故 $\sigma^2$ 是场值²。

**第三步：换成变异函数。** 代入 $C_{ij}=C_{00}-\gamma_{ij}$、$C_{i0}=C_{00}-\gamma_{i0}$，并用 $(\sum_iw_i)^2=1$：

$$
\sigma^2=C_{00}\Big[\big(\textstyle\sum_iw_i\big)^2-2\sum_iw_i+1\Big]+2\sum_iw_i\gamma_{i0}-\sum_i\sum_jw_iw_j\gamma_{ij}
=2\sum_iw_i\gamma_{i0}-\sum_i\sum_jw_iw_j\gamma_{ij}.
$$

方括号里恰好是 $1-2+1=0$：这就是要用 $\gamma$ 而不用 $C$ 的理由——均值与 $C_{00}$ 一起被消掉，只剩可估的差值统计量。

**第四步：拉格朗日。** 约束写成 $1-\sum_iw_i=0$，构造 $L(\mathbf w,\mu)=2\sum_iw_i\gamma_{i0}-\sum_i\sum_jw_iw_j\gamma_{ij}+2\mu(1-\sum_iw_i)$，对 $w_k$ 求偏导并令零（用 $\gamma_{ij}=\gamma_{ji}$）得 $2\gamma_{k0}-2\sum_j\gamma_{kj}w_j-2\mu=0$，即 $\sum_j\gamma_{kj}w_j+\mu=\gamma_{k0}$，$k=1,\dots,n$。连同约束写成矩阵：

$$
\begin{bmatrix}\boldsymbol\Gamma&\mathbf 1\\ \mathbf 1^\top&0\end{bmatrix}
\begin{bmatrix}\mathbf w\\ \mu\end{bmatrix}
=\begin{bmatrix}\boldsymbol\gamma_0\\ 1\end{bmatrix},
\qquad \Gamma_{ij}=\gamma(\mathbf x_i-\mathbf x_j),\ \ \gamma_{0,i}=\gamma(\mathbf x_i-\mathbf x_0).
$$

这就是**普通 Kriging 方程组**。$\Gamma$ 的对角元是 $\gamma(\mathbf 0)=0$，所以 $\Gamma$ 本身不定、整个系数矩阵也不正定——它靠"最后一行恰好补上约束"才可解。

**第五步：化简方差（最容易漏的一步）。** 用 $\boldsymbol\Gamma\mathbf w+\mu\mathbf 1=\boldsymbol\gamma_0$ 左乘 $\mathbf w^\top$，得 $\mathbf w^\top\boldsymbol\Gamma\mathbf w+\mu=\sum_iw_i\gamma_{i0}$，即 $\sum_i\sum_jw_iw_j\gamma_{ij}=\sum_iw_i\gamma_{i0}-\mu$。代回第三步：

$$
\sigma_0^2=\sum_iw_i\gamma_{i0}+\mu .
$$

量纲：两项都是场值²。边界检查一：若 $\mathbf x_0$ 正落在测点 $k$ 上，则 $\gamma_{0,k}=0$，解得 $w_k=1$、其余为 0、$\mu=0$，于是 $\sigma_0^2=0$，预测精确等于 $z_k$——这就是"精确插值器（exact interpolator）"。边界检查二：目标远离所有测点、测点彼此靠近时，$\gamma_{i0}\to c_0+c$、$\gamma_{ij}\to c_0+c$，对称性给出 $w_i=1/n$、$\mu=(c_0+c)/n$，于是 $\sigma_0^2\to(c_0+c)(1+1/n)$：比基台还大一点，多出的 $1/n$ 是"连均值也要估"的代价（SK 此时正好趋于基台）。方差有上界、不随距离无限增大，这是 Kriging 与 IDW 的本质差别。

**第六步：半正定条件。** 要求 $\sigma^2(\mathbf w)\ge0$ 对一切 $\sum w_i=1$ 成立。把两组权重之差记为 $a_i$（$\sum a_i=0$），则 $\mathrm{Var}(\sum_ia_iZ_i)=\sum_i\sum_ja_ia_jC_{ij}=-\sum_i\sum_ja_ia_j\gamma_{ij}\ge0$（第二个等号用 $\sum a_i=0$ 消去 $C_{00}$）。所以"$\gamma$ 条件负定"与"存在合法协方差"是同一件事，也是 $\sigma^2\ge0$ 的充要条件。检验：指数协方差 $C(h)=ce^{-h/a}$ 的两点相关矩阵 $\begin{bmatrix}1&e^{-1}\\ e^{-1}&1\end{bmatrix}$ 特征值为 $1\pm e^{-1}=1.632$ 与 $0.368$，皆正，合法；球状模型只在 $d\le3$ 维正定；高斯模型 $C(h)=ce^{-(h/a)^2}$ 也正定，但过于光滑，矩阵常近奇异。

**第七步：模型族与参数。** 球状 $\gamma(h)=c_0+c(\tfrac{3h}{2a}-\tfrac{h^3}{2a^3})$（$h\le a$），在 $h=a$ 处斜率归零、相关性精确截止，物理上最像"有限影响范围"；指数 $\gamma(h)=c_0+c(1-e^{-h/a})$ 要到约 $3a$ 才接近基台；高斯 $\gamma(h)=c_0+c(1-e^{-h^2/a^2})$ 在原点斜率为 0，对应极光滑的场，易致数值病态。三者都满足 $\gamma(0)=0$，而 $\gamma(0^+)=c_0$ 的跳变就是块金。

#### 7. 求解方法

- **解析解**：只有 $n=1,2$ 时有闭式；一般只能解线性系统。**不要**用 `chol`（矩阵对角为 0，非正定），用 `\`（部分主元 LU），或先消去 $\mu$ 再解 $\Gamma\mathbf w=\boldsymbol\gamma_0-\mu\mathbf 1$。
- **规模**：$n\lesssim2000$ 直接解，$O(n^3)$；更大用移动邻域（每点取最近 30—50 点，局部相容性换速度），或固定秩/低秩、SPDE 近似。
- **病态来源**：两点距离远小于变程时 $\gamma_{ij}\approx0$，矩阵两行几乎相同，条件数爆炸；高斯模型最严重。对策：合并近重合点、加很小块金、改用 Matérn（$\nu=3/2$）。
- **变异函数拟合**：先分箱算 $\hat\gamma(h)=\frac{1}{2N(h)}\sum_{(i,j)\in N(h)}(z_i-z_j)^2$，每箱至少 30 对点；再用加权最小二乘（权重 $N(h)/\gamma^2$）或 REML 在合法族内拟合 $(c_0,c,a)$。
- **必报三个诊断**：变异函数拟合图（模型是否跟住前几个箱）、留一交叉验证 RMSE（预测精度）、预测方差图（哪里不可信）。
- **MATLAB**：核心只需 `pdist2` 造距离矩阵加 `\` 解方程；没有工具箱也能十行写完。

#### 8. 适用场景

**该用**：测点 20—500 个、场在空间上连续（温度、浓度、磨损深度、高程、土壤水分）；实验变异函数有明显变程（前几个箱上升后走平）；需要"哪里不可信"的不确定度图；只做凸包内插值或极短距离外推。

**不该用**：测点少于 8 个（三个参数不可辨识）；场有强趋势或突变界面（断层、海岸线、城市边界）——先去掉趋势或升 UK，否则 OK 把界面抹成斜坡；目标是**极值**而非均值（Kriging 是最小方差平滑器，必然削峰）；观测是有无或类别（用指示 Kriging）；测点是"方便采样"而非随机布设（有偏采样让变异函数失真）。

#### 9. 局限

- **平滑效应**：最小方差必然削峰填谷，极值低估 10—40%，数据越稀越严重；要极值必须用条件模拟。
- **外推不可靠**：超出凸包后 $\gamma_{i0}$ 饱和，预测趋于局部均值，而方差被低估，真实误差远大于报出的 $\sigma_0$。
- **各向异性/非平稳**：单一变程把长轴方向的真实变化抹平，误差可达一个变程量级。
- **计算**：$O(n^3)$，$n=10^4$ 时单点求解已不现实。
- **小样本不可辨识**：$n<30$ 时 $(c_0,c,a)$ 互相补偿，不同参数组给出几乎相同的拟合图，但预测方差可差一倍。
- **只保证线性最优**：非线性、非高斯场（有检出限的浓度、计数数据）需变换或指示方法。

#### 10. 改进

1. **加各向异性**：先做坐标线性变换 $\mathbf x'=\mathbf A\mathbf x$（旋转+拉伸），再按各向同性拟合。解决"变程被平均"，代价是多 3—4 个参数与分方向的变异函数图。
2. **泛 Kriging**：均值写成 $\sum_k\beta_kf_k(\mathbf x)$（常数、$x$、$y$、$x^2$），无偏约束变为 $\sum_iw_if_k(\mathbf x_i)=f_k(\mathbf x_0)$。解决趋势导致的偏差，代价是外推更不稳。
3. **回归 Kriging/外部漂移**：用协变量（高程、NDVI、遥感指数）解释均值，Kriging 只拟合残差。解决"稀疏 + 有强相关协变量"，代价是要有全覆盖的协变量图层。
4. **指示 Kriging / 序贯高斯模拟（SGS）**：要极值、要多个等概率实现时用。解决平滑效应，代价是失去单一最优图，输出变成集合。
5. **贝叶斯/低秩空间模型**：$n$ 大时用固定秩加 SPDE 近似。代价是实现复杂度与近似误差。

#### 11. 与其他模型的组合

- **Kriging + 机理模型（PDE）**：机理模型给趋势 $m(\mathbf x)$（如磨损率的物理预测），Kriging 拟合残差场，最终预测 = 机理 + 残差校正。接口是"残差场"，适用于"有物理模型但参数空间变异"的题（2025 磨损深度场正是这种结构）。
- **Kriging + 蒙特卡洛/敏感性（4.4、3.10）**：把 $\sigma_0^2$ 当输入不确定度传播，得到"参数→场→指标"的区间；接口是方差图。
- **Kriging + 优化（3.9）**：用 $\hat z$ 当目标、$\sigma_0^2$ 当采集函数（expected improvement），决定"下一个测点打在哪"——贝叶斯优化的标准接口。

#### 12. 手算算例

原文无算例，这里给一个两点、可手算的例子。取指数变异函数 $c_0=0$、$c=1$、$a=1$（距离单位任意），则 $\gamma(h)=1-e^{-h}$。测点 $x_1=0$、$z_1=0$；$x_2=1$、$z_2=1$；预测 $x_0=0.25$。

先算三个值：$\gamma_{01}=1-e^{-0.25}=0.2212$，$\gamma_{02}=1-e^{-0.75}=0.5276$，$\gamma_{12}=1-e^{-1}=0.6321$。方程组为

$$
\begin{bmatrix}0&0.6321&1\\ 0.6321&0&1\\ 1&1&0\end{bmatrix}
\begin{bmatrix}w_1\\ w_2\\ \mu\end{bmatrix}=\begin{bmatrix}0.2212\\ 0.5276\\ 1\end{bmatrix}.
$$

第 1、2 行相减消去 $\mu$：$0.6321(w_2-w_1)=0.2212-0.5276=-0.3064$，故 $w_2-w_1=-0.4848$；与 $w_1+w_2=1$ 联立得 $w_1=0.7424$、$w_2=0.2576$。回代第 1 行：$\mu=0.2212-0.6321\times0.2576=0.0584$。

预测 $\hat z_0=0.7424\times0+0.2576\times1=0.2576$，几乎等于线性插值 $0.25$：因为 $a=1$ 与跨度同量级，场非常光滑，远处那个点也有可观权重。方差 $\sigma_0^2=0.7424\times0.2212+0.2576\times0.5276+0.0584=0.3585$，标准差 $0.599$（相对基台 1 很大——只有两个点）。

两个检查。其一，取 $x_0\to0$：$\gamma_{01}\to0$，解出 $w_1\to1$、$w_2\to0$、$\mu\to0$、$\sigma^2\to0$，回到精确插值。其二，同样两点若用 IDW（$p=2$）：$w_2=(1/0.75^2)/[(1/0.25^2)+(1/0.75^2)]=0.100$，预测 $0.100$，与 Kriging 的 $0.2576$ 差一倍多——IDW 只看距离，Kriging 还看"这个距离上数据实际差多少"。

#### 13. 最小可运行 MATLAB 示例

```matlab
%% 普通 Kriging（指数变异函数）+ 留一交叉验证
rng(0);
xy = [0 0; 1 0; 0 1; 1 1; 0.5 0.5];  z = [1 2 1.5 2.5 2.2]';  % 5 个测点
x0 = [0.3 0.4];                                   % 待预测位置

g  = @(h) 1 - exp(-h/0.8);                        % 指数模型: c0=0, c=1, a=0.8
D  = pdist2(xy,xy);  G = g(D);  n = numel(z);     % Gamma_ij = gamma(||xi-xj||)
d0 = pdist2(xy,x0)'; g0 = g(d0);                  % gamma_0i
K   = [G ones(n,1); ones(1,n) 0];                 % [Gamma 1; 1' 0]
rhs = [g0; 1];
sol = K\rhs;  w = sol(1:n);  mu = sol(end);       % 解 Kriging 方程组
zhat = w'*z;  s2 = w'*g0 + mu;                    % 预测值与方差（第 6 小节公式）
fprintf('预测 z(0.3,0.4)=%.4f, 预测标准差=%.4f\n', zhat, sqrt(max(s2,0)));

% 留一交叉验证：每次删一个点，用其余点预测它
err = zeros(n,1);
for i = 1:n
    idx = setdiff(1:n,i);
    Ki  = [G(idx,idx) ones(n-1,1); ones(1,n-1) 0];
    si  = Ki\[G(idx,i); 1];
    err(i) = z(i) - si(1:n-1)'*z(idx);
end
fprintf('LOOCV RMSE = %.4f\n', sqrt(mean(err.^2)));
```

输入是坐标矩阵 `xy`、观测列向量 `z`、预测点 `x0` 与变异函数参数 $(c_0,c,a)$；输出是预测值、预测标准差与 LOOCV RMSE。
**判对准则一**：把 `x0` 设成某个测点坐标，应得 `zhat` 等于该点观测值、`s2` 在 $10^{-15}$ 量级。
**判对准则二**：把 `a` 调到 1e4，`zhat` 应趋于 `mean(z)`、`s2` 应趋于 $c_0+c$——这两个极限正是第 6 小节的"精确插值"与"远点极限"检查。

#### 14. 真实美赛 A 题应用

- **2025**：由有限测点重建磨损深度场。模型手要做的不是插一张彩色图，而是：先用实验变异函数说明场确实有空间相关性，再报 LOOCV RMSE 支撑"这个场可信"，最后用预测方差图标出哪些区域是外推、不可信。
- **2023**：土壤水分空间分布。典型陷阱是测点沿样带分布（各向异性加有偏采样），需分方向看变异函数，并说明采样设计对结论的限制。
- 通用打法：先判断"趋势 vs 残差"——若场沿某方向单调漂移，直接 OK 会把漂移当相关性，此时应升 UK 或先回归掉趋势再对残差做 Kriging。

#### 15. 常见错误

1. **不拟合变异函数直接插值**：用默认参数或直接套 IDW，等于放弃 Kriging 的全部优势，却仍以"Kriging"名义报结果。
2. **对 $\Gamma$ 用 `chol`**：对角为 0、非正定，必然报错或给出无意义结果；应解一般线性系统。
3. **近重合测点不处理**：矩阵两行几乎相同，权重振荡到 ±10，预测跑到观测范围之外；应先合并或加小块金。
4. **把 $\sigma_0^2$ 当"预测误差"**：它是模型方差，不含变异函数误设造成的偏差，通常比真实误差小；两者应并报。
5. **用样本内残差报精度**：Kriging 在测点处误差恒为 0（精确插值），样本内 RMSE 必然好看；必须 LOOCV 或留出独立测点。
6. **忽略各向异性**：变程被方向平均，长轴方向系统性偏差，诊断图上表现为残差沿某方向成带。
7. **变异函数用任意曲线拟合**：$R^2$ 很高但模型非条件负定，出现负方差。
8. **不报预测方差图**：评审无法判断外推区域，等于隐藏了最大的不确定性来源。
9. **经纬度当米用**：变程失去物理意义（1° 纬度约 111 km，经度随纬度收缩），参数无法解释也无法复用。

#### 16. 学习要求（验收）

- 能独立完成"无偏约束 → 方差表达式 → 拉格朗日 → 方程组 → $\sigma_0^2=\sum_iw_i\gamma_{i0}+\mu$"这条链，并说清每一步用了哪个恒等式。
- 能解释为什么系数矩阵不能 Cholesky、为什么变异函数必须条件负定、负方差意味着什么。
- 给定变异函数与两点数据，能手算权重、预测值与方差，并用"目标落在测点上"和"变程趋于无穷"两个极限自检。
- 能对一份真实插值结果报出三个诊断（变异函数拟合图、LOOCV RMSE、方差图），并据此说出"哪里不可信、为什么"。

### 4.8 多目标优化与 Pareto

**定位**：P1｜难度 ★★★｜4 h｜前置：3.9｜成本—效益、生态—经济冲突。

#### 1. 定义

多目标优化（multi-objective optimization, MOO）的标准形式为
$$
\min_{\mathbf x\in X}\ \mathbf f(\mathbf x)=\big(f_1(\mathbf x),\dots,f_m(\mathbf x)\big),\qquad m\ge 2,
$$
其中可行域 $X\subset\mathbb R^n$ 由约束划定（已知），$f_k$ 是第 $k$ 个目标的取值（已知表达式或可计算的仿真）。
未知量不是单点，而是整条 **Pareto 前沿**（Pareto front）：$m\ge2$ 时目标此消彼长，"最优"只能定义为"不被别人全面压过"。
称 $\mathbf x$ **支配**（dominate）$\mathbf y$，若 $f_k(\mathbf x)\le f_k(\mathbf y)$ 对所有 $k$ 成立且至少一个严格，记 $\mathbf x\prec\mathbf y$；$X$ 中不被任何点支配的解称 Pareto 最优（Pareto optimal），其像点集合 $\mathcal P=\{\mathbf f(\mathbf x^*)\}$ 即前沿。求解因此分两层：先把前沿（或其一段）算准，再按明确规则从上面挑一个方案。

#### 2. 核心思想与大白话

笨办法是把多个目标拍成一个总分：收益减成本、加权打分表，或者只优化一个、其余"顺便看看"。问题不在算，而在**权重是拍的**——权重换个数量级，推荐方案就换一个；更糟的是有些方案（哪怕真的无可替代）永远进不了候选。
Pareto 的思路是先不谈偏好：把所有"没人能在各方面都压过它"的方案摊在桌上，连成一条曲线，再由决策者在曲线上谈条件。
生活类比一：买房比"总价"和"通勤时间"。没有哪套既更便宜又更近，剩下那几套就是前沿；挑哪套取决于你更在乎钱还是时间。
生活类比二：五个人分蛋糕，切法很多；如果一种切法做不到"每人都不比另一种少、且有人更多"，那这些切法就是有效的——剩下的分歧是公平观，不是数学。
所以前沿不是"最优解"，而是"没法再用数学淘汰的候选集"。
把前沿交给评审时必须同时给出挑选规则（权重来源、红线约束、膝点理由），否则等于把决策偷偷留给自己。

#### 3. 数学基础

- **支配关系是偏序**：用来在完全不引入权重的情况下定义"更好"，只需各目标同向（都最小化或都最大化）。
- **凸集与凸包**：加权和只能触到前沿的凸包边界，这是它漏掉非凸段的几何原因。
- **支撑超平面**：若 $\mathbf f(\mathbf x^*)$ 是某个正权重的最小值，则以 $\mathbf w$ 为法向的超平面在该点支撑像集；反过来，不被任何超平面支撑的点就取不到。
- **KKT 条件**：ε-约束问题的最优解满足 $\nabla f_1+\lambda\nabla f_2=\mathbf 0$、$\lambda\ge0$，$\lambda$ 是"放松 ε 一单位能省多少 $f_1$"的影子价格。
- **包络定理**：$\lambda=-\mathrm df_1^*/\mathrm d\varepsilon$，用来判断还该不该继续压第二个目标。
- **非支配排序**：按"被支配层数"给种群分层，$O(mN^2)$，是 NSGA-II 选择的基础。
- **拥挤距离**：同层内按相邻点间距排序，保证前沿铺得开而不是挤成一坨。
- **超体积（hypervolume）与 IGD**：把"前沿好坏"压成一个数，便于比较算法与随机种子。
- **归一化**：加权和之前必须除以参考量级，否则权重的含义被量纲吞掉。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
| --- | --- | --- | --- | --- |
| 各目标已统一方向 | 全部转成最小化（越大越好者乘 $-1$） | 支配关系要求逐项可比 | 漏改一个符号，前沿整体错位，推荐方案恰好取反 | 对每个 $f_k$ 做单变量扰动，看符号是否与"更好"一致 |
| 目标之间无隐藏硬约束 | 权衡可以连续交换 | Pareto 框架本身不设红线 | 若其实有排放/预算上限，前沿上一半方案不可行，结论被打回 | 把红线写成约束重算，比较前沿收缩幅度 |
| 决策变量维数可控 | $n$ 在数十以内 | 进化算法可搜索 | $n$ 上千时 200 代仍欠收敛，前沿点间距偏大 30% 量级 | 加倍代数看前沿是否仍移动 |
| 目标函数确定、可重复 | 同一 $\mathbf x$ 得同一 $\mathbf f$ | 支配关系可判定 | 目标噪声 5% 时接近点的支配关系约三分之一会翻转，前沿长毛刺 | 同一 $\mathbf x$ 重复评估，看方差 |
| 目标数 $m$ 不大 | $m=2$ 或 $3$ | 前沿可视化与采样 | $m\ge4$ 时同样分辨率需指数级样本，报出的"前沿"极稀疏 | 报告点数与目标维数之比 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
| --- | --- | --- | --- | --- |
| $\mathbf x$ | 决策变量向量 | 视题而定（艘、亿元、MW） | $n=2\sim10$ | 可调的方案参数 |
| $f_k$ | 第 $k$ 个目标 | 与题一致（亿元、吨 CO$_2$、人） | $m=2\sim3$ | 要压低的指标 |
| $w_k$ | 第 $k$ 个目标的权重 | —（无量纲） | $0\sim1$，$\sum_kw_k=1$ | 偏好强度 |
| $\varepsilon$ | 第二目标的允许上限 | 同 $f_2$ | 由决策者给 | 把软目标变成硬约束 |
| $N$ | 种群规模 | 个 | $50\sim200$ | 每代保留的候选解数 |
| $P$ | 非支配解个数 | 个 | $10\sim100$ | 前沿点密度 |
| $d_i$ | 拥挤距离 | 目标空间长度（已归一化） | $0\sim\infty$ | 相邻前沿点的间距 |
| $HV$ | 超体积 | 目标空间体积（已归一化） | 越大越好 | 收敛性与分布的综合度量 |

#### 6. 从基本原理推导核心公式

**第一步：加权和的解必是 Pareto 最优（反证）。** 设 $\mathbf x^*$ 最小化 $\sum_kw_kf_k(\mathbf x)$，$w_k>0$。若存在 $\mathbf y$ 支配 $\mathbf x^*$，则 $f_k(\mathbf y)\le f_k(\mathbf x^*)$ 且至少一个严格，乘以正权重求和得 $\sum_kw_kf_k(\mathbf y)<\sum_kw_kf_k(\mathbf x^*)$，与 $\mathbf x^*$ 最优矛盾。故加权和给出的解一定在前沿上，且这一步**不需要** $f$ 凸。
量纲检查：证明只用了 $w_k>0$ 与逐项比较，对任何单位都成立；但数值上 $\sum_kw_kf_k$ 只有各 $f_k$ 同量纲时才有意义（见第五步）。

**第二步：反过来不成立（非凸反例）。** 取三方案 A$(0,3)$、B$(2,2)$、C$(3,0)$，两目标都最小化、$w_1+w_2=1$。A 的值 $3w_2$，B 的值 $2$，C 的值 $3w_1$。B 要赢 A 需 $2\le3w_2$，即 $w_1\le1/3$；B 要赢 C 需 $2\le3w_1$，即 $w_1\ge2/3$，无解。所以 **B 在任何权重下都不是最优**，可它不被 A 或 C 支配（A 的 $f_1$ 更好但 $f_2$ 更差，C 反之），是货真价实的 Pareto 点。
几何：B 落在弦 AC 的上方（弦中点 $(1.5,1.5)$ 比 B 两个分量都小），任何正法向的支撑线都碰不到它。

**第三步：为什么非凸段取不到（支撑超平面）。** 对任意 $w_k>0$，$\min\sum_kw_kf_k$ 的解都是像集凸包的支撑点：超平面 $\{\mathbf f:\mathbf w^\top\mathbf f=c\}$ 在该点与凸包相切。落在凸包内凹一侧的点没有这样的切平面，于是任何权重都选不到它。推论：加权和扫描得到的点集 ⊆ 前沿 ∩ 凸包边界；前沿向内凹的地方就是漏掉的方案。

**第四步：ε-约束与 KKT。** 把第二目标变成上限：$\min f_1(\mathbf x)$ s.t. $f_2(\mathbf x)\le\varepsilon$，$\mathbf x\in X$。拉格朗日函数 $f_1+\lambda(f_2-\varepsilon)$ 的 KKT 条件为
$$
\nabla f_1+\lambda\nabla f_2=\mathbf 0,\qquad \lambda\ge0,\qquad \lambda\,(f_2-\varepsilon)=0 .
$$
$\lambda\ge0$ 的道理：放松上限只能让 $f_1^*$ 变小，故 $\lambda=-\mathrm df_1^*/\mathrm d\varepsilon\ge0$；它能读成"$\varepsilon$ 放宽一个单位，第一目标省下 $\lambda$ 个单位"。
关键性质：**每个 Pareto 点都是某个 ε 的最优解**。取 $\varepsilon=f_2(\mathbf x^*)$，若存在 $\mathbf x$ 使 $f_1(\mathbf x)<f_1(\mathbf x^*)$ 且 $f_2(\mathbf x)\le f_2(\mathbf x^*)$，则 $\mathbf x$ 支配 $\mathbf x^*$，矛盾。所以扫 ε 不漏点，代价是解一串单目标问题、分辨率由 ε 步长决定。
边界检查：$\varepsilon\to+\infty$ 退化为 $\min f_1$（前沿一端）；$\varepsilon<\min_Xf_2$ 不可行（越过前沿另一端）。

**第五步：量纲检查与归一化。** 加权和 $\sum_kw_kf_k$ 要求各项可加。若 $f_1$ 是亿元、$f_2$ 是吨 CO$_2$，则 $w_1f_1$ 与 $w_2f_2$ 单位不同，相加没有物理意义；数值上等于把"1 亿元"与"1 吨"当同一把尺子，等效权重会随单位选择变化（亿元换成万元，相对权重差 $10^4$ 倍）。
正确做法：先归一化 $\tilde f_k=f_k/f_k^{\text{ref}}$（$f_k^{\text{ref}}$ 取该目标的量级或"最差减理想"的跨度），再对 $\tilde f_k$ 加权，并在报告里注明"权重作用在归一化目标上"。

**第六步：NSGA-II 的两个核心算子。** 非支配排序：对每对个体比较 $m$ 个目标，共 $O(mN^2)$，把种群分成层 $F_1,F_2,\dots$，$F_1$ 即当前非支配集。拥挤距离：同层内按每个目标排序后
$$
d_i=\sum_{k=1}^{m}\frac{f_k^{(i+1)}-f_k^{(i-1)}}{f_k^{\max}-f_k^{\min}},
$$
两端点取 $\infty$ 以免丢掉边界解。选择规则：先比层号（小者优），同层比拥挤距离（大者优）——既向前沿收敛，又不让点挤在一处。

**第七步：极限与退化。** $m=1$ 时所有点可比，前沿退化为单点，整套机制回落到 3.9 的单目标优化；$w=(1,0)$ 时加权和退化为只优化 $f_1$，非凸段仍取不到；若两目标完全正相关（$f_2=af_1+b$，$a>0$），前沿退化为一点，问题其实只有一个目标。

#### 7. 求解方法

1. 统一方向：把"越大越好"的目标乘 $-1$，全部转成最小化。
2. 归一化：各目标除以参考量级。
3. 离散候选集（方案表几十个以内）：直接两两比较求非支配集，$O(mN^2)$，最稳、最该先做。
4. 加权和扫描：$w$ 取 $0\sim1$、步长 $0.02\sim0.05$，记录每个权重下的最优点；得到的是凸包段。
5. ε-约束扫描：对 $f_2$ 上限取 10—20 个水平逐个解 $\min f_1$，补上非凸段，用上一个解暖启动加速。
6. 连续问题：`gamultiobj` 或 `paretosearch`，种群 50—200、代数 200—500，交叉概率约 0.8、变异率约 $1/n$。
7. 验证：换 10 个随机种子重跑，比较前沿点集是否稳定；画超体积随代数的曲线。
8. 落表：每个前沿点必须能回溯到决策变量 $\mathbf x$，并标注由哪种方法得到。
9. 选择：给出权重/阈值/膝点的明确来源，做 $\pm20\%$ 敏感性。

#### 8. 适用场景

**该用**：目标两个以上且确实冲突（经济 vs 生态、成本 vs 精度、覆盖 vs 预算）；决策者尚未表态，需要先展示"可能的菜单"。
也适用于要论证"我的方案没被明显更优的方案压过"——多目标语言比单目标打分更难被挑刺。
**不该用**：目标其实能合成一个（都能折成钱且折算率公认）——单目标更省事；存在硬红线（先把红线写成约束再谈前沿）；一次评估要数小时而你需要上千次评估（先用 4.7 的代理模型）；只有一个目标却硬拆成两个以显方法高级。

#### 9. 局限

- 前沿点数随目标数指数膨胀：$m=2$ 时 50 点够用，$m=5$ 要同等分辨率需 $50^4$ 量级样本，实际只能看到稀疏点。
- 加权和漏非凸段：本题的 C$(8,8)$ 是前沿点，却在整个 $w_1\in[0,1]$ 上从不被选中——这属于结构性盲区，不是精度问题。
- NSGA-II 无收敛保证：200 代、种群 100 时与真前沿的超体积差常在 5%—15%，换种子前沿会移动。
- 主观性不会消失，只是从"建模"搬到"决策"；报告里不写权重/ε 的来源等于没做。
- 噪声下支配关系不稳：目标噪声 5% 时，相邻很近的两点谁支配谁可能翻转，前沿出现毛刺。
- 目标强相关时前沿退化成近乎一条线，"多目标"只是装饰。

#### 10. 改进

1. 归一化 + 权重敏感性表（几乎零成本，先做）。
2. ε-约束扫描补齐非凸段，与加权和结果合并去重（代价：多解 10—20 个单目标）。
3. 膝点（knee）识别：在前沿上最大化到两端弦的距离，给出"性价比拐点"作为推荐并写明理由。
4. 代理模型：用 Kriging 拟合 $\mathbf x\mapsto\mathbf f$，在代理上跑进化算法，只对候选点做真实评估（代价：几十次真实仿真）。
5. 参考点方法或噪声鲁棒前沿（NSGA-III、MOEA/D、期望目标），代价是实现与调参。

#### 11. 与其他模型的组合

- 与 3.9 单目标优化：ε-约束的每一次求解都是标准单目标问题，求解器直接复用。
- 与 4.7 Kriging：昂贵仿真下用代理模型加速前沿搜索；注意代理误差会让"支配"误判，候选点必须在真实模型上复核。
- 与 4.9 系统动力学：政策杠杆（努力量、税率）作决策变量，"长期均衡存量"与"累计收益"作两个目标，前沿给出稳捕捞与高产出的取舍。
- 与 4.10 时间序列：预测值（需求、价格）作为目标输入；预测区间宽时在区间两端各算一次前沿，看推荐方案是否翻转。
- 与第五章敏感性分析：权重与 ε 的敏感性是多目标特有的必做项。

#### 12. 手算算例

五个方案，两目标都最小化：A$(10,5)$、B$(9,6)$、C$(8,8)$、D$(7,9)$、E$(9,7)$。
支配判定：B 与 E 的 $f_1$ 同为 9，而 B 的 $f_2=6<7$，故 B 支配 E，E 出局。
其余两两比较：A 与 B（$10>9$ 但 $5<6$）、B 与 C（$9>8$ 但 $6<8$）、C 与 D（$8>7$ 但 $8<9$）都互不支配，故前沿 $=\{A,B,C,D\}$。
加权和扫描（$w$ 给 $f_1$，$1-w$ 给 $f_2$）：$w=0$ 只压 $f_2$，取 A（$f_2=5$）；$w=0.5$ 时 A 与 B 同为 $7.5$；$w=0.6$ 时 B 与 D 同为 $7.8$；$w=1$ 只压 $f_1$，取 D（$f_1=7$）。
解不等式得支撑区间：A 在 $w\in[0,0.5]$ 最优，B 在 $[0.5,0.6]$，D 在 $[0.6,1]$；C 无解——赢 A 需 $w\ge0.6$，赢 B 需 $w\ge2/3$，赢 D 需 $w\le0.5$，三条件矛盾。
所以 C$(8,8)$ 是合法 Pareto 点，却是**非支撑点**（non-supported）：弦 BD 的中点 $(8,7.5)$ 与 C 的 $f_1$ 相同而 $f_2$ 更小，说明 C 在凸包上方；若允许在 B、D 之间随机混合，混合方案的期望目标就更优。
决策演示：若决策者给出红线"生态损失不超过 8"，即 $\varepsilon=8$，则可行集为 A、B、C、E，$\min f_1$ 取 C$(8,8)$——这正是加权和给不出、ε-约束给得出的方案。

#### 13. 最小可运行 MATLAB 示例

```matlab
% 方案表：每行一个方案，两列目标（均已统一为最小化）
F = [10 5; 9 6; 8 8; 7 9; 9 7];
n = size(F,1); dom = false(n,1);
for i = 1:n
    for j = 1:n
        if all(F(j,:) <= F(i,:)) && any(F(j,:) < F(i,:))
            dom(i) = true;          % i 被 j 支配
        end
    end
end
front = find(~dom);                 % 非支配集下标
fprintf('被支配标记: %s\n', mat2str(dom'));
fprintf('Pareto 前沿: %s\n', mat2str(front'));
% 加权和扫描：看哪些方案在某个权重下最优
w = (0:0.02:1)';
V = w*F(:,1)' + (1-w)*F(:,2)';      % 行=权重，列=方案
[~,pick] = min(V, [], 2);
fprintf('被权重选中的方案: %s\n', mat2str(unique(pick)'));
```

输出应为：`front = [1 2 3 4]`（E 被 B 支配），`unique(pick) = [1 2 4]`（C 从不被选中）。
改 `w` 的步长或给 `F` 换一组数，可直接观察"支撑点"与"非支撑点"的差别。
若方案是连续变量，把内层循环换成 `paretosearch` 或 `gamultiobj`，其余判断逻辑不变。

#### 14. 真实美赛 A 题应用

2020 年 A 题（渔业）把"经济收益"与"资源可持续"作为两个方向相反的目标：短期提高努力量能增收益，却压低存量，正是 Pareto 框架的典型场景，前沿上每个点对应一套捕捞与政策参数。
2019 年 A 题（龙）的"生态影响"与"社区需求"同理：保护强度与当地生计此消彼长，只报一个"最优强度"无法回答"代价是什么"。
写法要点：这类题要的是给决策者的建议，因此必须交前沿 + 挑选规则 + 敏感性，而不是一个孤立的最优点。

#### 15. 常见错误

1. 不做方向统一，把"越大越好"的目标直接塞进最小化。
2. 不归一化就加权：亿元与吨相加，权重被单位吞掉。
3. 只跑加权和就宣布"这就是全部前沿"——非凸段被系统性丢掉。
4. 把前沿当"最优解"直接交卷，不写挑选规则。
5. 支配判定漏掉"至少一个严格"，把所有并列点都判成被支配。
6. 遗传算法只跑一次、不换种子、不报超体积，前沿其实是随机产物。
7. ε 扫描步长太粗漏掉窄的非凸段；ε 取到 $\min f_2$ 以下导致不可行却写成"无解"。
8. 目标里塞进互相等价的指标（收益与利润），造成虚假的多目标。
9. 忘记报告决策变量，前沿上的点没有对应 $\mathbf x$，方案无法执行。

#### 16. 学习要求（验收）

- 能对 5 个方案的离散例子手判支配关系，指出谁被谁支配、前沿是哪几个。
- 能说清加权和为什么取不到非凸段，并写出一个具体反例（点或权重不等式）。
- 能写出 ε-约束的 KKT 条件并解释 $\lambda$ 的含义与符号。
- 给定一条前沿，能按"权重来源、红线约束、膝点"三种规则各给一个推荐并说明代价。

### 4.9 系统动力学

**定位**：P1｜难度 ★★☆｜3 h｜前置：3.2、反馈｜政策、资源管理、存量流量。

#### 1. 定义

系统动力学（system dynamics, SD）用**存量**（stock）与**流量**（flow）描述系统，核心是守恒方程
$$
\frac{\mathrm d\mathbf S}{\mathrm dt}=\mathbf F_{\text{in}}(\mathbf S,\mathbf u,t)-\mathbf F_{\text{out}}(\mathbf S,\mathbf u,t),\qquad \mathbf S(t_0)=\mathbf S_0,
$$
已知：存量清单、每条流量依赖存量与参数的方式（机制来自题面或数据拟合）、控制变量 $\mathbf u$；未知：轨迹 $\mathbf S(t)$、平衡点、稳定性、临界参数值。
四个要素：存量、流量、反馈回路（feedback loop）、延迟（delay）。
它不是"画框图讲故事"：可检验的输出是 ODE 的解、平衡点位置、稳定性判据与临界值。

#### 2. 核心思想与大白话

笨办法是把年度数据画成折线，"看出"拐点再外推；没有机制，政策一改，曲线立刻失效。
系统动力学要求先写清"水怎么流"：浴缸里的水是存量，水龙头与下水口是流量。存量只能靠流量累积，任何"凭空变化"的假设都会被一眼识破。
类比一（浴缸）：水位不会自己跳，只会被进出水累积；把下水口接一个"水越深流得越快"的阀，就得到负反馈与自稳定——这正是 logistic 方程。
类比二（空调）：室温偏高，制冷量加大，温度回落，偏差被压小；绕回路一圈符号为负，负反馈（调节型）；若增益过大或测温滞后，温度就会来回过冲。
类比三（淋浴）：拧了龙头要等几秒水温才变，你嫌慢就再拧，结果忽冷忽热——这就是延迟 τ 引起的振荡。
所以判据很硬：把模型写成方程，找出平衡点，算它稳不稳、临界值在哪，再谈政策。

#### 3. 数学基础

- **守恒与积分形式**：$S(t)=S(0)+\int_0^t(F_{\text{in}}-F_{\text{out}})\,\mathrm dt$，是离散仿真的直接依据。
- **欧拉与 Runge–Kutta 离散**：$\Delta S\approx(F_{\text{in}}-F_{\text{out}})\Delta t$，$\Delta t$ 太大时数值解会自己振荡。
- **平衡点**：$f(\mathbf S^*)=\mathbf 0$ 的物理可行根（$S^*\ge0$）。
- **线性化与雅可比**：$\mathbf J=\partial\mathbf f/\partial\mathbf S\vert_{\mathbf S^*}$，特征值实部决定稳定性。
- **回路极性**：绕回路一周的负号个数；符号增益为正 → 正反馈（放大），为负 → 负反馈（调节）。
- **时间常数**：$\tau=1/\vert\lambda\vert$，偏离平衡后衰减到 $1/e$ 所需时间。
- **Routh–Hurwitz**：二阶系统 $\lambda^2-\mathrm{tr}\mathbf J\,\lambda+\det\mathbf J=0$ 两根实部为负 ⟺ $\mathrm{tr}\mathbf J<0$ 且 $\det\mathbf J>0$。
- **延迟方程**：$\dot S=f(S(t-\tau))$，特征方程含 $e^{-\lambda\tau}$，可发生 Hopf 分岔（稳定平衡 → 周期振荡）。
- **量纲**：$\dot S$ 的单位必须是"$S$ 的单位 / 时间"；单位对不上通常意味着漏了参数或写错了机制。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
| --- | --- | --- | --- | --- |
| 边界内守恒 | 存量只通过已列出的流量改变 | 物质/能量守恒 | 漏一条流（偷捕、迁出），长期平衡偏高 20%—50% | 对每个存量算"流入−流出−$\mathrm dS/\mathrm dt$"残差，应近 0 |
| 流量只依赖存量与已知参数 | 无未建模的外部驱动 | 机制可解释 | 政策期外部条件突变，轨迹系统性偏离 | 分段拟合，看参数是否漂移 |
| 参数在预测期内恒定 | $r,K,q$ 不随时间变 | 稳态假设 | 参数漂移 30% 时平衡点位置可差一倍量级 | 滚动估计参数并画置信带 |
| 空间均匀（well-mixed） | 不区分位置 | 单室模型 | 局部过捕被空间平均掩盖，低估崩溃风险 | 拆成子区域建模对比 |
| 时间连续、无随机 | 忽略噪声与灾变 | 确定性 ODE | 小种群（<50）的随机灭绝风险被完全忽略 | 加随机项做 Monte Carlo |
| 延迟可忽略或已知 | 反馈即时生效 | 简化 | 忽略 τ 会把振荡型系统误判成单调稳定 | 对 τ 扫敏感性并求临界值 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
| --- | --- | --- | --- | --- |
| $S$（或 $N$） | 存量：种群、水量、资本 | 个、m$^3$、亿元 | 由题给，如 $N\sim10^3$ | 系统的"记忆"，只能被流量改变 |
| $F_{\text{in}},F_{\text{out}}$ | 流入、流出 | 存量的单位/时间 | 与 $S$ 同量级/时间 | 改变存量的唯一途径 |
| $\mathbf u$（或 $E$） | 控制/政策变量 | 努力量（艘·年）、税率 | $E=100$ | 决策者可调的手段 |
| $r$ | 内禀增长率 | 1/年 | 0.1—0.5（鱼类） | 无竞争时的增长速度 |
| $K$ | 环境容纳量 | 与 $S$ 同 | 1000（本题） | 资源上限 |
| $q$ | 可捕系数 | 1/(努力量·年) | $10^{-3}$ | 单位努力量的捕获效率 |
| $\tau$ | 反馈延迟 | 年 | 0.5—10 | 政策或生理过程的滞后 |
| $\lambda$ | 特征值 | 1/年 | 实部可正可负 | 稳定性与振荡频率的载体 |

#### 6. 从基本原理推导核心公式

**第一步：存量与流量的积分关系。** 对任一存量，$\mathrm dS/\mathrm dt=F_{\text{in}}-F_{\text{out}}$ 积分得
$$
S(t)=S(0)+\int_0^t\big(F_{\text{in}}(s)-F_{\text{out}}(s)\big)\,\mathrm ds .
$$
量纲检查：右边两项都是"$S$ 单位"，左边也是；被积函数是"$S$ 单位/时间"乘 $\mathrm ds$，一致。

**第二步：离散化与数值稳定性。** 前向欧拉 $S_{n+1}=S_n+\Delta t\,f(S_n)$，在平衡点附近令 $e_n=S_n-S^*$，得 $e_{n+1}=(1+\Delta tf')e_n$。要收敛须
$$
\vert 1+\Delta t f'\vert<1\iff \Delta t<\frac{2}{\vert f'\vert},
$$
本题 $f'(N^*)=-0.2$ 给出 $\Delta t<10$ 年；取 $\Delta t=1$ 年余量充足，若取 $\Delta t=20$ 年则会算出虚假振荡。量纲：$1/\vert f'\vert$ 是时间，与 $\Delta t$ 同量纲。

**第三步：一维线性化与时间常数。** 令 $S=S^*+\delta$，泰勒展开保留一阶：$\dot\delta=f'(S^*)\delta$，故 $\delta(t)=\delta_0\mathrm e^{f'\cdot t}$。$f'<0$ 稳定，衰减时间常数 $1/\vert f'\vert$；$f'>0$ 指数发散（正反馈失控）。

**第四步：二维的 Routh–Hurwitz 判据。** 线性化 $\dot{\boldsymbol\delta}=\mathbf J\boldsymbol\delta$，特征方程 $\lambda^2-\mathrm{tr}\mathbf J\,\lambda+\det\mathbf J=0$。两根实部为负 ⟺ 系数同号 ⟺
$$
\mathrm{tr}\mathbf J<0,\qquad \det\mathbf J>0 .
$$
若判别式 $(\mathrm{tr}\mathbf J)^2-4\det\mathbf J<0$，根为共轭复数，轨迹是阻尼螺旋，角频率 $\omega=\sqrt{\det\mathbf J-(\mathrm{tr}\mathbf J)^2/4}$。量纲：$\mathrm{tr}\mathbf J$ 与 $\det\mathbf J$ 的不同次幂出现在同一多项式里，说明变量必须先无量纲化，否则方程本身不成立。

**第五步：回路极性与增益。** 把回路写成乘法链 $S\to x_1\to\cdots\to S$，回路增益 $g=\prod_i\partial(\text{下一环节})/\partial(\text{上一环节})$。$g>0$ 时偏离被同向放大（正反馈），$g<0$ 时被反向压制（负反馈）；一维情形 $g=f'(S^*)$，与第三步一致。

**第六步：延迟如何制造振荡。** 取 $\dot N=rN\big(1-N(t-\tau)/K\big)-qEN$，在正平衡 $N^*$ 处线性化：由 $r(1-N^*/K)=qE$ 得
$$
\dot\delta=-a\,\delta-b\,\delta(t-\tau),\qquad a=qE,\quad b=r-qE=\frac{rN^*}{K}.
$$
令 $\delta=\mathrm e^{\lambda t}$，得特征方程 $\lambda+a+b\mathrm e^{-\lambda\tau}=0$。设 $\lambda=\mathrm i\omega$ 代入并分开实虚部：$a+b\cos\omega\tau=0$、$\omega=b\sin\omega\tau$，于是
$$
\cos\omega\tau=-\frac{a}{b},\qquad \omega=\sqrt{b^2-a^2},\qquad \tau_c=\frac{\arccos(-a/b)}{\sqrt{b^2-a^2}},
$$
存在条件 $b>a$，即 $qE<r/2$（$E<E_{\text{MSY}}$）。纯 logistic（$E=0$）时 $a=0$、$b=r$，得 $\tau_c=\pi/(2r)=5.236$ 年、周期 $2\pi/r\approx20.9$ 年；本题 $E=100$ 时 $a=0.1$、$b=0.2$，$\tau_c\approx12.09$ 年、周期约 36.3 年——捕捞本身提供阻尼，把临界延迟推后了。
量纲：$r\tau$ 无量纲，$\tau_c$ 的量纲是时间，两处都对。直觉：响应越快（$r$ 大）、反馈越迟（$\tau$ 大），越容易冲过头；$\tau>\tau_c$ 后平衡失稳，出现 Hopf 振荡。

**第七步：极限检查。** $\tau\to0$ 时 $\mathrm e^{-\lambda\tau}\to1$，回到无延迟 logistic，单调趋近 $K$；$r\to0$（系统极慢）时 $\tau_c\to\infty$，永不振荡；$a\ge b$（过度捕捞到 $E\ge E_{\text{MSY}}$）时不存在 $\tau_c$，系统要么单调趋稳要么直接崩溃。

#### 7. 求解方法

1. 画存量-流量图，给每条流标注依赖关系与单位。
2. 写 ODE，逐项核对量纲——这是最省时间的除错手段。
3. 求平衡点：解 $f(\mathbf S)=0$，只留 $S\ge0$ 的根，并检查"正平衡不存在"的情形。
4. 线性化求 $\mathbf J$，用 $f'$（一维）或 tr/det（二维）判稳定性。
5. 数值积分：`ode45` 用于非刚性，`ode15s` 用于刚性/长时间；用 $\Delta t<2/\vert f'\vert$ 检查步长。
6. 延迟项用 `dde23`，扫描 $\tau$ 找临界值。
7. 参数扫描与敏感性（单参数 $\pm50\%$、Latin 超立方），报告结论是否翻转。
8. 校准与验证：与历史数据比 RMSE，检查守恒残差，检查极端参数下是否出现负存量。

#### 8. 适用场景

**该用**：政策与资源管理（捕捞配额、碳税、水资源调度）；有明确反馈与积累的题。
也适合回答"短期有效、长期失效"或"临界点在哪"这类问题；需要情景对比而非单点预测时同样是首选。
**不该用**：纯数据预测（用 4.10）；空间异质性是主角（用 4.7 或 PDE）；个体异质性关键（用 ABM）；一次性事件、无反馈链的静态分配问题。

#### 9. 局限

- 参数多且难辨识：5 个参数配 20 个数据点时存在大量等效组合，预测区间很宽。
- 结构不确定性：漏一条流就系统性偏移；漏掉自然死亡会让 $K$ 的估计偏大 20% 以上。
- 无空间、无个体差异、无随机性；小种群（<50）的随机灭绝风险被完全忽略。
- 容易滑向"画图讲故事"：框图漂亮但没有校准、没有稳定性分析。
- 长期外推不可靠：政策参数在十年尺度上很少恒定。
- 验证难：历史数据往往只覆盖一种政策，模型在未观测政策下的行为难以证伪。

#### 10. 改进

1. 量纲检查 + 守恒残差表（零成本，先做）。
2. 参数辨识与置信带（最小二乘 + bootstrap），报告参数不确定性而不是只报点估计。
3. 加入延迟与随机项（`dde23` / SDE），用 Monte Carlo 给分布而非单条轨迹。
4. 空间或年龄结构扩展（多斑块、Leslie 矩阵），代价是参数翻倍。
5. 与优化耦合：政策参数作决策变量，做最优控制或多目标前沿（见 4.8）。

#### 11. 与其他模型的组合

- 与 3.2 微分方程：系统动力学就是"带语义、带反馈分析"的 ODE 建模法，两者共用数值解法。
- 与 4.8 多目标优化：政策杠杆作决策变量，"长期存量"与"累计收益"作两个目标。
- 与 4.10 时间序列：用 ARIMA 拟合外生驱动或检验残差自相关，残差有自相关说明漏了状态变量。
- 与 4.7 空间插值：初始场或局部参数的空间分布可由 Kriging 补齐。
- 与第五章敏感性分析：参数扫描是系统动力学报告的必备一节。

#### 12. 手算算例

模型 $\dot N=rN(1-N/K)-qEN$，$r=0.3$/年，$K=1000$，$q=10^{-3}$，$E=100$。
求平衡：$rN(1-N/K)=qEN$，$N\neq0$ 时 $1-N/K=qE/r=0.1/0.3=1/3$，故 $N^*=K(1-qE/r)=666.7$。
稳定性：$f(N)=rN(1-N/K)-qEN$，$f'(N)=r(1-2N/K)-qE$，代入 $N^*$ 得 $f'(N^*)=-rN^*/K=-0.3\times666.7/1000=-0.2<0$，稳定，时间常数 $1/0.2=5$ 年。
临界努力量：正平衡存在要求 $qE<r$，即 $E<E_c=r/q=300$；$E=400$ 时 $qE/r=4/3>1$，$N^*=1000(1-4/3)<0$，生物上不可行，种群崩溃。
最大可持续产量：$E_{\text{MSY}}=r/(2q)=150$，$N_{\text{MSY}}=K/2=500$，$\text{MSY}=rK/4=75$（单位：个体/年）。
延迟检查：$E=0$ 时 $\tau_c=5.24$ 年；$E=100$ 时 $\tau_c\approx12.1$ 年。若政策反馈滞后 3—5 年，仍在稳定侧但已有阻尼振荡，因此实际配额应比 $E_c=300$ 更保守——这正是原文"延迟要求留余量"的量化版本。

#### 13. 最小可运行 MATLAB 示例

```matlab
r = 0.3; K = 1000; q = 1e-3;
for E = [0 100 150 300 400]
    f = @(N) r*N.*(1 - N/K) - q*E*N;
    Ns = roots([-r/K, r - q*E, 0]);      % 平衡方程 r*N*(1-N/K)-q*E*N = 0
    Ns = sort(Ns(Ns > 0));
    if isempty(Ns)
        fprintf('E = %3d : 无生物可行正平衡\n', E);
    else
        fprintf('E = %3d : 正平衡 N* = %s\n', E, mat2str(Ns', 4));
    end
end
E = 100; f = @(N) r*N*(1 - N/K) - q*E*N;
Nst = 666.667; dp = 1e-6;
fp = (f(Nst+dp) - f(Nst-dp)) / (2*dp);   % 数值求 f'(N*)
fprintf('f''(N*) = %.3f, 时间常数 = %.2f 年\n', fp, -1/fp);
[t, N] = ode45(@(t, N) f(N), [0 120], 100);
plot(t, N); xlabel('t / 年'); ylabel('N'); ylim([0 800]);
```

输出应为：$E=0$ 时 $N^*=1000$；$E=100$ 时 $666.7$；$E=150$ 时 $500$；$E=300$ 时只剩 0；$E=400$ 时无正平衡。
`f'(N*) = -0.200`、时间常数 5 年，与手算一致。
轨迹从 $N=100$ 单调上升到 666.7，可用来对照"稳定但慢"的直觉。
把 `-q*E*N` 换成含延迟的项需改用 `dde23`，平衡与稳定性判据不变。

#### 14. 真实美赛 A 题应用

2020 年 A 题（渔业）：把捕捞努力量/配额作政策变量，用存量-流量结构预测长期存量与收益，核心结论是可持续努力量存在上界。
2023 年 A 题（群落长期存续）：多种群相互作用的长期存续，本质是平衡点是否存在、是否稳定，正是系统动力学的主场。
写法要点：先用存量-流量讲清机制与反馈，再给平衡点、稳定性与临界值，最后做政策情景对比；只有框图没有方程与稳定性分析会被判"讲故事"。

#### 15. 常见错误

1. 把流量直接写成存量（"今年产量=今年需求"），漏掉累积。
2. 回路极性数错符号，把正反馈标成负反馈，结论从"崩溃"翻成"稳定"。
3. 忽略延迟，把振荡型系统报成单调收敛。
4. 时间步长过大（$\Delta t>2/\vert f'\vert$），数值振荡被误当成真实振荡。
5. 只取正根，忘了检查 $E>E_c$ 时正平衡根本不存在。
6. 参数没有来源（"取 $r=0.5$ 因为好算"），也没做敏感性。
7. 单位混用（年与月、个与千个），量纲检查本可以发现。
8. 只交框图、不交方程、不做数据校准。
9. 把短期最优政策当长期最优：$E_{\text{MSY}}$ 在有延迟时未必稳定。

#### 16. 学习要求（验收）

- 能对一个存量-流量图写出 ODE 并逐项核对量纲。
- 能给出一维模型的平衡点、稳定性判据，并解释时间常数的物理意义。
- 能用 tr/det 判二阶系统稳定性；能说明延迟为什么导致振荡并估算临界 $\tau$。
- 能对渔业例子手算 $N^*$、$E_c$、MSY，并写出对应的政策建议。

### 4.10 时间序列

**定位**：P1｜难度 ★★★｜5 h｜前置：3.8｜有历史序列时的预测与输入。

#### 1. 定义

时间序列（time series）是按等间隔 $t=1,\dots,T$ 采样的一列观测 $y_t$。ARIMA$(p,d,q)$ 的标准形式为
$$
\Phi(B)(1-B)^d y_t=c+\Theta(B)\varepsilon_t,\qquad \varepsilon_t\sim\text{WN}(0,\sigma^2),
$$
其中 $B$ 是滞后算子（$By_t=y_{t-1}$），$\Phi(B)=1-\phi_1B-\cdots-\phi_pB^p$ 是自回归（autoregressive, AR）多项式，$\Theta(B)=1+\theta_1B+\cdots+\theta_qB^q$ 是滑动平均（moving average, MA）多项式，$c$ 是常数/漂移项。
已知：观测序列（可能含趋势、季节、异常点）；未知：阶数 $(p,d,q)$、参数 $(\phi,\theta,c,\sigma)$，以及 $h$ 步预测与预测区间。

#### 2. 核心思想与大白话

笨办法是拿最后两个点连线外推，或者拟合一条多项式曲线穿过历史点：前者对噪声过敏，后者在任何转折处都失效，两者都给不出"有多不确定"。
ARIMA 的出发点是：序列里可预测的部分来自它自己的过去（自相关）与过去的意外（残差）；先用差分把趋势与季节剥掉，剩下的平稳部分用少数几个参数刻画。
类比一（惯性）：今天气温与昨天高度相关，最好的短预报常常是"昨天 + 最近的平均变化"。
类比二（回味）：MA 部分像"昨天的意外今天还留一点余味"，一次冲击不会立刻消失，而是按 $\theta$ 衰减。
类比三（剥洋葱）：差分把"一直在涨"的趋势剥掉，让序列回到围绕常数波动；剥多了（过差分）会把有用信号也削掉，区间反而变宽。
定阶靠两条"尾巴"：ACF 拖尾而 PACF 截尾 → AR；ACF 截尾而 PACF 拖尾 → MA。
但看图只是初判，最后要靠残差白噪声检验与回测误差说话。

#### 3. 数学基础

- **滞后算子 $B$**：把差分写成多项式，方程才能因式分解、判根。
- **平稳性（stationarity）**：均值、方差、自协方差不随时间变；AR 部分的根全在单位圆外。
- **可逆性（invertibility）**：MA 部分的根全在单位圆外，保证残差能由观测反推。
- **自协方差与 ACF**：$\rho_k=\gamma_k/\gamma_0$，识别 MA 阶与季节周期。
- **PACF**：剔除中间滞后后的净相关，识别 AR 阶。
- **差分算子 $(1-B)^d$**：$d=1$ 去线性趋势，$d=2$ 去"增长在变"的趋势；$d$ 过大会放大噪声。
- **Yule–Walker 方程**：矩估计，$\boldsymbol\rho=\mathbf R\boldsymbol\phi$，解出 $\boldsymbol\phi$。
- **条件似然（CSS）与极大似然（MLE）**：前者等价于最小二乘，后者给出更小的方差估计偏差。
- **ψ 权重**：把冲击传播到 $h$ 步之后，决定预测误差方差与区间宽度。
- **AIC/BIC**：定阶时在拟合优度与复杂度之间取舍。

#### 4. 基本假设

| 假设 | 含义 | 依据 | 被违反的后果 | 如何检验 |
| --- | --- | --- | --- | --- |
| 等间隔采样、无缺失 | 时间步长一致 | 滞后与差分按等间隔定义 | 时间戳错位时 ACF 估计偏差可达 30%，定阶全错 | 检查时间戳间隔与缺失率 |
| 差分后平稳 | 均值方差为常数 | ARMA 理论的前提 | 仍含单位根时 t 统计量失真，预测系统性漂移 | ADF/KPSS 检验，看差分后 ACF 是否快速衰减 |
| 残差为白噪声 | 无剩余自相关 | MLE 有效性与区间正确性 | 残差有自相关说明漏了状态/季节，区间偏窄 | Ljung–Box 检验，$p>0.05$ 才算过 |
| 参数在预测期内不变 | 结构稳定 | 估计的一致性 | 政策或结构突变后预测系统性偏移 | 滚动回测，看误差是否随时间上升 |
| 线性、可加 | ARIMA 是线性模型 | 计算可行 | 阈值、饱和、爆发式增长被抹平 | 残差对拟合值作图，找非线性结构 |
| 季节周期 $m$ 已知且稳定 | SARIMA 用固定 $m$ | 季节机制 | $m$ 设错时季节差分无效并引入虚假周期 | 看季节 ACF 峰值是否落在 $m$ 的倍数 |

#### 5. 变量与参数

| 符号 | 含义 | 单位 | 典型值/量级 | 物理意义 |
| --- | --- | --- | --- | --- |
| $y_t$ | 观测值 | 与题一致（°C、人次、亿元） | 50—500 期 | 被建模的序列 |
| $p$ | AR 阶 | —（无量纲） | 0—3 | 用过去多少期解释现在 |
| $d$ | 差分阶 | — | 0—2 | 去趋势的力度 |
| $q$ | MA 阶 | — | 0—3 | 冲击"余味"的长度 |
| $\phi$ | AR 系数 | — | 0.184（本题） | 惯性强度，越接近 1 越持久 |
| $c$ | 常数/漂移 | 每期增量（与 $y$ 同） | 0.241 | 长期平均每期变化 |
| $\sigma$ | 残差标准差 | 与 $y$ 同 | 由数据估得 | 预测区间的基准宽度 |
| $h$ | 预测步长 | 期 | 1—12 | 外推距离 |
| $m$ | 季节周期 | 期 | 12（月）、24（时） | 季节重复长度 |

#### 6. 从基本原理推导核心公式

**第一步：AR(1) 的平稳性与 ACF。** 对 $y_t=c+\phi y_{t-1}+\varepsilon_t$ 反复代入得 $y_t=\mu+\sum_{j\ge0}\phi^j\varepsilon_{t-j}$，$\mu=c/(1-\phi)$。几何级数收敛要求 $\vert\phi\vert<1$（$\phi=1$ 是单位根，方差随时间发散）。于是
$$
\gamma_0=\frac{\sigma^2}{1-\phi^2},\qquad \gamma_k=\phi^k\gamma_0,\qquad \rho_k=\phi^k ,
$$
即 ACF 指数拖尾、永不截尾——这是"AR 看 PACF 截尾"的根源。

**第二步：差分算子与过差分代价。** $(1-B)y_t=y_t-y_{t-1}$，$(1-B)^2y_t=y_t-2y_{t-1}+y_{t-2}$。若 $y_t=a+bt+u_t$，则 $(1-B)y_t=b+(1-B)u_t$ 已平稳，再差一次会把 $u$ 也差掉：$\mathrm{Var}((1-B)u)=2\gamma_0(1-\rho_1)$，$\rho_1$ 接近 1 时该值很小，但若 $u$ 已白噪声（$\rho_1=0$）方差直接翻倍——这就是"过差分使区间变宽"的量化解释。量纲：差分不改变单位，$\Delta y$ 与 $y$ 同单位。

**第三步：Yule–Walker 矩估计。** 用 $y_{t-k}$ 乘方程两边取期望，得 $\gamma_k=\phi_1\gamma_{k-1}+\cdots+\phi_p\gamma_{k-p}$（$k\ge1$），除以 $\gamma_0$：
$$
\boldsymbol\rho=\mathbf R\boldsymbol\phi,\qquad \mathbf R_{ij}=\rho_{\vert i-j\vert}\ \text{（Toeplitz）},
$$
故 $\hat{\boldsymbol\phi}=\mathbf R^{-1}\boldsymbol\rho$。$\mathbf R$ 正定（自协方差阵半正定且 $\gamma_0>0$），所以总能求逆。$\phi$ 无量纲，与量纲无关。

**第四步：条件似然与 MLE。** 给定初值，$\varepsilon_t=y_t-c-\sum_i\phi_i y_{t-i}$，条件似然
$$
L\propto\sigma^{-n}\exp\!\Big(-\frac{1}{2\sigma^2}\sum_{t}\varepsilon_t^2\Big),
$$
故对 $(\phi,c)$ 最大化 $L$ 等价于最小化残差平方和（CSS）；对 $\sigma^2$ 求导得 $\hat\sigma^2=\sum\hat\varepsilon_t^2/n$。

**第五步：$h$ 步预测与误差方差。** 用 MA$(\infty)$ 表示 $y_{T+h}=\mu+\sum_{j\ge0}\psi_j\varepsilon_{T+h-j}$，取条件期望时未来冲击的期望为 0：
$$
\hat y_{T+h}=\mu+\sum_{j\ge h}\psi_j\varepsilon_{T+h-j},\qquad e_h=\sum_{j=0}^{h-1}\psi_j\varepsilon_{T+h-j},
$$
$$
\mathrm{Var}(e_h)=\sigma^2\sum_{j=0}^{h-1}\psi_j^2 .
$$
方差随 $h$ 单调增，所以区间必然变宽；AR(1) 时 $\psi_j=\phi^j$，$\mathrm{Var}(e_h)=\sigma^2(1-\phi^{2h})/(1-\phi^2)\to\sigma^2/(1-\phi^2)$，即收敛到无条件方差。

**第六步：ARIMA(1,1,0) 的具体化。** 模型即 $\Delta y_t=c+\phi\Delta y_{t-1}+\varepsilon_t$，等价于对一阶差分建 AR(1)。代入本题 $\hat\phi=0.184$、$\hat c=0.241$：单步预测 $\hat y_{T+1}=y_T+\hat c+\hat\phi\Delta y_T$；增量长期收敛到 $c/(1-\phi)=0.241/0.816=0.295$/期，收敛半衰期 $\ln0.5/\ln0.184\approx0.41$ 期（几乎两步到位）。
把 $\Delta y$ 的 MA 权重累加得 $y$ 的 $\psi_j=(1-\phi^{j+1})/(1-\phi)$，即 $\psi=(1,1.184,1.218,1.224,\dots)$，故 $\sum_{j<h}\psi_j^2$ 在 $h=1,3,12$ 时为 $1$、$3.89$、$17.4$。

**第七步：ACF/PACF 判据的由来。** MA$(q)$ 的 $\rho_k=0$（$k>q$），ACF 在 $q$ 处截尾；AR$(p)$ 的 $\rho_k$ 拖尾，但 PACF 在 $p$ 处截尾——PACF 的第 $k$ 个值是把 $y_t$ 对前 $k$ 期回归的最后一个系数，AR$(p)$ 只需 $p$ 期，多出来的滞后系数为 0。实际数据常常两者都拖尾，此时只能靠 AIC/BIC 加残差检验定阶。

**第八步：极限检查。** $\phi\to0$：退化为白噪声加常数，预测就是均值；$\phi\to1$：接近随机游走，$h$ 步方差 $\approx h\sigma^2$ 线性增长；$h\to\infty$：ARMA 的预测收敛到 $\mu$，区间宽度收敛到无条件标准差。这三条可直接用来检验软件输出是否合理。

#### 7. 求解方法

1. 画图与分解：先看趋势、季节、异常点；异常点不处理会污染整条 ACF。
2. 平稳化：ADF/KPSS 定 $d$；季节用 $D$，$m$ 从季节 ACF 峰位读。
3. 定阶：ACF/PACF 初判，再对 $(p,q)\in[0,3]^2$ 网格搜 AIC/BIC。
4. 估计：CSS-ML；检查 $\vert\hat\phi\vert<1$、特征根在单位圆外。
5. 诊断：残差 ACF、Ljung–Box（滞后 $m$ 阶）、正态性。
6. 回测：留出最后 12—24 期做滚动预测，报 RMSE/MAPE，不要报样本内 $R^2$。
7. 预测：点预测 + 80%/95% 区间，并画在同一张图上。
8. 报告：$(p,d,q)$、系数与标准误、诊断 $p$ 值、回测误差四件套。

#### 8. 适用场景

**该用**：有 50 期以上等间隔历史；变化主要由自相关与惯性驱动。
需要短期预测、或要把预测当输入喂给别的模型时也是首选；还必须给不确定区间，这时 ARIMA 的区间是现成的。
**不该用**：样本不足 30 期；机制清楚可写方程（用 4.9）；强外生变量是主角（先回归再用 ARIMA 处理残差，即 ARIMAX）；长期外推；政策断点之后。

#### 9. 局限

- 线性可加：阈值、饱和、爆发式增长被抹平；遇到结构断点，误差可翻倍。
- 长程预测迅速退化：AR(1) 的 $h$ 步方差趋于无条件方差，$h>12$ 基本只能当"回到均值"。
- 定阶主观：ACF/PACF 都拖尾时不同准则给出不同 $(p,q)$，预测差 5%—10% 很常见。
- 单变量：外生信息要靠 ARIMAX 或先回归。
- 对异常值敏感：一个尖峰会抬高 $\hat\sigma$ 并扭曲 ACF。
- 纯曲线拟合不足以支撑 A 题结论，2026 年题目已明确排斥。

#### 10. 改进

1. 异常值处理 + 稳健估计（零成本，收益最大）。
2. 季节模型 SARIMA$(p,d,q)(P,D,Q)_m$，或 STL 分解后分别建模。
3. ARIMAX：把外生变量（政策、价格、气温）放进回归项。
4. 状态空间/ETS/结构时间序列，支持缺失值与多季节。
5. 组合预测：机制模型（4.9）给结构、ARIMA 给残差，按回测误差加权。

#### 11. 与其他模型的组合

- 与 3.8 回归与统计基础：ARIMA 可看作"误差带自相关"的回归，检验与估计思路一致。
- 与 4.9 系统动力学：机制给趋势与政策响应，ARIMA 处理残差；残差仍显著说明漏了状态变量。
- 与 4.8 多目标优化：预测的需求/价格作为目标输入，在区间两端各算一次前沿看结论是否翻转。
- 与 4.7 空间插值：每个测点先建时间模型，再对同一时刻做空间插值。
- 与第六章模型检验：回测与残差诊断是必写项，不能只报拟合图。

#### 12. 手算算例

合成"趋势 + 季节"序列拟 ARIMA(1,1,0)，得 $\hat\phi=0.184$、$\hat c=0.241$，12 步预测首值 $56.98$。
模型：$\Delta y_t=0.241+0.184\,\Delta y_{t-1}+\varepsilon_t$，即对一阶差分建 AR(1)。
若最后一期 $\Delta y_T=0.5$，则 $\Delta\hat y_{T+1}=0.241+0.184\times0.5=0.333$，$\hat y_{T+1}=y_T+0.333$；本题给出的首步预测为 $56.98$。
长期：每期增量按 $0.184$ 衰减并收敛到 $c/(1-\phi)=0.241/0.816=0.295$；半衰期 $\ln0.5/\ln0.184\approx0.41$ 期，两期后基本到位。
12 步：$\hat y_{T+12}\approx56.98+11\times0.295\approx60.2$（实际略高，因为首步增量通常高于长期漂移）。
区间：$\psi=(1,1.184,1.218,1.224,\dots)$ 给出 $\sum_{j<h}\psi_j^2$ 在 $h=1,3,12$ 时为 $1$、$3.89$、$17.4$，故 12 步预测标准差约为 1 步的 $\sqrt{17.4}\approx4.2$ 倍——区间必须画出来，否则读图的人会高估 12 步的可靠性。

#### 13. 最小可运行 MATLAB 示例

```matlab
y  = load('series.txt');            % 一列观测值，等间隔
Mdl = arima(1, 1, 0);               % ARIMA(1,1,0)：p=1, d=1, q=0
Est = estimate(Mdl, y);             % CSS-ML 估计
phi = Est.AR{1};  c = Est.Constant;  s = sqrt(Est.Variance);
fprintf('phi=%.3f  c=%.3f  sigma=%.3f\n', phi, c, s);
[yF, yCI] = forecast(Est, 12, 'Y0', y);
disp([yF yCI]);                     % 点预测 + 95% 区间
res = infer(Est, y);                % 残差
if exist('lbqtest', 'file')
    [~, p] = lbqtest(res, 'Lags', 12);
    fprintf('Ljung-Box p = %.3f\n', p);
end
plot(1:numel(y), y, 'k'); hold on;
plot(numel(y) + (1:12), yF, 'r-o');
```

前三行是"定阶—估计"，`forecast` 给区间，`infer` 取残差做白噪声检验（`lbqtest` 属于 Econometrics Toolbox，没有就手算残差 ACF）。
判对的标准：$\hat\phi$ 落在 $(-1,1)$ 内、残差 ACF 基本落在 $\pm2/\sqrt{T}$ 带内、区间随 $h$ 单调变宽。
回测时把 `'Y0'` 换成截断序列，逐点滚动，统计 RMSE 与 MAPE。

#### 14. 真实美赛 A 题应用

2020 年 A 题：海温序列的趋势与周期是鱼群分布变化的驱动，ARIMA 用来回答"未来几年海温大概在哪"，并给出区间。
2026 年题：用户负载的日内模式，SARIMA 可刻画"每天同一时刻"的重复结构，但题目明确不欢迎纯曲线拟合，必须与机制结合。
正确用法有三种：作输入（把预测喂给 4.9 的机制模型）、作基线（证明你的机制模型优于纯统计外推）、作诊断（残差仍有自相关说明机制漏了状态）。
写法：报回测误差与区间，不要只报拟合图。

#### 15. 常见错误

1. 不检验平稳性就套 ARMA，$\hat\phi$ 逼近 1 却当成稳定。
2. 过差分（$d=2$ 甚至 3），方差被人为放大、区间虚胖。
3. 用样本内 $R^2$ 报拟合优度，不给回测误差。
4. 只给点预测，不画预测区间。
5. 定阶只看 ACF/PACF 图，不做残差白噪声检验。
6. 忽略异常点，整条 ACF 被一个尖峰带偏。
7. 季节周期设错（月度数据用 $m=4$），季节差分无效。
8. 拿 ARIMA 外推 50 期还宣称"可靠"。
9. 把预测值当确定值塞进后续优化，不传区间与不确定性。

#### 16. 学习要求（验收）

- 能写出 ARIMA$(p,d,q)$ 的算子形式，并说明平稳性与可逆性条件各作用在哪个多项式上。
- 能由 ACF/PACF 的截尾/拖尾判断该用 AR 还是 MA，并说明理由。
- 能推导 AR(1) 的 ACF 与 $h$ 步预测方差，并解释区间为什么变宽。
- 能对一条真实序列完成"定阶—估计—诊断—回测"四步，并报告 RMSE/MAPE 与区间图。
