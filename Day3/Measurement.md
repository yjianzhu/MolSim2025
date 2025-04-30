# Monte-carlo simulation for more complex problems

## Event-chain Monte-carlo simulation
Non-Reversibility (in some variants): Some advanced ECMC methods (like Forward ECMC or Irreversible ECMC) intentionally break detailed balance but maintain global balance. This non-reversibility can lead to even faster convergence by suppressing random walk behavior and promoting directed exploration of the state space.

## Measurement
好的，我们来推导一下格林-久保关系（Green-Kubo relations）。这个关系是统计力学中非常重要的结果，它将宏观的输运系数（如电导率、热导率、粘度等）与微观粒子运动产生的流（如电流、热流、动量流等）在平衡态下的时间自关联函数联系起来。推导主要基于线性响应理论。

**核心思想：**
系统在外加微扰场的作用下偏离平衡态，产生宏观响应（例如产生电流）。线性响应理论假设微扰足够小，系统的响应与微扰成线性关系。格林-久保关系表明，描述这种线性响应的比例系数（即输运系数）可以通过计算平衡态下系统内部自发涨落的关联性质得到。

**推导步骤（以量子力学体系为例，经典体系类似）：**

1.  **系统哈密顿量与微扰：**
    考虑一个系统，其未受扰动的哈密顿量为 $H_0$。系统最初处于温度 $T$（$\beta = 1/(k_B T)$）的平衡态，其密度算符为：
    $$\rho_{eq} = \frac{e^{-\beta H_0}}{Z_0}$$
    其中 $Z_0 = \text{Tr}(e^{-\beta H_0})$ 是配分函数。

    从 $t = -\infty$ 开始，施加一个随时间变化的微扰。总哈密顿量变为：
    $$H(t) = H_0 + H'(t)$$
    通常，微扰可以写成 $H'(t) = -A F(t)$ 的形式，其中 $A$ 是系统的一个厄米算符（例如电偶极矩、粒子数密度等），$F(t)$ 是一个外部标量场（例如电场、化学势梯度等）。为了推导方便，我们假设在 $t \to -\infty$ 时 $F(t) \to 0$。

2.  **密度算符的演化（线性响应）：**
    系统的密度算符 $\rho(t)$ 遵循刘维尔方程（Liouville equation）：
    $$i\hbar \frac{d\rho(t)}{dt} = [H(t), \rho(t)] = [H_0 + H'(t), \rho(t)]$$
    为了求解 $\rho(t)$，我们通常转换到相互作用绘景（Interaction Picture）。定义相互作用绘景中的算符：
    $$O_I(t) = e^{i H_0 t / \hbar} O e^{-i H_0 t / \hbar}$$
    其中 $O$ 是薛定谔绘景中的算符。
    相互作用绘景中的密度算符 $\rho_I(t)$ 和微扰哈密顿量 $H'_I(t)$ 分别为：
    $$\rho_I(t) = e^{i H_0 t / \hbar} \rho(t) e^{-i H_0 t / \hbar}$$
    $$H'_I(t) = e^{i H_0 t / \hbar} H'(t) e^{-i H_0 t / \hbar} = -A_I(t) F(t)$$
    刘维尔方程在相互作用绘景中变为：
    $$i\hbar \frac{d\rho_I(t)}{dt} = [H'_I(t), \rho_I(t)]$$
    将此方程积分，得到：
    $$\rho_I(t) = \rho_I(-\infty) + \frac{1}{i\hbar} \int_{-\infty}^{t} [H'_I(t'), \rho_I(t')] dt'$$
    由于 $t \to -\infty$ 时系统处于平衡态 $\rho_{eq}$，且 $\rho_I(-\infty) = \rho_{eq}$ （因为 $H_0$ 与 $\rho_{eq}$ 对易），我们得到：
    $$\rho_I(t) = \rho_{eq} + \frac{1}{i\hbar} \int_{-\infty}^{t} [H'_I(t'), \rho_I(t')] dt'$$
    在线性响应近似下，我们只保留 $H'_I$ 的一阶项。这意味着在积分号内，可以用 $\rho_{eq}$ 替换 $\rho_I(t')$：
    $$\rho_I(t) \approx \rho_{eq} + \frac{1}{i\hbar} \int_{-\infty}^{t} [H'_I(t'), \rho_{eq}] dt'$$
    $$\rho_I(t) \approx \rho_{eq} - \frac{1}{i\hbar} \int_{-\infty}^{t} [A_I(t'), \rho_{eq}] F(t') dt'$$

3.  **可观测量的期望值：**
    我们关心另一个可观测量 $B$ 的期望值随时间的演化。$B$ 的期望值为：
    $$\langle B(t) \rangle = \text{Tr}(\rho(t) B) = \text{Tr}(\rho_I(t) B_I(t))$$
    代入 $\rho_I(t)$ 的近似表达式：
    $$\langle B(t) \rangle \approx \text{Tr}(\rho_{eq} B_I(t)) + \frac{1}{i\hbar} \int_{-\infty}^{t} \text{Tr}([H'_I(t'), \rho_{eq}] B_I(t)) dt'$$
    第一项 $\text{Tr}(\rho_{eq} B_I(t)) = \text{Tr}(\rho_{eq} B) = \langle B \rangle_{eq}$ 是 $B$ 的平衡期望值（假设 $B$ 本身不显含时间）。
    第二项是响应部分 $\Delta \langle B(t) \rangle = \langle B(t) \rangle - \langle B \rangle_{eq}$。利用迹的循环不变性 $\text{Tr}(XYZ) = \text{Tr}(ZXY)$：
    $$\Delta \langle B(t) \rangle \approx \frac{1}{i\hbar} \int_{-\infty}^{t} \text{Tr}(\rho_{eq} [B_I(t), H'_I(t')]) dt'$$
    $$\Delta \langle B(t) \rangle \approx -\frac{1}{i\hbar} \int_{-\infty}^{t} \text{Tr}(\rho_{eq} [B_I(t), A_I(t')]) F(t') dt'$$
    $$\Delta \langle B(t) \rangle = \int_{-\infty}^{t} \Phi_{BA}(t, t') F(t') dt'$$
    其中，响应函数 $\Phi_{BA}(t, t')$ 定义为：
    $$\Phi_{BA}(t, t') = -\frac{1}{i\hbar} \text{Tr}(\rho_{eq} [B_I(t), A_I(t')]) = -\frac{1}{i\hbar} \langle [B_I(t), A_I(t')] \rangle_{eq}$$
    由于平衡态是时间平移不变的，响应函数只依赖于时间差 $\tau = t - t'$：
    $$\Phi_{BA}(t - t') = -\frac{1}{i\hbar} \langle [B_I(t - t'), A_I(0)] \rangle_{eq}$$
    其中 $\langle \dots \rangle_{eq}$ 表示在平衡态 $\rho_{eq}$ 下的系综平均。
    注意：$A_I(0) = A$，$B_I(t-t')$ 是在相互作用绘景（或等效地，对于平衡态关联函数，在海森堡绘景）中演化 $t-t'$ 时间的算符 $B$。
    因此，响应可以写成卷积形式：
    $$\Delta \langle B(t) \rangle = \int_{-\infty}^{t} \Phi_{BA}(t - t') F(t') dt' = \int_{0}^{\infty} \Phi_{BA}(\tau) F(t - \tau) d\tau$$
    （令 $\tau = t - t'$）

4.  **响应函数与关联函数的关系：**
    格林-久保关系的关键在于将响应函数 $\Phi_{BA}(\tau)$ 与某个流 $J$ 的平衡态时间自关联函数联系起来。通常，我们感兴趣的响应 $B$ 本身就是一个流（例如电流密度 $J_e$），而微扰 $A$ 是对应的势场（例如矢量势 $A$ 或标量势 $\phi$ 的空间导数，导致电场 $E = -\nabla \phi - \partial A/\partial t$）。
    考虑一种常见情况，响应 $B$ 正是与微扰 $A$ 共轭的流 $J_A$，即 $J_A = \dot{A}_I(0) = \frac{1}{i\hbar}[A_I(0), H_0]$。（这里假设 $A$ 不显含时间）。
    我们考虑一种特殊但重要的响应函数，它与关联函数的时间导数有关。使用 Kubo 恒等式：
    $$[A, e^{-\beta H_0}] = -e^{-\beta H_0} \int_0^\beta d\lambda e^{\lambda H_0} [A, H_0] e^{-\lambda H_0}$$
    可以推导出 (细节略)：
    $$\Phi_{BA}(\tau) = \beta \langle \dot{A}_I(0) ; B_I(\tau) \rangle_{eq} = \beta \int_0^1 d\lambda \langle e^{\lambda \beta H_0} \dot{A}_I(0) e^{-\lambda \beta H_0} B_I(\tau) \rangle_{eq}$$
    其中 $\dot{A}_I(0) = \frac{1}{i\hbar}[A_I(0), H_0]$。这个形式称为 Kubo 变换或正则关联。
    在经典极限下 ($\hbar \to 0$, $[,]/(i\hbar) \to \{,\}_{PB}$ 泊松括号)，响应函数变为：
    $$\Phi_{BA}(\tau) = -\beta \langle \{A(0), B(\tau)\}_{PB} \rangle_{eq}$$
    并且可以证明 (通过分部积分等技巧)：
    $$\Phi_{BA}(\tau) = -\beta \frac{d}{d\tau} \langle A(0) B(\tau) \rangle_{eq}$$
    这个关系（或其量子版本）是连接响应函数和平衡态关联函数的桥梁。

5.  **输运系数：**
    输运系数通常定义为稳态（或低频）极限下的响应。例如，考虑直流电导率 $\sigma$。此时外场 $F(t)$ 是一个恒定的电场 $E$（$F(t)=E$ for $t>0$），响应 $B$ 是电流密度 $J_e$。相应的算符 $A$ 是总偶极矩 $P$（使得 $H' = -P \cdot E$）。那么电流 $J_e = \dot{P}$。
    稳态电流为 $t \to \infty$ 时的响应：
    $$\Delta \langle J_e \rangle_{steady} = \lim_{t\to\infty} \int_0^t \Phi_{J_e P}(t-t') E dt' = E \int_0^\infty \Phi_{J_e P}(\tau) d\tau$$
    电导率张量 $\sigma_{\alpha\beta}$ 定义为 $J_{e,\alpha} = \sum_\beta \sigma_{\alpha\beta} E_\beta$。这里假设各向同性， $J_e = \sigma E$。
    $$\sigma = \int_0^\infty \Phi_{J_e P}(\tau) d\tau$$
    现在使用第 4 步的关系。如果 $B = J_A = \dot{A}$, 那么响应函数 $\Phi_{J_A A}(\tau)$ 与 $\langle A(0) J_A(\tau) \rangle_{eq}$ 的导数相关。但更常用的是将 $J_A$ 视作基本流。
    假设 $B = J_B$ (某个流)，$A$ 是与 $J_B$ 共轭的量，使得 $H' = -A F$。
    输运系数 $L$ (例如 $\sigma, \kappa, \eta$) 通常与响应函数的积分有关：
    $$L \propto \int_0^\infty \Phi_{J_B A}(\tau) d\tau$$
    利用 $\Phi_{BA}(\tau) = -\beta \frac{d}{d\tau} \langle A(0) B(\tau) \rangle_{eq}$ （经典形式，量子形式类似但更复杂，涉及 Kubo 变换），或者直接从 $\Phi_{BA}(\tau) = -\frac{1}{i\hbar} \langle [B_I(\tau), A_I(0)] \rangle_{eq}$ 出发，经过一些推导（例如对 Onsager 回归假设的运用），可以得到最终的 Green-Kubo 公式形式：
    $$L \propto \int_0^\infty \langle J_B(0) J_B(\tau) \rangle_{eq} d\tau$$
    这里 $J_B$ 是与输运系数 $L$ 对应的微观流。这个形式表明输运系数由相应微观流在平衡态下的时间自关联函数的积分决定。因子 $\beta = 1/(k_B T)$ 和体积 $V$ 等通常包含在比例常数中。

**具体例子：**

* **电导率 (Electrical Conductivity, $\sigma$)：**
    $J$ 是电荷流密度。
    $$\sigma = \frac{1}{V k_B T} \int_0^\infty \langle J_e(0) \cdot J_e(t) \rangle_{eq} dt$$
    (严格来说是电导率张量，这里假设各向同性系统)
* **热导率 (Thermal Conductivity, $\kappa$)：**
    $J_q$ 是热流密度。
    $$\kappa = \frac{1}{V k_B T^2} \int_0^\infty \langle J_q(0) \cdot J_q(t) \rangle_{eq} dt$$
* **剪切粘度 (Shear Viscosity, $\eta$)：**
    $P_{xy}$ 是应力张量的非对角元（动量流）。
    $$\eta = \frac{1}{V k_B T} \int_0^\infty \langle P_{xy}(0) P_{xy}(t) \rangle_{eq} dt$$
* **扩散系数 (Diffusion Coefficient, D)：**
    $v_i(t)$ 是标记粒子 $i$ 的速度。
    $$D = \frac{1}{3} \int_0^\infty \langle \mathbf{v}_i(0) \cdot \mathbf{v}_i(t) \rangle_{eq} dt$$
    (这是通过速度自关联函数，也可以通过其他流推导)

**总结：**
格林-久保关系的推导依赖于线性响应理论，核心是将宏观系统对微小外扰的响应（由响应函数描述）与系统在没有外扰时（平衡态）内部微观量的自发涨落（由时间关联函数描述）联系起来。最终结果表明，输运系数正比于相应微观流的时间自关联函数的积分。这为从微观模拟（如分子动力学）计算宏观输运性质提供了理论基础。

## example2
g(r) 
用green-kubo关系来计算g(r)
radial distribution function
$$g(r) = \frac{1}{\rho} \left( \frac{N}{V} \right) \int_0^r \langle J(r) J(0) \rangle_{eq} dr$$
$$\rho = \frac{N}{V}$$

## example3
stucture factor

## Sampling observable quantities

### pressure
