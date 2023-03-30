# 以下由陆星欢（2021213053023）进行编写

## 2.3 最大似然估计

&emsp; &emsp;  $DeGroot$ 和 $Goel(1980)$ 考虑了一个基于一个破碎样本的双变量正态分布的相关系数的 $MLE$ 。换句话说，从分布中抽取 $n$ 对随机样本，但观察到的数据只是 $n$ 对的第一分量和另外的 $n$ 对的第二分量的未知排列，称为分布中的破碎随机样本。研究了两种方法，其中排列向量要么被视为固定的未知参数，要么是需要从观察到的可能性中积分的随机变量。\
&emsp; &emsp;他们发现在第一种方法下， $MLE$ 是在所有可能的样本相关性中可以计算出的最大或最小样本相关性，这作为一般方法是不合理的。类似地，将我们自己限制在 $A$ 和 $B$ 之间的一对一匹配的情况下，这种方法除了将 $\omega$ 作为基于链接数据的分析参数外，还作为需要估计的参数。例如，如果分析中涉及的变量独立于用于记录链接的关键变量的误差，那么 $\omega$ 的可能性可以从增益参数的可能性中进行因式分解。那么，以上对协议分区的讨论表明， $\omega$ 的 $MLE$ 在目前的情况下不太可能做得很好。例如，在 $|\Omega_0| = 1$ 的情况下，一个合理的关键变量误差模型可能很好地暗示了链接矩阵 $\omega_L$ 作为真实匹配矩阵的MLE，并且不会对链接误差进行有用的调整。或者，在 $|\Omega_0| > 1$ 的情况下，很可能存在多个 $MLE$ 的情况。\
&emsp; &emsp;对于第二种方法， $DeGroot$ 和 $Goel(1980)$ 考虑了对排列向量进行积分后观察到的似然性，以及由此推导出的轮廓似然。基于模拟在 $n = 5$ 的情况下，其中显式积分是可行的，他们发现轮廓的可能性普遍较差。此外，虽然综合似然看起来与完整样本似然非常相似，但对于其他样本，两者看起来非常不同。这仍然并不一定意味着破碎的样本 $MLE$ 是渐近不一致的。然而，要研究这一点，需要一种基于大破碎样本的 $MLE$ ，在这方面显式集成是不切实际的。作者认为，在这种情况下， $MLE$ 可以通过 $EM$ 算法获得，通过将置换向量视为缺失的数据。但是他们并没有检查这个 $MLE$ 的性能。\
&emsp; &emsp;下面我们用一个基于链接数据的回归分析的例子来研究这一点。不幸的是，基于支持 $EM$ 算法的缺失信息原理（ $Orchard$ 和 $Woodbury$ ， $1972$ ）， $MLE$ 可能是有偏差的和不一致的。\
&emsp; &emsp;考虑一对一匹配的设置，其中 $U_A = U_B = \phi$ 和 $n_A = n_B = n$，以及 $\omega$ 是一个排列矩阵，即 $\omega\omega^T = \omega^T\omega = I$ ，其中 $I$ 是单位矩阵。考虑线性回归模型:\
 $$ y_i = x_i^T\beta + \epsilon_i $$ \
&emsp; &emsp;其中， $\epsilon_i$ 是从 $N(0，\sigma^2)$得到的 $IID$ ， $\beta$ 是增益的回归系数的 $p$ 向量。让  $X_A=X
_nA \times k$ 相关的协变量，和 $(X_B，y_B)=(X_{_nB\times(p−k)}，y_{_nB\times1})$ 剩下的协变量和因变量与 $B$ ，这样给定真正的匹配矩阵 $\omega$ ，回归是基于 $y_M = \omega y_B$ 和分区设计矩阵 $X_M=\lbrack X_A:\omega X_B\rbrack$ 。\
&emsp; &emsp;观测数据由 $z=(G、X
_A、X_B、y_B)$ 组成，其中G由可用的比较数据组成。例如，在最小的情况下，我们有 $G =\omega_L$ 提供了一对一的链接。完整的数据为 $(\omega，z)$ ，另外还包含了未观察到的 $\omega$ 。基于 $z$ 的 $MLE$ 需要从完整数据的分布中积分出 $\omega$ 得到。通过 $EM$ 算法，可以将观测数据的期望完全数据对数似然条件最大化。完整数据的对数似然值和分数由\
 $$\ell m(\beta,\sigma^2;\omega,z) = -\frac n2 \log\sigma^2 - \frac12\sigma^2(y_M - X_M^T\beta)^T( y_M - X_M^T\beta)$$ \
 $$u_M(\beta;\omega,z) = \sigma^{-2}( X_M^Ty_M - X_M^TX_M \beta)\\
= \sigma^{-2}( {X_A^T\omega\brack X_B^T}y_B - \begin{bmatrix}
   X_A^TX_A^T & X_A^T\omega X_B \\
   X_B^T\omega^TX_A & X_B^TX_B
\end{bmatrix}\beta ) $$
$$u_M(\sigma^2;\omega,z) = -\frac n{\sigma^2} + \frac{( y_M - X_M\beta)^T(y_M - X_M\beta)} {\sigma^4}$$ \
&emsp; &emsp;为了目前的目的，我们现在进一步简化了情况，并假定已知给定$z$的条件分布，因此只需要估计增益的 $\beta$ 。可以很容易地验证，最大值 $E(l_M|z)$ 的 $MLE$ 与 $E[u(\beta;\omega,z)|z]$ 的解相同，以上的解由\
 $$\^\beta=E(X_M^TX_M|z)^{-1}E(X_M^Ty_M|z) \\
=\begin{bmatrix}
   X_A^TX_A^T & X_A^TE(\omega|z) X_B \\
   X_B^TE(\omega^T|z) X_A & X_B^TX_B
\end{bmatrix}^{-1}{X_A^TE(\omega|z)\brack X_B^T}y_B$$ \
 $$\^{\beta}_M = ({X_M^TX_M })^{-1}(X_M^Ty_M) = \begin{bmatrix}
   X_A^TX_A^T & X_A^T\omega X_B \\
   X_B^T\omega^TX_A & X_B^TX_B
\end{bmatrix}^{-1}{X_A^T\omega\brack X_B^T}y_B$$ \
其中 $\beta$ 的真实完整数据 $MLE$ 是由 $\^\beta$ 和 $\^\beta_M$ 提供的，除非在没有链接错误的情况下，否则 $\^\beta$ 和 $\^\beta_M$ 肯定会大概率有所不同。接下来，让 $z(Y)=(G，X_A，X_B)$ ，即不包括 $y_B$ 。当 $E(y_M|\omega,z_{(Y)})= X_M\beta$ ，我们有\
 $$E( y_M|z_{(Y)})=E(\omega^TE( y_M|\omega,z_{(Y)})|z_{(Y)})=E(\omega^TX_M\beta|z_{(Y)})\\
=E(\omega^T[X_A:\omega X_B]\beta|z_{(Y)})=[E(\omega^T|z_{(Y)})X_A:X_B]\beta$$ 

 $$E(\^\beta|z_{(Y)})=\begin{bmatrix}
   X_A^TX_A^T & X_A^TE(\omega|z_{(Y)})X_B \\
   X_B^TE(\omega^T|z_{(Y)}) X_A & X_B^TX_B
\end{bmatrix}^{-1}\\
\begin{bmatrix}
   X_A^TE(\omega|z_{(Y)})E(\omega^T|z_{(Y)})X_A^T & X_A^TE(\omega|z_{(Y)})X_B \\
   X_B^TE(\omega^T|z_{(Y)})X_A & X_B^TX_B
\end{bmatrix}\beta\\$$ \
因此， $MLE$ 是偏向条件 $z(Y)$ 的，除非只有在唯一确定给$z(Y)$的情况下才会发生，就像在没有链接错误的情况下一样。此外，由于偏见不是渐近消失的，就如 $n\to\infty$ ， $MLE$ 的 $\beta$ 也都不是一致的。\
 $$E(\omega|z_{(Y)})E(\omega|z_{(Y)})=I,$$ \
&emsp; &emsp;以上 $MLE$ 的困难在于，由于匹配矩阵空间 $\Omega$ 大，只要建模方法明确考虑到链接数据结构，直接评估观察到的似然是不可行的。这并不意味着 $MLE$ 不能是一些统计模型下的一个合理的估计器，这些统计模型提供了对观测数据 $(X_A，Y_B，\omega_L)$ 变化的合理描述，即使它并不严格依赖于链接数据结构。在这个阶段似乎是一个开放式的问题是关于：是否或如何有可能开发一个普遍可行的基于似然的方法，并且注重链接数据结构。
## 2.4 比较数据模型下的分析

### 2.4.1 链接模型下的线性回归
&emsp; &emsp;线性回归一直是许多现有的频率主义分链接数据方法的重点。这里的一个中心概念是链接矩阵。为了集中精力于想法，假设 $A$ 和 $B$ 之间是一对一的匹配，根本没有不匹配的实体。假设 $A$ 包含所有协变量 $X_{n \times p}$ ，并修复 $A$ 中记录的排序，使得 $X_M\equiv X$ 。让 $y_L$ 是通过完成一个记录链接从 $B$ 获得的依赖变量的向量，其中 $y_L$ 的第 $i$ 个分量和 $X$ 的第 $i$ 行形成一个链接。另 $P$ 是未知的链接矩阵，这样 $y_L=Py_M$ ，其中 $y_M$ 的第 $i$ 个分量和 $X$ 的第 $i$ 行相匹配，即是 $\omega=I_{n\times n}$ 的构造。注意到 $P$ 不是观察到的链接矩阵 $\omega_L$ 。提供非信息链接，并根据定义， $X$ 作为条件，我们有，\
 $$E(y_L|X)=E[E(Py_M|G^*)|X]=E[E(Py_M|G^*)E(y_M|X)|X]\\
=E[E(P|G^*)|X]E(y_M|X)=E(P|X)X\beta=QX\beta$$ \
其中 $G^*$ 是数据链接器可访问的比较数据，数据链接器可能包含 $X$ 的一部分。当 $(X^TX)^{-1}(X^TQX)\ne I$ 时，面值最小二乘估计器 $(X^TX)^{-1}X^Ty_L$ 是偏置的，除非 $Q=I$ 。 $Lahiri$ 和 $Larsen(2005)$ 提出了调整的无偏置估计器\
 $$\^\beta_{LL}=(X^TQ^TQX)^{-1}X^TQ^Ty_L$$ \
&emsp; &emsp; $Kim$ 、 $Chambers(2012a，2012b)$ 、 $Hof$ 和 $Zwinderman(2012)$ 在特定假设下提出了允许两个以上数据集和样本寄存器链接的扩展。 $Chambers(2009)$ 指出，尽管原则上可以使用最优 $(X^TQ^T\sum_l^{-1}QX)^{-1}X^TQ^T\sum_L^{-1}y_L$ ，其中 $\sum_L=Cov(y_L，y_L|X)$ ，但增益可能较小，因为需要估计协方差矩阵 $\sum_L$ ，并且有避免误判的风险。\
&emsp; &emsp;假设上面的非信息链接意味着 $X$ 作为条件，关键变量的失真与依赖变量 $y$ 无关。它在许多应用中似乎是合理的 $(Lahiri和Larsen，2005)$ 。对于 $C-PR$ 链接，如果其他关键变量（例如，姓名，出生日期）的分布与研究变量（例如，健康）无关，则会出现这种情况，例如，考虑到上例的变量X中，可能包括性别和年龄。许多大型事业采用类似的记录链接方法；参见，例如，《澳大利亚人口普查纵向数据集》 $(Zhang$ 和 $Campbell$ ， $2012$ 年 $)$ 。信息链接模型目前仍然是一个开放式的问题  $(Chambers$ 和 $Kim$ ，$2015$ 年 $)$ 。\
&emsp; &emsp;矩阵 $Q=E(P|X)$ 包含成对真 $(qii)$ 和假 $(qij)$ 链路概率，跨越 $G^*$ 的失真过程和记录链路过程。然而，分析通常需要的是分布 $f(P|X）)$ ,例如，为了评估 $\^\beta_{LL}$ 的方差。除以下特殊情况外，此分发不能从 $Q$ 派生。让 $A=\{a_1，a_2\}$ 和 $B=\{b_1，b_2\}$ 。在完整的一对一链接下，只有两个可能的链接矩阵，即 $P=I_{2\times2}$ 或者 $P=J_{2\times2}-I$ ，其中 $J$ 是单位矩阵，因此 $q_{11}=q_{22}$ 必须是 $P=I=\omega$ 的概率，并且 $q_{12}=q_{21}$ 必须是 $P=J_{2\times2}-I$ 的概率。否则二级分析师无法访问 $f(P|X)$ 。\
&emsp; &emsp; $Chambers(2009)$ 提出了最简单链接模型下的二次分析，其中 $q_{ii}=\lambda$ 和 $q_{ij}=(1−\lambda)/(n−1)$ 用于 $j\ne i$ ，这被称为可交换链接误差 $(ELE)$ 模型。这具有允许二级分析师基于与链接误差相关的单个参数获得矩阵 $Q$ 的优点。尽管其具有简单性和二次分析中的适用性，但在 $ELE$ 模型下并不清楚 $f(P|X;\lambda)$ 是什么，除非二次分析的范围受到一定的限制。
### 2.4.2 比较数据模型下的线性回归

&emsp; &emsp;我们现在考虑比较数据模型 $f(G|\omega)$ 下的方法。下面的两节将讨论 $f(G|\omega)$ 本身的建模。如上所述，让 $X$ 来自 $A$ ， $y$ 来自 $B$ ，假设完成一对一匹配，这样 $X=X_M$ ， $y_M=\omega y$ 和 $y=\omega^Ty_M$ 。我们假设 $\omega$ 在 $\Omega$ 是均匀分布的，所提供的记录是随机排列在 $A$ 和 $B$ 中的。 $DeGroot$ 和 $Goel(1980)$ 也使用了同样的假设。非资料性比较数据，即：\
 $$f(G|\omega,X,Y;\zeta)=f(G|\omega;\zeta)$$ \
其中 $\zeta$ 是比较数据模型的参数向量。然后，
 $$E(y|G,X)=E(\omega^T|G)E(y_M|X)=Q_GX\beta$$ \
其中当 $f(\omega)$ 在数字器和面额器中被取消时， $Q_G=E(\omega^T|G)$ 是通过以下公式计算而得到的
 $$f(\omega|G;\zeta)=f(G|\omega;\zeta)/\sum_{\mathclap{\omega \in\Omega}}f(G|\omega;\zeta)\tag{2.1}$$ \
以这种方式，条件分布 $f(\omega|G)$ 可以从比较数据模型 $f(G|\omega;\zeta)$ 中得到，这样就可以对记录在给定情况下如何实际比较进行条件推理。\
&emsp; &emsp;在线性回归分析的情况下，类似于上面的 $\^\beta_{LL}$ ，可以给出 $\beta$ 给定 $（G，X）$ 的无条件无偏差估计数 $(X^TQ_G^TQ_GX)^{−1}$ ，其中 $H=（X^TQ_G^TQ_GX）^{−1}X^TQ_G^T$ 。\
 $$\^\beta_G = (X^TQ_G^TQ_GX)^{−1}X^TQ_G^Ty=Hy\tag{2.2}$$ \
对于方差协方差矩阵 $V(\^\beta_G)$ ，我们有
 $$V(\^\beta_G|X)=H(\triangle_G+\sigma^2Q_GQ_G^T)H^T$$ \
其中 $V[E(y|y_M)|G，X]=\sigma^2Q_GQ_G^T$ ，和 $\triangle_G=E[V(y|y_M)|G，X]$ 具有 $ij$ 个元素：  $\triangle_{G,ij}=\beta^TX^T\tau_{ij}X\beta+\sigma^2Trace(\tau_{ij})$ ，对于 $\tau
_{ij}=E
(\omega_i\omega_j^T|G)-E(\omega_i|G)E(\omega_j^T|G)$ 和 $\omega_i$ 是 $\omega$ 和 $\omega_j$ 的第 $i$ 列。为了估计 $\sigma^2$ ，考虑\
 $$\^e=y-\^y_M=\omega^Ty_M-\^y_M=\omega^Ty_M-X\^\beta_G=(I-XH)\omega^Ty_M$$ \
 $$E(\^e^T\^e|G)=\sigma^2Trace(D)+\beta^TX^TDX\beta$$ \
其中 $D=E(\omega R\omega^T|G,X)$ 和 $R=(I-XH)^T(I-XH)$ 。可以用 $\beta$ 替代 $\^\beta_G$ ，从而计算出 $\sigma^2$ 的插入估计。\
&emsp; &emsp;从二次分析的角度来看，比较数据模型下 $f(\omega|G)$ 的使用是与之前链接模型下方法的关键区别。值得强调的是，这两种模型都能提供有效的分析，但适用于不同的情况。如果数据链接器与链接数据集一起发布管理链接矩阵分布 $f(P|X)$ ，那么二级分析人员当然可以在链接模型下开发增益分析。然而，在现实中，除了简单的 $ELE$ 模型之外，数据提供商到目前为止甚至无法提供 $Q$ 矩阵，更别说是 $f(P|X)$ 。到如此程度上仅提供链接比较数据 $G$ ，同时比较数据建模提供了一个原则选项。
