&emsp;&emsp;**示例 4. 列联表‑** 假设给定一个具有 $K$个类别的离散$r.v.X$，$Y$和$Z$是分别具有$\psi$和$\phi$类别的离散$r.v.s$,不一定是有序的。在这种情况下，到目前为止开发的整个理论仍然可以通过简单地用概率函数替换d.f.s来应用。不失一般性，符号$k = 1,...,K,\psi = 1,...,\Psi$和$\phi=1,...,\Phi$分别表示$X、 Y$和$Z$所采用的类别。

&emsp;&emsp;设$\theta_{\psi\phi|k}$为$(Y=\psi,Z=\phi|X=k)$的概率，边缘$\theta_{\psi.|k}$表示事件的概率$(Y=\psi|X=k)$,$\theta_{.\phi|k}$表示事件的概率$(Z=\phi|X=k)$。根据

$$max(\theta,\theta_{\psi.|k}+\theta_{.\phi |k}−1)\le\theta_{\psi,\phi|k}\le min(\theta_{\psi.|k},\theta_{ .\phi |k})  \tag{4.37}$$

条件不确定性测度结果等于
 
$$\Delta^{x=k}=\sum_{\psi=1}^\Psi\sum_{\phi=1}^\Phi\{min(\theta_{\psi.|k}\theta_{.\phi |k})-max(\theta,\theta_{\psi.|k}+\theta_{.\phi |k}-1)\}\theta_{\psi.|k}\theta_{.\phi |k}$$

 
总体不确定性度量是
$$\Delta = \sum_{k=1}^K\Delta^{x=k}\theta_k $$

其中$\theta_k$表示事件的概率$(X=k)$。当对$(X,Y,Z)$所采用的类别进行排序时，可以获得更清晰的结果[<font color=Blue>17</font>]。为了简单起见，我们使用自然数的习惯顺序。在这种情况下，累积 $d.f.s$ 是
$$H_{\psi,\phi|k}=\sum_{\psi=1}^\Psi\sum_{\phi=1}^\Phi\theta_{\psi\phi|k},\tag{ $\psi=1,...,\Psi,\phi=1,...,\Phi$}$$

$$F_{\psi|k}=\sum_{\psi=1}^\psi\theta_{y.|k},\tag{ $\psi=1,...,\Psi,k=1,...,K$}$$

$$G_{\phi|k}=\sum_{\phi=1}^\phi\theta_{z.|k},\tag{$\Psi,k=1,...,K$}$$

那么,下面的不等式

$$max(0,F_{\psi|k}+G{\phi|k}-1)\le H_{\psi,\phi|k}) \le min(F_{\psi|k},G{\phi|k})\tag{4.38}$$

请注意，不等式 $(4.38)$ 意味着
$$\theta^-_{\psi\phi|k}\le \theta_{\psi\phi|k}\le \theta^+_{\psi\phi|k}\tag{4.39}$$

此时
$$\theta^-_{\psi\phi|k}=L(F_{\psi|k},G_{\phi|k})-L(F_{\psi-1|k},G_{\phi|k})-L(F_{\psi|k},G_{\phi-1|k})+L(F_{\psi-1|k},G_{\phi-1|k})$$

$$\theta^+_{\psi\phi|k}=U(F_{\psi|k},G_{\phi|k})-U(F_{\psi-1|k},G_{\phi|k})-U(F_{\psi|k},G_{\phi-1|k})+U(F_{\psi-1|k},G_{\phi-1|k})$$

给定两个实数a、 b，两个量L、U定义如下
$$L(a,b)=max(0,a+b-1),U(a,b)=min(a,b)\tag{4.40}$$

不难证明
$$\theta^-_{\psi\phi|k}\ge L(\theta_{\psi.|k},\theta_{.\phi|k})$$

$$\theta^+_{\psi\phi|k}\le U(\theta_{\psi.|k},\theta_{.\phi|k})$$

因此不等式 (4.39) 比 (4.37)更尖锐。无论如何，条件不确定性度量是

$$\Delta^{x=k}=\sum_{\psi=1}^\Psi\sum_{\phi=1}^\Phi\{U(F_{\psi|k},G_{\phi |k})-L(F_{\psi|k},G_{\phi |k})\}\theta_{\psi.|k}\theta_{.\phi |k}\tag{4.41}$$

无条件不确定性测度是
$$\Delta = \sum_{k=1}^K\Delta^{x=k}\theta_k \tag{4.42}$$

&emsp;&emsp;与之前发现的相比，不确定性度量不等于$1/6$，因为不确定性度量取决于$Y|X$和$Z|X$的边际概率。正如之前所强调的，这并不成立。
当$Y$和$Z$连续时。有关详细信息，请参阅[<font color=Blue>17</font>]，其中调查了有序分类变量的统计匹配的不确定性。

-------------------------------------

**4.4	复杂样本调查的统计匹配**
为了估计合理的匹配分布以及对不确定性度量$(4.31)$和$(4.33)$进行推断，有必要对样本$s_A$、$s_B$抽取所依据的抽样设计进行假设。文献中关于统计匹配的大多数论文都基于样本$s_A$、$s_B$均由独立同分布观察组成的假设，例如参见[<font color=Blue>35</font>]及其中的参考文献。

在[<font color=Blue>16</font>]、   [<font color=Blue>17</font>]和[<font color=Blue>19</font>]中介绍和研究了$i.i.d$观察的非参数设置中统计匹配的不确定性。具体地， 在[<font color=Blue>16</font>]和[<font color=Blue>16</font>]中，讨论了统计匹配中的不确定性，并引入了基于Frechet类的不确定性度量$（4.31）$，$（4.33）$,以衡量引入对模型不确定性的影响逻辑约束$(4.17)$。


然而，独立同分布假设几乎对调查数据无效，因为存在分层、不同级别的聚类以及基于规模   度量的不同包含概率，样本设计通常很复杂。


[<font color=Blue>42</font>]和[<font color=Blue>39</font>]研究了复杂样本调查中的统计匹配。 [<font color=Blue>42</font>]提出的方法包括通过计算相对于人工设计 $p(s)$ 的新采样权重来连接两个文件$s_A$和$s_B$ ，其中$s = s_A\cup s_B$，由作用于$s_A$和$s_B$的采样设计的并集给出，分别。这种方法在实践中很少使用，因为它需要知道一个样本中的单位在另一个样     本的抽样设计下的包含概率。    [<font color=Blue>39</font>]中提出的方法包括调整两个不同样本sA和sB的实际调查权 重，以便通过重复应用校准程序在$X、 Y |X$和$Z|X$上具有均匀的公共分布。





D'Orazio 等人。[<font color=Blue>24</font>]比较它们的效率，包括基于最大伪似然的额外估计，如[<font color=Blue>48</font>]中的那些。无论如何，在有限的样本量下，这些估计器中没有任何突出的估计器。

当根据复杂的调查设计抽取样本时，[<font color=Blue>18</font>]给出了一种可能性。他们将注意力限制在粗略地说⾼熵的抽样设计上（实际假设在第 $4.4.1$节中）。在此设置中，可以通过迭代比例拟合(IPF) 算法（第$4.4.2$节）获得感兴趣变量的合理联合分布的存在。孔蒂等人。 [<font color=Blue>18</font>]通过根据渐近特性（第  $4.4.3$节）分析估计器的可靠性，克服了为有限样本量建立统计匹配的最佳方法的问题允许使用一些额外的工具，作为不确定宽度测试的定义。

------------------------------------
 **4.4.1  复杂抽样调查的统计匹配**	
样本设计的技术假设抽取样本$s_A$、$s_B$所依据的抽样设计假设参考了以下框架。
&emsp;&emsp;对于总体 $U$ 的每个单元$i$ ,设$D_{i,\alpha}(α = A,B)$为伯努利随机变量 $(r.v.)$，如果$\mathbf D_{i,\alpha}= 1$，则$i$在样本$s_\alpha$中，而%i%不在$s_\alpha$中如果$D_{i,\alpha}= 0$。让$\mathbf D_{i,\alpha}= (D_{1,\alpha} ... D_{N,\alpha})$和$(\alpha = A,B)$。 $(\alpha = A, B)$。在 $A$ 抽样设计中， $P_\alpha$是特定$D_{N,\alpha}$的概率分布$\pi_{i,\alpha}=E_{P_\alpha}[D_{i,\alpha}]$是单元i在抽样设计$P_\alpha (\alpha = A,B)$下的一阶包含概率。$s_\alpha$的（有效）大小.接下来，我们将只考虑固定的$n_{s,\alpha}= D_{1,\alpha} + ··· + D_{N,\alpha}$ （有效）尺寸抽样设计，使得$n_{s,\alpha} \equiv  n_\alpha (\alpha = A,B)$。

&emsp;&emsp;对于每个单元$i$，令$p_{i,\alpha}$为正数，其中$p_{i,\alpha}+ ··· + p_{N,\alpha}=n_\alpha(\alpha=A,B)$。具有参数$p_{i,\alpha},...,p_{N,\alpha}$的泊松抽样设计的特征在于假设$r.v.s$D_{i,\alpha}独立于$P_{r_{p_{o_\alpha}}}(D_{i,\alpha}=1)=p_{i,\alpha}(\alpha=A,B)$

&emsp;&emsp;拒绝采样或归一化条件泊松采样$P_{R,\alpha}$（参见[<font color=Blue>27</font>]、 [<font color=Blue>47</font>])是通过条件$w.r.t.n_{s,\alpha}= n_\alpha(\alpha= A,B)$从泊松采样获得的。在符号中：

$$P_{r_{P_{R_\alpha}}}(\mathbf D_{N,\alpha}) = P_{r_{P_{o_\alpha}}}(\mathbf D_{N,\alpha}|n_{s,\alpha}=n_\alpha),$$

其中$\alpha=A,B$
 
抽样设计$P_\alpha$与拒绝之间的Hellinger距离符号$P_{R,\alpha}(\alpha = A,B)$定义为

$$d_H(P_\alpha,P_{R,\alpha})=\sum(\sqrt{P_{r_{P_\alpha}}(\mathbf D_{N,\alpha})}-\sqrt{P_{r_{P_{R,\alpha}}}(\mathbf D_{N,\alpha})})^2,\alpha=A,B \tag{4.43}$$

其中和在$D_{1,\alpha},...,D_{N ,\alpha}$ 上扩展。抽样设计的基本假设如下：

$A1.(\mathcal{U}_N;N\ge1)$是$N$递增的有限种群序列。

$A2.$对于每个$N, (x_1, y_1, z_1), ..., (x_N , y_N , z_N )$是超总体&(X_1,Y_1,Z_1), ..., (X_N , Y_N , Z_N)&的实现$i.i.d. r.v.s. (X_i , Y_i , Z_i )$与公共 $d.f. H(x, y, z)$。从现在开始，符号 $P$ 用于表示 $r.v.s. (X_i,Y_i,Z_i)s,$的（超种群）概率分布，$\Bbb{E}$和$\Bbb{V}$分别表示相应的均值和方差算子。此外，假设Xi是离散的,$r.v.s.$，取值为正概率$p(x^1), ..., p(x^K) $。

$A3$.样本$s_A、s_B$是从总体$U$中独立选择的。此外，根据固定大小的样本设计$P_\alpha$选择$s_\alpha$,统计匹配中的不确定性和估计,具有正一阶包含概率$\pi_{1,\alpha}, ..., \pi_{N,\alpha}$和样本大小$nα = \pi_{1,\alpha} + ··· + \pi_{N,\alpha }(\alpha = A,B)$。尽管样本量、包含概率和 $r.v.s.  D_{i,A^s},D_{i,B^s}$取决于$N$，但为了简化相应的符号系统，将省略后缀$N$。还假设

$$d_{N,\alpha}=\sum_{i=1}^N\pi_{i,\alpha}(1-\pi_{i,\alpha})\rightarrow\infty and  \cfrac 1Nd_{N,\alpha}\rightarrow d_\alpha,\alpha=A,B$$

其中$0\lt d_\alpha<\infty(\alpha=A,B)$

$A4$.样本量$n_A、 n_B$随着总体规模$N$的增加而增加，其中
$$\lim_{N\to\infty}\cfrac {n_\alpha} N=\int_\alpha\in(0,1), \alpha=A,B $$

$A5$.对于每个总体$\mathcal{U}_N ; N \ge 1)$，令$P_{R,\alpha}$为具有包含概率$\pi_{1,\alpha}, ..., \pi_{N,\alpha}$ 的拒绝抽样设计，并令$P_\alpha$为实际抽样设计，具有相同的包含概率作为$P_{R,\alpha} (\alpha= A,B)$。然后

$$d_H (P_\alpha, P_{R,\alpha})\rightarrow  0 as N\rightarrow\infty,\alpha = A, B。$$

$A6$.存在正实数$\varsigma_A,\varsigma_B,\zeta_A,\zeta_B$使得

$$\lim_{N\to\infty}\cfrac 1 N\sum_{i=1}^N\cfrac 1 {\pi_{i,\alpha}}=\varsigma_A\lt\infty,   \lim_{N\to\infty}\sum_{i=1}^N\cfrac 1{(i\pi_{i,\alpha})^2}=\zeta_A\lt\infty,\alpha=A,B$$

&emsp;&emsp;假设 A1‑A6 类似于[<font color=Blue>12</font>]和[<font color=Blue>13</font>] 中使用的假设。它们主要用于确保后续提出的估计量的渐近正态性。特别地，假设A5意味着抽样设计$P_A、P_B$具有最大渐近熵。此类抽样设计的示例有[<font color=Blue>6</font>]和[<font color=Blue>5</font>]中所示的简单随机抽样、连续抽样、Sampford设计、Chao设计、分层设计等。总之，如[12]和[13]中所述，它们甚至可能过于严格   并被用于简化符号和证明：例如，  $（X_i ，Y_i，Z_i）$的独立性或同分布对于它们的有效性并不是真正必要的。

&emsp;&emsp;鉴于[<font
color=Blue>18</font>]中的主要发现是渐近的,如果人口规模N趋于无穷大，还要考虑相应的序列
$$\mathbf x_\infty=(x_1,x_2,...),\mathbf y_\infty=(y_1,y_2,...),\mathbf z_\infty=(z_1,z_2,...) \tag{4.44}$$

分别为$x_is、 y_is、 z_is$值。实际的$\mathbf x_N , \mathbf y_N , \mathbf z_N$ (4.1)分别是序列$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty $中前N 个$x_is、 y_is、 z_is$的“段” 。

序列$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty $由 A2 中的超种群模型产生。因此$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty $存在于概率空间中$((\Bbb R^3)^\infty,\mathcal{B}(\Bbb R^3)^\infty,\Bbb P^\infty)$其中 $\mathcal{B}(\Bbb R^3)^\infty$是$(\Bbb R^3)^\infty$ 上的 Borel σ 场的乘积\ , $\Bbb P$上的乘积测度$(\Bbb R^3)^\infty$ 由 $\Bbb P^\infty$生成。我们考虑的概率陈述的形式为$P_{r_P}(·|\mathbf x_N , \mathbf y_N , \mathbf z_N$)，其中N趋于无穷大。后缀P表示概率指的是抽样设计$P_A、P_B$。我们将获得的结果适用于A2中的超种群模型可以产生的“几乎所有”序列$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty$，即对于具有P∞‑概率 1 的一组序列。 我们将使用表达式“with P‑probability 1”。

&emsp;&emsp;从现在开始，符号

$$F(x,y) = H(x, y, +\infty); G(x, z) = H(x, +\infty, z); \tag{4.45}$$
$$Q(x) = H(x,+\infty, +\infty); p(x) = Q(x) − Q(x^-)  \tag{4.46} $$

将表示$(X_i , Y_i ),(X_i,Z_i)$的联合超种群 $d.f.s. ， X_i$的边缘超种群 $d.f.$ ，
和$X_i$的概率函数,（注意$p(x) \gt0 $当且仅当$x\in \{x^1,...,x^K\}$,最后，符号

$$F(y|x)=\cfrac {F(x,y)}{p(x)},F(z|x)=\cfrac {G(x,z)}{p(x)},H(y,z|x)=\cfrac {H(x,y,z)}{p(x)}\tag{4.47}$$

将表示$Y_i, Z_i,(Y_i;Z_i)$各自给出$X_i=x$.

**4.4.2 选择匹配分布的建议**

统计匹配的第一个目标是为给定$X$的$Y$和$Z$选择匹配分布，即在类(4.27)中，分别由样本$s_A$和$s_B$估计裕度。在约束条件下(4.17)，CIA假设$H_N(y,z|x)=F_N(y|x)G_N(z|x)$是不允许的，并且创建一个满足约束的发行版似乎不是一件容易的任务。

&emsp;&emsp;在样本水平上，$d.f.s F_N,G_N$是未知的，必须在样本数据的基础上估计。Conti等人[<font color=Blue>18</font>]建议使用H´ajek估计量:

$$\widehat {F_H}(y|x)=\cfrac{\sum_{i=1}^N\cfrac{D_{i,A}}{\pi_{i,A}}I_{(y_i\le y)}I_{(x_i= x)}}{\sum_{i=1}^N\cfrac{D_{i,A}}{\pi_{i,A}}I_{(x_i= x)}} , \widehat {G_H}(z|x)=\cfrac{\sum_{i=1}^N\cfrac{D_{i,B}}{\pi_{i,B}}I_{(z_i\le z)}I_{(x_i= x)}}{\sum_{i=1}^N\cfrac{D_{i,B}}{\pi_{i,B}}I_{(x_i= x)}} \tag{4.48}$$

接下来定义$d.f. \widetilde{H_H}(y,z|x)$

$$\int_Ld\widetilde{H_H}(y,z|x)=\widehat C\int_LI_{(a_x\le f_x(y,z)\le b_x)d[\widehat F_H(y|x)\widehat G_H(z|x)]} \tag{4.49}$$

其中$L$是$\Bbb R^2$中的任何$Borel  \sigma‑field$，并且

$$\widehat C^{-1}=\int_{\Bbb R^2}I(a_x\le f_x(y,z)\le b_x)d[\widehat F_H(y|x)\widehat G_H(z|x)]=\cfrac{\sum_{i\in s_A}\sum_{j\in s_B}I(a_x\le f_x(y,z)\le b_x)\pi ^{-1}_{i,A}\pi ^{-1}_{j,B}I_{(x_i=x)}I_{(x_j=x)}}{\sum_{i\in s_A}\sum_{j\in s_B}\pi ^{-1}_{i,A}\pi ^{-1}_{j,B}I_{(x_i=x)}I_{(x_j=x)}}$$













 




