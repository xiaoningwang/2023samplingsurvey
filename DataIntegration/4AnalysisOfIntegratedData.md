&emsp;&emsp;根据前面的例子，可以说，通常情况下，参数设置中缺乏可识别性的主要后果是，模型的一些参数无法根据可用的样本信息进行估计。例如，人们只能合理地构建合理的估计集，而不是参数点估计，与实际估计的一致。这些集合（通常是区间）正式提供了模型参数不确定性的表示。[16]中给出了分析统计匹配不确定性的参数和非参数方法的统一框架。  
&emsp;&emsp;当额外的样本信息可用时，可以改进Frechet边界（4.12）。在统计实践中，一种经常可用的额外样本信息由逻辑约束组成，即对给定X的变量（Y，Z）的支持的限制。一个足够通用的约束具有以下形式

$$a_{x}\leq f_{x}(y,z)\leq b_{x}\tag{4.17}$$

给定 $X=x$ ，其中 $f_{x}(y,z)$ 是固定 $z$ 的 $y$ 的递增函数，是固定 $y$ 的 $z$ 的递减函数[19]。  
&emsp;&emsp;在实践中，约束（4.17）意味着，对于人口 $U_{N}$ 的每个单位 $i$ ，值 $(x_{i},y_{i},z_{i})$ 必须满足一对不等式 $a_{x_{i}}\leq f_{x_{i}}(y_{i},z_{i})\leq b_{x_{i}}$ 。粗略地说，约束（4.1 7）修改了条件概率分布函数 $H_{N}(y,z|x)$ 的支持，即

$$\\{(y_{i},z_{i}):a_{x_{i}}\leq f_{x_{i}}(y_{i},z_{i})\leq b_{x_{i}}andx_{i}=x\\}.$$

&emsp;&emsp;在约束条件（4.17）下，Frechet边界（4.12）减小到

$$K_{N-}^{x}(y,z)\leq H_{N}(y,z|x)\leq K_{N+}^{x}(y,z)\tag{4.18}$$

其中（使用 $\wedge $ 表示两个数字之间的最小值）

$$K_{N-}^{x}(y,z)=max(0,G_{N}(z|x)\wedge G_{N}(\gamma _{y}(a_{x})|x)+F_{N}(y|x)\wedge F_{N}(\delta_{z}(b_{x})|x)-1,F_{N}(y|x)+G_{N}(z|x)-1)\tag{4.19}$$

$$K_{N+}^{x}(y,z)=min(G_{N}(z|x),G_{N}(\gamma _{y}(a_{x})|x),F_{N}(y|x),F_{N}(\delta_{z}(b_{x})|x))\tag{4.20}$$

且对于固定的 $y$ 和 $z$ 来说， $\gamma_{y}(\cdot )$ , $\delta_{z}(\cdot )$ 就分别是 $f_{x}(y,z)$ 的反函数。当

$$K_{N-}^{x}(y,z)\geq max(0,F_{N}(y|x)+G_{N}(z|x)-1),$$

$$K_{N+}^{x}(y,z)\leq min(F_{N}(y|x),G_{N}(z|x)),$$

边界（4.18）实际上是对无约束Frechet边界（4.12）的改进。

&emsp;&emsp;例2.假设存在常数 $a_{x}$ 、 $b_{x}$ ，设 $f_{x}(Y,Z)=Y/Z$ ，则对常数的约束形式为 $a_{x}\leq Y/Z\leq b_{x}$ 。例如，在家庭调查中， $X$ 是家庭的社会经济特征（即家庭成员的数量）， $Y$ 是家庭消费， $Z$ 是家庭收入。然后，使用国民核算技术，可以为每个家庭规模产生相当合理的平均消费倾向的下限和上限，即消费支出和收入之间的比率。因此，对于每个家庭规模，可以计算给定收入的消费支出的下限和上界。[14]分析了这种约束，其中对意大利银行的家庭收入和财富调查（简称SHIW）和意大利国家统计局的家庭预算调查（简称HBS）进行了统计匹配。  
&emsp;&emsp;上面引入的约束意味着，对于每个给定的 $x$ ，条件概率分布函数 $H_{N}(y,z|x)$ 的支持度在两条直线 $z=y/b_{x}$ 和 $z=y/a_{z}$ 之间，即，是圆锥体的（适当或不适当的子集）

$$\\{(y_{i},z_{i}):y_{i}\geq 0,y_{i}/b_{x}\leq z_{i}\leq y_{i}/a_{x}\\}.\tag{4.21}$$ 

&emsp;&emsp;由（4.19）、（4.20）的结果，令 $\gamma_{y}(a_{x})=y/a_{x}$ ， $\delta_{z}(b_{x})=b_{x}z$ ，我们得到以下的Frechet界

$$K_{N-}^{x}(y,z)=max(0,G_{N}(z|x)\wedge G_{N}(\frac{y}{a_{x}} |x)+F_{N}(y|x)\wedge F_{N}(b_{x}z|x)-1,F_{N}(y|x)+G_{N}(z|x)-1)\tag{4.22}$$

$$K_{N+}^{x}(y,z)=min(G_{N}(z|x),G_{N}(\frac{y}{a_{x}}|x),F_{N}(y|x),F_{N}(b_{x}z|x))\tag{4.23}$$

&emsp;&emsp;约束条件 $a_{x}\leq Y/Z\leq b_{x}$ 的另一个例子发生在商业调查中， $X$ 可以是活动类型， $Y$ 是总销售额，而 $Z$ 是员工人数。

&emsp;&emsp;例3.在实践中经常出现的另一种约束是 $Y\geq Z$ ，给定 $X$ 。例如，在[34]的情况下， $Y$ 为总收入， $Z$ 则为应纳税收入。约束 $Y\geq Z$ 意味着，对于每个给定的 $X$ ，条件概率分布函数 $H_{N}(y,z|x)$ 的支持度变成

$$\\{(y_{i},z_{i}):y_{i}\geq z_{i}\\}.\tag{4.24}$$

则它是直线 $z=y$ 下方的半平面（的子集）。

&emsp;&emsp;设置 $f_{x}(y,z)=y/z$ ，其中 $a_{x}=1$ ， $b_{x}\rightarrow \infty$ ，根据结果（4.19）、（4.20），我们得到了以下的新Frechet界

$$K_{N-}^{x}(y,z)=max(0,G_{N}(z|x)\wedge G_{N}(y|x)+F_{N}(y|x)-1,F_{N}(y|x)+G_{N}(z|x)-1)\tag{4.25}$$

$$K_{N+}^{x}(y,z)=min(F_{N}(y|x),G_{N}(z|x),G_{N}(y|x))\tag{4.26}$$

&emsp;&emsp;请注意，不确定性的概念可以用于不同的任务。D’Orazio等人在[22]使用不确定性来从两个数据集中可用的匹配变量中选择最佳的匹配变量集。其他人，如[33]，试图通过最小推理，给出在不可量化问题中选择参数值的替代方法。在本文中，我们将主要使用不确定性来表示不存在联合观测值的一对变量的联合分布，以及当明确做出一个选择时匹配误差的表示（见4.3.1）。

4.3.1 匹配误差  
由于只有 $(X,Y)$ 和 $(X,Z)$ 的联合观测可用，因此只能根据样本数据估计 $Q_{N}(x)$ ， $F_{N}(y|x)$ ， $G_{N}(z|x)$ 。即使条件概率分布函数 $F_{N}(y|x)$ 和 $G_{N}(z|x)$ 完全已知，变量 $(X,Y,Z)$ 缺乏联合观测也是 $H_{N}(y,z|x)$ 不确定性的原因。也就是说，可用信息不能在 $(X,Y,Z)$ 的一组合理的联合分布之间进行区分。  
&emsp;&emsp;就 $H_{N}(y,z|x)$ 而言，除非做出特殊假设，否则我们唯一能确定的是，它位于所有合理的联合概率分布函数中，定义如下

$$H_{N}^{x}=\\{H_{N}(y,z|x):H_{N}(y,\infty |x)=F_{N}(y|x),H_{N}(\infty ,z|x)=G_{N}(z|x),a_{x}\leq f_{x}(y,z)\leq b_{x}\\}\tag{4.27}$$

集合 $H_{N}^{x}$ 包含具有边缘 $F_{N}(y|x)$ 和 $G_{N}(z|x)$ 并且满足所施加的逻辑约束 $(Y,Z)|X$ 的所有联合概率分布函数。  
&emsp;&emsp;于是，类（4.27）中的每个二元分布都是给定 $X$ 的 $Y$ 和 $Z$ 的联合分布的候选者。从现在起，类（4.27）中的每一个二元分布函数都将被称为给定 $X$ 的 $Y$ 和 $Z$ 的匹配分布。  
&emsp;&emsp;统计匹配问题本质上在于选择匹配分布，即 $H_{N}^{x}$ 类中的分布函数。  
&emsp;&emsp;匹配过程本质上是在类（4.27）中选择一种分布函数的过程，或者更好地，在（4.27）所定义的类中选择分布函数，但边缘分布函数是基于样本数据估计的。  
&emsp;&emsp;假设选择 $H_{N}^{\ast }(y,z|x)\in H_{N}^{x}$ 作为给定 $X$ 的 $(Y,Z)$ 的匹配分布，但“真实”的分布是 $H_{N}(y,z|x)$ （当然也是在 $H_{N}^{x}$ 中）。所选 $H_{N}^{\ast }(y,z|x)$ 和实际 $H_{N}(y,z|x)$ 之间的差异为匹配误差。匹配误差的概念在评估匹配过程中至关重要，因为匹配误差越小，匹配过程越好。  
&emsp;&emsp;匹配误差最有利的情况发生在CIA下，因为类（4.27）折叠成单个分布函数并且从（4.11）为了构建 $H_{N}(y,z|x)$ 的一致估计，构建 $F_{N}(y|x)$ 和 $G_{N}(z|x)$ 的一致估算就足够了。在这种情况下，随着 $s_{A}$ 和 $s_{B}$ 的大小增加，匹配误差变为负值。  
&emsp;&emsp;[30]中评估了CIA简化背景下非参数插补程序和独立同分布观测结果产生的匹配误差。该类基于 $s_{B}$ 中 $Z$ 对  $X$ 回归函数的 $k$ 最近邻（kNN）非参数估计，包括固定和可变数量的供体 $k$ ，并包括一些最流行的非参数插补程序，如距离和随机热甲板。对插补程序的渐近财产进行了初步分析，然后通过仿真进行了研究。在[15]中，我们进一步引入了基于局部线性回归的新的非参数匹配技术，其性能通过匹配误差来衡量。  
&emsp;&emsp;不幸的是，在许多具有实际重要性的情况下，CIA是不合理的。在这种情况下，类（4.27）不会减少到单个点，并且即使 $s_{A}$ 和 $s_{B}$ 的大小增加，匹配误差也不能变得可以忽略。  
&emsp;&emsp;尽管存在这一缺陷，但对匹配误差的研究仍然具有重要意义，因为“小”的匹配误差意味着所选择的匹配分布 $H_{N}^{\ast }(y,z|x)$ “接近”实际的 $H_{N}(y,z|x)$ ，因此用 $H_{N}^{\ast }(y,z|x)$ 代替 $H_{N}(y,z|x)$ 不会产生“大”的误差。  
&emsp;&emsp;从现在起，为了评估匹配分布 $H_{N}^{\ast }(y,z|x)$ 的准确性，作为真实概率分布函数 $H_{N}(y,z|x)$ 的估计器，作为匹配误差度量（有条件地基于 $X$ ），我们将使用以下内容：

$$ME_{x}(H_{N}^{\ast },H_{N})={\int _{R^2}}|H_{N}^{\ast }(y,z|x)-H_{N}(y,z|x)|d\lbrack F_{N}(y|x)G_{N}(z|x)\rbrack .\tag{4.28}$$

&emsp;&emsp;通过类似的推理，作为匹配误差的无条件度量，我们可以考虑以下内容：

$$ME(H_{N}^{\ast },H_{N})=\int _{R}ME_{x}(H_{N}^{\ast },H_{N})dQ_{N}(x)=\int _{R}\left\\{\int _{R^2}|H_{N}^{\ast }(y,z|x)-H_{N}(y,z|x)|\times d\lbrack F_{N}(y|x)G_{N}(z|x)\rbrack \right\\}dQ_{N}(x).\tag{4.29}$$

4.3.2 通过不确定性测量来限制匹配误差  
本节专门讨论统计匹配框架中的不确定性概念。这种不确定性是由于缺乏有关感兴趣变量的联合信息，当联合概率分布函数的统计模型不可识别时，这种不确定性也是相关的。  
&emsp;&emsp;统计匹配中的不确定性度量很重要，因为它与匹配误差有关，因此可以提供关于合理匹配分布估计的可靠性的信息。有关技术细节，请参阅[19]。  
&emsp;&emsp;请注意，如果一方面匹配误差的概念对统计匹配过程的评估至关重要，另一方面 $s_{A}$ 、 $s_{B}$ 中的样本数据不包含足够的信息来估计匹配误差（4.28）。然而，可以使用4.3.1.中的分析来克服这一缺陷

$$|H_{N}^{\ast }(y,z|x)-H_{N}(y,z|x)|\leq K_{N+}^{x}(y,z)-K_{N-}^{x}(y,z)\forall x,y,z$$

易得出

$$ME_{x}(H_{N}^{\ast },H_{N})\leq \bigtriangleup ^{x}(F_{N},G_{N})\forall x\tag{4.30}$$

其中

$$\bigtriangleup ^{x}(F_{N},G_{N})=\int _{R^2}(K_{N+}^{x}(y,z)-K_{N-}^{x}(y,z))d[F_{N}(y|x)G_{N}(z|x)]$$

$$=\frac{1}{N^2p_{N}(x)^2}\sum _{i=1}^{N}\sum _{j=1}^{N}(K_{N+}^{x}(y_{i},z_{j})-K_{N-}^{x}(y_{i},z_{j}))I_{(x_{i}=x)}I_{(x_{j}=x)}.\tag{4.31}$$

&emsp;&emsp;从某种意义上说，（4.31）中的数量 $\bigtriangleup ^{x}(F_{N},G_{N})$ 代表了类（4.27）“规模”的衡量标准。  
&emsp;&emsp; $\bigtriangleup ^{x}(F_{N},G_{N})$ 越小，匹配误差越小，匹配分布越接近给定 $X$ 的 $Y$ 和 $Z$ 的实际分布。因此， $\bigtriangleup ^{x}(F_{N},G_{N})$ 可以合理地提出作为一种衡量标准，以评估使用 $(Y,Z)|X$ 的匹配分布作为实际分布函数的替代品的可靠性。也就是说，作为类（4.27）可靠性的衡量标准。  
&emsp;&emsp;在（4.30）中， $\bigtriangleup ^{x}(F_{N},G_{N})$ 可以定义为当真实的总体分布函数 $H_{N}(y,z|x)$ 被可信联合概率密度函数的类（4.27）中的匹配分布 $\tilde{H_{N}}(y,z|x)$ 替换。  
&emsp;&emsp;那么，不确定性度量越小，实际使用的匹配分布就越可靠。不确定度测量（4.31）最有趣的方面是，它可以根据观察到的样本数据进行估计。  
&emsp;&emsp;可以为无条件匹配错误（4.29）构造类似于（4.30）的边界。事实上，立刻就能得出

$$ME(H_{N}^{\ast },H_{N})\leq \bigtriangleup (F_{N},G_{N})\tag{4.32}$$

其中

$$\bigtriangleup (F_{N},G_{N})=\sum _{x}\bigtriangleup ^{x}(F_{N},G_{N})p_{N}(x).\tag{4.33}$$

是对 $（X，Y，Z）$ 联合总体分布的无条件不确定性的度量。显然，无条件不确定度量（4.33）是不确定性条件度量（4.31）的平均值，也就是 $X$ 的边际总体分布。

备注1.不确定度测量（4.31）和（4.33）的一个有趣的特性是可以计算它们的最大值。正如处理多元分布函数时经常发生的情况一样，copula函数的使用简化了问题，例如，如[33]中所示。设 $H_{N}$ 是一个分别具有裕度 $F_{N}$ 和 $G_{N}$ 的概率分布函数。如[19]中所述，存在一个copula 函数 $C^{x}$ ，使得对于所有 $(z,y)$ 

$$H_{N}(y,z|x)=C^{x}(F_{N}(y|x),G_{N}(z|x)).\tag{4.34}$$

根据Sklar定理，如果 $F_{N}$ 和 $G_{N}$ 是连续的，那么 $C^{x}(\cdot ,\cdot )$ 是唯一的，并且它等于 $H_{N}(F_{N}^{-1}(y|x),G_{N}^{-1}(z|x)).$ 。copula函数（4.34）表示 $Y$ 和 $Z$ 之间的“内在”关联，而不考虑它们的边际概率分布函数。在匹配问题中，可以根据样本数据估计边际概率分布函数 $F_{N}(y|x)$ ， $G_{N}(z|x)$ ，“实际不确定性”仅涉及copula函数（4.34）表示的关联。Frechet界（4.12）的copula函数递推由下式给出

$$\\{C^{x}(u,v):W^{x}(u,v)\leq C^{x}(u,v)\leq M^{x}(u,v),\forall u,v\in I\\}.\tag {4.35}$$

在条件 $x$ 上，无论 $F_{N}$ 和 $G_{N}$ 的形状如何， $U=F_{N}(Y|x)$ 和 $V=G_{N}(Z|x)$ 在（0,1）中都具有均匀分布。此外， $W^{x}(u,v)=max(0,u+v-1)$ ， $M^{x}(u,v)=min(u,v)$ 本身就是copula函数；无论是消极的还是积极的，它们都代表着完美的依赖。因此，不确定度测量的最大值（4.31）变为表面 $M^{x}(u,v)$ 和 $W^{x}(u,v)$ 之间的体积：

$$\int _{0}^{1}\int _{0}^{1}\\{min(u,v)-max(0,u+v-1)\\}dudv=\frac{1}{6}.\tag{4.36}$$

值1/6表示当没有超出裕度 $F_{N}$ 和 $G_{N}$ 的知识范围的外部辅助信息可用时所实现的不确定性。不确定性的无条件测量也是如此（4.33）。粗略地说，当唯一可用的信息是数据时，不确定性测量等于1/6，与边际分布函数 $F_{N}$ ， $G_{N}$ 无关。这与直觉一致，因为一方面不确定性只取决于copula函数 $C^{x}(u,v)$ 的最大值和最小值，而不取决于边缘 $F_{N}$ ， $G_{N}$ ，另一方面，数据不提供任何关于 $C^{x}(u,v)$ 的信息。当没有关于 $C^{x}(u,v)$ 的额外样本辅助信息可用时， $C^{x}(u,v)$ 的最小值和最大值分别为 $max(0,u+v-1)$ 和 $min(u,v)$ ，与可用数据无关。从某种意义上说，这是最大不确定性的情况。

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

序列$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty $由 A2 中的超种群模型产生。因此$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty $存在于概率空间中$((\Bbb R^3)^\infty,\mathcal{B}(\Bbb R^3)^\infty,\Bbb P^\infty)$其中 $\mathcal{B}(\Bbb R^3)^\infty$是$(\Bbb R^3)^\infty$ 上的 $Borel \sigma$场的乘积\ , $\Bbb P$上的乘积测度$(\Bbb R^3)^\infty$ 由 $\Bbb P^\infty$生成。我们考虑的概率陈述的形式为$P_{r_P}(·|\mathbf x_N , \mathbf y_N , \mathbf z_N$)，其中N趋于无穷大。后缀P表示概率指的是抽样设计$P_A、P_B$。我们将获得的结果适用于A2中的超种群模型可以产生的“几乎所有”序列$\mathbf x_\infty , \mathbf y_\infty , \mathbf z_\infty$，即对于具有P∞‑概率 1 的一组序列。 我们将使用表达式“with P‑probability 1”。

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













 







