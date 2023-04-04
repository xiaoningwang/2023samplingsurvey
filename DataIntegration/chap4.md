#### 以下内容是黄张欢所写
# 4

## *统计匹配中的不确定性和估计概述*

##### Pier Luigi Conti

罗马第一大学，罗马，意大利

##### Daniela Marella

罗马TRE大学，罗马，意大利

##### Mauro Scanu

意大利国家统计局，意大利国家统计研究所，罗马，意大利

------



### 目录

##### 4.1  简介

##### 4.2  统计匹配问题：符号和技术细节

##### 4.3  未共同观测变量的联合分布：估计和不确定性

##### 4.3.1  匹配误差

##### 4.3.2  通过不确定性度量来缩减匹配误差

##### 4.4  复杂抽样调查的统计匹配

##### 4.4.1  关于样本设计的技术假设

##### 4.4.2  关于选择匹配分布的建议

##### 4.4.3  匹配分布的可靠性

##### 4.4.4  作为一个假设问题的匹配可靠性的评估

##### 4.5  结论和有待解决的问题：统计匹配问题与生态推断之间的关系

##### 参考文献

------



### 4.1  简介

&emsp;&emsp;对于当今社会日益增长的信息需求提供快速答案的必要性，可以通过联合使用从同一人群中抽取的两个或多个样本调查来实现。在它们的联合使用中，根据推论来看有趣的是，对于大多数样本设计，两个样本观察到一个公共单位子集的概率是可以忽略的；因此，不可能重建成对变量的联合观测，其中一个变量只在一个样本中观测到，而在另一个样本不观测到，另一个变量则相反。这个问题通常被称为统计匹配，人们对它越来越感兴趣（参见[18]）。

&emsp;&emsp;长期讨论的统计匹配应用之一是在家庭层面的收入和支出数据框架内，决策者广泛用于提供对一些领域的见解（参见[34]，[43]；同时可以参见[14]及其参考文献）。首先，分析财政政策提案的影响通常需要有关家庭支出和收入的信息。这通常通过微观模拟模型来实现，该模型模拟家庭在税收改革或养老金制度变化的情况下的短期或长期支出/储蓄行为。

&emsp;&emsp;其次，大多数旨在提高生活水平的政策举措倾向于使用收入或消费来衡量贫困。然而，如果单独考虑，收入和消费都不是衡量贫困的充分标准。更好的方法是同时使用两者，特别是如果从达到的生活水平的角度考虑贫困。因此，商品和服务的消费可以说是经济福祉的更重要的决定因素，而不仅仅是收入。所有贫困和不平等措施通常都是收入和支出分配的分位估计的功能。

&emsp;&emsp;一个主要缺点是，大多数国家没有单一的微观数据来源，包括高质量的收入和支出分类联合数据，原因有二：避免回复负担和确保对问卷的高质量回答。因此，由于在同一项调查中询问有关收入和消费的详细问题可能是问题，因此调查往往只专注于两个主题中的一个。

&emsp;&emsp;在意大利，通过使用EU-SILC（欧盟收入和生活条件统计）调查提供的关于家庭收入的可靠信息来解决这个问题，该调查每年由国家统计局（意大利国家统计局）进行。另一方面，每年进行的ISTAT抽样调查（家庭预算调查，简称HBS）提供了有关消费支出的可靠和准确的信息。两个样本都是根据多级复杂采样设计进行选择的。ISTAT数据库分别提供收入和消费支出的相关信息。然而，它们不包含对*家庭收入和消费支出*的任何联合观察。然后一个自然的问题出现了：*在某种程度上，是否有可能在上述两个数据库的基础上估计收入和消费支出的共同分配？*简而言之，这就是统计匹配的问题。上述问题本质上等同于此：*在某种程度上，是否有可能（重新）构建一个包含收入和消费支出的独特数据库，而不需要进行特别的、昂贵的调查？*

&emsp;&emsp;这个问题通常可以通过对同一目标人群进行独立抽样调查的统计匹配来解决（参见[23]；同时可以参见[21]及其参考文献）。

&emsp;&emsp;尽管大多数申请都是关于收入和支出的，但这个问题将在官方统计中变得无处不在。削减成本和及时取得成果的需求日益增长，导致许多国家统计研究所面临现代化问题，其中包括来自不同来源的数据在一个独特的登记册中共存。这些寄存器增加了将统计匹配的使用扩展到所有那些变量对的可能性，对于这些变量对，不存在直接观察到的联合数据，并且无法重新创建属于同一统计单元的链接数据。此外，还出现了其他问题，比如在大数据可用的情况下，是否可能使用统计匹配技术参见[20]，[40]）。

&emsp;&emsp;关于这一主题，有大量的参考书目，可以追溯到20世纪60年代末。在这种情况下，两个突出的问题似乎是缺乏对感兴趣的变量对的联合观察，这导致了一种特定形式的不确定性；如何联合使用两个样本进行统计匹配，这两个样本可以根据不同的复杂调查设计绘制。本章分别在4.2和4.4中回顾了这两个问题的可用解决方案。最后，在4.5节中讨论了统计匹配问题和独立开发的研究问题（生态推理）之间的联系。

------



### 4.2  统计匹配问题：符号和技术细节

&emsp;&emsp;设 U 是 N 个单位的由整数 1,...,N 标记的有限总体，并且设 (X，Y，Z) 是感兴趣的三个变量。进一步用 $x _ { i } , y _ { i } , z _ { i }$ 分别表示  X， Y，Z  的值；对于单元 i （=1,...,N），并且依其次序构成向量 $x _ { N } , y _ { N } , z _ { N }$ ：

$$x _ { N } = ( x _ { 1 } \cdots , x _ { N } ) , y _ { N } = ( y _ { 1 } \cdots , y _ { N } ) , z _ { N } = ( z _ { 1 } \cdots , z _ { N } ).\tag{4.1}$$

&emsp;&emsp;(X，Y，Z) 在整个总体上的概括描述能力本质上等同于联合总体分布函数（简称p.d.f.）的概括描述能力 

 $$H _ { N } ( x , y , z ) = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } I _ { ( x _ { i }\leq x ) }{ I } _ { ( y _ { i } \leq y ) }{ I } _ { ( z_ { i } \leq z ) }, x,y,z\in \mathbb{R} \tag{4.2}$$ 

其中 $I _ { ( \phi ) }$ 是所有单元总体 $\phi$  的标志。设 

$$F _ { N } ( x , y ) = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } I _ { ( x _ { i \leq x } ) } I _ { ( y _ { i } \leq y ) },\tag{4.3}$$

$$G _ { N } ( x , z ) = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } I _ { ( x _ { i } \leq x ) } I _ { ( z _ { i } \leq z ) },\tag{4.4}$$

$$Q _ { N } ( x ) = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } I _ { ( x _ { i } \leq x ) }\tag{4.5}$$

分别为 (X，Y) , (X，Z) 的联合 p.d.f.s 和X的边际 p.d.f. ;又设

$$p _ { N } ( x ) = Q _ { N } ( x ) - Q _ { N } ( x ^ { - } ) = \frac { 1 } { N } \sum _ { i = 1 } ^ { N } I _ { ( x _ { i } = x ) }\tag{4.6}$$

为总体所有单元的分配比例，使得 X = x 。再进一步定义条件p.d.f.s

$$F _ { N } ( y | x ) = \frac { F _ { N } ( x , y ) } { p _ { N } ( x ) },\tag{4.7}$$

$$G _ { N } ( z | x ) = \frac { G _ { N } ( x , z ) } { p _ { N } ( x ) },\tag{4.8}$$

$$H _ { N } ( y , z | x ) = \frac { H _ { N } ( x , y , z ) } { p _ { N } ( x ) }\tag{4.9}$$

无论何时只要 $p _ { N } ( x ) \neq 0$ 则可以任意定义。

&emsp;&emsp;例如，在国家统计局定期收集进行的抽样调查中，经常可以获得使用关于 X, Y, Z  三个特征量的部分抽样信息。我们假设：

1.通过复杂抽样样本设计从有限总体 U 中选择 $n _ { A }$ 个单元的样本 $s _ { A }$ ，并观察值 $( x _ { i } , y _ { i } )，i\in s _ { A }$ ；

2.通过复杂抽样样本设计，从有限总体 U 中独立于 $s _ { A }$ 的部分选择 $n _ { B }$ 个单元的样本 $s _ { B }$ ，并观察值 $( x _ { i } , z _ { i } )，i\in s _ { B }$ ；

3.两个样本 $s _ { A }$ ,  $s _ { B }$ 是独立的。

&emsp;&emsp;由于至少在真实的社会调查中，在 $s _ { A }$ 和 $s _ { B }$ 中选择共有单元的概率基本上可以忽略不计（参见[42]），因此可以假设两个样本 $s _ { A }$ , $s _ { B }$ 不重叠。粗略地说，观测机制是这样的：（i）在 $s _ { A }$ 中只观察变量（X，Y），（ii）在 $s _ { B }$ 中只观察变量（X，Z）。变量 X 对样本 $s _ { A }$ , $s _ { B }$ 是共有的，并起到匹配变量的作用。因此，变量 X, Y ,Z 并不是共同观察到的。



------





### 4.3 未共同观测变量的联合分布：估计和不确定性

&emsp;&emsp;统计匹配旨在结合 $s _ { A }$ 和 $s _ { B }$ 中的信息。更准确地说，在“宏观”层面上，统计匹配的目标在于使用样本数据 

$$ \{ ( x _ { i } , y _ { i } ) ; i \in s _ { A } \} , \{ ( x _ { i } , z _ { i } ) ; i \in  s _ { B } \} \tag{4.10}$$

来估计联合p.d.f.(4.2)。

&emsp;&emsp;由于X，Y，Z不是联合观测的，除非做出特殊假设，否则（X，Y，Z）联合分布的统计模型通常是不可识别的。也就是说，不能基于可用的样本信息来估计联合p.d.f.(4.2)。从形式上讲，（X，Y，Z）的联合p.d.f.的估计本质上等价于（i） X 的边际分布的估计；（ii）（Y，Z）在 X 上有条件的联合分布。

&emsp;&emsp;由于没有对 Y 和 Z 的联合观测（既没有边际观测，也没对 X 有条件观测），唯一无法从样本数据中估计的数量是（4.9）中的 $H _ { N } ( y , z | x )$ 。这种可识别性的缺乏是（X，Y，Z）统计模型“内在”不确定性的原因。即使样本量 $n_ { A }$ , $n _ { B }$ 很大，也只能期望分别对（X，Y）和（X，Z）的二元分布进行合理精确的估计，而不能期望对（X，Y，Z）整个分布进行合理准确的估计。也就是说，统计模型仍然无法识别。

&emsp;&emsp;手头上的数据的一个可识别模型是假设给定X的Y和Z的独立性的模型。这个假设是众所周知的条件独立性假设（简称CIA）。根据条件独立性假设的规定，联合p.d.f.（4.2）可以按如下公式分解

$$H _ { N } ( y , z | x ) = F _ { N } ( y | x ) G _ { N } ( z | x )\tag{4.11}$$

有几篇论文讨论了条件独立性假设的恰当性。除此外，我们引用了[43]和[41]。

&emsp;&emsp;然后一个自然的问题出现了：当条件独立性假设不成立时，我们能对 $H _ { N } ( y , z | x )$ 提出些什么？

&emsp;&emsp;在没有任何关于Y和Z之间相关性的假设或先前信息的情况下，如果只知道p.d.f.s（4.3）--（4.5），那么只能说

$$m a x ( 0 , F _ { N } ( y | x ) + G _ { N } ( z | x ) - 1 ) \leq H _ { N } ( y , z | x ) \leq \min ( F _ { N } ( y | x ) , G _ { N } ( z | x )).\tag{4.12}$$

（4.12）中的界限是众所周知的Fréchet界限。它们分别表示 $H _ { N } ( y , z | x )$ 可以取的最小和最大逐点值。我们强调，在给定 X 的情况下，Fréchet下界对应于 Y 和 Z 之间的最大负关联；当且仅当给定 X时，（iff）Y 是 Z 的严格递减函数（反之亦然） 。类似地，给定 X 时，Fréchet上界对应于Y和Z之间的最大正关联；这是真的，当给定 X 时，Y是Z的严格递增函数时（反之亦然）。

&emsp;&emsp; **示例1--多正态分布--** 自[2]以来，这个例子已经得到了广泛的研究。设 X，Y，Z 为联合多正态分布，平均向量和协方差矩阵分别等于

$$
\mu=\begin{bmatrix}
\mu_{x}\\
\mu_{y}\\
\mu_{z}\\
\end{bmatrix},
\Sigma=\begin{bmatrix}
\sigma _{x}^2&\sigma _{x y}&\sigma _{x z}\\
\sigma _{x y}&\sigma _{y}^2&\sigma _{y z}\\
\\sigma _{x z}&\sigma _{y z}&\sigma _{z}^2\\
\end{bmatrix}\tag{4.13}
$$

 
在此条件下，X,（Y，Z）确实具有联合二元正态分布，其均值向量和协方差矩阵很容易从（4.13）中获得，由下式给出

$$
\mu  _ { y z | x }=\begin{bmatrix}
{ \mu _ { y } + \beta _ { y / x } ( x - \mu _ { x } ) }\\
{ \mu _ { z } + \beta _ { z / x } ( x - \mu _ { x } ) }\\
\end{bmatrix},\tag{4.14}
$$

$$
\\Sigma_ { y z | x }=\begin{bmatrix}
{ \sigma _ { y } ^ { 2 }(1-\rho _ { x y } ^ { 2 })}&{ \sigma _ { y } \sigma _ { z }(\rho _ { y x }-\rho _ { x y } \rho _ { x z })}\\
{ \sigma _ { y } \sigma _ { z }(\rho _ { y x }-\rho _ { x y } \rho _ { x z }) }&{ \sigma _ { z }^ { 2 }(1-\rho _ { x z } ^ { 2 }) }\\
\end{bmatrix}\tag{4.15}
$$

&emsp;&emsp;设 $\rho _ { xy }$ , $\rho _ { xz }$ , $\rho _ { yz }$ 分别是（X，Y）,（X，Z）,（Y，Z）之间的相关系数。唯一未确定的参数是 $\rho _ { yz }$ ，即 Y 和 Z 之间的相关系数。

&emsp;&emsp;从（4.15）中可以立即看出，给定 $\rho _ { xy }$ 和 $\rho _ { xz }$ , $\rho _ { yz }$ 在区间中的范围

$$[ \rho _ { xy }\rho _ { xz } - \sqrt { ( 1 - \rho  _ { xy }^ { 2 } ) ( 1 - \rho _ { xz }^ { 2 } ) },\rho _ { xy }\rho _ { xz }  + \sqrt { ( 1 - \rho  _ { xy }^ { 2 } ) ( 1 - \rho _ { xz }^ { 2 } ) } ].\tag{4.16}$$.

（4.16）中的所有值都同样合理。请注意，在条件独立性假设下，参数 $\rho _ { yz }$ 位于区间（4.16）的中点，如[32]所述。

&emsp;&emsp;根据Slepian不等式[45]，在给定 X 的情况下，Y，Z 的联合 d.f. 是给定 X 的 Y 和 Z 之间的相关系数的递增函数：

$$\rho  _ { y z | x } = \frac { \rho _ { y z } - \rho _ { x y } \rho _ { x z } } { \sqrt { ( 1 - \rho _ { x y }^ { 2 } ) ( 1 - \rho _ { x z }^ { 2 } ) } }$$

即它被证明是 $\rho _ { y z }$ 的单调函数。换句话说，条件 d.f.s   $H _ { N } ( y , z | x )$ 是在 $\rho _ { y z }$ 的基础上完全有序的。  $\rho  _ { y z } = \rho  _ { x y } \rho  _ { x z } - \sqrt { ( 1 - \rho  _ { x y } ^ { 2 } ) ( 1 - \rho  _ { x z } ^ { 2 } ) }$ 对应于 $\rho _ { y z | x } = - 1$ ，因此在条件限定 X 下， Y 是 Z 的线性递减函数（反之亦然）。因此，当唯一可用的信息是数据时，Fréchet下界是 $\max ( 0 , F _ { N } ( y | x ) + G _ { N } ( z | x ) - 1 )$ 。类似地， $\rho _ { y z } = \rho _ { x y } \rho _ { x z } + \sqrt { ( 1 - \rho _ { x y } ^ { 2 } ) ( 1 - \rho _ { x z } ^ { 2 } ) }$ 的情况对应于 $\rho _ { y z | x } = 1$ ，即对应于Fréchet上界 $\min  ( F _ { N } ( y | x ) , G _ { N } ( z | x ) )$ 。


P85——P90 程佳悦2021211143010
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
