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


#### 以下部分为史雯杉所写（P79-P84）

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

且对于固定的 $y$ 和 $z$ 来说， $\gamma_{y}(\cdot )$ ， $\delta_{z}(\cdot )$ 就分别是 $f_{x}(y,z)$ 的反函数。当

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

#### 4.3.1 匹配误差  
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
&emsp;&emsp;从现在起，为了评估匹配分布 $H_{N}^{\ast }(y,z|x)$ 的准确性，作为真实概率分布函数 $H_{N}(y,z|x)$ 的估计量，作为匹配误差度量（有条件地基于 $X$ ），我们将使用以下内容：

$$ME_{x}(H_{N}^{\ast },H_{N})={\int _{R^2}}|H_{N}^{\ast }(y,z|x)-H_{N}(y,z|x)|d\lbrack F_{N}(y|x)G_{N}(z|x)\rbrack .\tag{4.28}$$

&emsp;&emsp;通过类似的推理，作为匹配误差的无条件度量，我们可以考虑以下内容：

$$ME(H_{N}^{\ast },H_{N})=\int _{R}ME_{x}(H_{N}^{\ast },H_{N})dQ_{N}(x)=\int _{R}\left\\{\int _{R^2}|H_{N}^{\ast }(y,z|x)-H_{N}(y,z|x)|\times d\lbrack F_{N}(y|x)G_{N}(z|x)\rbrack \right\\}dQ_{N}(x).\tag{4.29}$$

#### 4.3.2 通过不确定性测量来限制匹配误差  
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


#### 以下部分为程佳悦2021211143010所写（P85——P90）

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





D'Orazio 等人。[<font color=Blue>24</font>]比较它们的效率，包括基于最大伪似然的额外估计，如[<font color=Blue>48</font>]中的那些。无论如何，在有限的样本量下，这些估计量中没有任何突出的估计量。

当根据复杂的调查设计抽取样本时，[<font color=Blue>18</font>]给出了一种可能性。他们将注意力限制在粗略地说⾼熵的抽样设计上（实际假设在第 $4.4.1$节中）。在此设置中，可以通过迭代比例拟合(IPF) 算法（第$4.4.2$节）获得感兴趣变量的合理联合分布的存在。孔蒂等人。 [<font color=Blue>18</font>]通过根据渐近特性（第  $4.4.3$节）分析估计量的可靠性，克服了为有限样本量建立统计匹配的最佳方法的问题允许使用一些额外的工具，作为不确定宽度测试的定义。

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



------------------------------------
**以下内容是郭慧慧所写（P91-100）**

统计匹配中的不确定性与估计

也就是说，$H _ { H } ( y , z | x )$可以改写为:

$$ \tilde{H}_{H}(y,z|x)=\frac{\sum _ { i \in s _ { A } } \sum _ { j \in s_{B} } I _ { ( a _ { x } \leq f _ { x } ( y _ { i } , z _ { j } )\leq b_{x}) } \pi _{i,A}^{-1}\pi _{j,B}^{-1}I_{(x_{i}=x)}I_{(x_{j}=x)}I_{(y_{i}\leq y)}I_{(z_{j}\leq z)} }{\sum _ { i \in s _ { A } } \sum _ { j \in s_{B} } I _ { ( a _ { x } \leq f _ { x } ( y _ { i } , z _ { j } \leq b _ { x } )} \pi _{i,A}^{-1}\pi _{j,B}^{-1} I _ { ( x _ { i } = x ) }I _ { ( x _ { j } = x ) }}$$

其基本思想是应用IPF算法，以$ \tilde{H}_{H}(y,z|x)$为“起点”，并通过交替地重新比例其边缘w.r.t.(4.48)。注意起始分布$ \tilde{H}_{H}(y,z|x)$一个模拟分布条件独立假设(由于选择与相对于$ \widehat{F}_{H}(y|x)$和$ \widehat{G}_{H}(z|x)$分布的乘积)，而不是表示CIA，由于$ \tilde{H}_{H}(y,z|x)$在约束区域，而边际分布在相应的域中不一定为零。的边际分布也要注意$ \tilde{H}_{H}(y,z|x)$分别不是$ \widehat{F}_{H}(y|x)$和$ \widehat{G}_{H}(z|x)$。因此，指规数算法沿$ \widehat{F}_{H}(y|x)$和$ \widehat{G}_{H}(z|x)$方向应用，以求得解估计匹配分布$ \widehat{H}_{H}^{*}(y,z|x)$，即，一个d.f.有同样的时间边缘(4.48)和满足约束(4.17)。
Conti等人[18]建议使用估计量(4.48)，因为在假设A1-A6，它们倾向于在概率上有相同的极限估计量是基于泊松抽样设计，在最后的情况下估计量的渐近正态性可以很容易地得到。

**4.4.3匹配分布的可靠性**

为了评估匹配分布的准确性，我们继续估计不确定性度量（4.31）和（4.33）。
关于不确定性的条件测度（4.31），一个公平同步的方法包括首先估计（4.9）给出的条件d.f.s$F_{N}(y|x)$，$G_{N}(z|x)$，然后将这些估计插入（4.31）中。利用估计量（4.48），得到了条件测度数确定性$ \Delta ^{x}(F_{N},G_{N})$的以下估计量

$$ \widehat{\Delta}_{H}^{x}= \int _{R^{2}}(\widehat{K}_{H+}^{x}(y,z)- \widehat{K}_{H-}^{x}(y,z))d \left[ \widehat{F}_{H}(y|x)\widehat{G}_{H}(z|x)\right] 
= \frac{1}{\widehat{N}_{A}(x)\widehat{N}_{B}(x)}\sum _{i,j=1}^{N}(\widehat{K}_{H+}^{x}(y_{i},z_{j})- \widehat{K}_{H-}^{x}(y_{i},z_{j}))\frac{D_{i,A}}{\pi _{i,A}}\frac{D_{i,B}}{\pi _{i,B}}I_{(x_{i}}=x){I}_(x_{j}=x)\tag{4.50}$$

其中$ \widehat{K}_{H-}^{x}$，$ \widehat{K}_{H+}^{x}$分别精确地定义为(4.19)和(4.20)，但用$ \widehat{F}_{H}, \widehat{G}_{H}.$代替$F_{N},G_{N}$。

Conti等人[18]证明了定理1，其中得到了估计量$ \widehat{\Delta}_{H}^{x}$的大样本分布。让我们定义一下量

$$ \widehat{p}_{H,A}(x)= \frac{\sum _{i=1}^{N}\frac{D_{i.A}}{n_{iA}A}(x_{i}=x)}{\sum _{i=1}^{N}\frac{D_{i,A}}{\pi _{i,A}A}}, \widehat{p}_{H,B}(x)= \frac{\sum _{i=1}^{N}\frac{D_{i,B}}{\pi _{i,B}}I_{i_{i}-x})}\tag{4.51})$$

$$ \widehat{n}_{H,A}(x)= \frac{N \widehat{p}_{H_{i}}(x)}{\frac{1}{N}\sum _{i=1}^{N}\pi _{i,A}^{-1}-1}, \widehat{n}_{H,B}(x)= \frac{N \widehat{p}_{H.B}(x)}{\frac{1}{N}\sum _{i=1}^{N}\pi _{i,B}-1}\tag{4.52}$$

其中的$ \widehat{p}_{H,A}(x)$代表分别表示从两个样本$S_{A},S_{B}$中得到的$p_{N}(x)$的Hajek估计量。
因为，如[18]所示

$$ \widehat{n}_{H,A}(x)(\frac{n_{A}p(x)}{f_{A}(\zeta _{A}-1)})^{-1}\rightarrow 1, \widehat{n}_{H,B}(x)(\frac{n_{B}p(x)}{f_{B}(\zeta _{B}-1)})^{-1}\rightarrow 1 \tag{4.53}$$

随着N的增加，可以在定理1中建立出$ \widehat{\Delta}_{H}^{x}$的渐近正态性。
定理1：让$x \in \left\{ x^{1}, \ldots ,x^{K}\right\}$，并假设它

$$\frac{ \frac{n_{A}p(x)}{f_{A}(\varsigma_{A}-1)}\frac{n_{B}p(x)}{f_{B}(\varsigma_{B}-1)} }{\frac{n_{A}p(x)}{f_{A}(\varsigma_{A}-1)}+ \frac{n_{B}p(x)}{f_{B}(\varsigma_{B}-1)} } \rightarrow \alpha   as   N \rightarrow \infty \tag{4.54}$$

然后，对于几乎所有的$(x_{i},y_{i},z_{i})s$的值，有条件地在$x_{N},y_{N},z_{N}$上，随着N的增加：

$$\sqrt{\frac{\widehat{n}_{H,B}(x)\widehat{n}_{H,B}(x) }{\widehat{n}_{H,A}(x)+ \widehat{n}_{H,B}(x) }}  (\widehat{\Delta}_{H}^{x}- \Delta ^{x}(F_{N},G_{N}))\rightarrow N(0,V(F,G;x)) \tag{4.55}$$

$N(0,V)$表示在[18]中给出均值0和方差 为$V(F,G;x)$的正态分布r.v.  对于不确定性的无条件度量(4.33) ，其思想是使用形式为 $p_{N}(x)$的估计量：

$$ \widehat{p}_{\tau}(x)= \tau \widehat{p}_{H,A}(x)+(1- \tau)\widehat{p}_{H,B}(x)\tag{4.56}$$

对于$0 \leq \tau \leq 1$，根据两个观察到的样本已经被定义为（4.51）。$\tau$的一个相当自然的选择是使（4.56）的渐近方差最小化的选择。
因此，得到了所提出的无条件度量（4.33）的估计量，并对X的所有的类别：$x^{k}k=1, \ldots ,K of X$求和：

$$\widehat{\Delta}_{H} = \sum _ { k = 1 } ^ { K } \Delta _ { H } ^ { x ^ { k } } \widehat{p}_{H}AB(x^{k})=  \overline{p}_{H,AB}(x)^{\prime}\widehat{\Delta}_{H}^{x}\tag{4.57}$$

这从本质上模仿了（4.33）的结构。
统计匹配中的不确定性与估计

定理2：设

$$ \widehat{n}_{A}\quad = \frac{n_{A}}{\frac{n_{A}}{N}(N^{-1}\sum _{i=1}^{N}\pi _{i,A}^{-1}-1)}= \frac{N}{N^{-1}\sum _{i=1}^{N}\pi _{i,A}^{-1}-1}\tag{4.58}$$

$$ \widehat{n}_{B}\quad = \frac{n_{A}}{\frac{n_{B}}{N}(N^{-1}\sum _{i=1}^{N}\pi _{i,B}^{-1}-1)}= \frac{N}{N^{-1}\sum _{i=1}^{N}\pi _{i,B}^{-1}-1}\tag{4.59}$$

对于几乎所有的$( x _ { i } , y _ { i } , z _ { i } ) s$的值，有条件地在$x_{N},y_{N},z_{N}$上，随着N的增加，以下结果成立：

$$ \sqrt{\frac{\widehat {n}_{A}\widehat {n}_{B}}{\widehat{n}_{A}+ \widehat{n}_{B}}} (\widehat{\Delta}_{H}- \Delta(F_{N},G_{N}))^{d}\rightarrow N(0,V(F,G))\quad asN \rightarrow \infty\tag{4.61}$$

以及

$$V(F,G)= \sum _{k=1}^{K}p(x^{k})V(F,G;x^{k})+ \frac{(\zeta _{A}-1)(\zeta_{B}-1)}{(\zeta_{A}+\zeta_{B}-2)^{2}}\Delta ^{x}(F,G)^{\prime}\sum \Delta ^{x}(F,G).(4.61)$$$$V(F,G)= \sum _{k=1}^{K}p(x^{k})V(F,G;x^{k})+ \frac{(\zeta _{A}-1)(\zeta_{B}-1)}{(\zeta_{A}+\zeta_{B}-2)^{2}}\Delta ^{x}(F,G)^{\prime}\sum \Delta ^{x}(F,G) \tag{4.61}$$

在[14]中证明，估计量$A _ { H } ^ { x }$和估计量$ \widehat{\Delta}_{H}$是渐近一致的（在Brewer定理中）。

**4.4.4匹配可靠性的评估作为一个假设问题**

渐近结果（4.55）和（4.60）允许我们分别对不确定性度量（4.31）和（4.33）定义检验假设步骤。考虑到条件不确定性度量$ \Delta ^{x}(F_{N},G_{N})$可以解释为当真实的总体分布函数被$H _ { N } ^ { x }$4.27）中的匹配分布取代时发生的误差的上界，这些测试程序是有效的。例如，是否不确定性实际上可以视为零，或限制可忽略的数量，或者如果数据不支持这一假设。设ex besuch正数；例如，它可以是$ \Delta ^{x}(F_{N},G_{N})$的最大值的一个百分比，即$1/6 \approx 0.17$（见[16]）。当$ \Delta ^{x}(F_{N},G_{N})$小于$ \epsilon ^{x}$时，保守的方法会认为匹配分布是“可靠的”分布。
对匹配分布的可靠性的评估可以通过检验以下假设来处理：

$$\left\{ \begin{matrix} H_{0}: \quad \Delta ^{x}(F_{N},G_{N})\leq \epsilon ^{x}\\ H_{1}: \quad \Delta ^{x}(F_{N},G_{N})> \epsilon ^{x}\\ \end{matrix} \right.\tag{4.62}$$

给定(渐近)显著水平a，如果满足以下条件，则零假设$H_{0}$被接受

$$ \widehat{\Delta}^{x}\leq \epsilon ^{x}+z_{\alpha}\sqrt{\overline{V}_{x}}(\frac{\widehat{n}_{H,A}(x)\widehat{n}_{H,B}(x)}{\widehat n_{H,A}(x)+ \widehat{n}_{H,B}(x)})^{-1/2}$$

其中，$z_{\alpha}$是标准正态分布的路径分位数，$ \widehat{V}_{x}$是方差的$\widehat{V}(F,G;x)$的适当估计量。
注意，前面公式中涉及的渐近方差是很尴尬的。一种建议是通过通常的插件规则$F_{N}(y|x)$分别用$G_{N}(z|x)$以及$\widehat{F}_{H}(y|x)$和$\widehat{G}_{H}(z|x),$，或者用复杂设计的重采样方法（cfr.[3]，[13]）来构造估计量$\widehat{V}_{x}$。显然，对于无条件的不确定性测量（4.33），使用（4.60）也有类似的考虑。

注2.在假设的情况下，也就是说，Bootstrap的经典版本可以用来估计方差，如在[18]中所应用的，它包括：(I)产生大小为$n_{A}$的样本$n_{A}$；(ii)产生大小为$n_{B}$的样本$n_{B}$；显然，经典的Bootstrap方法不适用于复杂的调查抽样。Conti等人[18]还详细介绍了来自有限总体的重抽样方法的参考文献，并提到了两种主要方法：即席方法([31]，[37]，[44]，[4]，[11]，[13])和插入式方法([26]，[10]，[28]，[36])。

### 4.5 结论和待决问题：统计匹配问题与生态学推论之间的关系

随着不同的数据源开始共享一个共同的元数据框架，统计匹配将是一个无处不在的问题，并且在一个独特的地方，可以看到既没有联合操作服务的数据，也不能通过单元标识符链接的数据。虽然第一次应用是在20世纪60年代末，但由于重要问题，统计匹配应用不是很复杂：除非引入特定的假设，否则数据不能为描述从未共同观察过的变量之间的联合关系的参数提供点估计；数据可以通过不同的复杂调查设计获得，而且不能直接调和它们。 正如在第4.2节中所介绍的，统计匹配解决了推断两个从不联合服务的变量的联合分布的问题。这个问题已经在一个不同的框架中进行了独立的分析：生态学推论。使这个框架特定的经济学推理：无论如何，这两个领域之间的接触点是如此之多，以至于一起分析它们并可能交换方法和解决方案似乎很重要。“生态学推论”利用人口中群体聚集形式的数据，得出关于个体层面关系的结论。这些群体通常在地理上被定义为“[46]”。生态学这个术语似乎是由于这些群体的存在，以及他们通常是根据地理位置设计的事实。事实上，统计匹配并不涉及任何类型的分组，而是使用公共变量X来证实任何类型的对从未共同观察到的变量的推断。尽管有这些相似之处，但在这两种情况下，解决方案并没有被共享。这也由典型的生态学推论例子所代表：在选举时，两个政党的选举人票是已知的，以及一些选区根据种族（白人、非白人）划分的选民人数。 因此，需要建立关于种族和选举投票率扫描的双向列联表，其中边界已知，表内的单元格未知。这正是用X表示区域的统计匹配情况。正如在统计匹配环境中已经指出的，给定符合统计匹配符号的precincts，i.e.（Y，ZX）的种族和选票的联合分布是不可识别的分布。这也是生态学推论的兴趣对象。

给定区域的列联表中的百分比，在一无所知的情况下，是介于0和1之间的值。当每个区域的边界已知时，列联表中百分比的可容许值的空间减少到一条线，即 King [29]称为断层扫描线的那条线。这对应于减少统计匹配的不确定性对二项随机变量给定他们的条件边界。一个细微的差别是，生态学推论的重点是在群体(分区)水平上已经汇总的数据。因此，King [29]假设模型在每个区域内生成条件百分比，并在其中工作。相反，统计匹配使用直接模拟任何表格中的任何单个观察([46]将这种方法形式化，区分个人和总体水平的抽样)。此外，给出了群 X，而在统计匹配中，它们被一起建模，感兴趣的变量和有时感兴趣的对象的统计匹配只是成对分布(Y，Z)。请注意，关于 Y 和 Z 的联合分布的不确定性的概念与 X 与 Y 和 Z 之间的关联严格相关。可能，这种轻微的观点变化在两个领域产生了独立的解决方案，因此是时候在同一框架下进行两个领域的研究了。例如，在生态学推论中，King [29]对每个分区的表格内的百分比做出了假设，可以这样总结:

1.每个区域的YZ的两个条件百分比共同遵循特定的分布（例如，截断的正态分布）

2.根据观察到的种族百分比，观察到的在不同选区的选民的每个比例是平均独立的（没有空间自相关） 

3.观察到的种族边际比例独立于点1中建模的两个条件百分比。 

在这种情况下，可以对给定种族的选民的条件概率进行最大似然分析或贝叶斯分析。相反，统计匹配已经尽可能地发展为无假设框架（或者更好的是，关于给定X的Y和Z的连续独立性的初始假设显然不令人满意）。统计匹配框架也可以建立在合理的假设之上了。在生态学推论或其他领域（例如，小区域模型）中，当数据不能支持分析的目标时，假设并不一定是邪恶的，即使在大多数没有假设的官方统计框架中也是如此。因此，生态学推论中发展的方法也应应用于统计匹配目的。一些方法上的挑战仍然存在，这只是一个初步的清单：

来自生态学推论的方法表明，可以重新考虑两个样本调查中共同变量的作用：从共同变量来看，一方面可以考虑协变量（即匹配变量），以及第二阶段的组（阶层）。第二组是建立聚合数据建模的可能性的基础； 似乎生态学推论已经发展为二元Y和Z，而统计匹配不限制可变变量匹配的性质。 我们认为，这是在这两个领域中交换和重用方法应该遵循的路线。

### 参考文献

[1] A. Agresti and M.C. Yang. An empirical investigation of some effects of sparseness in contingency tables. Computational Statistics & Data Analysis, 5:9–21, 1987.

[2] T.W. Anderson. Maximum likelihood estimates for a multivariate normal distribution when some observations are missing. Journal of theAmerican Statistical Association, 52:200–203, 1957.  Uncertainty and Estimation in Statistical Matching 97  

[3] E. Antal and Y. Tille. A direct bootstrap method for complex sampling ´designs from a finite population. Journal of the American Statistical Association, 106:534–543, 2011.

[4] J-F. Beaumont and Z. Patak. On the Generalized Bootstrap for Sample Surveys with Special Attention to Poisson Sampling. International Statistical Review, 80:127–148, 2012.

[5] Y.G. Berger. Asymptotic consistency under large entropy sampling designs with unequal probabilities. Pakistan Journal of Statistics, 27:407– 426, 2011.

[6] Y.G. Berger. Rate of convergence to normal distribution for the HorvitzThompson estimator. Journal of Statistical Planning and Inference, 67:
209–226, 1998.

[7] W. Bergsma. A bias-correction for Cramer’s V and Tschuprow’s. Journal of the Korean Statistical Society, 42:9–21, 2013.

[8] L. Breiman. Random forests. Machine Learning, 45:5–32, 2001.

[9] L. Breiman, J. H. Friedman, R. A. Olshen, and C. J. Stone. Classification and Regression Trees. Wadsworth, 1984.

[10] M.-T. Chao and S.-H. Lo. A bootstrap method for finite population. Sankhya¯, pages 399–405, 1985.

[11] A. Chatterjee. Asymptotic properties of sample quantiles from a finite population. Annals of the Institute of Statistical Mathematics, 63:157–179, 2011.

[12] P.L. Conti. On the estimation of the distribution function of a finite population under high entropy sampling designs, with applications. Sankhya¯ B. DOI: 10.1007/s13571-014-0083-x, 2014.

[13] P.L. Conti and D. Marella. Inference for quantiles of a finite population: asymptotic vs. resampling results. Sandinavian Journal of Statistics,
2014.

[14] P.L. Conti, D. Marella, and A. Neri. Statistical matching and uncertainty analysis in combining household income and expenditure data. Statistical Methods & Applications, 26:485–505, 2017.

[15] P.L. Conti, D. Marella, and M. Scanu. Evaluation of matching noise for imputation techniques based on nonparametric local linear regression estimators. Computational Statistics and Data Analysis, 53:354–365, 2008.

[16] P.L. Conti, D. Marella, and M. Scanu. Uncertainty analysis in statistical matching. Journal of Official Statistics, 28:69–88, 2012. 98 Analysis of Integrated Data

[17] P.L. Conti, D. Marella, and M. Scanu. Uncertainty analysis for statistical matching of ordered categorical variables. Computational Statistics and
Data Analysis, 68:311–325, 2013.

[18] P.L. Conti, D. Marella, and M. Scanu. Statistical matching analysis for complex survey data with applications. Journal of the American Statistical Association, 111:1715–1725, 2016.

[19] P.L. Conti, D. Marella, and M. Scanu. How far from identifiability? a systematic overview of the statistical matching problem in a non parametric framework. Commmunications in Statistics - Theory and Methods, 46:967–994, 2017.

[20] T. De Waal. Statistical matching: experimental results and future research questions. Discussion paper n. 19, CBS, 2015.

[21] G. Donatiello, M. D’Orazio, D. Frattarola, A. Rizzi, Scanu M., and M. Spaziani. The role of the conditional independence assumption in
statistically matching income and consumption. Statistical Journal of the IAOS, 77:667–675, 2016.

[22] M. D’Orazio, M. Di Zio, and M. Scanu. Auxiliary variable selection in a statistical matching problem. In L.-C. Zhang and R. L. Chambers, editors, Ecological Inference: new methodological strategies. CRC: Chapman and Hall, London.

[23] M. D’Orazio, M. Di Zio, and M. Scanu. Statistical matching: theory and practice. Wiley, Chichester, 2006.

[24] M. D’Orazio, M. Di Zio, and M. Scanu. Uncertainty intervals for nonidentifiable parameters in statistical matching. Proceedings of the 57th
Session of the International Statistical Institute World Congress, Durban - South Africa, 2009.

[25] B. Efron. Bootstrap methods: another look at the jackknife. The Annals of Statistics, 7:1–26, 1979.

[26] S.T. Gross. Median estimation in sample surveys in: Proceedings of the section on survey reasearch methods. Proceedings of the Section on Survey Research Methods, American Statistical Association, pages 181–184, 1980.

[27] J. Hajek. Asymptotic theory of rejective sampling with varying proba- ´ bilities from a finite population. The Annals of Mathematical Statistics, 35:1491–1523, 1964.

[28] A. Holmberg. A bootstrap approach to probability proportional-to-size sampling. Proceedings of the ASA Section on Survey Research Methods, pages 378–383, 1998.
Uncertainty and Estimation in Statistical Matching 99

[29] G. King. A solution to the ecological inference problem: reconstructing individual behavior from aggregate data. Princeton University Press, Princeton, 1997.

[30] D. Marella, M. Scanu, and P.L. Conti. On the matching noise of some nonparametric imputation procedures. Statistics and Probability Letters, 78:1593–1600, 2008.

[31] P J. McCarthy and C B. Snowden. The bootstrap and finite population sampling. In Vital and health statistics, pages 1–23. Public Heath Service Publication, U.S. Government Pronting, Washington, DC, 1985.

[32] C. Moriarity and F. Scheuren. Statistical matching: A paradigm for assessing the uncertainty in the procedure. Journal of Official Statistics, pages 407–422, 2001.

[33] R.B. Nelsen. An introduction to copulas. Springer, New York, 1999.

[34] B.A. Okner. Constructing a new data base from existing microdata sets: the 1966 merge file. Annals of Economic and Social Measurement, 1:325– 342, 1972.

[35] S. Raessler. Statistical matching: A Frequentist Theory, Practical Applications, and Alternative Bayesian Approaches. Springer, New York, 2002.

[36] M. G. Ranalli and F. Mecatti. Comparing recent approaches for bootstrapping sample survey data: a first step towards a unified approach. In Proceedings of the ASA Section on Survey Research Methods, pages 4088–4099, 2012.

[37] J N K. Rao and C F J. Wu. Resampling inference with complex survey data. Journal of the American Statistical Association, 83:231–241, 1988.

[38] J. Reiter. Using multiple imputation to integrate and disseminate confidential microdata. International Statistical Review, 77:179–195, 2009.

[39] R.H. Renssen. Use of statistical matching techniques in calibration es. Survey Methodology, 24:171–183, 1998.

[40] D. Rivers. Sampling fom web surveys. In Proceedings of the ASA Section on Survey Research Methods, 2007.

[41] W. L. Rodgers. An evaluation of statistical matching. Journal of Business and Economic Statistics, pages 91–102, 1984.

[42] D B. Rubin. Statistical matching using file concatenation with adjusted weights and multiple imputations. Journal of Business and Economic Statistics, 4:87–94, 1986. 100 Analysis of Integrated Data

[43] C.A. Sims. Comments on: “Constructing a new data base from existing microdata sets: the 1966 merge file”, by B.A. Okner. Annals of Economic and Social Measurements, 1:343–345, 1972.

[44] R R. Sitter. A resampling procedure for complex data. Journal of the American Statistical Association, 87:755–765, 1992.

[45] D. Slepian. The one-sided barrier problem for gaussian noise. Bell System Technical Journal, 41:463–501, 1962.

[46] D.G. Steel, E.J. Beh, and R. L. Chambers. The information in aggregate data. In G. King, O. Rosen, and M.A. Tanner, editors, Ecological Inference: new methodological strategies. Cambridge University Press, Cambridge, 2004.

[47] Y. Tille.´ Sampling Algorithms. Springer, New York, 2006.

[48] C. Wu. Combining information from multiple surveys through the empirical likelihood method. Canadian Journal of Statistics, 32:15–26, 2004.

[49] L.-C. Zhang and R. L. Chambers. Minimal inference from incomplete 2 × 2-tables. In L.-C. Zhang and R. L. Chambers, editors, Analysis of Integrated Data. CRC: Chapman and Hall, London.

















