***以下部分为王若琳同学的作业：

​      &emsp;我们顺便注意到每个清单变量与X的相关性，可以用一个简单的形式来表示覆盖不足和覆盖过多的比率：

$$\left. \begin{array}  { l  }  { C o v [ Y , X ] = \pi _ { 11 } ^ { Y|X } - \pi _ { 1 } ^ { X } \pi _ { 1 } ^ { Y }\\ { = \pi _ { 1 } ^ { X } \pi _ { 0 } ^ { X } ( 1 - \pi _ { 0|1 } ^ { Y | X }- \pi_  { 1|0 } ^ { Y|X } } ) } \end{array} \right.$$ 

​        &emsp;因此，如果源Y的覆盖率不足率和覆盖率过高率的总和接近于1，那么Y通常会与X弱相关，并且它对模型的贡献会很差。如果错误率之和超过1，则相关性将是负的，并且Y的作用可能混杂在潜在变量的识别和解释中。因此，如果我们知道一个源的质量很差，有很高的错误率，在某些情况下，从模型中排除它可能是最好的选择。

​        &emsp;可以涵盖这两个方面(捕获中未观察到的异构性概率和真/假捕获)在同一模型中使用多个潜在变量或多个潜在类。在潜在类别模型中使用多个潜在变量早在[14]中就有暗示(另见[16]，[4])。包含两个(或更多)潜在变量允许我们指定更复杂的变量依赖结构在我们的模型和调整在一个更精细的交互道路。然而，人们应该识别潜在变量，它们各自的模态的数量，以及它们在图g中的依赖性。所以，除非我们有一个附加潜在变量所反映的主题的特定知识一个关于产生数据机制的具体假设，一个更简单的解决方案是否考虑一个潜在变量与许多潜在类别更大超过两个。这样一来，我们就必须选择最好数量的类，然后通过标记它们代表真捕获或假捕获，将潜在类别分为两组。设g为所选的潜数类 $x_1,\cdots,x_g$ ,设 $w_1$ 和 $w_2$ 是表示true的类集合和假案例。然后，属于后验概率剖面为 $\vec{y}$ 的单位的目标人口数为:

$$ \frac{\sum_{i\in W _{1}}\pi_{x_{i},y}}{\sum_{i=1}^g\pi_{x_{i},y}} $$ 

​       &emsp; 这两种方法的应用示例请参见[10]中的医学诊断测试，[18]为一个记录联动。

# 8.5 贝叶斯方法

​        &emsp;在本节中，我们将介绍我们的模型的贝叶斯方法。我们将只关注可分解模型，因为这些模型的贝叶斯方法很简单，而一般情况下会出现一些计算困难。

​        &emsp;贝叶斯方法的一个明显优势是可以以一种简单的方式包含手头数据的先验知识。我们可以在对参数值的信念建模之前设置一个信息，而不是像8.2.4节中所见的那样固定它，我们可以平滑地调整我们的置信度。此外，信息先验可以帮助我们确定所需的潜在亚群体。

​        &emsp;在贝叶斯方法中，我们可以很容易地计算群体大小的区间估计，这在捕获中一直很难于重新研究。在非贝叶斯背景中，我们被限制于自举方法，这是计算密集型的，因为潜在变量是的。两个数学家[31]扩展了[7]的结果，提出了一种基于对数线性模型中使用剖面对数似然的方法来有条件地估计少计数的置信区间协变量的值。然而，他们的方法不能轻易地扩展到计算总少计的置信区间。相反，在贝叶斯方法中，当我们检查其后验分布时，自然地获得了总体大小的区间估计。

​        &emsp;另一个有利的观点是，可以使用多种工具进行解释模型的不确定性。特别地，我们在模型平均中得到了一些结果以及专门为可分解对数线性模型设计的模型选择(见[20]，[21])。

​         &emsp;对于可分解模型的先验分布，我们有两个建议。第一个是一类先验，称为超狄利克雷它是共轭的模型(8.4)，具有在边缘化下封闭的性质(见[9])。对于每一个极大团 $\mathcal C$ ,称为 $\mathcal h$ 变量 ；对于每一个分布 $Π^c$ ,设置一个有参数的(h−1)维单形设狄利克雷函数 $\alpha_{y_c}$ 为每个可能的值组合定义 $\vec{y_c}\in\lbrace{0,1}\rbrace^h$ 在极大团 $\mathcal C$ 中。这些先验不是独立的，因此，为了保证任何交点的边缘分布是一致的，我们对参数施加如下限制:
 
 $$ \sum\alpha_{y_{C_{i-1}}}=\alpha_{S_{{i}}}=\sum\alpha_{y_{c_ {i}}}\quad i = 2 ,\cdots,g, $$
 
​         &emsp;其目的在于对变量 $\vec{y_{S_{i}}}$ 进行一致的组合。这类先验的应用可以在[21]中找到。   

​        &emsp;在第二种方法中，我们引用表达式(8.5)。对于每一个 $Π^{C_i|S_i}$ 和每一个 $S_i$ 的固定值，我们设置一个狄利克雷分布限制。这些分布是相互独立的，就像 $Π^{C_i|S_i}$ 独立于结构之中，不难看出这类先验与(8.5)共轭。这种方法被用于局部独立下潜在类模型的贝叶斯分析(8.2)。(见[41])

​        &emsp;示例如下，思考模型 ${[ABX][CDX]}$ , $C_1 = \{A,B,X\},C_2=\{C,D,X\},S_2=\{X\}$ 第一种方法将设置两个Dirichlet分布:

$$\pi ^ { A B X } \sim D i r ( \alpha _ { a b x } ^ { A B X } ) \quad  \quad \pi ^ { C D X } \sim  ( \alpha _ { c d x } ^ { C D X } )$$

 &emsp;限制如下: $\alpha _ { x } = \sum _ { a , b } \alpha _ { a b x } = \sum _ { c , d } \alpha _ { c d x }$ ，
 
 &emsp;目的是为了确保 $Π^X$ 会 在这两种情况中保持一致。从 $\sum _ { a , b } n _ { a b x } = \sum _ { c , d } n _ { c d x }$ 尾部的公式会一直持续。

​        &emsp;在第二种方法中，我们将设置四个狄利克雷分布:

$$\pi ^ { A B | X } \sim D i r ( \alpha _ { a b | x } ) ( x = 0 , 1 ) , \quad \pi ^ { C D | X } \sim Dir (\alpha_ { cd|x } ^ {C D | X }  ) ( x = 0 , 1 )$$

和一个Beta分布：

$$\pi ^ { X } \sim B e t a ( \alpha _ { 0 } ^ { X } , \alpha _ { 1 } ^ { X } )$$ 

​        &emsp;Dirichlet的非信息先验可以设置所有参数(参数)为1/2(遵循多项抽样的Jeffreys先验)，或为1(遵循均匀分布)。相反，如果我们想包含一些先验知识，我们可以将狄利克雷参数的通常解释作为“伪计数”添加到实际计数 $n_{\vec{y}}$ 中。例如，如果我们认为源Y几乎没有过度覆盖，我们会设置 $\alpha_{0|0}^{Y|X}$ 比 $\alpha _ {1|0}^{Y|X}$ 大.
​        &emsp;仍然需要在 $n_ {\vec{0}}$ 上设置先验分布，或者等价地，在N上设置先验分布。以下是文献中关于贝叶斯捕获-再捕获的常见选项:

•不恰当的平坦先验:P (N)∝1/N;
•泊松分布，最终与其pa参数上的超先验:N∼Poi (λ)， λ∼Gamma(α， β);
•里萨宁分布([29])，它总是恰当的，由p (N)∝2−log∗(N)给出，其中log∗(N)是数列{ $\log_2(N)，log_2(log_2(N))$ ，…}中正数项的和。

​        &emsp;我们进一步假设N和 $\theta$  的先验分布是相互独立的。

## 8.5.1 蒙洛卡特算法

​        &emsp;在本节中，我们详细介绍了基于Gibbs的蒙洛卡特算法从 $N_1$ 的后验分布进行抽样的步骤。让我们将模型的参数表示为 $\theta$ ，无论它们是否为 $\{Π^{C_i|S_i}\}$ 或 $\{Π^{C_i}\}$ 然后，在迭代(t+1)时，

1、狄利克雷分布从后验条件分布中对所有参数Θ(t+1)进行抽样；

2、

$$P ( N | \theta , T ) = P ( N | \theta , n _ { 0 , b s ) } = \frac { P ( N ) } { P ( n _ { 0  b s |\theta } ) } P ( n_ { obs } | N,\theta ) \propto  P(N)\dbinom{N}{n_{obs}} \pi _{\vec{0}}^{N-n_{0bs}}(1-\pi_{\vec{0}})^n_{obs}$$ 

​         &emsp;如果我们选择非正常先验 $P(N)\propto 1/N$ ,会得到 $N ^ { ( t + 1 ) } \sim \operatorname { NegBin } ( n _ { 0 b s } , 1 - \pi _ { \vec{0} } ^ { ( t + 1 ) } )$ 
3、
$$N _ { x , \vec{y} } ^ { ( t + 1 ) } \sim B i n ( n _ { \vec{y} } , \pi _{x|\vec{y}}^{(t+1)}\qquad n _ { \vec{0} }) = N ^ { ( t + 1 ) } - n _ { 0 b s }$$

 &emsp; 如果我们参考模型 $\{[ABX][CDX}$ 的例子，并参考参数化(8.5)，第一步包括以下三个步骤:

1、从后验的 $P ( \pi ^ { X } | N _ { x , y } )$ 得到的样本 $\pi _ { x } ^ { ( t + 1 ) }$ ，是Beta分布 $\operatorname { Beta } ( n _ { x } ^ { ( t ) } + \alpha _ { x } ) \qquad n _ { x } ^ { ( t ) } = \sum _ { \vec{y} } ^ { n } x _ { \vec{y} }$ ；

2、从 $P( \pi ^ { A B | X } | N _ { x , \vec{y} } )$ 中得到的样本 $\pi ^ { ( t + 1 ) }$ ，是狄利克雷分布 $Dir( n _ { a b x }^{(t)} + \alpha _ { a b| x } ) \qquad n_{abx} ^ { ( t ) } = \sum _ {\vec{y}:(A=a,B=b)}n_{x,\vec{y}}$ ；

3、从 $P ( \pi ^ { C D | X } | N _ { x , \vec{y} } )$ 中得到样本 $\pi _ { C d|x } ^ { ( t + 1 ) }$ ，是狄利克雷分布 $D i r ( n _ { c d x }^{(t)} + \alpha _ { c d |x } ) \qquad n _{cdx}^ { ( t ) } = \sum _ {\vec{y}:(C=c,D=d) } n _ { x ,\vec{y} }$ 并且，在第三步中，

$$
\pi _ { x|\vec{y} } = \frac { \pi _{x}^{( t + 1 )} \pi _{ab|x}^{( t + 1 )}\pi_{cd|x}^{(t+1)} } { \sum_x\pi_{x}^{(t+1)} \pi_{ab|x} ^ { ( t + 1 ) } \pi _{cd|x}^{( t + 1 )} }
$$

​         &emsp;如果对N采用泊松或里萨宁先验而不是1/N，则上述算法会稍微复杂一些。大都会-黑斯廷斯台阶上的第2步(大都会-黑斯廷斯-吉布斯内算法)，从 $\pi ( N | \theta , n _ { 0 b s } )$ 中采样一个值 $N ^ { ( t + 1 ) }$ 。 

## 8.5.2 仿真结果

  &emsp;&emsp;在本节中，我们报告模拟的结果，以经验地评估估计算法。

​         &emsp;我们考虑了两种情况:在第一种情况(见图8.2)中，用于生成数据的模型和估计模型在 $[ A X ] , [ B X ] , [ C D X ]$ 中具有相同的用法。在第二个场景中(见图8.3)，我们通过从模型 $[ A B X ] , [ C D X ]$ 生成数据并使用模型 $[ A X ] , [ B X ] , [ C D X ]$ 进行估计来测试模型对错误规范的鲁棒性。对于每个场景，我们生成了两个具有两种不同总体规模的数据集:一个是N = 500，一个是N = 1,000,000。在每个场景的生成模型中，我们设 $\pi _ { 0 } ^ { X } = 0.4$ ，即范围外单位(包括捕获和未捕获)的比例为40%，未观测单位(包括范围内和范围外)的比例为23%。在情景2中唯一的区别是A和B之间存在一种关系式的相互作用，特别是，在X = 1和X = 0下，A和B的相关性约为0.6，而在estimating配对模型中，它们是条件独立的。

![](D:\bigdata21\tupian.png)

**图8.2**

​         &emsp;情景1中 $N_1$ 的后验分布。N=500(左图)和N=1,000,000(右图)，非信息先验(上图)和信息先验(下图)的结果。实线表示的真实值 $N_1$ ,虚线为最大似然估计值 $N_1$ .

​         &emsp;为了评估先验分布对估计的影响，我们设置了两个案例。在第一种方法中，我们为所有参数(所有狄利克雷参数等于1且P (N) = 1/N)设置非信息先验。在第二种方法中，我们模拟来自审计样本的信息场景:我们从生成的完整总体 $[ X A B C D ]$ 中选取5%的样本，并将狄利克雷参数(parameter)设置为与该样本中观察到的计数相等。此外,为了模拟N的信息先验，我们设置泊松分布为情景2中 $N_1$ 的后验分布。N=500(左图)和N=1,000,000(右图)的结果，非信息先验(上图)，以及信息先验(底部图表)。实线表示 $N_1$ 的真实值(当N=1,000,000时不可见)，虚线表示最大似然估计( $N_1$ )。

![](D:\bigdata21\微信图片_20230329204824.png)

**图8.3**

​        &emsp;在图8.3中我们可以看到，当总体规模不是太大(N = 500)时，正确的信息先验可以有效地缓解估计中缺少参数所带来的误差模型。当总体规模较大(N=1,000,000)时，即使信息先验在正确的方向上影响了后验，它们的贡献似乎不足以弥补缺失的参数(至少在我们介绍的情况下是这样)。事实上 $N_1$ (600,000)的真实值甚至没有包含在99%最高后验密度区间中。