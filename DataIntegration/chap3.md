$\quad\quad$与Fellegi和Sunter [21]的方法不同，贝叶斯方法并不将不同记录对的匹配视为彼此不相关的，这是合乎逻辑的。在实践中，需要在初始链接步骤（[33]）之后施加一个线性约束，以避免在Fellegi和Sunter [21]的决策规则方法下可能存在的多个链接。贝叶斯公式更自然地考虑了通过矩阵C的约束。最后，根据贝叶斯的观点进行推理，该方法允许人们将链接过程的不确定性传播到后续的链接数据分析中，这将在第3.4节中所示。

$\quad\quad$贝叶斯方法的一个实际困难是缺乏对大型数据集的可扩展性，就像在类似人口普查的人口规模估计的情况下一样。从理论上讲，目标参数N的先验分布的选择似乎有些随意。在[59]的激励2列表例子中，我们有$\left(n_{1+}, n_{+1}\right)=(34,45)$，其中$n_11^* = 25&对有完全匹配的链接变量。na ̈ıve DSE是61，这可能是由于缺少匹配而高估的;而45将是N的下界，假设没有错误的枚举。由于没有信息，建议的先验不考虑这种考虑。结果N的后验中位数为55，97.5%分位数为65，后者似乎na ̈ıve DSE相当高。无论如何，尚不清楚结果对先验分布、信息性或非信息性有多敏感。

$\quad\quad$Steorts等人[55]提出了另一种贝叶斯方法，它允许同时链接和删除来自多个列表的记录。这个想法（类似于[44]）是将记录链接看作是一个识别单独文件背后的真正潜在"实体"的过程。因此，每个记录指一个目标人口单位（即潜在实体），由一个整数从1到Nmax表示，称为链接结构，其中Nmax是所有列表的记录总数，所以在极端情况下列表可以没有共同的实体。所有引用同一实体的记录都是"链接的"。假设连杆结构的先验分布是均匀。建模方法在其他方面与上述概述的方法类似。采用混合MCMC算法生成连杆结构的后验分布。引入了最可能最大匹配集的规则，保证了多个列表匹配的传递性。详情请参阅[57]。

$\quad\quad$对于基于捕获-再捕获数据的种群大小估计，[57]的方法原则上允许根据它们的后验分布采样细胞计数，如表3.1和3.2中的样本。例如，附加到一个记录中的潜在实体只在该列表中捕获，附加到两个记录中的潜在实体在这两个列表中捕获，以此类推。如果每个记录都附加到至少一个潜在实体，细胞计数集可以形成捕获-再捕获数据，以拟合对数线性模型，从而得出种群大小估计。

$\quad\quad$[60]在相同的贝叶斯设置中扩展了[57]，通过将种群大小作为模型参数，允许在多列表记录链接和重复数据删除问题中进行种群大小估计。
他们观察到，在统一之前的单位（即独立随机抽样替换人口N），样本标签的分布N诱发分布分区空间（即不同的潜在个体的数量）取决于N，从而估计标签的分区将允许同时产生推理N和估计连锁结构。该模型具有完整的先验性$$P{N}=\frac{1}{Z{g}N^g} \qquad （1<N<+∞)$$
其中Z （g）是黎曼泽塔函数。同样，连锁结构的先验（均匀）分布的选择并不适应捕获-再捕获概率的不同情况。如果能够建立并合并对数线性模型的条件单元概率的信息先验分布，它可能会更有用。
## 3.4存在链接错误时的DSE
### 3.4.1Ding和Fienderg估计器
[32]研究了连锁误差对捕获-再捕获估计的影响。在存在连锁错误的情况下，覆盖率，如（3.2）和（3.3），以及人口规模的估计，如（3.1），可能会有偏差，需要进行调整。在2个列表的情况下，[19]提出了一种修正彼得森估计器的方法。[32]$n=n_11^* + n_10^* + n_01^*$根据表3.1给出，但基于连锁数据，其中$n_{1+}^* = n_{1+}$和$n_{+1}^* = n{+1}$。使用N∗11而不是n11的面值DSE通常是有偏倚的。给定n的（n∗11，n∗10，n∗01）的条件似然为$$L{\pi_11^*，\pi_10^*，\pi_01^*}= \frac{n！}{{n_11^*}！{n_10^*}！{n_01^*}！}\frac{{\pi_11^*}^{n_11^*}{\pi_10^*}^{n_10^*}{\pi_01^*}
{n_01^*}}{{{\pi_11^*}+{\pi_10^*}+{\pi_01^*}}^n}$$ 其中，$\pi_10^* = \pi_1^* - \pi_11^*$和$\pi_01^* = \pi_2^* - \pi_11^* = \pi_2 - \pi_11^*$
为了将这些参数与连锁误差联系起来，Ding和finberg（1994）做出了以下假设：

1.有一个假设的链接方向，其中L1链接到L2，即，对于L1中的每一个记录，如果可能的话，我们可以在L2中找到一个链接，而不是相反;

2.L1和L2之间的真实链接的概率为α;

3.涉及L1中匹配记录的错误链接可以忽略不计，涉及L1中不匹配记录的错误链接发生的概率为β。
 

$\quad\quad$请注意，除了链接相互依赖关系外，真实匹配率与（3.8）相同，但错误链接率与（3.7）不相同。在DSE的这些附加假设下，我们有$\pi_{11}^{*}=\alpha \pi_{1} \pi_{2}+\beta \pi_{1}\left(1-\pi_{2}\right)$其中这两项分别来自n11和n10。最大化关于π1和π2的条件似然（3.10），对于给定的β和α值，第一个列表的估计覆盖范围由$$\hat{\pi}_{1, D F}=\frac{-n_{11}^{*}+\beta\left(n_{11}^{*}+n_{10}^{*}\right)}{(\beta-\alpha)\left(n_{11}^{*}+n_{01}^{*}\right)}=\left(\frac{1}{n_{+1}}\right)\left(\frac{n_{11}^{*}-\beta n_{1+}}{\alpha-\beta}\right)$$
第二个列表是由$$\hat{\pi}_{2, D F}=\frac{-n_{11}^{*}+\beta\left(n_{11}^{*}+n_{10}^{*}\right)}{(\beta-\alpha)\left(n_{11}^{*}+n_{10}^{*}\right)}=\left(\frac{1}{n_{1+}}\right)\left(\frac{n_{11}^{*}-\beta n_{1+}}{\alpha-\beta}\right)$$

$\quad\quad$将（3.11）和（3.12）与（3.2）和（3.3）进行比较，可以确定（3.11）和（3.12）的通用项为真实匹配数的链接误差调整估计。由此
可见，N的条件MLE，也是调整后的彼得森估计量，由$$\tilde{N}_{D F}=(\alpha-\beta) \frac{n_{1+} n_{+1}}{n_{11}^{*}-\beta n_{1+}}=\frac{(\alpha-\beta) n_{11}^{*}}{n_{11}^{*}-\beta n_{1+}} N_{p}^{*}$$
其中，$\hat{N_p^*}$是直接使用n 的面值DSE;参见[19]
### 3.4.2修改后的Ding和Fienberg估计器
Di Consiglio和Tuoto[16]提出了一种推广，放宽了单向链接的限制。管理来源的链接激发了这种需求。请注意，丁和芬伯格。
[19]处理的是传统的人口普查覆盖率不足评估，即人口普查计数和计数后调查之间的联系是在一个方向上工作的。当不是单向链接时，L1和L2中不匹配的记录可能会出现错误链接，因此由于后者，我们使用了一个附加项的$\pi_{11}^{*}=\alpha \pi_{1} \pi_{2}+\beta \pi_{1}\left(1-\pi_{2}\right)+\beta \pi_{2}\left(1-\pi_{1}\right) \\$

$\quad\quad$保留上述关于连锁误差的其他假设，最大化条件似然（3.10），对于给定的β和α值，修正的覆盖率的Ding和feng（MDF）估计分别由，$$\hat{\pi}_{1, M D F}=\frac{2 \beta n_{11}^{*}+\beta x_{10}^{*}+\beta n_{01}^{*}-n_{11}^{*}}{(2 \beta-\alpha)\left(n_{11}^{*}+n_{01}^{*}\right)}=\left(\frac{1}{n_{+1}}\right)\left(\frac{n_{11}^{*}-\beta\left(n_{1+}+n_{+1}\right)}{\alpha-2 \beta}\right)$$ $$\hat{\pi}_{2, M D F}=\frac{2 \beta n_{11}^{*}+\beta n_{10}^{*}+\beta n_{01}^{*}-n_{11}^{*}}{(2 \beta-\alpha)\left(n_{11}^{*}+n_{10}^{*}\right)}=\left(\frac{1}{n_{1+}}\right)\left(\frac{n_{11}^{*}-\beta\left(n_{1+}+n_{+1}\right)}{\alpha-2 \beta}\right)$$同样，（3.14）和（3.15）的共同项是n11的连锁误差调整估计，因此MLE，也是n的调整彼得森估计，由$$\tilde{N}_{M D F}=(\alpha-2 \beta) \frac{n_{1+} n_{+1}}{n_{11}^{*}-\beta\left(n_{1+}+n_{+1}\right)}=\frac{(\alpha-2 \beta) n_{11}^{*}}{n_{11}^{*}-\beta\left(n_{1+}+n_{+1}\right)} \tilde{N}_{p}^{*}$$

$\quad\quad$如上所述，DF和MDF估计器都是基于链接误差总体上是恒定的假设。如果这一假设至少在子群中成立，那么估计量可以应用于连锁误差概率（和捕获概率）比在整个种群中更均匀的地层中。

$\quad\quad$当假定错误率已知时，可以考虑真实值n11，n10，n01是由$n_{11}^*$，$n_{10}^*$，$n_{01}^*$的代数确定地得到的。那么，MDF的方差估计量与标准的DSE方差估计量（见[63]）相同，它是由$$\hat{V}\left(\tilde{N}_{M D F}\right)=N \frac{\left(1-\hat{\pi}_{1}\right)\left(1-\hat{\pi}_{2}\right)}{\hat{\pi}_{1} \hat{\pi}_{2}}$$其中，$\hat{\pi}_{1}=\hat{n}_{11} / n_{+1}$和$\hat{\pi}_{2}=\hat{n}_{11} / n_{1+}$在[58]中也采用了同样的方法来校正捕获-重捕获双样本丰度估计器，以解释识别中的假阴性误差。此外，他们还引入了一种自举方法来估计修正估计量的方差。自助法允许合并连锁错误率的估计不确定性。   
### 3.4.3一些评论
下面我们在存在连锁误差的情况下对DSE进行一些评论，并在没有连锁误差的情况下对其与DSE的偏差和方差进行比较。

$\quad\quad$首先，在DSE估计器的标准动机下，两个列表都被认为是随机的，独立枚举，每个都具有恒定的覆盖率，分别用τ1和τ2表示。我们已经$$\left\{\begin{array}{l}
E\left(n_{1+}\right)=N \tau_{1} \\
E\left(n_{+1}\right)=N \tau_{2} \\
E\left(n_{11}\right)=N \tau_{1} \tau_{2}
\end{array} \quad \Rightarrow \quad \frac{E\left(n_{1+}\right) E\left(n_{+1}\right)}{E\left(n_{11}\right)}=N \Rightarrow \hat{N}=\frac{n_{1+} n_{+1}}{n_{11}}\right.$$ 正如[65]和[66]所示的那样，一个更简单的动机是对待其中一个列表，比如列表1是固定的，它需要更少的假设，更适合来自管理来源的列表。假设随机n+1，与在整个目标种群中的恒定捕获率，我们以n1+为条件，$$\left\{\begin{array}{c}
E\left(n_{1+} \mid n_{1+}\right)=n_{1+} \\
E\left(n_{+1} \mid n_{1+}\right)=E\left(n_{+1}\right)=N \tau \\
E\left(n_{11} \mid n_{1+}\right)=n_{1+} \tau
\end{array} \quad \Rightarrow \quad \frac{n_{1+} E\left(n_{+1}\right)}{E\left(n_{11} \mid n_{1+}\right)}=N \quad \Rightarrow \hat{N}=\frac{n_{1+} n_{+1}}{n_{11}}\right.$$我们将在下文中采用这种简化的条件观点。

$\quad\quad$设α为列表1和列表2之间的匹配包含在链接的关节子集中的概率。设β为列表1或列表2中不匹配的记录包含在链接的关节子集中的概率。我们有$$E\left(n_{11}^{*} \mid n_{1+}\right)=\alpha E\left(n_{11} \mid n_{1+}\right)+\left[n_{1+}-E\left(n_{11} \mid n_{1+}\right)\right] \beta+\left[E\left(n_{+1}\right)-E\left(n_{11} \mid n_{1+}\right)\right] \beta \\
=(\alpha-2 \beta) E\left(n_{11} \mid n_{1+}\right)+\left[n_{1+}+E\left(n_{+1}\right)\right] \beta$$给定（α，β），我们可以得到以下估计量$$\hat{n}_{11}=\hat{n}_{11}(\alpha, \beta)=\frac{n_{11}^{*}-\beta\left(n_{1+}+n_{+1}\right)}{(\alpha-2 \beta)}$$ $$\hat{N}(\alpha, \beta)=\frac{n_{1+} n_{+1}}{\hat{n}_{11}(\alpha, \beta)}$$这与（3.16）中的MDF估计器相同。β= 0的特殊情况是特别有趣的，例如，在动物丰度估计下，我们有$$\hat{N}(\alpha)=\hat{N}(\alpha, 0)=\frac{\alpha n_{1+}{ }_{+1}}{n_{11}^{*}}$$

$\quad\quad$下面我们对$n_{+1}$的调整$\hat{N}$DSE(α)和标准DSE进行了明确的比较。这是应用DSE的通常情况，前提是观察并考虑每个单一源的计数。有和没有连锁误差之间的差异是观察$n_{11}^*$或$n_{11}$之间。对于标准的DSE$\hat{N}$，$$\hat{N}=\frac{n_{1+} n_{+1}}{n_{11}} \approx \frac{n_{1+} n_{+1}}{E\left(n_{11} \mid n_{1+}, n_{+1}\right)}-\frac{n_{1+} n_{+1}}{E\left(n_{11} \mid n_{1+}, n_{+1}\right)^{2}}\left[n_{11}-E\left(n_{11} \mid n_{1+}, n_{+1}\right)\right] \\
+\frac{1}{2} \frac{2 n_{1+} n_{+1}}{E\left(n_{11} \mid n_{1+}, n_{+1}\right)^{3}}\left[n_{11}-E\left(n_{11} \mid n_{1+}, n_{+1}\right)\right]^{2}$$ $$E\left(\hat{N} \mid n_{1+}, n_{+1}\right) \approx \frac{n_{1+} n_{+1}}{E\left(n_{11} \mid n_{1+}, n_{+1}\right)}-\frac{n_{1+} n_{+1}}{E\left(n_{11} \mid n_{1+}, n_{+1}\right)^{3}} V\left(n_{11} \mid n_{1+}, n_{+1}\right)$$ $$V\left(\hat{N} \mid n_{1+,} n_{+1}\right) \approx \frac{\left(n_{1+} n_{+1}\right)^{2}}{E\left(n_{11} \mid n_{1+,} n_{+1}\right)^{4}} V\left(n_{11} \mid n_{1+,} n_{+1}\right)$$

$\quad\quad$我们有n11个∼Bin（$n_{+1}$，$n_{1+}$/N）条件为（$n_{+1}$，$n_{1+}$），因为在列表2中枚举列表1内或列表外的记录的概率是相同的。因此，$E\left(n_{11} \mid n_{1+} n_{+1}\right)=n_{1+} n_{+1} / N$和$V\left(n_{11} \mid n_{1+}, n_{+1}\right)=n_{+1}\left(\frac{n_{1+}}{N}\right)\left(1-\frac{n_{1+}}{N}\right)$，从而$$E\left(\hat{N} \mid n_{1+}, n_{+1}\right) \approx N+\frac{N}{\hat{N}^{2}} n_{+1} \frac{n_{1+}}{N}\left(1-\frac{n_{1+}}{N}\right)=N+O(1)$$ $$V\left(\hat{N} \mid n_{1+}, n_{+1}\right) \approx \frac{N^{2}}{\hat{N}^{2}} n_{+1} \frac{n_{1+}}{N}\left(1-\frac{n_{1+}}{N}\right)$$其中，我们假设$\frac{n_{1+}}{N}=\frac{n_{+1}}{N}=\frac{N}{N}=O(1)$为N→∞。接下来，对于(α)，我们也有$$E\left(\hat{N}(\alpha) \mid n_{1+}, n_{+1}\right) \approx \frac{\alpha n_{1+} n_{+1}}{E\left(n_{11}^{*} \mid n_{1+} n_{+1}\right)}-\frac{\alpha n_{1+} n_{+1}}{E\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)^{3}} V\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)$$ $$V\left(\hat{N}(\alpha) \mid n_{1+}, n_{+1}\right) \approx \frac{\left(\alpha n_{1+} n_{+1}\right)^{2}}{E\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)^{4}} V\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)$$由于$n_{11}^*$∼Bin（$n_{11}$，α）条件是$n_{11}$和$n_{11}$∼Bin（$n_{+1}$，$n_{1+}$/N）条件（$n_{+1}$，$n_{1+}$），我们有$n_{11}^*$∼Bin（$n_{+1}$，α$n_{1+}$/N）条件（$n_{+1}$，$n_{1+}$）。因此，$E\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)=\alpha n_{1+} n_{+1} / N$和$V\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)=n_{+1}\left(\alpha n_{1+} / N\right)\left(1-\alpha n_{1+} / N\right)$，这样$$E\left(\hat{N}(\alpha) \mid n_{1+} n_{+1}\right) \approx N+\frac{N}{\hat{N}(\alpha)^{2}} n_{+1} \frac{\alpha n_{1+}}{N}\left(1-\frac{\alpha n_{1+}}{N}\right)=N+O(1)$$ $$V\left(\hat{N}(\alpha) \mid n_{1+}, n_{+1}\right) \approx \frac{N^{2}}{\hat{N}(\alpha)^{2}} n_{+1} \frac{\alpha n_{1+}}{N}\left(1-\frac{\alpha n_{1+}}{N}\right)$$

$\quad\quad$由此可见，$\hat{N}$(α)的条件偏差和方差，即存在连锁误差的方差，与$\hat{N}$的条件偏差，即没有连锁误差的顺序相同。然而，如果$n_{1+}$/N和α都接近于1，我们就有了$$V\left(n_{11}^{*} \mid n_{1+}, n_{+1}\right)>V\left(n_{11} \mid n_{1+}, n_{+1}\right)$$也就是说，当$\hat{N}$（1−α）$\hat{n}_{11}$=（$n_{1+}$+$n_{+1}$−2$\hat{n}_{11}$）β时，我们有$\hat{n}_{11}=$n_{11}^*$，也就是说，如果两种类型的链接错误发生“相互抵消”。为了衡量$\hat{n}_{11}$的期望和方差，让$$n_{11}^{*}=n_{11}^{*}(1,1)+n_{11}^{*}(0,1)+n_{11}^{*}(1,0)$$依次来自列表1和列表2之间的匹配，列表2中不匹配的记录，以及列表1中列表2中不匹配的记录。有条件的（$n_{+1}$，$n_{1+}$），我们有$$n_{11}^{*}(1,1) \sim \operatorname{Bin}\left(n_{+1}, \alpha n_{1+} / N\right)$$ $$n_{11}^{*}(0,1) \sim \operatorname{Bin}\left(n_{+1}, \beta\left(1-n_{1+} / N\right)\right)$$
同样的$$n_{11}^{*}(1,0) \sim \operatorname{Bin}\left(n_{1+},\left(1-n_{+1} / N\right) \beta\right)$$ 由此得出结论$$E\left(n_{11} \mid n_{1+}, n_{+1}\right)=(\alpha-2 \beta) \frac{n_{1+} n_{+1}}{N}+\left(n_{1+}+n_{+1}\right) \beta$$ $$E\left(\hat{n}_{11} \mid n_{1+}, n_{+1}\right)=\frac{n_{1+} n_{+1}}{N}$$ $$E\left(\hat{N}(\alpha, \beta) \mid n_{1+}, n_{+1}\right)=N+\frac{N}{\hat{N}(\alpha, \beta)^{2}} V\left(\hat{n}_{11} \mid n_{1+}, n_{+1}\right)$$

$\quad\quad$观察$n_{11}^*$（1、1）、$n_{11}^*$（0、1）和$n_{11}^*$（1、0）不相互独立，例如，较小n11偶然使较小$n_{11}^*$（1、1）和较大$n_{11}^*$（0、1）和$n_{11}^*$（1、0）同时更可能。

$\quad\quad$虽然V（$\hat{n}_{11}$|$n_{1+}$，$n_{+1}$）在分析上似乎难以处理，但可以通过蒙特卡罗模拟来计算它。我们推测\mid V\left(\hat{n}_{11}^{1} \mid n_{1+}, n_{+1}\right) / N=O(1)渐近为N→∞，，在这种情况下，$\hat{N}$（α，β）的偏差和方差与标准DSE$\hat{N}$.的偏差和方差的顺序相同。
### 3.4.4例子
下面我们将给出一些存在连锁错误的两个例子。

$\quad\quad$例子1。考虑[23]表3.3和3.4中的数据，有三个来源：1990年美国人口普查，相应的后普查调查（PES），行政清单补充，圣路易斯的PES抽样层（ALS）。表3.3给出了Census-PES捕获-再捕获数据，而匹配误差研究（[42]）的结果见表3.4。  
   表3.3  
1990年美国人口普查，第11层，圣路易斯，美国人口普查<h4></h4>   
<table border="1" width="500px" cellspacing="10">
<tr>
<td colspan="2" align="center"></td>
  <td colspan="2" align="center">人口普查</td>
</tr>
<tr>
  <td></td>
  <td></td>
  <td>目前</td>
  <td>缺失</td>
</tr>
<tr>
  <td> PES</td>
  <td>目前</td>
  <td>487</td>
  <td>129</td>
</tr>
<tr>
  <td></td>
</td>
  <td>缺失</td>
  <td>217</td>
  <td>-</td>
</tr>
</tr>
</table>
