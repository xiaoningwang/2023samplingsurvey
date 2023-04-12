##以下是邹欣写的
### 2.1简介

​       在记录联系中（Fellegi和Sunter，1969年;Herzog等人，2007年），也称为实体或共同参考决议（Christen，2012年），其目的是将记录（与相关观察结果）识别并合并到单独的文件中，它们对应于相同的实体（或单元）。除非存在允许精确匹配的唯一标识符，当链接的记录实际上并未引用同一实体，或者未能链接引用同一实体的记录时，链接错误是不可避免的。
$$
表2.1
\space \space 2011年英国人口普查汇总和患者等级链接\\
\begin{array}
{cll}
\hline
	类型  & 通过  & 链接数 & 错误链接率 \\
\hline
	确定性的 &1&30780660   & 0.00011 \\
	        &2&11733197   &0.00389 \\
	        &3&1513471    & 0.00561 \\
	        &4&2444838    & 0.00375\\
	        &5&1346432    &0.00748\\
	        &6&121483     &0.00886\\
	        &7&1007293    &0.00100\\
	        &8&825069     &0.01485\\
	        &9&35432      &0.00100\\
\hline
概率性的&1&511239&0.02948\\
       &2&298645&0.07165\\
\hline
全部的& & &50617759\\
\hline
\end{array}
$$
​       表2.1提供了国家统计局（ONS）开展的2011年英格兰人口普查和患者登记（C-PR)联系的汇总。详见Owen等。（2015）关于基本方法的一些相关方面。是一个过程中的每一个都使用匹配（或关键）变量的不同组合，包括姓名、出生日期、性别、邮政编码等。最后，大约300万人口普查记录没有与患者登记册链接，而所有链接的记录都不是真是匹配的，如在不同通过处的估计错误链接率（FLR）所指示的。

​      二级分析员如何根据C-PR关联数据进行有效分析？这是我们在这里考虑的问题。有两个条件有助于澄清二级分析师的立场。一方面，可以有把握地假设二级分析师不能直接访问所有关键变量和单独的数据文件，也没有复制实际链接程序的详细知识或工具。例如，姓名、出生日期和邮政编码都用于C-PR链接，但出于隐私保护，这些数据不太可能被发布给二级分析师。另一方面，还必须假设可以提供一些关于记录如何实际相互比较的非公开链接比较数据。这对于允许任何其他可能性而不是按原样分析链接数据集是必要的；例如，参见Gilbert等人（2017）关于用于数据链接的信息交换协议。例如，可以为一对链接记录分配一个分数，以反映它们在关键变量方面的一致性。

#### 2.1.1 相关工作

​       当涉及到由于链接误差传播不确定性的记录链接数据分析时，可以在文献中辨别出不同的方法。不幸的是，现有的方法都不能充分解决上述二次分析的问题。

​        几乎所有的频繁方法都是基于一个数据集中的记录的概率的链接模型。Scheuren和Winkler(1993,1997)、Lahir和larsen(2005)、Chambers(2009)、Chipperfield和Chambers(2015)使用Monte Carlo模拟两个数据集之间的一致模式和后续记录链接过程来估计链接模型。该方法仅限于能够访问所有关键变量和实际链接过程的人员。因此，二级分析师将不得不采用大大简化的假设，如分庭的可交换链接错误模型（2009年），或采用特别办法。例如，Hof和Zwinderman（2015）假设要链接的两个数据集中的每对记录独立贡献成对为伪似然。但目前还不清楚这种伪概率如何产生有效的推论。

​        对于贝叶斯框架下的记录链接和分析，可以参考Tancredi和Liseo（2011），Gutman等人（2013）、Sadinle（2014）和 Stoerts等人（2016）等。推理基于匹配数据集的后验分布。观察到的关键变量的建模方法，如Stoerts等人所述。（2016），建立在Copas和Hilton（1990）的撞中模式。然而，失真模型适用于次要分析人员无法访问的关键变量，除非在特殊情况下。将大量后验匹配数据集分发给二级分析师可能并不总是切实可行的，也不能像C-PR案例中那样扩展到大型数据集。

​      Goldstein等人(2012)应用多重插补，其中被接受为正确链接的单元具有完整的观察结果，而其他记录对应的单元被认为具有缺失的观察结果。Gutman等人（2015）采用不同的方法对链接的记录对进行乘法插补，无论它们是否为真匹配，其中假定未链接的真匹配随机丢失。应该注意的是，这些估算方法不受底层链接数据结构的限制。例如，链接的记录是否为真匹配将影响其他记录对是否可能为真匹配。这些插补方法中没有内置类似这样的限制。

####  2.1.2勘察情况

​      我们将在续集中研究三种不同的方法。我们在第2.3节考虑了感兴趣参数的最大似然估计（MLE）。这基于这样一个概念，即比较数据是根据概率模型产生的，包括可能导致链接误差的关键变量真实值的失真。然而，它将会被证明，基于"链接数据结构"（详见第2.2节）模型的MLE一般不简单。讨论涉及到从二元正态变量的破碎随机样本对相关系数的MLE分析（DeGroot和Goel,1980).

​      接下来，使用线性回归作为案例，我们在第2.4节中考虑无偏回归估计，条件是可用的比较数据，并将其与现有的频繁方法（例如，Lahiri和Larsen, 2005）。本文概述了一种适用于常见的概率性和确定性记录链接过程的比较数据建模方法，并讨论了C-PR案例中的可扩展性问题。我们发现，对于这种重量级的数据集，旨在使用所有记录的有效分析是极其不切实际的。

​      第三种方法包括仅利用具有几乎零 FLR的链路，以及基于所产生的链路子集针对剩余的错误链路和缺失的真匹配来调整所述分析。虽然这显然会导致效率的损失，但在诸如C-PR链接，其中仍然可以获得大的精确数据集，例如在估计的FLR大约为一万分之一的情况下，从第一次通过获得3000万条记录。当所有链路都是真匹配时，调整问题就变成了数据丢失问题。虽然对缺失数据的分析本身就是一个具有挑战性的问题，但使用可扩展到大型数据集的重权重和估算方法的传统由来已久。我们在第2.5节中说明了无偏链接子集回归的条件是如何发展的，然后是对面值回归分析进行简单的渐进调整，与允许在存在缺失数据的情况下基于观察数据进行分析的假设类似。

### 2.2联动数据结构

​       用于记录链接的独立数据集的链接数据结构将该问题与数据丢失或测量错误的问题区分开来。在后一种情况下，分析单位周围没有歧义，但是与单位相关的值可以观察到有误差或没有误差完全被观察到了。通常，可以合理地假设随机测量误差或缺失指标在一个单元与另一个单元之间是独立的，或者在一些已识别的单元簇之间是独立的。然而，在记录链接中，记录是否链接到另一条记录，或链接到哪条记录，必要地影响其他记录的链接方式，但注意力被限制在每个可链接记录集块内时也是如此。记录的不同配对会导致结果分析的不确定性，即使所有相关结果值都被视为固定值。一个关键问题是如何适当地适应这种相互依赖单元模糊性，而不是控制与识别单元相关的输出变量的随机机制。由于链接数据结构，没有明确说明这种相互依赖性的任何建模方法可能难以从一般意义上证明其合理性。此外，与实现相关的一个问题是如何处理两个或多个记录集之间的这种相互依赖所必需的组合学。

#### 2.2.1 定义

现在，我们提供两个独立数据集基础上的链接是数据结构正式描述。该定义可以很容易地推广到涉及三个或更多数据集的情况。分别用大小为
$$
n_A和n_B
$$
的A和B表示单独的数据集。用
$$
a\in A
$$
表示A中的任意记录，类似于
$$
b\in B
$$


如果是任意固定整数，则表示
$$
a\in A
$$


所指的真实实体。

如果
$$
i_a = i_b
$$
则对于
$$
a \in A
$$
和
$$
b \in B
$$
则记录对（ab)被称为匹配，在这种情况下，a和b都被称为A和B之间的匹配记录。如果对于所有的
$$
b\in B有i_a\neq i_b
$$
则a是不匹配的记录。类似地，对于B中的b。我们应该参考
$$
{\Lambda}_{AB} = \{i_a;a\in A\}\cup\{i_b;b\in B\}
$$
作为A和B之间的链接数据结构，除非存在允许精确匹配的唯一标识符，否则链接数据结构是未知的。A和B之间记录链接的结构可以看作是一个估计的链接数据结构，表示为
$$
{\widehat\Lambda}_{AB}=\{\widehat i_a;a\in A\}\cup\{\widehat i_b;b\in B\}
$$
其中
$$
\widehat i_a = \widehat i_b
$$
提供（ab)是一对链接记录。

​      当只有两个数据集A和B时，链接数据结构可以由
$$
n_A \times n_B
$$
匹配矩阵
$$
\omega
$$
表示

其中
$$
{\omega}_{ij}=1
$$
如果在A中的第i个记录和B中的第j个记录是一对匹配，否则

***因为下文中有很多字母和字符，所以从下文开始中文和字母字符都在$$里面写，除标题外其他默认居中对齐***


$$
{\omega}_{ij}=
$$

$$
我们假设已删除重复记录，因为{\sum}_j{\omega}_{ij}\leq 1 和{\sum}_i{\omega}_{ij}\leq 1，用\Omega表示匹配矩阵空间，这个空间包含了所有可能的\omega。由于每个不同的\omega\\对应于不同的链路数据结构{\Lambda}_{AB},如果i_a和i_b的值不相关，则空间\Omega有效地定义了链接数据结构的空间。
$$

$$
设\omega 是任意记录对（ab）的匹配指示符。设M(\omega)=\{(ab);{\omega}_{ab}=1,a\in A，b\in B\}是A和B之间的\\匹配集合。设\cup_A(\omega)是A中不匹配记录的集合，\cup_B(\omega)是B中不匹配记录的集合。在没有重复记录的情况下，对应于\omega的子矩阵\\是置换矩阵,其可以由单位矩阵的列（或行）的置换给出。这是匹配矩阵表示的一个有用属性。设 \widehat\omega=\omega_L是链接矩阵，\\其中如果(ab)是一对链接记录，则{\omega}_{L,ab}=1,否则\omega = 0。设A_L由来自A的链接记录组成，B_L由来自B的链接记录组成。\\记录链接可视为M，\cup_A和\cup_B的产量估计，其中\widehat M=M(\omega_L)=\{(ab);\omega_{L,ab}=1,a\in A_L,b\in B_L\},\\并且各自有\widehat \cup_A=A或者A_L和\widehat \cup_B = B或者B_L。\\
$$

$$
我们通常观察到的链接和匹配之间的术语上的明显区别：匹配（ab)意味着\omega_{ab}=1并且事实上i_a=i_b\\一个链接(ab)意味着\widehat \omega_{ab}=\omega_{L,ab}=1,可能匹配也可能不匹配。
$$

$$
用\delta_{ab}表示a\in A和b\in B 之间单独比较数据，供二级分析员使用。设G为可用\delta_{ab}的集合。为了区分，\\让\delta_{a}^*和\delta_{b}^*是执行记录链接的数据提供者可访问的完全比较数据。在C-PR链接中，\delta_a^*包含\\年龄、性别、姓名、出生日期、邮政编码等。不是所有的这些都会被释放。事实上，即使其中一些因素，\\例如年龄和性别，也很可能不会被指定为关键变量。
$$

$$
作为可用比较数据的可能示例，\delta_{ab}可以是总结a如何与b进行绝对或相对比较的值。特别的，如果a和b\\被认为彼此一致，则一致指示符由\delta_{ab}=1给出，否则由\delta_ab=0给出。链路指示器，其中\delta_{ab}=1,当且仅当\omega_{L,ab}=1,可视为特殊协议指示符，\\不仅涉及A和B之间的比较，还涉及是否将它们联系起来的决定。例如，从表2.1的确定性通过，如果a和b在所有相关关键变量方面\\彼此相等，则可以设置一致指示符\delta_{ab}=1。其中\delta_{ab}=1和\delta_{ab'}=1也许会出现，尽管b\neq b'\in B。然而，如果\delta_{ab}是来自相同通道的指示符，\\则通常不允许b=b'的\delta_{ab}=\delta_{ab'}=1,因此\delta_{ab}=1现在表示唯一一致。
$$

$$
作为最小可用比较数据的场景，让辅助分析仅给定链接数据集，但是没有任何无链接记录，G将会是\\链接矩阵\omega_L的数量，而不是其他。
$$

#### 2.2.2 匹配空间的协议划分

$$
随着数据集大小的增加，匹配空间\Omega的完整枚举很快失去控制。为了便于有效的部分枚举，可以引入\Omega 的合理划分，\\并且在某些假设的统计模型下仅处理其中\omega 最可能的子空间。下面我们描述\Omega的一致分割，并考虑它的一些含义。
$$

$$
在不失一般性的情况下，设D为由第二分析员基于可用比较数据G定义的一致矩阵，其中D_{ab}=1,\\
如果记录a\in A与b\in B'一致',否则D=0。例如，我们有D_{n\times n}=\omega_L提供的最小信息，\\其中n=|A_L|=|B_L|。或者，如果另外提供了A和B未链接的记录，则可以使\\D_{n\times n}=\omega_L。如果Fellegi-Sunter(FS)范式（Fellegi和Sunter，1969)下\\提供来自概率记录链接的权重\lambda_{ab}，则可以将一致性定义为如果\lambda_{ab}超过某个阈值的情况。
$$

$$
设\kappa = {\sum}_i{\sum}_jD_{ij}是D中的总协议数。对任意的\omega \in \Omega,设d = {\sum}_a{\sum}_bD_{ab}\omega_{ab}是匹配协议数，\\其中(\omega_{ab},D_{ab})=(1,1),即，根据\omega匹配且根据D彼此一致的记录的编号对。例如，假设D_{n\times n}和一一匹配\\,因此，每个\omega必然是一个置换矩阵，我们还需要附加的(\omega_{ab},D_{ab})=(1,0)的n-d记录对，\\(\omega_{ab},D_{ab})=(0,1)的k-d记录对，以及(\omega_{ab},D_{ab})=(0,0)的n^2-n-(k-d)记录对。
$$

$$
设\Omega_0是达到最大值的d的集合，用d_0表示。\Omega_0中的矩阵被称为最一致匹配矩阵（MAMM_s)。\\例如，在一一链接和一一匹配下，其中D也是置换矩阵，并且\kappa=n,只有一个d_0=n的MAMM，\\即，\omega=D。更一般地，在一一链接下，D可以是一个具有\kappa\geq n。只要D允许，一个一一链接的\\MAMM可以有(\omega_{ab},D_{ab})=(1,1)中的n个元素,(\omega_{ab},D_{ab})=(1,0)中的0个元素，\\(\omega_{ab},D_{ab})=(0,1)中的\kappa-n个元素，以及(\omega_{ab},D_{ab})=(0,0)中的n^2-\kappa个元素。
$$

$$
只要这样的矩阵存在，设\Omega_1包含被d_0-1给出的所有带有协同协定的\omega s。这不一定总是可能。\\例如，在一一链接和意义匹配的情况下，D和a之间的下一个最高可能的一致数是n-2,而不是n-1，\\因此，我们有\Omega_1 = \phi。设\Omega_k是d_0-k给出的符合数的匹配矩阵集。\Omega的协议分割由下式给出：
$$

$$
\Omega(D)=\cup_{k=0}^{d_0}\Omega_k
$$

$$
如果对记录链接程序由一定的信心，那么实践中的任何建模方法都会使得\omega的可能性在\Omega的一致性\\分割上随k而减小，这似乎是很自然的。事实上，在最小比较数据的情况下，其中D=G=\omega_L,\\二级分析员将不会有信息对来自相同\Omega_k的\omega放置不同的可能性。此外，设\omega\in \Omega_k\\和\omega'\in \Omega_{k+1},使得它们仅在一对(ab)中不同，其中D_{ab}=1,\omega_{ab}=1，\\但是\omega_{ab}’=0。两者之间基于D的似然比由下式给出：
$$

$$
LR(\omega',\omega;D)=\frac{L(\omega';D)}{L(\omega;D)} = \frac {Pr(D|D_ab=1;\omega')} {Pr(D|D_ab=1;\omega)}·\frac{Pr(D_{ab}=1;\omega')}{Pr(D_{ab}=1;\omega)}
$$

$$
既然现在\omega和\omega'在除(ab)之外的任何地方都是一致的，那么最简单的建模方法可以假设
$$

$$
LR(\omega',\omega;D)=\frac{Pr(D_{ab}=1;\omega_{ab}'=0)}{Pr(D_{ab}=1;\omega_{ab}=1)}
$$

显然，作为二次分析的最低可行性要求，必须提供某些链接准数据，已便能够对该比率建模，从而建立尊重链接数据结构的方法。

## 以下由陆星欢进行编写
## 2.3 最大似然估计

&emsp; &emsp;  $DeGroot$ 和 $Goel(1980)$ 考虑了一个基于一个破碎样本的双变量正态分布的相关系数的 $MLE$ 。换句话说，从分布中抽取 $n$ 对随机样本，但观察到的数据只是 $n$ 对的第一分量和另外的 $n$ 对的第二分量的未知排列，称为分布中的破碎随机样本。研究了两种方法，其中排列向量要么被视为固定的未知参数，要么是需要从观察到的可能性中积分的随机变量。\
&emsp; &emsp;他们发现在第一种方法下， $MLE$ 是在所有可能的样本相关性中可以计算出的最大或最小样本相关性，这作为一般方法是不合理的。类似地，将我们自己限制在 $A$ 和 $B$ 之间的一对一匹配的情况下，这种方法除了将 $\omega$ 作为基于链接数据的分析参数外，还作为需要估计的参数。例如，如果分析中涉及的变量独立于用于记录链接的关键变量的误差，那么 $\omega$ 的可能性可以从增益参数的可能性中进行因式分解。那么，以上对协议分区的讨论表明， $\omega$ 的 $MLE$ 在目前的情况下不太可能做得很好。例如，在 $|\Omega_0| = 1$ 的情况下，一个合理的关键变量误差模型可能很好地暗示了链接矩阵 $\omega_L$ 作为真实匹配矩阵的MLE，并且不会对链接误差进行有用的调整。或者，在 $|\Omega_0| > 1$ 的情况下，很可能存在多个 $MLE$ 的情况。\
&emsp; &emsp;对于第二种方法， $DeGroot$ 和 $Goel(1980)$ 考虑了对排列向量进行积分后观察到的似然性，以及由此推导出的轮廓似然。基于模拟在 $n = 5$ 的情况下，其中显式积分是可行的，他们发现轮廓的可能性普遍较差。此外，虽然综合似然看起来与完整样本似然非常相似，但对于其他样本，两者看起来非常不同。这仍然并不一定意味着破碎的样本 $MLE$ 是渐近不一致的。然而，要研究这一点，需要一种基于大破碎样本的 $MLE$ ，在这方面显式集成是不切实际的。作者认为，在这种情况下， $MLE$ 可以通过 $EM$ 算法获得，通过将置换向量视为缺失的数据。但是他们并没有检查这个 $MLE$ 的性能。\
&emsp; &emsp;下面我们用一个基于链接数据的回归分析的例子来研究这一点。不幸的是，基于支持 $EM$ 算法的缺失信息原理（ $Orchard$ 和 $Woodbury$ ， $1972$ ）， $MLE$ 可能是有偏差的和不一致的。\
&emsp; &emsp;考虑一对一匹配的设置，其中 $U_A = U_B = \phi$ 和 $n_A = n_B = n$，以及 $\omega$ 是一个排列矩阵，即 $\omega\omega^T = \omega^T\omega = I$ ，其中 $I$ 是单位矩阵。考虑线性回归模型:\
 $$ y_i = x_i^T\beta + \epsilon_i$$ \
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

## 以下内容刘奕彤编写
### 2.4.3 比较数据建模(Ⅰ)

首先考虑一种比较数据建模的方法，它可以与FS范式下的概率记录联系起来。在保持普遍性的同时，让$g_a^*=h$作为$a \in A$的完整比较数据，同时$g_b^*=j$作为$b \in B$的完整比较数据，对于h,j=1,…,K 。正如Stoerts et al.（2015）等人，我们假定比较数据的潜在失真在实体之间是独立的，因此
$f(G^*|ω)= \prod\limits_{(ab) \in M}f(g_a^*,g_b^*|ω_ab=1) \prod\limits_{a \in U_B}f(g_a^*|a \in U_A) \prod\limits_{b \in U_B}f(g_b^*|b \in U_B)$

请注意，原则上，数据链接者可以提供这些概率作为比较数据，而不用披露关键变量。
当只有比较数据$G=\{g_{ab};a \in A_L,b \in B_L\}$ \\可用时，也就是实践中最有可能出现的情况，此时，可以从整合出$g_a^*$和$g_b^*$在完整的比较数据模型中的$g_a^*$和$g_b^*$进行积分，从而得到:

$f(G|\omega )= \prod\limits_{(ab):\omega_ab=1}f(g_ab|\omega_ab=1)$(2.3)

为了建立相关的概率模型，考虑Copas和Hilton（1990）提出的命中-失误失真模型模型，该模型由Copas和Hilton（1990）提出。设想一个二进制试验导致概率为$\alpha$的失误，随后观察到的值h是以概率$p_h$随机分布的，其中$p_h$是真实值h的频率。让$K=n_A +n_B-m$是A和B中实体的数量，其中$m=|M|$。假设所有实体的完整比较数据是不同的，因此以便$p_h=1/K$。对于一个简单的比较数据模型，让为协议指标。我们有:

$\mu_{ab}=f(g_{ab}=1|\omega_{ab}=1)=(1-\alpha)^2+2\alpha(1-\alpha)/K+\alpha^2/K$

$\mu_{ab}^c=f(g_{ab}=0|\omega_{ab}=1)=1-\mu_{ab}=(1-(1-\alpha)^2)(1-1/K)$

这可能被称为可交换协议（EA）模式。次要的分析师只需要提供K（或m）的估计值。

似乎很直观的是，在适当的条件下，EA模型可以产生ELE模型。适当的条件下，EA模型可以产生ELE模型。例如，假设一对一完全匹配，其中综合数据的分析h,j = 1,……,n 和$U_A = U_B =∅ $。在 EA 模型下，让$ \mu_ab=\mu(\alpha)$，其中$ g_{ab}$是协议指标。假设记录连接进行如下进程：如果它们是唯一一致的，则链接(ab)；如果没有唯一一致的记录，则随机链接记录中随机链接，但要遵守一对一的链接约束。这将导致ELE模型，其参数为$\lambda(\alpha,\phi)=q_{ii}=f(P_{ii}=1|\omega_{ii}=1)$, 其中$\omega=I_{n \times n}$的定义。

| $\alpha$ | n=3


| n=10 | n=50        | n=300       | n=1000      |             |             |
| ---- | ----------- | ----------- | ----------- | ----------- | ----------- |
| 0.02 | 0.974/0.987 | 0.965/0.981 | 0.961/0.975 | 0.961/0.964 | 0.960/0.962 |
| 0.05 | 0.935/0.966 | 0.912/0.950 | 0.904/0.925 | 0.903/0.910 | 0.903/0.907 |
| 0.10 | 0.873/0.937 | 0.829/0.895 | 0.814/0.846 | 0.811/0.827 | 0.810/0.824 |
| 0.20 | 0.760/0.878 | 0.676/0.790 | 0.647/0.712 | 0.641/0.691 | 0.640/0.688 |

虽然$\mu(\alpha)$有一个封闭的表达式，但$\lambda(\alpha,\phi)$只能通过G的蒙特卡洛模拟和随后的联系程序获得$\phi$。表2.2提供了$\mu$和$\lambda$给定的$\alpha$和n之间的比率，其中每个都是基于5000次模拟的结果。显然，在所有情况下$\mu /\lambda<1$。在一端，观察在任何固定的情况下，$\mu /\lambda→1$为$n\longrightarrow\infty$，因为在失误之后达成一致的机会趋于零，同时也趋于零。错过的机会趋向于零，而通过随机连接匹配的机会也是如此。在另一端，对于较小的n和较大的 ，$\mu /\lambda$进一步远离1，其中在这里，一对一联系限制的效果直观上是最强的。

表2.3

对于比较模型下的线性回归的一个非常简单的说明模型（2.3）下的线性回归，假设n=3的一对一完全匹配。矩阵$Q_G$为$\beta_G$由（2.2）取决于观察到的G。表2.3给出了两个例子。
$\alpha=0.02$时，观察到$G_1$的概率$f（G_1|\omega）$是0.923$G_1$，鉴于此
$\beta_G$的效率几乎与已知的$\omega$一样，因为$Q_G \approx I$ 。当$\alpha =0.1$观察到$G_1$的概率下降到0.666。回归的有效性低于$G=G_2$，从相应的$Q_G$可以看出。差异根植于$f(G|\omega)=\mu^{d}(1-\mu)^{n-d}$ 其中是观察到的$M(\mu)$的吻合度。让$d_0(G)$为给定G下最大可能的d，当$d_0(G_1)=3$且$d_0(G_2)=2$事，我们有：

$\sum_{\omega\in\Omega}f(G|\omega)=\sum_{d=0}^{d_{0}(G)}m_{d}f_{d}$

其中$f_d=f(G|\omega)$对于任何具有相同d的$\omega$来说是相同的，$m_d$是这样的数量$\omega's$. 似然比为$f_{d+1}/fd=\mu/(1-\mu)$.由$G_1$引起的协议分区是$m_d=(2,3,1)$，对于d=(0,1,3)，并且$\omega=G_1$的可能性是$f_3$,是次高的$f_1$的$\mu^2/(1-\mu)^2$倍。而鉴于$G_2$，我们有$m_d=(2,2,2)$，因为$d=(0,1,2)$，其中有两个$\omega's$最高的可能性$f_2$，是次高$f_1$的$\mu/(1-\mu)$倍。观察到$G=G_2$的不确定性显然大于$G=G_1$。

### 2.4.4 比较数据建模(Ⅱ)

在实践中，多通道的确定性联系是很常见的。让$L_k = A_k\times B_k$包含
k=1,…,K的第k次传递的$n_k$个链接。这些是在第k次传递中彼此完全一致的记录。没有链接的的记录留待以后的检验。让$A_0$包含A和$B_0$中的记录。联系的结果是分区:

$A\times B=\bigcup_{k=0}^{K}\bigcup_{j=0}^{K}A_{k}\times B_{j}$

让$g_{ab} = (k, j)$表示记录对(ab)属于分区的这些域的集合构成了比较数据G。

对于任何$（ab）\in  L_k$，让$\pi_{k}=\,f\bigl(\omega_{a b}\,=\,1\vert g_{a b}\,=\,(k,k)\bigr)$，这样，$1-\pi_{k}$   是$L_k$中的FLR（如表2.1）。我们有

$\mu_{a b}=f{\big(}g_{a b}=(k,k)|\omega_{a b}=1{\big)}=\pi_{k}f{\big(}g_{a b}=(k,k){\big)}/f(\omega_{a b}=1)$

人们可以通过$\pi_k n_k/m$来估计$\mu_{ab}$。理想情况下，我们希望数据提供者
能够提供$\mu_{a b}=\,f\bigl(g_{a b}=(k,j)\bigr)\bigl\vert\omega_{a b}=\,1\,\bigr)$对所有其他的组合

(k, j)也是如此。然而，这在实践中可能是很苛刻的。作为一个可能的
简化，让$L=\bigcup_{k=1}^{K}L_{k}$包含所有的链接对，并让

$\mu_{a b}=f\Bigl((a b)\ll L|\omega_{a b}=1\Bigr)=1-\sum_{k=1}^{K}f\Bigl(g_{a b}=(k,k)|\omega_{a b}=1\Bigr)$

其中假定$g_{ab}$在连接对之外的均匀分布。

当然，上述比较模型是对现实的极大简化。
因为匹配不可能均匀地分布在
$R = (A\times B) / L$。事实上，在所有未链接的记录对中，我们可以期待
预计$m_{R}=m-\sum_{k=1}^{K}\pi_{k}n_{k}$的匹配。我们可以探索$\mu_{ab}$的界限.例如，让$R_k$成为R的子集，其中涉及$a \in A_k$或$b \in B_k$，在这里我们最多可以期待$2n_{k} 1-\pi_{k})$匹配。因此，

$ f{\bigl(}(a b)\in R_{k}{\vert}\omega_{a b}=1{\bigr)}\leq2(1-\pi_{k})n_{k}/mf$，这就提供了一个$\mu_{ab}$在$R_k$上的上限.此外，让$R_0=A_0 \times B_0$，其中可以预期至少有

$m_{R}-2\sum_{k=1}^{K}(1-\pi_{k})n_{k}$k的匹配，因此$\mu_{ab}$在$R_0$上的下限是

$R_{0}|\omega_{a b}=1\rangle\geq1-\sum_{k=1}^{K}(2-\pi_{k})n_{k}/m$

让我们考虑C-PR的情况来说明。让$n_A$是人口普查数据集的大小，$n_B$是病人登记册的大小。让$n_L$是链接的数量（表2.1），m为匹配实体的数量。为了便于理解，假设m=51，nA=53，nB=55，都是以百万计。

* 对于第k个确定性通道，其中k=1,…9，FLR $1 -\pi_k$给出如表2.1我们设定  $\mu_{ab} = \pi_k n_k/m$ 为$（ab）\in  L_k$。

- 除了表2.1之外，我们没有其他关于两个概率性，分别用$L_10$和$L_11$表示。因此我们把它们就像从确定的联系中产生的一样。
- 我们还计算了$R_k$中的上界概率$\mu_{ab}$，和$R_0$中的下限概率 ˜$\mu_{ab}$。

表2.4

当涉及到大型人口时，可扩展性是一个重要的问题，比如在C-PR案例中。考虑一下评估的一般任务：

$E(t_{\omega}|G)=\sum_{\Omega}t_{\omega}f(\omega|G)=\sum_{\Omega}t_{\omega}f(G|\omega)/\sum_{\Omega}f(G|\omega)$

此时$t_\omega$ 是$\omega$的一个函数，正如在2.2中作为$Q_G$的$t_\omega=\omega^T$

完整列举$\omega$在大多数情况下是不可行的。让我们从$\omega$
中去一个样本$ s=\{\omega_{1},...,\omega_{n_{\circ}}\}\subset\Omega$.如果对$\omega's$进行随机抽样，那么H´ajek类型的估计是$E(t_\omega|G)$是$\bar{t}_{H}:=\sum\omega{\in}s\,t_{\omega}f(G|_{c o})/{\sum}\omega{\in}s\,f(G|_{c o})$而如果$\omega\sim f(\omega|G)$,蒙特卡洛估计是$\bar{t}=\sum_{\omega\in s}t_{\omega}/n_{s}$ 。Metropolis-Hastings算法提供了一种解决后者的一般方法：在每一步中从当前的中提出一个新的$\omega's$，其中的提议是对称的，接受率为$min|\{1,\;f(\omega^{\prime}|G)/f(\omega|G)\}.$

在上面的比较数据模型下，我们可以直接分析C-PR联系数据的接受率。让$\omega=\omega_L$ 为与链接矩阵完全一致的匹配矩阵，其中当且仅当$(ab) \in L $时$\omega_{ab}=1$。在以下子矩阵中得到了一般性的描述$\{a,a^{\prime},a^{\prime\prime}\}\times\{b,b^{\prime},b^{\prime\prime}\}$:

$[\omega_{L}]=\left[\begin{array}{c c c}{{1}}&{{0}}&{{0}}\\ {{0}}&{{1}}&{{0}}\\ {{0}}&{{0}}&{{0}}\end{array}\right]$

$[\omega_{1}]=\left[\begin{array}{c c c}{{1}}&{{0}}&{{0}}\\ {{0}}&{{1}}&{{0}}\\ {{0}}&{{0}}&{{1}}\end{array}\right]$

$[\omega_{2}]=\left[\begin{array}{c c c}{{0}}&{{0}}&{{1}}\\ {{0}}&{{1}}&{{0}}\\ {{0}}&{{0}}&{{0}}\end{array}\right]$

$[\omega_{3}]=\left[\begin{array}{c c c}{{0}}&{{1}}&{{0}}\\ {{1}}&{{0}}&{{0}}\\ {{0}}&{{0}}&{{0}}\end{array}\right]$

其中$1\leq k\neq k^{\prime}\leq11$ $(ab) \in L_K$

每一矩阵中的第一行是子向量$(\omega_{a b},\omega_{a b^{\prime}},\alpha_{a b^{\prime\prime}})$，第二行和第三行也是如此.

* $\omega_L$和$\omega_1$之间的唯一区别是，$\omega_1$包含一个额外的匹配 (a'',b'')。根据表2.4,我们有$f(\omega_{1}|G)/f(\omega_{L}|G)=\mu_{a b;0}=0.0100,$，或者至少是0.0055.
* 根据$\omega_2$，$R_k$中的一对实际上是一个匹配，而不是$L_k$中的一对。有$f(\omega_{2}|G)/f(\omega_{L}|G)=f{\big(}(a b^{\prime\prime})\in R_{k}|\omega_{a b^{\prime}}=1{\big)}/\mu_{a b;k}\leq{\tilde{\mu}}_{a b;k}/\mu_{a b;k}$，这对于表2.4中的所有1＜k＜ 11接近于零。
* 对于$\omega_3$，$L_k$中的(ab)和$L_{k'}$中的(a',b')是错误的链接，而$R_k$中的(ab')和(a'b)匹配未能被链接。而$R_k$中的(ab')和$R_{k'}$中的(a'b)不能被链接。我们有$f(\omega_{3}|G)/f(\omega_{L}|G)\ \leq\tilde{\mu}_{a b;k}\tilde{\mu}_{a b;k^{\prime}}/\mu_{a b;k}\mu_{a b;k^{\prime}}$这对所有$1\leq k\mp k^{\prime}\leq1^{\circ}$来说几乎是零。

因此，Metropolis-Hastings算法在数值上是不够有效的。考虑到$m>n_L$，匹配矩阵有$\omega_L$ 作为链接记录的子矩阵，再加上$m-n_L$在$A_0 \times B_0$之间的额外匹配,都有相同的可能性。这种矩阵的数量是$\left((n_{A}-n_{L})!/(n_{A}-m)!\right)\Bigl((n_{B}-n_{L})!/(n_{B}-m)!\Bigr)/(m-n_{L})!$

在结果空间如此之大的情况下，任何对$f(\omega|G)$的再抽样方法都是非常不切实际的。



## 以下内容侯宇欣编写
# 2.5关于链接子集的分析
## 2.5.1非信息性平衡选择
&ensp;&ensp; 出现上述可伸缩性问题是因为我们的目标是利用存在于单独数据集之间的所有匹配实体。然而，由于错过了比赛，ω的可能性可能在一个非常大的空间上是平坦的。因此，在C-PR案例中，似乎很自然地问，是否有可能将分析限制在5000万条链接记录，甚至是它的一个子集。再次使用线性回归作为点情况，我们研究下面有效链接子集分析的条件。  
&ensp;&ensp;首先，假设A和B之间有一个匹配。令$\left(A_{1}, A_{2}\right)$是A和$\left(B_{1}, B_{2}\right)$是B的一个双分区，其中$\left(A_{1}, B_{1}\right)$形成一个完整的匹配空间（CMS），这意味着记录$A_{1}$没有一个与记录$b < B_{1}$匹配。此外，$\left(A_{2}, B_{2}\right)$形成一个CMS，提供$\left(A_{1}, B_{1}\right)$形成一个CMS。设$X_{1}$和$X_{2}$分别为$a_{1}$与和$a_{2}$相关的回归设计矩阵。设$y_{1}$和$y_{1}$分别为与$b_{1}$和$b_{2}$相关的因变量的向量。对于给定的$\left(X_{1}, X_{2}\right)$，分别用y_{M 1}=\omega_{1} y_{1}和$y_{M 2}=\omega_{2} y_{2}$表示因变量的真实向量。A和B之间的真正匹配矩阵是块对角线矩阵，其中$\omega_{1}$和$\omega_{2}$是对角线上的块。设$g_{1}$和$g_{1}$分别为比较数据。我们说，如果$\left(A_{1}, B_{1}\right)$形成一个CMS，并且比较数据G1是非信息的，那么$\left(A_{1}, B_{1}\right)$形成一个$\left(A, B\right)$的非信息完全选择（NICS）。NICS $\left(A_{2}, B_{2}\right)$也是如此。  
引理1  &ensp;设$\left(A_{1}, A_{2}\right)$是A和$\left(B_{1}, B_{2}\right)$与B的双分法，其中$\left(A_{1}, B_{1}\right)$和$\left(A_{2}, B_{2}\right)$构成$\left(A, B\right)$的NICS。让$\hat{\beta}_{G}$、$\hat{\beta}_{G 1}$和$\hat{\beta}_{G 2}$分别由基于$\left(X, y\right)$、$\left(X_{1}, y_{1}\right)$和$\left(X_{1}, y_{2}\right)$的（2.2）给出。让$T_{1}=X_{1}^{T} Q_{1}^{T} Q_{1} X_{1}$和$T_{2}=X_{2}^{T} Q_{2}^{T} Q_{2} X_{2}$使用。让$r=T_{1}^{-1} T_{2}$和$\phi=(I+r)^{-1} r $。我们有  

$$
\hat{\beta}_{G}=(I-\phi) \hat{\beta}_{G 1}+(I-\phi) r \hat{\beta}_{G 2}=\hat{\beta}_{G 1}+\left[(I-\phi) r \hat{\beta}_{G 2}-\phi \hat{\beta}_{G 2}\right] \\
$$
证明：通过NICS，我们得到了$T=(Q X)^{T}(Q X)=T_{1}+T_{2}$, 当$Q_{1}=E\left(\omega_{T} \mid G\right)$ 为 $y_{M}=\omega_{y}$,$Q_{1}=E\left(\omega_{1}^{T} \mid G\right)$ for yM1 = ω1y1 和 $Q_{2}=E\left(\omega_{2}^{T} \mid G\right)for yM2 = ω2y2.同样地，我们也有$(Q X)^{T} y=X_{1}^{T} Q_{1}^{T} y_{1}+X_{2}^{T} Q_{2}^{T} y_{2}$。使用$T_{1}^{-1}-(I+r)^{-1} r T_{1}^{-1}=(I-\phi) T_{1}^{-1}$，我们得到  

$$
\hat{\beta}_{G}=(I-\phi) T_{1}^{-1}\left(X_{1}^{T} Q_{1}^{T} y_{1}+T_{2} T_{2}^{-1} X_{2}^{T} Q_{2}^{T} y_{2}\right)=(I-\phi) \hat{\beta}_{G 1}+(I-\phi) r \hat{\beta}_{G 2} \\
$$
请注意，虽然$\hat{\beta}_{G}$不仅仅是$\hat{\beta}_{G 1}$和$\hat{\beta}_{G 2}$的凸组合，但我们确实有$(I-\phi)+(I-\phi) r=I$，因为  

$$
\left.(I-\phi) r-\phi=\left[I-\phi-(I+r)^{-1}\right)\right] r=\left[I-(I+r)^{-1}(r+I)\right] r=0 \\
$$
因此，引理1表明，基于一个链接子集的回归是有效的，只要它是所有记录的某些NICS-分区的一部分。当然，除非FLR为零，否则我们不能确定用$\left(A_{L}, B_{L}\right)$表示的链接记录是否形成了$\left(A, B\right)$的NICS。但是，也可以考虑将记录链接作为识别一个可能的NICS的一种手段。设$\pi_{L}$为$\left(A_{L}, B_{L}\right)$形成NICS的概率。然后我们有$E\left(\hat{\beta}_{G L} \mid G, X\right)=\beta$，其中$\hat{\beta}_{G_L}$、由（2.2）给出，如果$\left(A_{L}, B_{L}\right)$形成一个NICS。此外，当$\left(A_{L}, B_{L}\right)$不形成NICS时，设ω(L)为真匹配矩阵$\omega_{n_{A}} \times n_{B}$的$n_{A} \times n_{L}$子矩阵，由对应于BL的ω列组成。提供$\omega_{(L)}^{T} \mathbf{1}=\mathbf{1}$，也就是说，所有的记录匹配即使他们不是所有，真正的回归设计矩阵的$X_{M L}=\omega_{(L)}^{T} X $.引理2下面显示$\hat{\beta}_{G}$的偏见是有限的概率$\left(A_{L}, B_{L}\right)$不形成CMS，即$1-\pi_{L}$ 。  
引理2  让$\hat{\beta}_{G L}=T_{L}^{-1}\left(X_{L} Q_{L}^{T} y_{L}\right)$，其中$T_{L}=X_{L}^{T} Q_{L}^{T} Q_{L} X_{L}$和QL计算为$\left(A_{L}, B_{L}\right)$形成一个NICS。设$f\left(\Omega_{L} \mid G\right)$为$\left(A_{L}, B_{L}\right)$形成NICS的概率。假设$\omega_{(L)}^{T} \mathbf{1}=\mathbf{1}$，其中$\omega_{(L)}$是$B_{L}$的真正匹配子矩阵，我们有  

$$
E\left(\hat{\beta}_{G L} \mid G, X\right)=\pi_{L} \beta+\left(1-\pi_{L}\right) T_{L}^{-1} X_{L}^{T} Q_{L}^{T} E\left(\omega_{(L)}^{T} \mid G\right) X \beta \\
$$
&ensp;&ensp;请注意，在$B_{L}$中不是所有的记录都在A中有匹配的情况下，如果不对BL中不匹配的y值做一些额外的假设，就不能给出$\hat{\beta}_{G L}$的期望。然而，很明显，在创建链接集$\left(A_{L}, B_{L}\right)$时，它有助于增加$\left(A_{L}, B_{L}\right)$形成CMS的概率。虽然很难评估大型数据集的概率πL，但限制FLR似乎有助于该课程，因为如果FLR为零，对于任何链接集，我们都必须有$\pi_{L}=1$。  
&ensp;&ensp;上述结果涉及调整后的$\hat{\beta}_{G L}$（2.2）。对于大的人口数据集，如果能够使用面值估计器，将会以一种更直接的方式更加实用。   

$$
\hat{\beta}_{L}=\left(X_{L}^{T} X_{L}\right)^{-1}\left(X_{L}^{T} y_{L}\right)
$$
基于引理1和引理2背后的推理，我们现在开发了一个基于大数据集的βˆL的渐近结果。  
&ensp;&ensp;设$\left(A_{ML}, A_{cL}\right)$为AL和设$\left(B_{ML}, B_{cL}\right)$与BL的双分区，其中$A_{ML}$包含正确链接到$B_{ML}$的$n_{M}$记录，其中$n_{M} \leq n_{L} $和$A_{cL}$的$n_{c} = n_{L}−n{M}$记录。让设$\left(X_{ML}, X_{cL}\right)$与设$\left(A_{ML}, A_{cL}\right)$和设$\left(y_{ML}, y_{cL}\right)$与$\left(B_{ML}, B_{cL}\right)$关联。让$T_{M L}=X_{M L}^{T} X_{M L}$，$T_{c L}=X_{c L}^{T} X_{c L}, r_{L}=T_{M L}^{-1} T_{c L} $和$\phi_{L}=\left(I+r_{L}\right)^{-1} r_{L}$。即使不知道哪些记录在$A_{ML}$和$B_{ML}$中，我们也有  

$$
\hat{\beta}_{L}=\left(X_{M L}^{T} X_{M L}+X_{c L}^{T} X_{c L}\right)^{-1}\left(X_{M L}^{T} y_{M L}+X_{c L}^{T} y_{c L}\right) \\
=\left(I-\phi_{L}\right) T_{M L}^{-1} X_{M L}^{T} y_{M L}+\left(I-\phi_{L}\right) r_{L} T_{c L}^{-1} X_{c L}^{T} y_{c L} \\
$$
假设ycL存在XMcL，其中$E\left(y_{c L} \mid X_{M c L}\right)=X_{M c L} \beta$，尽管BcL中的记录是错误链接的。考虑到大数据集和低FLR，我们提出了以下条件。渐近地，如nL→∞，设$\left(A_{L}, B_{L}\right)$被认为形成一个非信息平衡选择（NIBS）  

$
\text { (B1) } n_{M} / n_{L} \stackrel{P}{\rightarrow} \pi_{M}, \text { where } \pi_{M} \in(0,1] \\
$

$
\text { (B2) } n_{M}^{-1} X_{M L}^{T} X_{M L} \stackrel{P}{\rightarrow} \mu_{x} \mu_{x}^{T}+\Sigma_{x} \text { and } n_{c}^{-1} X_{c L}^{T} X_{c L} \stackrel{P}{\rightarrow} \mu_{x} \mu_{x}^{T}+\Sigma_{x} \\
$

$
\text { (B3) } n_{c}^{-1} X_{c L}^{T} X_{M c L} \stackrel{P}{\rightarrow} \mu_{x} \mu_{x}^{T} \\
$


例如，如果XML、XcL和XMcL是X的随机子集。  
定理1  提供了（B1）-（B3），渐近为$n_{L} \rightarrow \infty$，我们有  

$$
E\left(\hat{\beta}_{L}\right) \stackrel{P}{\rightarrow}\beta-\left(1-\pi_{M}\right)(I-\psi) \beta \\
$$

$$
n_{L} V\left(\hat{\beta}_{L}\right) \stackrel{P}{\rightarrow} \sigma^{2} \Delta_{x}^{-1} \\
$$

$
\text { where } \psi=\left(\mu_{x} \mu_{x}^{T}+\Sigma_{x}\right)^{-1} \mu_{x} \mu_{x}^{T} \text { and } \Delta_{x}=\mu_{x} \mu_{x}^{T}+\Sigma_{x} \\
$

证明：从（B1）和（B2）中，我们有$r_{L}\stackrel{P}{\rightarrow} I\left(1-\pi_{M}\right) / \pi_{M} $和$I-\phi_{L} \stackrel{P}{\rightarrow} \pi_{M} I$。对于（2.4）的第一项，我们可以得出这个结论  

$$
\left(I-\phi_{L}\right) T_{M L}^{-1} X_{M L} y_{M L} \stackrel{P}{\rightarrow} \pi_{M} \beta=\beta-\left(1-\pi_{M}\right) \beta \\
$$
对于（2.4）的第二项，从（B3）得出，$T_{c L}^{-1} X_{c L}^{T} X_{M c L} \stackrel{P}{\rightarrow} \psi$，即$\left(I-\phi_{L}\right) r_{L} T_{c L}^{-1} X_{c L}^{T} y_{c L} \stackrel{P}{\rightarrow}\left(1-\pi_{M}\right) \psi \beta$。结合这两项，可以得到第一个结果。对于ˆβL的方差，我们可以观察到  

$$
V\left(\hat{\beta}_{L}\right)=\left(I-\phi_{L}\right) V\left(\hat{\beta}_{M L}\right)\left(I-\phi_{L}\right)^{T}+\left(I-\phi_{L}\right) r_{L} T_{c L}^{-1} V\left(X_{c L}^{T} y_{c L}\right) T_{c L}^{-1} r_{L}^{T}\left(I-\phi_{L}\right)^{T} \\
$$

$
\text { 当 } \hat{\beta}_{M L}=T_{M L}^{-1} X_{M L}^{T} y_{M L} \text {. } \\
$
请注意，V（βˆL）并没有用$\left(X_{M L}, X_{c L}, X_{M c L}\right) $来指定，因为这些都没有被观察到。因此，使用作为a的NIBS链路子集的渐近偏差以FLR 1−πM为界，并给出了β的渐近偏差调整估计量为  

$$
\hat{\beta}_{A}=\left(I-\left(1-\pi_{M}\right)(I-\hat{\psi})\right)^{-1} \hat{\beta}_{L} \\ 
$$
估计量（2.5）计算简单，只需要由（2.4）的面值βˆL和XL的估计均值和协方差矩阵。  
## 2.5.2 C-PR数据的说明
上述结果为C-PR数据提供了一种实用的方法（表2.1），即只制作具有几乎为零的错误链接概率的链接，并基于产生的链接子集调整分析，即所有本来可以创建的链接。当这个链接子集中的所有记录都是真匹配时，它们必然会形成一个CMS，而调整问题将成为缺失的数据之一。开发调整这种链接子集的一般方法超出了这里的范围。下面我们将考虑一些广泛的替代方案，以及数据提供商和二级分析师的隐含成本和收益。  
&ensp;&ensp;考虑线性回归$y_{i}=x_{i}^{T} \beta+\epsilon_{i} $，其中yi和xi分别驻留在两个数据集中，而$\beta_{p \times 1}$是回归系数的向量。假设$\quad \operatorname{Cov}\left(\epsilon_{i}, \epsilon_{j}\right)=\sigma^{2}$，如果$i=j$，0，如果$i \neq j $。为了我们这里的目的，假设xi是标准化的，其中$\mu_{x}=E\left(x_{i}\right)=0$和$\Sigma_{x}=V\left(x_{i}\right)$，因此定理1中的$\psi=0 $和$V(\hat{\beta})=O\left(1 / n_{L}\right) $提供了NIBS链接子集。  
&ensp;&ensp;从表2.1中可以看出，超过97%的链接记录的FLR低于0.01，其中的例外情况是确定性传递8和两个概率传递。请注意，这两个概率传递对数据提供者来说是最昂贵的，同时与最大的flr相关联。  
&ensp;&ensp;考虑L1，由从第一个确定性通道开始的3000万条链路组成，其中的FLR低至0.00011。二级分析人员可以(i)假设NICS并获得基于L1（2.2）给出的βˆGL1，或（ii）假设NIBS并根据定理1（2.5）获得$\hat{\beta}_{A 1}=\hat{\beta}_{L 1} / \pi_{1}$，其中ψ = 0和βˆL1是基于L1和$1-\pi_{1}=0.00011$的面值OLS（2.4）。显然，对于检测到的FLR，NICS假设并非如此，因此选项(i)并不完全合理，正如第2.4.4节所讨论的，为L1计算$\mathrm{Q}_{1}=E\left(\omega_{1}^{\mathrm{T}} \mid G\right)$是不切实际的。相比之下，选项（ii）是简单的，并且避开了与条件分布$f(\omega \mid G)$相关的计算。假设NIBS假设，$\hat{\beta}_{L 1} $和βˆA1的合理性似乎是直观的低FLR。  
&ensp;&ensp;类似的推理也分别适用于所有其他链接子集L2，...，L11。现在考虑结合子集估计数$\hat{\beta}_{A k}=\hat{\beta}_{L k} / \pi_{k}$。让$S \subseteq${1,2，...，11}是一个索引集。一种基于s中的链路子集的组合估计器  

$$
\hat{\beta}_{A s}=\sum_{k \in s} w_{k} \hat{\beta}_{A k} \quad \text{和 }V\left(\hat{\beta}_{A s}\right)=\sum_{k \in s} w_{k}^{2} V\left(\hat{\beta}_{A k}\right) \\
$$
其中，$1 / w_{k}=v_{k} W$、$W=\sum_{k \in s} 1 / v_{k }$和vk是可选择的常数。vk的一个自然选择是$v_{k}=V\left(\hat{\beta}_{A k}\right)$，前提是NIBS假设适用于$k \in S$的每个子集Lk。另一种可能性是简单的$v_{k}=1 / n_{k} $，所以工作=nk/ns，$n_{s}=\sum_{k \in s} n_{k}$。如果NIBS假设在所有子集中都成立，那么这比$v_{k}=V\left(\hat{\beta}_{A k}\right) $的效率要低一些。但是，如果NIBS假设对于子集$L_{s}=\bigcup_{k \in s} L_{k}$的合并成立，但在每个子集中不是单独存在，那么它仍然是可信的。   

$$
\begin{table}[!htbp]
\centering
\caption{人口普查和患者登记链接子集回归的说明}
\label{tab:pagenum}
\begin{tabular}{llll}
\toprule
& 链接子集 &组合FLR & 标准错误顺序\\
\midrule
{1} &1 &1\\
{1,2,3,4,5,6,7,9} &1 &1\\
{1,2,3,4,5,6,7,8,9,10,11} &1 &1\\
\bottomrule
\end{tabular}
\end{table}
$$

$$
表2.5
\space \space 人口普查和患者登记链接子集回归的说明 \\
\begin{array}
{cll}
\hline
	链接子集  & 组合FLR  & 标准错误顺序 \\
\hline
	\{1\} &1.1 \times 10^{-4}  &1.80 \times 10^{-4} \\
	\{1,2,3,4,5,6,7,9\} &16.1 \times 10^{-4} &1.43 \times 10^{-4} \\
	\{1,2,3,4,5,6,7,8,9,10,11\} &25.2 \times 10^{-4} &1.41\times 10^{-4} \\
\hline
\end{array}
$$

&ensp;&ensp;表2.5显示了s的三个备选设置的一些汇总统计数据。组合的FLR由$1-\sum_{k \in s} n_{k} \pi_{k} / n_{s}$给出，标准误差（SE）的阶数为$O\left(1 / \sqrt{n_{s}}\right) $。第一种替代方案只使用$L_{1}$，提供了NIBS假设。第二种选择是使用flr低于0.01的8个子集。合并的FLR是单独的$L_{1}$的15倍，而SE有可能降低约20%。如果NIBS假设合理地支持相关子集的并集，则该替代方案可能近似有效。它确实需要为数据提供者做额外的工作。所有链接的记录都用于最终的选择。联合FLR比以前的替代方案高出约50%，而SE的潜在降低是最小的，可能不值得增加偏倚风险。如前所述，与直接确定性传递相比，这两个概率记录链接过程更需要资源。   
&ensp;&ensp;简而言之，基于具有接近零flr的链接子集的分析在所需的资源、二级分析人员和数据提供者的资源方面具有实际优势。似乎对潜在的连锁误差偏差的必要调整可能不那么复杂，并且由于比较数据模型的错误说明，可能具有较低的风险。  
## 2.6 结语
我们考虑了三种广泛的方法，关于辅助用户如何对不能没有错误地链接的数据集进行分析，提供相关记录之间的非分解比较数据。  
&ensp;&ensp;关于在尊重链接数据结构的模型下的MLE，目前仍然是一个有待解决的问题。对于一个旨在在单独的数据集中使用所有匹配的实体的分析，有两个关键的方法挑战。首先，可能需要更准确的比较数据模型，例如，以适应关键变量的异构失真误差，特别是对于小到中等大小的数据集。接下来，与条件分布$f(\omega \mid G)$相关的计算目前似乎无法伸缩。  
&ensp;&ensp;对于大种群数据集的分析，链接子集分析可以提供一个实用的替代方法。以上得到的结果，提供了渐近的NIBS链接子集，需要扩展，以处理不统一的或潜在的信息选择的链接子集。  
&ensp;&ensp;重要的是要观察数据提供者的关键角色，从我们的调查，准备数据集和传播适当的链接比较数据，为了二级分析师有可能进行有效的统计分析，并占的不确定性，因为不知道真正的匹配实体。  
### 参考文献
[1] Chambers, R.C. (2009). Regression Analysis of Probability Linked Data.Statisphere, Vol. 4: http://ro.uow.edu.au/eispapers/762/  
[2] Chambers, R.C. and Kim, G. (2015). Secondary Analysis of LinkedData. In Methodological Developments in Data Linkage (eds. K. Harron,H. Goldstein and C. Dibben), Chapter 5.  
[3] Chipperfield, J.O., Bishop, G.R. and Campell, P. (2011). Maximum likelihood estimation for contingency tables and logistic regression with incorrectly linked data. Survey Methodology, 37, 13–24.  
[4] Chipperfield, J.O. and Chambers, R.C. (2015). Using bootstrap to account for linkage errors when analysing probabilistically linked categorical data. Journal of Official Statistics, 31, 397–414.  
[5] Christen, P. (2012). A survey of indexing techniques for scalable record linkage and deduplication. ISEE Transactions on Knowledge and Data Engineering, 24.  
[6] Copas, J.B. and Hilton, F.J. (1990). Record linkage: Statistical models from matching computer records. Journal of the Royal Statistical Society,Ser. A, 153, 287–320.  
[7] DeGroot, M.H. and Goel, P.K. (1980). Estimation of the correlation coefficient from a broken random sample. The Annals of Statistics, 8, 264–278.  
[8] Fellegi, I.P. and Sunter, A.B. (1969). A theory for record linkage. Journal of the American Statistical Association, 64, 1183–1210.  
[9] Gilbert, R., Lafferty, R., Hagger-Johnson, G., Harron, K., Zhang,L-C., Smith, P., ... Goldstein, H. (2017). GUILD: GUidance for Information about Linking Data sets. Journal of Public Health, DOI:10.1093/pubmed/fdx037  
[10] Goldstein, H., Harron, K. and Wade, A. (2012). The analysis of recordlinked data using multiple imputation with data value priors. Statistics in Medicine, 31, 3481–3493.  
[11] Gutman, R., Afendulis, C.C. and Zaslavsky, A.M. (2013). A Bayesian procedure for file linking to analyze end-of-life medical costs. Journal of the American Statistical Association, 108, 34–47.  
[12] Gutman, R., Sammartino, C.J., Green, T.C. and Montague, B.T. (2015).Error adjustments for file linking methods using encrypted unique client identifier (eUCI) with application to recently released prisoners who are HIV+. Statistics in Medicine, 35, 115–129.  
[13] Herzog, T.N., Scheuren, F.J. and Winkler, W.E. (2007). Data Quality and Record Linkage Techniques. Springer.  
[14] Hof, M.H.P. and Zwinderman, A.H. (2012). Methods for analysing data from probabilistic linkage strategies based on partially identifying variables. Statistics in Medicine,31, 4231–4242.  
[15] Hof, M.H.P. and Zwinderman, A.H. (2015). A mixture model for the analysis of data derived from record linkage. Statistics in Medicine, 34,74–92.  
[16] Kim, G. and Chambers, R.C. (2012a). Regression analysis under incomplete linkage. Computational Statistics and Data Analysis, 56, 2756–2770.  
[17] Kim, G. and Chambers, R.C. (2012b). Regression analysis under probabilistic multi-linkage. Statistica Neerlandica, 66, 64–79.  
[18] Lahiri, P. and Larsen, M.D. (2005). Regression analysis with linked data. Journal of the American Statistical Association, 100, 222–230.  
[19] Orchard, T. and Woodbury, M.A. (1972). A missing information principle: theory and applications. Proc. Sixth Berkeley Symp. on Math. Statist. and Prob., Vol. 1, 697–715.   
[20] Owen, A., Jones, P. and Ralphs, M. (2015). Large-scale Linkage for Total Populations in Official Statistics. In Methodological Developments in Data Linkage (eds. K. Harron, H. Goldstein and C. Dibben), Chapter 8.  
[21] Sadinle, M. (2014). Detecting duplicates in a homicide registry using a Bayesian partitioning approach. Annals of Applied Statistics, 8, 2404–2434.  
[22] Scheuren, F. and Winkler, W. E. (1993). Regression analysis of data files that are computer matched. Survey Methodology, 19, 39–58.  
[23] Scheuren, F. and Winkler, W. E. (1997). Regression analysis of data files that are computer matched – Part II. Survey Methodology, 23, 157–165.  
[24] Stoerts, R., Hall, R. and Fienberg, S. (2016). A Bayesian approach tographical record linkage and de-duplication. Journal of the American Statistical Association, 111, 1660–1672.  
[25] Tancredi, A. and Liseo, B. (2013). A hierarchical Bayesian approachto record linkage and population size problems. The Annals of Applied Statistics, 5, 1553–1585.  
[26] Zhang, G. and Campbell, P. (2012). Data Survey:Developing the Statistical Longitudinal Census Dataset and identifying its potential uses.Australian Economic Review, 45, 125–133.  
