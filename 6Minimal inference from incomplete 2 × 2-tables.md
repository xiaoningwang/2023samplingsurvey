2021216083006徐棱霄 （翻译第六章第二部分 P126-130）

$\left(\mathbf{M}_{2}\right)$ 参数 $\psi $ 的抽样分布是点可识别的，并且 MLE  $\widehat{\psi}$  使得 $\widehat{\psi} \stackrel{\operatorname{Pr}}{\rightarrow} \psi_{0}$ , 渐近为 n $\rightarrow \infty$ , 其中 $\psi_{0}$ 是真正的参数值。

在最小设置下，$\Theta(\psi)=[L(\psi), U(\psi)]$ , 其中 $L(\psi)$ 是 $\theta$  的下界，由  $\psi$ , 以及 $U(\psi)$ 的上限。标识区域为  $\Theta_{0}=\Theta\left(\psi_{0}\right)=\left[L_{0}, U_{0}\right]$ , 其中 $L_{0}=L\left(\psi_{0}\right)$ 以及 $U_{0}=U\left(\psi_{0}\right)$ 。因此，对于表 6.1 中的缺失数据设置，我们有  $\psi_{0}=\left(\lambda_{11}^{0}, \lambda_{01}^{0}, \lambda_{+0}^{0}\right)$ , 其中 $$\Theta_{0}=\left[L_{0}, U_{0}\right]=\left[\lambda_{11}^{0}, \lambda_{11}^{0}+\lambda_{+0}^{0}\right]$$对于匹配数据设置，我们有 $\psi_{0}=\left(\lambda_{1+}^{0}, \lambda_{+1}^{0}\right)$ , 和弗雷谢边界（Fréchet，1951）定义了识别区域  $$\Theta_{0}=\left[L_{0}, U_{0}\right]=\left[\max \left(\lambda_{1+}^{0}+\lambda_{+1}^{0}-1,0\right), \min \left(\lambda_{1+}^{0}, \lambda_{+1}^{0}\right)\right]$$使  $\widehat{L}=L(\widehat{\psi})$  and  $\widehat{U}=U(\widehat{\psi})$  分别是 $L_{0}$  和  $U_{0}$ 的MLEs, 并让  $\widehat{\Theta}=\Theta(\widehat{\psi})=[\widehat{L}, \widehat{U}]$  表示 $\theta$ 的最大配置文件似然估计量。 $\widehat{\Theta}$ 里面的点都可以被认为是同样最有可能的，即，根据观测数据模型下基于 $d_{n}$ 的可能性得到最佳支持。我们将 $\theta$ 的证实函数定义为 $\theta \in \Theta$，为 $$c(\theta ; \psi)=\operatorname{Pr}(\theta \in \widehat{\Theta} ; \psi)$$即, 给定值 $\theta$ 被 $\widehat{\Theta}$ 覆盖的概率，其中概率是相对于 $f\left(d_{n} ; \psi\right)$ 计算的。让实际的佐证是 $$c_{0}(\theta)=c\left(\theta ; \psi_{0}\right)$$即, 根据真实抽样分布进行评估。特别是，$c\left(\theta_{0} ; \psi_{0}\right)$ 是 $\widehat{\Theta}$ 的置信水平，作为 $\theta_{0}$ 的区间估计量。设观测到的佐证为 $\widehat{c}(\theta)=c(\theta ; \widehat{\psi})$ 由于 $\widehat{c}(\theta)$ 是 $c_{0}(\theta)$ 的 MLE，因此可以将观察到的佐证定义为给定观测数据的 $\theta$ 最可能的佐证级别。如图 6.1 所示的 OCBGT 数据，如果将观察到的佐证视为 $\theta$ 的函数，那么这个函数通常可以在 $\widehat{\Theta}$ 上变化，而不是在同一区域上平坦的剖面似然。请注意，在这种情况下，为了计算 $\widehat{c}(\theta)$, 其中 $\left(\widehat{\lambda}_{11}, \widehat{\lambda}_{+0}\right)=\left(n_{11} / n, n_{11} / n+n_{+0} / n\right)$ , 我们采用二元正态近似 $ \left(\widehat{\lambda}_{11}, \widehat{\lambda}_{+0}\right) \sim N_{2}(\mu, \Sigma)$ , 其中 $\mu=\left(\lambda_{11}, \lambda_{+0}\right)$  以及  $\sum$ 的独特元素是 $V\left(\widehat{\lambda}_{11}\right)=\lambda_{11}\left(1-\lambda_{11}\right) / n ,  V\left(\widehat{\lambda}_{+0}\right)=\lambda_{+0}\left(1-\lambda_{+0}\right) / n$  和 $\operatorname{Cov}\left(\widehat{\lambda}_{11}, \widehat{\lambda}_{+0}\right)=-\lambda_{11} \lambda_{+0} / n$ . 更一般地说，观察到的佐证可以通过模拟计算如下。

靴带法（Bootstrap）算 $\widehat{c}(\theta)$ 

对于给定的  $\theta$  和 MLE  $\widehat{\psi}$ , 重复  b=1, $\ldots B$  

- 从 $f\left(d_{n} ; \widehat{\psi}\right)$ 生成  $d_{n}^{(b)} $ 得到  $\widehat{\psi}^{(b)}$  和相应的 $\left[L\left(\widehat{\psi}^{(b)}\right), U\left(\widehat{\psi}^{(b)}\right)\right]$
- 如果 $\theta \in\left[L\left(\widehat{\psi}^{(b)}\right), U\left(\widehat{\psi}^{(b)}\right)\right]$ ，设 $\delta^{(b)}=1$，否则为0。

把 $\widehat{c}(\theta)=\sum_{b=1}^{B} \delta^{(b)} / B$  作为观察到的佐证的引导估计值，用于 $\theta $.

## 6.3 最大确证集

让level-$ \alpha$的确证集通过以下方式给出 $$A_{\alpha}(\psi)=\{\theta: c(\theta ; \psi) \geq \alpha\}$$假设这存在 $\theta \in A_{\alpha}(\psi)$  其中 $c(\theta ; \psi)=\alpha$ . 因此，根据定义，我们有 $c(\theta ; \psi)<\alpha$ , 对于任何 $\theta \notin A_{\alpha}(\psi)$ , 而我们不能有 $c(\theta ; \psi)>\alpha$  对所有的 $\theta \in A_{\alpha}(\psi) $.  $A_{\alpha}(\psi)$ 的一些特性给在下方. 注意我们用 $c(\theta)$  作为 $c(\theta ; \psi) $ 的简称以及 $A_{\alpha}$ 作为 $A_{\alpha}(\psi)$ 的简称, 这里没有必要强调它们对 $\psi$ 的依赖性 .

**定理 1**  假设适用最小推理设定，即，只要条件 $\left(M_{1}\right)$ 和 $\left(M_{2}\right)$ 成立。(i) 令  $A_{\alpha_{1}}=\left[L_{1}, U_{1}\right]$  及  $A_{\alpha_{2}}=\left[L_{2}, U_{2}\right]$ . 如果  $\alpha_{1}>\alpha_{2}$ , 那么  $\left[L_{1}, U_{1}\right] \subset\left[L_{2}, U_{2}\right]$ . (ii) 令  $\theta_{L}<\theta_{U}$ , 其中  $c\left(\theta_{L}\right)=c\left(\theta_{U}\right)=\alpha$ . 那么，对于任何 $\theta \in\left(\theta_{L}, \theta_{U}\right)$ ，$c(\theta) \geq \alpha$ 成立.

**证明** (i) 一方面，我们有 $A_{\alpha_{1}} \backslash A_{\alpha_{2}}=\emptyset$  因为，否则，一定存在着一些 $\theta \in A_{\alpha_{1}} \backslash A_{\alpha_{2}}$ 以致于 $c(\theta) \geq \alpha_{1}  (因为  \theta \in A_{\alpha_{1}}  )$ 及  $c(\theta)<\alpha_{2}  (因为  \theta \notin A_{\alpha_{2}}  )$ 同时存在, 与规定的 $\alpha_{1}>\alpha_{2}$ 相矛盾。另一方面，集合 $A_{\alpha_{2}} \backslash A_{\alpha_{1}}$ 是非空的，因为否则 $A_{\alpha_{2}}$中的每个 $\theta$ 必须属于$A_{\alpha_{1}}$ ，因此，$c(\theta) \geq \alpha_{1}$ ，因此，不存在$\theta \in A_{\alpha_{2}}$，使得$c(\theta)=alpha_{2}<\alpha_{1}$ ，这与$A_{\alpha_{2}}$的定义相矛盾。

(ii) 每个 $\widehat{\Theta}$ 可以分为4种不同的类型，分别表示为 (a)  $\widehat{\Theta}_{L} \bar{U}$  其中  $\theta_{L} \notin   \widehat{\Theta}$  以及 $\theta_{U} \notin \widehat{\Theta} , (b)  \widehat{\Theta}_{L U}$  其中  $\theta_{L} \in \widehat{\Theta}$  及  $\theta_{U} \in \widehat{\Theta}$  ，因此,  $\theta \in \widehat{\Theta}_{L U}$ , (c)  $\widehat{\Theta}_{L} $ 其中  $\theta_{L} \in \widehat{\Theta}$  及  $\theta_{U} \notin \widehat{\Theta}$ , (d)  $\widehat{\Theta}_{U}$  其中  $\theta_{L} \notin \widehat{\Theta}$  及 $ \theta_{U} \in \widehat{\Theta}$ . (c)型可进一步分为 (c.1)  $\widehat{\Theta}_{L 1}$  其中  $\theta \in \widehat{\Theta}_{L 1}$ 以及 (c.2)  $\widehat{\Theta}_{L 2}$  其中  $\theta \notin \widehat{\Theta}_{L 2}$ , 也就是说，取决于 $\theta$ 是否出现在  $\widehat{\Theta}$ . 类似地，(d)型进一步划分为(d.1)  $\widehat{\Theta}_{U 1}$ 其中 $\theta \in \widehat{\Theta}_{U 1}$ 及 (d.2)  $\widehat{\Theta}_{U 2}$ 其中  $\theta \notin \widehat{\Theta}_{U 2}$ . 我们有

$c\left(\theta_{L}\right)$ = $\operatorname{Pr}\left(\widehat{\Theta}_{L U}\right)$ + $\operatorname{Pr}\left(\widehat{\Theta}_{L}\right)$ = $\operatorname{Pr}\left(\widehat{\Theta}_{L U}\right)$+ $\operatorname{Pr}\left(\widehat{\Theta}_{L 1}\right)$ +$\operatorname{Pr}\left(\widehat{\Theta}_{L 2}\right) \\$

$c\left(\theta_{U}\right)$ = $\operatorname{Pr}\left(\widehat{\Theta}_{L U}\right)$+$\operatorname{Pr}\left(\widehat{\Theta}_{U}\right)$ =  $\operatorname{Pr}\left(\widehat{\Theta}_{L U}\right)$ + $\operatorname{Pr}\left(\widehat{\Theta}_{U 1}\right)$ + $\operatorname{Pr}\left(r \widehat{\Theta}_{U 2}\right) \\$
$c(\theta) \geq \operatorname{Pr}\left(\widehat{\Theta}_{L U}\right)$+$\operatorname{Pr}\left(\widehat{\Theta}_{L 1}\right)$+$\operatorname{Pr}\left(\widehat{\Theta}_{U 1}\right)$

因此，如果  $\operatorname{Pr}\left(\widehat{\Theta}_{U 1}\right) \geq \operatorname{Pr}\left(\widehat{\Theta}_{L 2}\right)$ ,那么  $c(\theta) \geq c\left(\theta_{L}\right)$ , 或者如果 $\operatorname{Pr}\left(\widehat{\Theta}_{U 1}\right) \leq \operatorname{Pr}\left(\widehat{\Theta}_{L 2}\right)$ ,，那么 $\operatorname{Pr}\left(\widehat{\Theta}_{L 1}\right) \geq \operatorname{Pr}\left(\widehat{\Theta}_{U 2}\right)$ 因为 $c\left(\theta_{L}\right)=c\left(\theta_{U}\right)$ , 这样 $c(\theta) \geq c\left(\theta_{U}\right)$ . 同样，在 $\operatorname{Pr}\left(\widehat{\Theta}_{L 1}\right)$ 和  $\operatorname{Pr}\left(\widehat{\Theta}_{U 2}\right)$ 的比较上.

**定理 2** 鉴于一个最小的推理设置，存在一个最大的确证值，用 $\theta^{\max }$ 表示，这样 $c\left(\theta^{max }\right) \geq c(\theta)$ 对任何 $\theta \neq \theta^{\max }$ 成立。

**证明** 取任何初始level-$\alpha_{1}$ 的佐证集 $ A_{\alpha_{1}}=\left[L_{\alpha_{1}}, U_{\alpha_{1}}\right]$ 。在不失一般性的情况下，根据定理1.i，其中一个端点必须有确凿的$\alpha_{1}$；假设 $ c\left(L_{\alpha_{1}}\right) \geq c\left(U_{\alpha_{1}}\right)=\alpha_{1}$ . 根据定义，对所有$\theta \in A_{\alpha_{1}}$ ， $c(\theta) \geq \alpha_{1}$ 成立。如果 $c(\theta)=c\left(L_{\alpha_{1}}\right)$ ，对所有 $L_{\alpha_{1}}<\theta<U_{\alpha_{1}}$ 成立, 那么  $\theta^{\max }=L_{\alpha_{1}}$ , 因为  $c(\theta)<\alpha_{1} \leq c\left(L_{\alpha_{1}}\right)$  对于任何  $\theta \notin A_{\alpha_{1}}$ 都成立。否则，出现 $L_{\alpha_{1}}<\theta<U_{\alpha_{1}}$ , 其中  $c(\theta)=\alpha_{2}>c\left(L_{\alpha_{1}}\right) \geq \alpha_{1}$ , 和相应的level-$\alpha_{2}$证实集，表示为 $A_{\alpha_{2}}=\left[L_{\alpha_{2}}, U_{\alpha_{2}}\right]$ . 根据定理1.i，我们有 $\left[L_{\alpha_{2}}, U_{\alpha_{2}}\right] \subset\left[L_{\alpha_{1}}, U_{\alpha_{1}}\right]$ . 因为  $\alpha \leq 1$ , 参数的迭代必须终止于某个最大的 level-$\alpha$ .

用 $A^{\max }=A^{\max }\left(\psi_{0}\right)$  表示最大佐证集，这样， $c_{0}(\theta)>c_{0}\left(\theta^{\prime}\right)$  对任何  $\theta \in A^{\max }$  和  $\theta^{\prime} \notin A^{\max }$ 成立，和  $c_{0}(\theta)=c_{0}\left(\theta^{\prime}\right)$  对任何 $\theta \neq \theta^{\prime} \in A^{\max }$ 成立。从（6.1）中可以看出，这些是$\widehat{\Theta}$意味着最高置信度的点，在这个意义上，我们可以认为这些是最难驳倒的参数值。用$\widehat{\psi}$代替$\psi_{0}$，我们得到$A^{\max }$的MLE或者观察到的最大确证集合。 $$\widehat{A}^{\max }=A^{\max }(\widehat{\psi})$$图6.2说明了在匹配数据环境下的确证性，其中 $\theta=\lambda_{11}$  。真实的抽样分布参数 $\left(\lambda_{1+}, \lambda_{+1}\right)$ 在左图中为(0.1,0.9)，右图为(0.3,0.3)。左边的样本大小为 $\left(n_{1}, n_{2}\right)=   (1000,500) $，右边为(200,300)。识别区域 $\Theta_{0}$ 是垂直虚线之间的区间，实心曲线显示实际确证度（图中表示cvalue）如何随$\theta$变化。一些 $\Theta_{0}$ 的内部点的佐证可以是1，而对于许多 $\theta \notin \Theta_{0}$ （非内部点），它可以是0。在左图中，$c_{0}\left(L_{0}\right)$ 和 $c_{0}\left(U_{0}\right)$都约为0.5；在右图中，我们有 $c_{0}\left(L_{0}\right)=1$ 和 $c_{0}\left(U_{0}\right) \approx 0.25$。

让 $\bar{c}(\theta ; \psi)=\lim _{n} c(\theta ; \psi)=\lim _{n} \operatorname{Pr}\left(\theta \in \widehat{\Theta}_{n} ; \psi\right)$  是$\theta$在$\psi$处评价的渐进确证 , 其中 $\lim _{n}$ 代表  $\lim _{n \rightarrow \infty}$ 和 $\widehat{\Theta}_{n}$ 明确了对样本量的依赖性。表6.2总结了渐进的实际确证率 $\bar{c}_{0}(\theta)=\bar{c}\left(\theta ; \psi_{0}\right)$ 对于这两个数据设置。让 $\bar{A}^{\max }$ 为渐进的最大实际确证集，基于  $\bar{c}_{0}(\theta)$ . Lemma 1指出，

**图 6.2**
匹配数据设置中的确证说明。左边： $\left(\lambda_{1+}, n_{1}\right)=   (0.1,1000)$ 和  $\left(\lambda_{+1}, n_{2}\right)=(0.9,500)$ . 右边：$\left(\lambda_{1+}, n_{1}\right)=(0.3,200)$ 和 $\left(\lambda_{+1}, n_{2}\right)=(0.3,300)$ 

**表 6.2**
在缺失和匹配数据的情况下，$\bar{c}_{0}(\theta)$渐进的实际确证

|数据设置|$\theta \notin\left[L_{0}, U_{0}\right]$|$\theta=L_{0}$|$\theta \in\left(L_{0}, U_{0}\right)$|$\theta=U_{0} \\$|
|--|--|--|--|--|
|**缺失**|0|0.5 if  $L_{0}>0$|1|0.5 if $U_{0}<1$|
|**匹配**|0|0.5 if  $\lambda_{1+}+\lambda_{+1} \geq 1$|1|0.5 if $\lambda_{1+} \neq \lambda_{+1}$|
|||1 if $ \lambda_{1+}+\lambda_{+1}<1$||0.25 if  $\lambda_{1+}=\lambda_{+1}$|

除了$L_{0}$和$U_{0}$的界限外，$\bar{A}^{\max }$与$\Theta_{0}$无法区分，$\bar{c}_{0}(\theta)$是$\Theta_{0}$的指标函数。定理3指出，观察到的最大确证集$\widehat{A}_{n}^{\max }$的内部在概率上收敛于$\Theta_{0}$的内部。

**定理 1** 给定一个最小推理设定,  $\theta \in \bar{A}^{\max }$  及  $\bar{c}_{0}(\theta)=1$  如果  $\theta \in  Int  \left(\Theta_{0}\right)=\left(L_{0}, U_{0}\right)$ , 也就是说，如果$\theta$属于$ \Theta_{0}$的内部，那么 $\theta \notin \bar{A}^{\max }$ 和 $\bar{c}_{0}(\theta)=0$ 对任何  $\theta \notin\left[L_{0}, U_{0}\right]$ 都成立.

**证明**  令 $\delta\left(\theta ; \widehat{\psi}_{n}\right)=1$，如果  $\theta \in \operatorname{Int}(\widehat{\Theta})=\left(\widehat{L}_{n}, \widehat{U}_{n}\right)$ , 否则为0，其中$\widehat{\psi}_{n} $为MLE。在不丧失一般性的情况下，对于任何 $ \theta=U_{0}-\epsilon$ , 其中 0 < 2$ \epsilon<U_{0}-L_{0}$ ,我们有 $\delta\left(\theta ; \widehat{\psi}_{n}\right)=1$  if  $\left|\widehat{U}_{n}-U_{0}\right|<\epsilon$ 和 $\left|\widehat{L}_{n}-L_{0}\right|<\epsilon$ , 的概率趋向于1，因为 $\widehat{\psi}_{n} \stackrel{P r}{\rightarrow} \psi_{0}$ . 因此,  $\delta\left(\theta ; \widehat{\psi}_{n}\right) \stackrel{P r}{\rightarrow}  1$, 即,  $\bar{c}_{0}(\theta)=1  和  \theta \in \bar{A}^{\max } $. 同样地，可以证  $\bar{c}_{0}(\theta)=1$ , 对于  $\theta \notin \Theta_{0}$ ,即  $\theta \notin \bar{A}^{\max }$ .

**定理 3** 鉴于最小推理的设定，我们有 $\operatorname{Int}\left(\widehat{A}^{\max }\right) \stackrel{\operatorname{Pr}}{\rightarrow} \operatorname{Int}\left(\Theta_{0}\right)$ ; 即，如果 $\theta \in \operatorname{Int}\left(\Theta_{0}\right)$，$\lim _{n} \operatorname{Pr}\left(\theta \in \widehat{A}_{n}^{\max }\right)=1$ 以及如果 $\theta \notin \Theta_{0}$，$\lim _{n} \operatorname{Pr}\left(\theta \in \widehat{A}_{n}^{\max }\right)=0$ .
**证明** 根据Slutsky定理的一般形式（例如定理7.1，Kapadia等人，2005），我们有 $\bar{c}\left(\theta ; \widehat{\psi}_{n}\right) \stackrel{\operatorname{Pr}}{\rightarrow} \bar{c}\left(\theta ; \psi_{0}\right)$ , 因为 $\widehat{\psi}_{n} \stackrel{P r}{\rightarrow} \psi_{0}$  以及对于所有的$\psi$ ，$\bar{c}(\theta ; \psi)$ 是一个有界的。因此，如果 $\theta \in\left(L_{0}, U_{0}\right)$ , 根据定理1，使得 $\bar{c}\left(\theta ; \psi_{0}\right)=1$，我们有  $\bar{c}\left(\theta ; \widehat{\psi}_{n}\right) \stackrel{\operatorname{Pr}}{\rightarrow} \bar{c}\left(\theta ; \psi_{0}\right)=1$，意味着  $\lim _{n} \operatorname{Pr}\left(\theta \in \widehat{A}_{n}^{\max }\right)=1$同样地，可以证明  $\lim _{n} \operatorname{Pr}\left(\theta \in \widehat{A}_{n}^{\max }\right)=0$ , 对于 $\theta \notin \Theta_{0}$ .

## 6.4 对 $\Theta_{0}$ 的高保证估计  

在最小推理环境下，$ \Theta_{0}$ 的置信区$C_{n}$（是一个区间）的置信度为$\operatorname{Pr}\left(\Theta_{0} \subseteq C_{n}\right)$；例如，见Chernozhukov等人（2007）。在高置信度下，$ C_{n}$包含不属于$\Theta_{0}$的点的概率也一定很高，这是由于抽样的可变性，因此$C_{n}$从 "外部 "渐进地向$\Theta_{0}$收缩。相比之下，$\Theta_{0}$中的任何一点都是不可辩驳的，而$\widehat{A}^{\max }$则确定了那些在观测数据下最难辩驳的参数值。因此，我们将$\widehat{A}^{\max }$的保证定义为 $$\tau_{0}=\operatorname{Pr}\left(\widehat{A}^{\max } \subseteq \Theta_{0}\right),$$其中概率是相对于$f\left(d_{n}; \psi_{0}\right)$评估的。也就是说，这是观察到的$\widehat{A}^{\max }$中的点确实都是不可辩驳的概率。如果$\widehat{A}^{\max }$有很高的保证，那么它包含$\Theta_{0}$以外的点的概率就会很低。随着样本量的增加，$\Theta_{0}$的高保证估计值应该从它的 "内部 "向$\Theta_{0}$增长。根据定理1，对于一些小常数 $h\geq 0$，$\Theta_{0}$的高保证估计值可以定义为

$$\widehat{A}_{h}=\left\{\theta: c(\theta ; \widehat{\psi}) \geq \max _{\theta} c(\theta ; \widehat{\psi})-h\right\}$$以下引导法可用于估计  $\widehat{A_{h}}$ , 包括 $\widehat{A}_{0}=\widehat{A}^{\text {max }}$ .

对 $\widehat{A}_{h}$ 进行引导分析 

给出 MLE  $\widehat{\psi}$  和相应的  $[\widehat{L}, \widehat{U}]$ , 重复 b=1, $\ldots B$  :

1. 从 $f\left(d_{n} ; \widehat{\psi}\right)$生成 $d_{n}^{(b)}$ , 并得到 $\widehat{\psi}^{(b)}$ ;
2. 对于任何给定的h，在 $0 \leq h<1$的情况下 ,  在 $\widehat{\psi}^{(b)}$ 处获得$\widehat{A}_{h}^{(b)} $ ，与在 $\widehat{\psi}$ 获得$\widehat{A}_{h}$ 的方法相同, 并且相应的  $L^{(b)}=L\left(\widehat{A}_{h}^{(b)}\right)$ 及 $U^{(b)}=U\left(\widehat{A}_{h}^{(b)}\right)$ ;
3. 设 $\delta^{(b)}=1$ 如果  $\widehat{L} \leq L^{(b)}<U^{(b)} \leq \widehat{U}$ ,否则  $\delta^{(b)}=0$。
