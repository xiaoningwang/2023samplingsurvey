## **以下内容由张馨引填写**
## **内容如下**
## **错误列表数据的对数-线性模型203**
9.3.2最大似然估计

最大似然估计 (MLE)  {$$\theta _ { ( a ) + 1 }$$}在模型 (9.3)  或 (9.4)  下，概念是简单的，只要有二项变量$$r _ { ( v ) }$$~二项式($$x _ { ( v ) }$$,  $$\theta _ { ( v ) }$$），  在每个交叉分类的列表域![img](file:///C:\Users\Silver\AppData\Local\Temp\ksohtml5976\wps8.jpg)中，(如果是 从a中随机抽样，  则只需要替换{ $$x _ { ( v ) }$$ }和{$$\theta _ { ( v ) }$$} 根据相应的观察域大小和误差 。然后给对数似然值为

 $$e(\alpha ;r)= \sum _{\omega \in \Omega _{K}}r_{\omega}\log \theta _{\omega}+y_{\omega}\log(1- \theta _{\omega})$$

此时$$y _ { ( v ) } = x _ { ( v ) } - r _ { ( v ) }$$。交叉分类的误差概率向量$$\theta$$可以从边缘指定的$$\theta$$ 中得到$$\theta _ { + }$$通过线性关系

$$E ( r _ { \alpha o t } ) = \sum _ { v _ { j \leq v } } E ( r _ { v } )$$ 

其中，$$\theta _ { + }$$ 由 (9.3)  或 (9.4)  指定 。设这两个向量为$$\theta$$ 和$$\theta _ { + }$$按相应的向量如下：

$$\left. \begin{array}  { l  }  { S _ { 2 } = \{ \{ 1 \} , \{ 2 \} , \{ 1 , 2 \} \} } \\ { C _ { 3 } = \{ ( 1 \} , \{ 2 \} , \{ 3 \} , \{ 1 , 2 \} , \{ 1 , 3 \} , \{ 2 , 3 \} , \{ 1 , 2 , 3 \} \} } \end{array} \right.$$

对于$$C _ { 2 }$$和$$ \Omega _{3}$$。排列$$x$$和$$x+$$，因此这样

$$ \theta =C \theta _{+}\quad for \quad C=Diag(x)^{-1}\Gamma Diag(x_{+})$$

其中，  对于K = 2，

$$T = \left( \begin{array}  { l l l  }  { 1 } & { 0 } & { - 1 } \\ { 0 } & { 1 } & { - 1 } \\ { 0 } & { 0 } & { 1 } \end{array} \right)$$

对于K =  3，

$$T = \left( \begin{array}  { l l l l l  }  { 1 } & { 0 } & { 0 } & { - 1 } & { - 1 } & { 0 } & { 1 } \\ { 0 } & { 1 } & { 0 } & { - 1 } & { 0 } & { 1 } \\ { 0 } & { 0 } & { 1 } & { 0 } & { - 1 } & { - 1 } & { 1 } \\ { 0 } & { 0 } & { 0 }$$

a的MLE现在可以用牛顿-拉夫森方法得到 。一阶二阶导数如下 。a的分数的相应分量，由

$$\frac { \partial l } { \partial \alpha _ { \alpha _ { \alpha } } } = d ^ { T } ( \frac { \partial \theta } { \partial \alpha _ { \alpha _ { \alpha } } } ) = d ^ { T } ( c \frac { \partial \theta _ { + } } { \partial n } \frac { \partial _ { \alpha } } { \partial _ { \alpha } } ) = d ^ { T } C  A _ { \alpha$$ 

 



204.综合数据分析

其中$$A _ { w }$$是∆的对应列，为$$A ^ { T } = T ^ { - 1 }$$和$$n = A \alpha$$在 (9.3)  和 (9.4)  下，以及衍生物的向量d都有分量

 $$d _ { \alpha } = \frac { \partial l } { \partial \theta _ { \alpha _ { \alpha } } } = \frac { r _ { \alpha } } { \theta _ { \alpha _ { \alpha } } } - \frac { y _ { \alpha _ { \alpha } } } { 1 - \theta _ { \alpha } }$$

和$$\partial \theta _ { + } / \partial n = D i a g ( w )$$，  向量 \omega 的分量为$$ \omega  _ {  \omega  } = \theta _ {  \omega  }$$，在对数模型 (9.3)  和$$w_{\omega}= \theta _{\omega}(1- \theta _{\omega})$$下的logit 模型  ( 9 . 4 )  。分数的向量是由

$$\frac { \partial l } { \partial \alpha } = ( \frac { \partial \theta } { \partial \alpha } ) ^ { T } d = ( \frac { \partial \theta _ { Q } } { \partial \alpha } ) ^ { T } C ^ { T } d = A ^ { T } D i a g (  \omega  ) C ^ { T } d$$

在进一步的推导中，  我们得到了二阶导数的矩阵

$$\frac { \partial ^ { 2 } l } { \partial \alpha \partial \alpha ^ { T } } = ( \frac { \partial \theta } { \partial \alpha } ) ^ { T }  ( \frac { \partial \theta } { \partial \alpha } ) + D$$

其中偏导数的向量b有分量

 $$b _ { ( v ) } = \frac { \partial d _ { c o } } { \partial \theta _ { c o } } = - \frac { r _ {  \omega  } } { \theta _ { c o } } + \frac { y _ { c o } } { ( 1 - \theta _ { c o } ) ^ { 2 } }$$

和D的列对应于$$W$$是由

$$D_{\omega}= \Delta ^{T}(\frac{\partial Diag( \omega )}{\partial \alpha _{\omega}})C^{T}d= \Delta ^{T}Diag(v)Diag(\Delta _{\omega})C^{T}d$$ 

而偏导数的向量$$V$$也有分量$$v _ { \alpha } = \partial  \omega  _ { \alpha } / \partial _ { T }$$和$$v _ { ( v ) } = \theta _ {  \omega  }$$，例如在 (9.3)  和$$v _ { \alpha } = \theta _ { \alpha } ( 1 - \theta _ { \alpha } ) ( 1 - 2 \theta _ { \alpha } )$$ 下 (9.4)  。

9.3.3基于列表调查数据的估计

在许多应用中观察到的人口单位可能来自于对目标人口的单独覆盖调查，但列举不足 。对于$$i E U$$， 让$$\operatorname { s s } ( i ) = 1$$ ，如果$$i E S$$，  否则为0 。让 $$n _ { \alpha } = | A _ { \alpha } n S |$$美国是列表域A中捕获的数量。列表命中和错误，例如 $$( y _ { \alpha } , r _ { \alpha } )$$本身并没有被直接观察到 。此外，名单外的子人口 $$U _ {  \theta } = U | A$$被部分枚举，产生$$n _ {  \theta }$$; 观察调查捕获和$$n _ {  \theta }$$;未观察到的调查错过了，在那里 $$N _ {  \theta } = n _ {  \theta } + m _ {  \theta }$$; 假设的完整数据包括两部分：该列表出现错误$$r _ { ( v ) }$$ (或者$$y _ { ( v ) }$$，在$$\alpha = Q _ { K }$$;

该调查捕获了$$n _ {  \omega  }$$，在A中，为$$w E Q _ { K }$$,同时$$n _ {  \theta }$$ 在 $$U _ {  \theta }$$，我们假设调查捕获到处遵循IID伯努利分布，捕获概率不变，$$\gamma = P r ( i  U _ { i }$$。在实践中，如果可行， 该假设可以与$$A U U$$的适当分层相结合 。给出$$x _ { ( v ) }$$，列表命中值和调查捕获值构成了从$$A _ {  \omega  }$$，有概率$$\sum _ { ( u ) }$$和$$\gamma$$，分别 。对于名单外的调查，捕获了$$n _ {  \theta }$$，我们将考虑以下两个选项。



205错误列表数据的对数-线性模型

**二项式模型** 让$$n _ {  \theta }$$遵循二项式分布$$( N _ {  \theta } , \gamma )$$, 其中N;被认为是一个固定的未知参数。

**负二项式模型**，让$$m _ {  \theta}$$遵循负二项分布$$( n _ { \theta} , 1 - \gamma )$$给定$$n _ {  \theta }$$预先确定 。 特别是$$N_{0}$$然后被视为一个随机变量 。与二项式模型相比，这将导致数据中损失一个自由度，同时同样减少了一个未知参数的数量 。最后，N∅的最佳预 测器 (BP)  由$$E ( N _ { \theta } | n _ {  \theta } ) = n _ {  \theta } / \gamma$$。

让$$y = ( y , \alpha )$$，且$$\alpha$$低于  (9. 3)  或  (9.4)  。在二项式模型下，基于观测数据$$( n , n _ { 0 } )$$的似然值为

$$L_{1}(\psi ,N_{\theta})\propto L \cdot L_{2}$$

$$L ( v ) \infty | \prod _ {  \omega  2 x } | _ { w } ^ { n _ { \omega } } ( 1 - p _ {  \omega  } ) ^ { x _ {  \omega  } - n _ {  \omega  } }$$

$$L _ { 2 } ( \gamma , N _ { \theta } ) \alpha \left( \begin{array}  { l  }  { N _ { \theta } } \\ { n _ { \theta } } \end{array} \right) \gamma ^ { n ^ { n } \theta } ( 1 - \gamma ) ^ { N _ { \theta } - n _ { \theta } }$$

(9.6)

 在$$p _ { w } = \gamma ( 1 - \theta _ { \alpha } ) = P r ( i  S | i E A _ { w } )$$情况下。

而在负二项模型下，基于观测数据的可能性$$( n , n _ { \theta } )$$是简单的$$L ( \phi )$$，在积分出m;从基于$$( n , n _ { \theta } , m _ { \theta } )$$，它只依赖于通过n得到的数据 。然而，众所周知，给定$$\gamma$$，$$L _ { 2 }$$最大值为$$N _ { \theta } = n _ { \theta } / \gamma _ { x }$$，不舍入为整数，即的BPN;在负二项式模型下 。最大限度地提高$$L _ { 2 }$$此外 ，在二项模型下等价于最大化$$L ( \phi )$$，并使用$$N_{\theta}$$的BP在负二项模型下 。换句话说没有实际的错误的模型假设。

应用EM算法来最大化 (9.6)  是很方便的，因为我们可以重用第9.3.2节中描述的MLE过程 。与$$L ( \phi )$$对应的算法，与$$L ( \phi )$$ 对应的全数据对数似然值是基于r和n的，  并由

$$l{C}(\psi ;r,n)= \sum _{\omega \in \Omega _{K}}r_{\omega}\log \theta _{\omega}+y_{\omega}\log(1- \theta _{\omega})+n_{A}\log \gamma +(y_{A}-n_{A})\log(1- \gamma)$$

其中$n_{A}= \sum _{\omega \in \Omega _{K}}^{n_{\omega}}$$和$$y _ { A } = \sum _ {  Q _ { K } } y _ { w }$$，EM算法如下。

**E-step** 给定当前的参数值，计算y和r的条件期望，给定$$n _ { ( w ) }$$,例如, 未观察到的列表命中和错误：

$$E(y_{\omega}|n_{\omega})=(x_{\omega}-n_{\omega})(1- \theta _{\omega})(1- \gamma)/(1-p_{\omega})$$

$$E(r_{\omega}|n_{\omega})=x_{\omega}-E(y_{\omega}|n_{\omega})$$

**M-step** 通过nA/E更新$$\gamma$$的估计值$$n _ { A } / E ( y _ { A } | n )$$，以及使用第9.3.2节中描述的基于$$E ( r | n )$$和$$E ( y | n )$$的程序。

在收敛时，我们得到了

 $$ \widehat{N}=n_{A}+n_{0}+n_{0}(1- \widehat{\gamma})/ \widehat{\gamma}=n_{A}+n_{0}/ \widehat{\gamma}$$



206.综合数据分析

 

![img](file:///C:\Users\Silver\AppData\Local\Temp\ksohtml5976\wps41.jpg) 

9.4具有零自由度的模型选择

9.4.1潜在似然比准则

在每一类日志或日志模型中，两个嵌套模型之间的模型选择可以同时利用拟合优 度和双参数数 。根据赤池国际银行的信息标准 。在它们之间存在一个问题，即从两个参数数量相同但具有不同链接函数的模型中选择一个，在这种情况下，形式上的选择只能取决于拟合优度 。然而，当可用的列表很少时，可能会发 生两种模型“耗尽”了数据中的所有自由度，而且都准确地符合数据。例如，在 K = 2的情况下，可能性 (9.6)  是基于3个自由变化的观察结果$$( n _ { 1 } , n _ { 2 } , n _ { 12 } )$$。调查捕获概率$$\gamma$$作为一个必要参数，通过设置a给出最大不饱和列表误差模型$$\alpha _ { 12 } = 0$$可以在 (9.3)  或 (9.4)  中使用。这使得正式测试拟合优度的自由度为零 。此外，这两种模型都可能精确地拟合数据，因此，即使是非正式地，也不可能知道哪一种拟合更好。下面我们阐述了一个潜在的似然比准则，最初在Zhang (1996)  中讨论过，它可能与这些情况下的模型选择有关。

其起点是基于似然比检验 (LRT)  的模型选择的标准理论 。例如Vuong (  1989)  基于独立和同分布 (IID)  的观测结果，解释了LRT的一般性质 。在目前的  背景下，完整的数据 $$( r , n )$$来自于的IID实现$$( 8 v , 8 s )$$，为所有的$$i \in A_{\omega}$$和$$w E Q _ { K }$$提供服务，其中r缺失，只观察到n 。观测数据n仍然是由IID伯努利试验有概率生成的$$P _ { omega }$$在每个$$A _ { omega }$$，可能性$$L _ { 1 }$$基于$$ n $$，由 (9.6)  给出，在工作模型(9.3)  或(9.4) 下。假设这两个模型有相同的$$\alpha$$项， 并且在链接函数中只有二项。

如果不假设一个工作模型一定是正确的，那么MLE$$\psi^{}$$在它收敛于所谓的伪真值 ，记为$$ \psi ^{**}=(\alpha ^{*}, \gamma ^{*})$$，  这是使偏差最小化的工作模型的参数值，例如, 散度 ( Vuong，  1989)  。伪真值可以等价地定义为由真期望推导出的MLE工作模型的对数似然，  在我们的情况下可以由

$$\psi ^ { * } = \operatorname { a r g m i n } E ^ { 0 } [ e ( \psi ; n ) ] = \arcsin ] ( ( \psi ; n ^ { 0 } )$$ 

 

其中$$n ^ { 0 } = E ^ { 0 } ( n )$$为nw的期望 。当对数似然在n中是线性的时， 数据的真实分布和最后一个等式都成立，这是指数分布族的任何模型的情况。

现在，为了区分两个模型 (9.3)  和 (9.4)  ，它们对a有相同的约束，我们将用来表示记录和logit。



错误列表数据的对数线性模型

LRT在两个模型之间的选择是不确定的

$$l ( \psi l _ { \log } ^ { * } n ^ { 0 } ) = l ( \psi l _ { \log i t } ^ { * } n ^ { 0 } ) = l ( n ^ { 0 } ; n ^ { 0 } )$$

207

 

 

(9.8)



其中$$l ( n ^ { 0 } ; n ^ { 0 } )$$表示根据定义使偏差最小化的饱和对数似然值 。提供 (9.8)  ，Vuong (1989)  将这两种模型作为等价模型 。然而，等价是不切实际的，只有两个模型给出N的不同估计 。此外，等价性不是真实的，而是由于r的缺失，因为两个模型可以被区分，否则两者记录和 对数几率是可识别的。

在 (9.8)  的情况下，我们可以考虑一个二次鲁棒性准则和一个基于模拟的模 型选择方法 。首先，给定缺失数据的情况和指数族的分布

完整数据 $$( r , n )$$，  由EM算法得到的MLE为

$$E(r|n; \widehat{\psi}_{\log})= \widehat{r}=E(r; \widehat{\alpha}_{\log})$$ 

其中第一个方程对应于在e步中的计算，第二个方程是在m步中求解的分数方程。提供了MLE的一致性 (登普斯特等人，1977；Wu，1983) 对于$$x _ { w } \rightarrow \infty$$在$$w E Q _ { K }$$中，我们有

$$x _ { \alpha } ^ { - 1 } E ( r _ { \alpha } | _ { \alpha } | _ { \alpha } | _ { \alpha } | _ { \alpha } | _ { \alpha } | _ { \alpha } | _ { \alpha } | _ { \alpha } | x _ { \alpha } | _ { \alpha } | _ { \alpha } | x _ { \alpha } | = x _ { \alpha } ^ { - 1 } | E ( r _ { \alpha } ^$$

换句话说，在模型$$\psi \log$$对数下，估计的条件期望收敛于假设的完全数据，这样就产生伪真实的结果记(⇤)录作为基于完全数据对数似然的MLE 。记录我们将$$E ( r ; \psi ) ^ { * } \log )$$ 

 称为模型$$\psi ^{*}\log$$下的伪真缺失数据记录.类似地，  设$$E ( r ; \psi ) ^ { * } \log )$$ 为模型下的 伪真缺失数据对数机率，同时，$$r _ { \omega }$$的真条件期望是这样的

$$x_{\omega}^{-1}E^{0}(r_{\omega}|n_{\omega})\rightarrow x_{\omega}^{-1}E^{0}(r_{\omega}|n_{\omega}^{0})=x_{\omega}^{-1}E^{0}(r_{\omega})$$

一般情况下，两个模型都不是真的，两种模型下的伪真数据都不相等，真数据也不相等$$E ^ { 0 } ( r _ { ( omega ) } )$$ ,  例如 ,

$$E(r_{\omega}|n_{\omega}; \psi _{\log}^{*})\neq E(r_{\omega}|n_{\omega}; \psi _{logit}^{*})\neq E^{0}(r_{\omega}|n_{\omega})$$

如果提供 (9.8)  ，  我们提出了一个潜在似然比准则 (LRC)  ，通过该$$ \psi _{\log}$$模型记录是否被认为比所提供的$$ \psi _{\logit}$$模型逻辑更健壮

 $$l_{C}(n,r_{logit}^{*};n,r_{\log it}^{*})-l{C}(\widehat{\psi}_{\log _{o}};n,r_{\log it}^{*}) \leq l_{C}(n,r_{\log}^{*};n,r_{\log}^{*})-\operatorname { l c } ( v \operatorname { l o g } i t ; n , r _ { l o g } )$$(9.9)
