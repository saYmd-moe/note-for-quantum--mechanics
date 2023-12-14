// 导入模板
#import "template.typ": *


// 文章信息设置
#show: conf.with(
  title: "微扰理论（量子力学）",
  shortitle: "Perturbation Theory in Quantum Mechanics",
  author: "saYmd"
)

*写在前面：*

谐振子和不考虑相对论效应的氢原子具有简单的哈密顿算符，可以求出精确解析解。但是，对于更复杂的原子或分子（多电子原子），一般很难得到解析解。这时，我们可以通过一些 *近似方法* 得到本征值方程的近似解。

微扰理论就是其中的一种近似方法
#footnote([微扰理论是一种数学方法，并不只是运用在量子力学中，经典力学同样有非常类似的近似方法。])
，根据体系是否显含时间，可以分为定态微扰理论和含时微扰理论，根据体系是否存在简并，可以分为非简并微扰理论和简并微扰理论。使用微扰理论遵循的思想是：我们一般优先处理体系中起到主导作用的效应，将这部分解决后再去考虑在近似中忽略掉的次要效应。

= 定态微扰理论
== 方法概论

首先我们需要对微扰理论这一方法有一个大概的印象，考虑体系的哈密顿算符可写作下列形式：
$ H=H_0+H' $

#h(2em) $H^0$是 *未微扰哈密顿算符* ，它的本征方程应是容易解的，而$H'$是 *微扰项* ，考虑$H'$与时间无关时的 *定态微扰* 。为了方便标记，我们假设未微扰 哈密顿 算符的本征值构成一个 *离散谱*
#footnote([这里离散谱的假设是为了简化书写，微扰理论同样可以用于处理连续谱。])
。使用整数指标 n 来标记： $E_n^0$，对应本征态 $ket(n^0_i)$
#footnote([$ket(n^0_i)$，即$ket(psi_n^0)$，下标 i 用于区分简并。])
。这个未微扰的 哈密顿 算符的本征方程是可以解析求解的：
$ H^0 ket(n^0_i) = E_n^0 ket(n^0_i) $<non-pert>
微扰理论的作用就是近似求解施加微扰后的本征方程：
$ H ket(n_i) = (H^0+H') ket(n_i) = E_n ket(n_i) $<perted>

#h(2em) $H'$是 *微扰项* ，在实际物理问题中，微扰项可以来自外界扰动、体系相互作用或者体系内部结构的影响，比如放置在电磁场中的氢原子，$H'$就可以用来表示氢原子与电磁场相互作用引起的势能。@fig:inf_with_pert 展示了一个带有扰动的无限深方势阱。
#figure(
  cetz.canvas({
    import cetz.draw: *

    // 坐标轴
    let x = (-1,3)
    let y = (-.5,2)

    set-style(mark: (fill: none))
    line((x.at(0),0),(x.at(1),0), mark: (end: ">"), stroke: (paint: luma(30%)))
    content((x.at(1)-.1,-.2),[$x$])
    line((0,y.at(0)),(0,y.at(1)), mark: (end: ">"), stroke: (paint: luma(30%)))
    content((0,y.at(1)+.3),[$V(x)$])

    // 无限深势阱
    on-layer(
      -1,{
        rect((x.at(0)/4,0),(0,0.7 * y.at(1)), stroke: (none), fill: luma(50%))

        rect((0.6 * x.at(1),0),(0.6 * x.at(1)-x.at(0)/4,0.7 * y.at(1)), stroke: (none), fill: luma(50%))
      }
    )

    // 势能曲线
    on-layer(
      1,{
        merge-path({
          line((0,0.7 * y.at(1)),(0,0))
          line((),(0.2 * x.at(1),0))
          bezier-through((),(0.3 * x.at(1),0.2 * y.at(1)),(0.4 * x.at(1),0))
          line((),(0.6 * x.at(1),0))
          line((),(0.6 * x.at(1),0.7 * y.at(1)))
        })
      }
    )
    
  }),
  caption: [带有扰动的无限深方势阱]
)<inf_with_pert>

下面引入一个连续实参量$lambda in [0,1]$：
$ H(lambda) = H^0 + lambda H' $
$lambda$用来表示我们引入微扰的强度，$lambda = 0$对应着无微扰的情况，$lambda = 1$对应着完全微扰的情况。对于某微扰体系，如果$lambda$从0连续变化到1时，$E_n^0$和$ket(n_i^0)$能够分别平滑过渡到$E_n$和$ket(n_i)$，微扰理论用在这个体系上就能够得到比较好的近似结果。

@eqt:non-pert[式] 说明$H^0$的本征函数集${ket(n^0_i)}$是完备的，进一步假设能量的谱不存在简并，将 @eqt:perted[式] 表示为$lambda$的函数形式：
$
  (H^0 + lambda H') ket(n (lambda)) = E_n (lambda) ket(n (lambda))
$<eigen-eqt>
大部分量子力学教材在这一步都会直接将$ket(n(lambda))$和$E_n (lambda)$对$lambda$展开，这里根据 sakurai 的《现代量子力学》
#footnote([Sakurai, Jun John, and Jim Napolitano. 
Modern Quantum Mechanics. 3rd ed, Cambridge University Press, 2021])
给出推导过程。

$lambda$从零递增意味着我们考虑的微扰强度越来越接近实际强度，定义第$n$个能级下微扰项提供的能量：
$
E'_n=E_n-E_n^0
$
将其带入 @eqt:eigen-eqt[式] 中，我们需要求解的方程就变为：
$
(E_n^0 - H^0) ket(n) = (lambda H' - E'_n) ket(n)
$<pert-lam-eqt>
其中$E_n^0$是已知项，通过该式解出微扰能量$E'_n$就能够得到能级的总能量$E_n$。看到这个式子不难想到需要取算符$(E_n^0 - H^0)$的逆算符，但由 @eqt:non-pert[式] 可知该算符作用于$ket(n^0)$上得到的结果是0，要保证算符可逆我们需要构造一个不包含$ket(n^0)$的子空间，引入 *投影算符* ：
$
phi_n equiv 1 - ketbra(n^0) = sum_(k eq.not n) ketbra(k^0)
$
在该子空间上，$(E_n^0 - H^0)$可逆，由于$ketbra(k^0)$同样是$H^0$的本征向量，所以：
$
1/(E_n^0 - H^0) phi_n = sum_(k eq.not n) 1/(E_n^0 - E_k^0) ketbra(k^0)
$<rev-op>
同时 @eqt:pert-lam-eqt[式] 右侧由于没有$ket(n^0)$方向上的分量，即$bra(n^0) (lambda H' - E'_n) ket(n) = 0$，因此我们有：
$
(lambda H' - E'_n) ket(n) = phi_n (lambda H' - E'_n) ket(n)
$

#h(2em) 现在对 @eqt:pert-lam-eqt[式] 左乘 @eqt:rev-op[式] ，得到：
$
ket(n) = 1/(E_n^0 - H^0) phi_n (lambda H' - E'_n) ket(n) + #text(red)[$display(c_n (lambda)ket(n^0))$]
$
