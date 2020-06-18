# 耦合相界面追踪

## VOF

| 第一相 |  $\alpha=1$  |
| :----: | :----------: |
| 第二相 |  $\alpha=0$  |
| 混合相 | $0<\alpha<1$ |

如果认为两相分别是水和空气，分别用$l$和$a$表示，则混合相密度和粘度表示为：


$$
\rho = \alpha\rho_l+(1-\alpha)\rho_a \\
\mu = \alpha\mu_l+(1-\alpha)\mu_a \tag{1}
$$

`界面的平移方程`写成
$$
\frac{\partial{\alpha}}{\partial t}+\frac{\partial{\alpha v_j}}{\partial x_j}-\alpha\frac{\partial v_j}{\partial x_j} = 0 \tag{2}
$$
结合不可压流体连续方程，有
$$
\frac{\partial{\alpha}}{\partial t}+\frac{\partial{\alpha v_j}}{\partial x_j}=0 \tag{3}
$$
使用*VoF*方法的主要挑战是保持相界面的清晰，同时也要保证质量守恒，**OpenFOAM**使用一个额外的逆梯度对流项来满足可压性同时保证有界性和守恒性[^1]。将公式(3)写成
$$
\frac{\partial\alpha}{\partial t}+\frac{\partial{\alpha v_j}}{\partial x_j}+\frac{\partial v_j^c \alpha \beta}{\partial x_j}=0 \\
v_j^c=v^l-v^a \\ \tag{4}
\beta = 1-\alpha
$$
其中的$v_j^c$引入了相界面的可压缩性，$\partial/\partial x_j$保证了守恒，$\alpha \beta$项确保有界性。这个逆梯度项在`alphaEqn.H`中实现，同时通量是被限制的，使用`MULES`算法修正[^2]。

### MULES算法

界面平移方程写成积分形式
$$
\int_{\Omega_i}\frac{\partial \alpha}{\partial t}dV+\int_{\partial \Omega_i}\alpha\mathbf{u}\cdot\mathbf{n}dS = 0 \tag{5}
$$
离散格式方面，时间离散可以是任意格式
$$
\frac{\alpha_i^{n+1}-\alpha_i^n}{\Delta t} = -\frac{1}{|\Omega_i|}\sum_{f\in \partial\Omega_i}(F_u+\lambda F_c)^n \\
F_u = \Phi_f \alpha_{f,upwind} \\
F_c = \Phi_f\alpha_f + \Phi_{rf}\alpha_{rf}(1-\alpha)_{rf}-F_u \tag{6}
$$
$\Phi_f$是体积面通量，$\lambda$在接触面取1，在别处取0。$\Phi_{rf}$的定义：

![山本卓也](https://image.slidesharecdn.com/s-clsvof-140610112254-phpapp02/95/openfoamsclsvoflaplace-7-638.jpg?cb=1402399482)

$C_{\alpha}$对应`cAlpha`，用来界定界面的压缩程度，通常取1.0。$S_f$是面积向量，即方向是面法向，由`owner`指向`neighbor`，模为面积。



Fluent采用的几何重构方法，理论上能够实现较高精度，但是对网格适应性不好。

## Level-Set

水平集方法最先由Osher和Sethian发明[^3]，后被Sussman[^4]用于多相流计算中。主要思想是利用一个带符号的相界面距离函数$\phi$来分辨不同的流体。界面的平移是通过求解带有人工时间的输运方程计算的，优点是相界面清晰，但是并不满足守恒律[^4]。

耦合算法CLSVOF(Coupled Level-Set Volume of Fluid)具体细节[^5]

首先使用*VoF*方法中的$\alpha$初始化水平集函数：
$$
\phi_0=(2\alpha -1 )\Gamma \\ 
\Gamma = 0.75\Delta x \tag{7}
$$
其中$\Delta x$是网格大小。此初始值是一个带符号距离方程，在液相中为正，气相中为负。

然后重新计算此距离函数，使用
$$
\frac{\partial \phi}{\partial \tau}=S(\phi_0)(1-\nabla\phi)\\
\phi(x,0)=\phi(x) \tag{8}
$$
$S$为符号函数，$\tau$是所谓的“人工时间”，选择为$0.1\Delta x$。重新计算后的距离函数应该是光滑的，并在$\phi_{corr} = \epsilon/\Delta \tau$次循环后收敛到$|\Delta \phi|=1$。$\epsilon=1.5\Delta x$是所需要的相界面厚度。

表面张力通过下式计算：
$$
F_{\sigma}=\sigma \kappa (\phi)\delta(\phi)\nabla\phi \tag{9}
$$
其中$\sigma$是表面张力系数，$\kappa(\phi)$为表面曲率，定义为
$$
\kappa = -\nabla\cdot \mathbf{n}_f = -\nabla \cdot (\frac{\nabla \phi_f}{\nabla \phi_f +\delta_N}) \tag{10}
$$
$\delta$为狄拉克函数，用来限制表面张力的影响仅存在相界面处，而在两种流体内取为0，定义为：
$$
\delta(\phi)=\left\{  
             \begin{array}{**lr**}  
             0   &|\phi|>\epsilon\\
			 \frac{1}{2\epsilon}[1+cos(\frac{\pi\phi}{\epsilon})]
			 &|\phi|\leq \epsilon  
             \end{array}  
\right. \tag{11}
$$
然后相分数可以使用一个阶跃函数来计算：
$$
H(\phi)=\left\{  
             \begin{array}{**lr**}  
             0   &\phi<-\epsilon\\
			 \frac{1}{2}[1+\frac{\phi}{\epsilon}+\frac{1}{\pi}sin(\frac{\pi\phi}{\epsilon})] &|\phi|\leq \epsilon  \\
			 1 &\phi>\epsilon
             \end{array}  
\right. \tag{12}
$$
使用$H(\phi)$代替(1)中的$\alpha$计算物性参数。使用通过$\phi$计算的表面张力替代原动量方程中的表面张力。



## 流程图[^6]

![flow](https://user-images.githubusercontent.com/39750942/85037302-c0cccc00-b1b7-11ea-8a24-f1631e1b2544.png)



## isoAdvector

2016年的这篇论文[^7]提出的这种算法已经在OpenFOAM-v1712及以后的版本中实现，经过我对`damBreak`算例的验证，它的液面明显分辨率更高，对小液滴的捕捉也更好，但是收敛较慢。



[^1]: H. G. Weller, “A new approach to vof-based interface capturing methods for incompressible andcompressible flows,” Tech. Rep. TR/HGW/04, 2008.

[^2]: S. T. Zalesak, “Fully multidimensional flux-corrected transport algorithms for fluids,”Journalof Computational Physics, vol. 31, no. 3, pp. 335 – 362, 1979.

[^3]: S. Osher and J. A. Sethian, “Fronts propagating with curvature-dependent speed:  Algorithmsbased on hamilton-jacobi formulations,”Journal of Computational Physics, vol. 79, no. 1, pp. 12– 49, 1988.

[^4]: M.  Sussman,  P.  Smereka,  and  S.  Osher,  “A  level  set  approach  for  computing  solutions  to  in-compressible two-phase flow,”Journal of Computational physics, vol. 114, no. 1, pp. 146–159,1994.

[^5]: A.  Albadawi,  D.  Donoghue,  A.  Robinson,  D.  Murray,  and  Y.  Delaur ́e,  “Influence  of  surfacetension implementation in volume of fluid and coupled volume of fluid with level set methodsfor bubble growth and detachment,”International Journal of Multiphase Flow, vol. 53, no. 0,pp. 11 – 28, 2013

[^6]: CLSVOFsf和CLSVOF的区别在后者在进入主循环前除了更新表面张力，还更新粘度和密度

[^7]: RoenbyJ,BredmoseH,Jasak H.2016Acomputationalmethodforsharp interfaceadvection.R.Soc.opensci.3:160405. http://dx.doi.org/10.1098/rsos.160405
