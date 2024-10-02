# 标准PCG公式1
https://flat2010.github.io/2018/10/26/%E5%85%B1%E8%BD%AD%E6%A2%AF%E5%BA%A6%E6%B3%95%E9%80%9A%E4%BF%97%E8%AE%B2%E4%B9%89/#12-%E9%A2%84%E5%A4%84%E7%90%86
$$
\begin{align}
  r_0&=b-Ax_0\\
  p_0&=M^{-1}r_0\\
  \alpha_k&=\frac{r_k^T M^{-1}r_k}{p_k^T A p_k}\\
  x_{k+1}&=x_k+\alpha_k p_k\\
  r_{k+1}&=r_k-\alpha_k A p_k\\
  \beta_{k+1}&=\frac{r_{k+1}^T M^{-1}r_{k+1}}{r_k^T M^{-1}r_k}\\
  p_{k+1}&=M^{-1}r_{k+1}+\beta_{k+1} p_k
\end{align}
$$

注意这个公式里面频繁出现了$M^{-1}r$，在写算法的时候频繁这么做是不经济的。我们记$M_{-1}r_i=z_i$，则公式转变为：

$$
\begin{align}
  r_0&=b-Ax_0\\
  p_0&=z_0=M^{-1}r_0\\
  \alpha_k&=\frac{r_k^T z_k}{p_k^T A p_k}\\
  x_{k+1}&=x_k+\alpha_k p_k\\
  r_{k+1}&=r_k-\alpha_k A p_k\\
  z_{k+1}&=M^{-1}r_{k+1}\\
  \beta_{k+1}&=\frac{r_{k+1}^T z_{k+1}}{r_k^T z_k}\\
  p_{k+1}&=z_{k+1}+\beta_{k+1} p_k
\end{align}
$$

这里补充一句，如果没有preconditioning，那么直接令$M=I$为单位阵就行了。

# 标准PCG公式2：
https://zhuanlan.zhihu.com/p/98642663
这个和上面那个公式的最大区别是，它取的残差是$Ax-b$而不是$b-Ax$.因此，这个式子里面的$r$和$z$和上面式子里面的$r$和$z$刚好差一个负号。但注意，由于$r$和$z$都差了一个负号，所以$\alpha$反而是相同的，二者的$x$和$p$的含义也一样。两个式子里面$\beta$是反的，我怀疑是打错了。
$$
\begin{align}
  r_0&=Ax_0-b\\
  z_0&=M^{-1}r_0\\
  p_0&=-z_0\\
  \alpha_k&=\frac{r_k^T z_k}{p_k^T A p_k}\\
  x_{k+1}&=x_k+\alpha_k p_k\\
  r_{k+1}&=r_k+\alpha_k A p_k\\
  z_{k+1}&=M^{-1}r_{k+1}\\
  \beta_{k+1}&=-\frac{r_{k+1}^T z_{k+1}}{r_k^T z_k}\\
  p_{k+1}&=-z_{k+1}+\beta_{k+1} p_k
\end{align}
$$

# 标准PCG公式1.1
所以我们还是采用公式1里面的形式。注意到$r^Tz$项也经常出现，计算这一项是一个点积，也比较麻烦，故新增一个变量$\gamma_i=r_i^T z_i$，则进一步改写成

$$
\begin{align}
  r_0&=b-Ax_0\\
  p_0&=z_0=M^{-1}r_0\\
  \gamma_0&=r_0^T z_0\\
  \alpha_k&=\frac{\gamma_k}{p_k^T A p_k}\\
  x_{k+1}&=x_k+\alpha_k p_k\\
  r_{k+1}&=r_k-\alpha_k A p_k\\
  z_{k+1}&=M^{-1}r_{k+1}\\
  \gamma_{k+1}&=r_{k+1}^T z_{k+1}\\
  \beta_{k+1}&=\frac{\gamma_{k+1}}{\gamma_k}\\
  p_{k+1}&=z_{k+1}+\beta_{k+1} p_k
\end{align}
$$


# Eigen使用的PCG公式

https://eigen.tuxfamily.org/dox/ConjugateGradient_8h_source.html

$$
\begin{align}
  r_0&=b-Ax_0\\
  p_0&=z_0=M^{-1}r_0\\
  \gamma_0&=r_0^T p_0 (=r_0^T z_0)\\
  Ap_k&=A\cdot p_k\\
  \alpha_k&=\frac{\gamma_k}{p_k^T A p_k}\\
  x_{k+1}&=x_k+\alpha_k p_k\\
  r_{k+1}&=r_k-\alpha_k A p_k\\
  z_{k+1}&=M^{-1}r_{k+1}\\
  \gamma_{k+1}&=r_{k+1}^T z_{k+1}\\
  \beta_{k+1}&=\frac{\gamma_{k+1}}{\gamma_k}\\
  p_{k+1}&=z_{k+1}+\beta_{k+1} p_k
\end{align}
$$

其终止判据有三个：第一是如果$|b|=0$则设$x=0$退出。设阈值$\eta=(\epsilon|b|)^2$，第二个判据是如果$|r_0|^2<\eta$则退出。注意如果不这么做，当初始猜测很准使得$|r_0|\approx 0$的时候，算$\alpha_0$会得到`nan`。第三个判据是在得到$r_{k+1}$之后，如果$|r_{k+1}|^2<\eta$则退出。

其中$A$叫做`mat`，$b$叫做`rhs`，$r$叫做`residual`，$\eta$叫做`threshold`，$Ap$叫做`tmp`，$\gamma$叫做`absNew`.

# CPX使用的PCG公式

这里采用零向量作为初始猜测。

$$
\begin{align}
  r_0&=b\\
  x_0&=0\\
  z_0&=M^{-1}r_0\\
  p_0&=z_0\\
  \gamma_0&=r_0^T z_0\\
  Ap_k&=A\cdot p_k\\
  \alpha_k&=\frac{\gamma_k}{p_k^T A p_k}\\
  x_{k+1}&=x_k+\alpha_k p_k\\
  r_{k+1}&=r_k-\alpha_k A p_k\\
  z_{k+1}&=M^{-1}r_{k+1}\\
  \gamma_{k+1}&=r_{k+1}^T z_{k+1}\\
  \beta_{k+1}&=\frac{\gamma_{k+1}}{\gamma_k}\\
  p_{k+1}&=z_{k+1}+\beta_{k+1} p_k
\end{align}
$$
