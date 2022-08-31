1. Euler方程

2. SW分裂方法

3. 计算域为：[0,1]\(x方向,Lx) &times; [0,2]\(y方向,Ly)；
 
   网格：160  &times; 321

4. 初始激波位置：y=3*Ly/40；初始界面中心位置：y=Ly/5

5. 界面的形状：余弦

6. 初始激波Mach数：1.5

7. 激波与界面之间的密度和压强设置为1

8. 激波前有来流速度：

   $$ v_d=-( 1-\frac{(\gamma-1)Ma^2+2}{(\gamma+1)Ma^2 } ); $$

   激波后初始速度设置为0

9. Atwood数为 -0.7575

10. x方向为周期性边界条件，y方向为对流边界条件：对y的空间导数为0

11. plt文件为初始场文件。文件中的变量依次是：密度，速度uv，压强，温度，组分浓度

12. 总时间：t_all=4.7

## References:

[^1]: Pooya Movahed & Eric Johnsen, A solution-adaptive method for efficient compressible multifluid simulations, with application to the Richtmyer–Meshkov instability, Journal of Computational Physics, 239: 166–186, 2013.
