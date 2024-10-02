|日期|机器|simplex sha|实验名称|实验方法|实验目的|结论|
|:---------|:---|:---------|:------|:-------|:------|:---|
|2021/12/30|王梦迪Sachem机器|13882b9bd4cd7fb25f9978a799d73f3706b62bc3|drop-2d|python easy_run.py drop-2d -driver 2 -test 2 -s 64 -lf 100 -d 2|表面张力基础测试：隐式表面张力，2D椭圆泡泡|行为正确，算上初始状态，100帧的时候达到第三次纵向更长的形态|
|2021/12/30|王梦迪Sachem机器|13882b9bd4cd7fb25f9978a799d73f3706b62bc3|drop-3d|python easy_run.py drop-3d -driver 2 -test 2 -s 64 -lf 100 -d 3|表面张力基础测试：隐式表面张力，3D椭圆泡泡|和2D相比，100帧时达到第四次纵向更长的形态，另外似乎“弹力”衰减得比2D也更快|
|2021/12/30|王梦迪Sachem机器|bc181a529eeaed2d8fbd8624e3c1ce8f11acb432|round-2d|python easy_run.py round-2d -driver 2 -test 2 -s 64 -lf 100 -d 2|把drop-2d的水滴形状改为正圆，研究算法表现|正圆水滴会逐步对称地向内收缩|
|2021/12/30|王梦迪Sachem机器|bc181a529eeaed2d8fbd8624e3c1ce8f11acb432|round-3d|python easy_run.py round-3d -driver 2 -test 2 -s 64 -lf 100 -d 3|把drop-3d的水滴形状改为正球，研究算法表现|正球水滴会逐步对称地向内收缩，但收缩的速度比round-2d实验更快一些|
