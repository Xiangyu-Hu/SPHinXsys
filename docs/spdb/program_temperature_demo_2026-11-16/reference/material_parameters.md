# 材料来源及适用范围

Ren F.-C., Chen J., Chen F. (2014). Constitutive modeling of hot deformation behavior of X20Cr13 martensitic stainless steel with strain effect. Trans. Nonferrous Met. Soc. China 24(5),1407–1413. DOI 10.1016/S1003-6326(14)63206-4。

依据前轮用户上传全文核实的式13—18；完整系数见config/ren2014_coefficients.json。多项式按c0到c5升序储存，使用Horner。alpha单位MPa^-1，n无量纲，Q已经是J/mol，lnA自然对数，速率单位s^-1。内部应力若用Pa，alpha除以1e6，不能乘以1e6。

rate=exp(lnA-Q/(R*T))*sinh(alpha*sigma_MPa)^n。
stress_MPa=asinh(exp((ln(rate)+Q/(R*T)-lnA)/n))/alpha。

原试验包含0.001s^-1，但该组多峰数据被排除；实际拟合0.01—10s^-1、真应变0.05—0.75、温度1173—1423K。该模型不是低应力长时蠕变试验拟合，低速用于demo属于外推。论文3.22%是其流动应力数据统计误差，不是本任务孔形变精度。

实际累积非弹性应变从0开始，仅系数评估取至少0.05。不要把内部应变设为0.05。不要用另文413kJ/mol替换Q多项式。温度低于拟合窗口的处理按主任务书，记录假设。

本包不含原论文PDF，服务器实现材料点只需这些已核实系数即可，不必重新下载。论文没有提供本任务完整E(T)、热膨胀、电热物性；配置中的fallback全是演示假设。标量参考程序不是SPH求解器，仅用于核对公式。
