# Hybrid-field Beam-split Pattern Detection-based Channel Estimation for THz XL-MIMO Systems

matlab主函数为nmse_pilot.m, channel_adjust_facor.m以及main_phi_spread.m。注意，matlab极化域码本生成方法相关代码参考[2]，支撑集偏移检测相关代码参考[3]。[2],[3]的代码见http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html。

This is the code for [1]. The main function of Matlab is nmse_pilot.m, channel_adjust_facor.m, and main_phi_spread.m. Note that the relevant codes of Matlab polar-domain codebook generation method refer to [2], and the relevant codes of support set offset detection refer to [3]. [2], [3] code see http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html.

[1] J. Lu, J. Zhang, H. Lei, H. Xiao, Z. Lu, and B. Ai,“Hybrid-field Beam-split Pattern Detection-based Channel Estimation for THz XL-MIMO Systems,” IEEE Trans. Veh. Technol., under review, 2024.
keywords: {XL-MIMO; hybrid-field communications; channel estimation; terahertz communications},

[2] M. Cui and L. Dai, "Channel Estimation for Extremely Large-Scale MIMO: Far-Field or Near-Field?," in IEEE Transactions on Communications, vol. 70, no. 4, pp. 2663-2677, April 2022, doi: 10.1109/TCOMM.2022.3146400. 
[3] M. Cui and L. Dai, “Near-field wideband channel estimation for extremely large-scale MIMO,” Sci. China Inf. Sci., vol. 66, p. 172303, Jun. 2023.

# Abstract
In this letter, we delve into the uplink channelestimation for terahertz extremely large-scale multiple-inputmultiple-output (XL-MIMO) communication systems, focusingon both the polar and the angle domains. Specifically, we explorethe pervasive implementation of hybrid-field communication inXL-MIMO scenarios, addressing challenges such as the beamsplit and power interference effects caused by bandwidth expansion and dense scatterer environments. As a solution, we present ahybrid-field channel estimation approach that utilizes the beamsplit pattern detection method. Furthermore, to enhance the estimation accuracy, we prioritize the near-field channel estimationand minimize its impact on the far-field path components as muchas possible. Simulation results validate that the proposed hybrid-field scheme outperforms benchmark algorithms in normalizedmean square error (NMSE) performance.


# License and Referencing
If you in any way use this code for research that results in publications, please cite our original article listed above ([1]).
