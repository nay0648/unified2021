# A UNIFIED SPEECH ENHANCEMENT FRONT-END FOR ONLINE DEREVERBERATION, ACOUSTIC ECHO CANCELLATION, AND SOURCE SEPARATION

Yueyue Na, Ziteng Wang, Zhang Liu, Yun Li, Gang Qiao, Biao Tian, Qiang Fu

Machine Intelligence Technology, Alibaba Group

{yueyue.nyy, ziteng.wzt, yinan.lz, yl.yy, songjiang.qg, tianbiao.tb, fq153277}@alibaba-inc.com

# ABSTRACT
Dereverberation (DR), acoustic echo cancellation (AEC), and blind source separation (BSS) are the three most important submodules in speech enhancement front-end. In traditional systems, the three submodules work independently in a sequential manner, each submodule has its own signal model, objective function, and optimization policy. Although this architecture has high flexibility, the speech enhancement performance is restricted, since each submodule's optimum cannot guarantee the entire system's global optimum. In this paper, a unified signal model is derived to combine DR, AEC, and BSS together, and the online auxiliary-function based independent component/vector analysis (Aux-ICA/IVA) technique is used to solve the problem. The proposed approach has unified objective function and optimization policy, the performance improvement is verified by simulated experiments.

# Links

[【技术揭秘】解决“鸡尾酒会问题”的利器-基于盲源分离的前端信号处理框架](https://mp.weixin.qq.com/s/wjVDgwHpLYM3vZWL5AJvDg)
