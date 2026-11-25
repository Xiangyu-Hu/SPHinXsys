# 给 Cursor 的启动说明

请先读CURSOR_DEVELOPMENT_PLAN.md，在用户当前SPDB分支新增独立case，直接编写、编译、运行和验证。复用旧分支已有材料、接触和准静态代码，不重新clone或切分支。

本轮：双槽试样对接，均匀程序温度20→300→500→700→950℃各段5分钟，最后暂按950℃保温10分钟，总1800秒；5MPa名义恒力。算热膨胀和受压非弹性变形，不算传热、电场和石墨。

geometry/lower_m.stl及upper_m.stl已按米转换并装配，不能再次缩放。保留geometry/specimen_original_mm.stl作为原始输入。tools/prepare_inputs.py可重现转换；tools/material_reference.py可独立核对材料点。

完整任务、系数、配置及几何均在包内。配置为待实现接口建议，不是已运行求解器。最后10分钟保温、材料fallback与低温非弹性门函数为明确演示假设。完成后提供真实计算结果，不仅回复计划。
