# Doc205：最小圆柱闭合环流

任务书：`docs/ophelie/205_OPHELIE_MINIMAL_CYLINDER_CLOSED_CIRCULATION_CURSOR_PLAN_2026-11-01.md`

## 目的

用与 C1 同源的水动力装配，在法国尺寸圆柱上验证**闭合回流**是否可出现。  
只改顶 BC（A 全无滑移 / B 顶滑移）。不做文献双涡、电磁、热方程。

## 运行矩阵

| ID | 顶 | beta | 命令骨架 |
|---|---|---:|---|
| A0 | no-slip | 0 | `--case=A0 --end-time=0.2` 再 5s |
| A1 | no-slip | 1e-5 | `--case=A1 --end-time=5` 再 20s |
| B0 | lid free-slip（开杯粘性 + AcousticWallFreeSlip） | 0 | `--case=B0 --dp=0.015 --end-time=0.2` |
| B1 | 同上 | 1e-5 | B0 通过后 `--case=B1 --end-time=5` 再 20s |

B 装配：`wall_closed`=声学/dens 全封闭；`wall_noslip`=侧+底粘性（无顶粒子）；顶盖粒子打 `AcousticWallFreeSlip=1`。

输出目录：`docs/ophelie/reference/doc205_runs/<run_tag>/`

## 禁止

复制法国主程序；quiet-reset-mass；用 U_rms/好看 Glyph 代替截面质量流量验收。
