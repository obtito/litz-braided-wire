"""赛题统一工况（表 1）——所有仿真共用，不得修改。"""
import math

SIGMA = 5.8e7            # S/m 无氧铜电导率
RHO = 1.0 / SIGMA        # Ω·m 电阻率
MU0 = 4e-7 * math.pi     # H/m 真空磁导率
MU_R = 1.0
F0 = 200e3               # Hz 工作频率
OMEGA = 2 * math.pi * F0
I_RMS = 20.0             # A 目标载流（有效值）
J_DC_MAX = 4e6           # A/m^2 直流电流密度上限（温升约束）
LENGTH = 1.0             # m 评估长度
A_CU_MIN = I_RMS / J_DC_MAX   # m^2 最小铜截面（= 5 mm²）


def skin_depth(f: float) -> float:
    """铜的趋肤深度 δ(f)，式 (A1)。"""
    return math.sqrt(RHO / (math.pi * f * MU0 * MU_R))
