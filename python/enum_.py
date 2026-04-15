from enum import Enum


# 定义一个枚举类
class BA(Enum):
    pba_lm = 1
    sba_lm = 2
    pba_gn = 3
    sba_gn = 4
    init = 5

class ObjectPointType(Enum):
    #零锚点
    xyz=0
    xy_inverse_z=1
    depth=2
    inverse_depth=3
    #单锚点
    archored_xyz = 4
    archored_xy_inverse_z = 5
    archored_depth=6
    archored_inverse_depth=7
    #双锚点
    parallax=8

    init=9
    ground_truth=10

class Rotation3DType(Enum):
    euler_angle=0
    axis_angle=1
    quaternion=2

class ImagePointType(Enum):
    uv=0
    light_cone=1
    gt=2
    init=3

class ParameterType(Enum):
    rotation_translation_landmark = 0,
    rotation_landmark = 1,
    translation_landmark = 2,
    rotation_translation = 3,
    landmark = 4,
    rotation = 5,
    translation = 6

class ManifoldType(Enum):
    lie = 0,
    quaternion_manifold = 1,
    euclidean_manifold = 2,
    sphere_manifold = 3,
    line_manifold = 4,
    none = 5


