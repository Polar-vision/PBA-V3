import open3d as o3d
import numpy as np
import camera
from enum_ import ObjectPointType


def normalization(color):
    a = np.min(color)
    b = np.max(color)
    r = (color[0] - a) / (b - a)
    g = (color[1] - a) / (b - a)
    b = (color[2] - a) / (b - a)
    return [r, g, b]

def draw_camera_point(optype, path):
    width = 500
    height = 500
    n_rows = 2
    n_cols = 2
    margin_x = 20
    margin_y = 20
    col = optype.value % n_cols
    row = optype.value // n_cols
    #
    left_ = col * (width + margin_x) + 20
    top_ = row * (height + margin_y) + 20

    # load a scene point cloud
    if optype == ObjectPointType.ground_truth:
        scene = o3d.io.read_point_cloud(path + 'G-XYZ.ply')
        pose = np.loadtxt(path + 'G-Cam.txt')
    elif optype == ObjectPointType.init:
        scene = o3d.io.read_point_cloud(path + 'XYZ.ply')
        pose = np.loadtxt(path + 'Cam.txt')
    elif optype == ObjectPointType.xyz:
        scene = o3d.io.read_point_cloud(path + 'SBA-LM-Final3D.ply')
        pose = np.loadtxt(path + 'SBA-LM-FinalPose.txt')
    elif optype == ObjectPointType.parallax:
        scene = o3d.io.read_point_cloud(path + 'PBA-LM-Final3D.ply')
        pose = np.loadtxt(path + 'PBA-LM-FinalPose.txt')

    cal = np.loadtxt(path + 'cal.txt')
    scene = scene.scale(1, center=scene.get_center())
    color = [8, 46, 84]  #蓝色
    scene.paint_uniform_color(normalization(color))
    rows, cols = pose.shape

    visualizer = o3d.visualization.VisualizerWithKeyCallback()  # .visualizer1WithEditing()
    visualizer.create_window(window_name=optype.name, width=width, height=height, left=50, top=50)
    render_option = visualizer.get_render_option()  # 设置点云渲染参数
    render_option.point_size = 1  # 设置渲染点的大小（范围并没有具体数据）
    for i in range(rows):
        r = camera.a2r(pose[i, :3])  # 转换为旋转矩阵
        c = np.array([[pose[i, 3]], [pose[i, 4]], [pose[i, 5]]])
        t = -r @ c
        # 创建4x4的单位矩阵
        T = np.identity(4)
        # 将R和t合并到T中
        T[:3, :3] = r  # 将旋转矩阵放入T的前三行前三列
        T[:3, 3:4] = t  # 将位移向量放入T的前三行最后一列
        intrinsics = cal
        wd = int(cal[0, 2] * 2)
        ht = int(cal[1, 2] * 2)


        camera_lines = o3d.geometry.LineSet.create_camera_visualization(view_width_px=wd,
                                                                        view_height_px=ht,
                                                                        intrinsic=intrinsics[:3, :3],
                                                                        extrinsic=T)
        # 缩放视锥
        scale_factor = 0.06  # 缩小为原来的一半
        camera_lines.scale(scale_factor, center=camera_lines.get_center())

        # 获取 LineSet 的 colors 属性
        colors = np.asarray(camera_lines.colors)

        # 视锥的线段在前 8 条（具体数量取决于实现）
        num_frustum_lines = 8  # 通常视锥有 8 条线段
        colors[:num_frustum_lines] = normalization([255, 0, 0])  # 将视锥颜色设置为红色

        # 更新 LineSet 的 colors
        camera_lines.colors = o3d.utility.Vector3dVector(colors)
        visualizer.add_geometry(camera_lines)

    # 可视化坐标轴. The x, y, z axis will be rendered as red, green, and blue arrows respectively.
    axis = o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.2, origin=[0, 0, 0])
    visualizer.add_geometry(axis)
    visualizer.add_geometry(scene)
    return visualizer, scene, axis


