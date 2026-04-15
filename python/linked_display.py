import os
from draw_camera_point import draw_camera_point
from enum_ import ObjectPointType
from argv_path import paths
vis_ctr = "p"


# 主函数
def linked_display(path1,path2):
    vis_gt, scene_gt, axis_gt = draw_camera_point(ObjectPointType.ground_truth, path2)
    vis_init, scene_init, axis_init = draw_camera_point(ObjectPointType.init, path1)
    vis_xyz, scene_xyz, axis_xyz = draw_camera_point(ObjectPointType.xyz, path1)
    vis_parallax, scene_parallax, axis_parallax = draw_camera_point(ObjectPointType.parallax, path1)

    # 获取初始视角控制对象
    ctr_gt = vis_gt.get_view_control()
    ctr_init = vis_init.get_view_control()
    ctr_xyz = vis_xyz.get_view_control()
    ctr_parallax = vis_parallax.get_view_control()
    # 设置初始视角
    ctr_gt.set_zoom(10)
    ctr_init.set_zoom(10)
    ctr_xyz.set_zoom(10)
    ctr_parallax.set_zoom(10)

    # 主循环，实现联动效果
    while True:
        param_gt = ctr_gt.convert_to_pinhole_camera_parameters()
        ctr_init.convert_from_pinhole_camera_parameters(param_gt)
        ctr_xyz.convert_from_pinhole_camera_parameters(param_gt)
        ctr_parallax.convert_from_pinhole_camera_parameters(param_gt)

        # 更新3个窗口
        vis_gt.update_geometry(scene_gt)
        vis_gt.update_geometry(axis_gt)
        vis_init.update_geometry(scene_init)
        vis_init.update_geometry(axis_init)
        vis_xyz.update_geometry(scene_xyz)
        vis_xyz.update_geometry(axis_xyz)
        vis_parallax.update_geometry(scene_parallax)
        vis_parallax.update_geometry(axis_parallax)

        vis_gt.poll_events()
        vis_gt.update_renderer()

        vis_init.poll_events()
        vis_init.update_renderer()

        vis_xyz.poll_events()
        vis_xyz.update_renderer()

        vis_parallax.poll_events()
        vis_parallax.update_renderer()

        vis_gt.capture_screen_image(path + 'G-XYZ.png', do_render=True)
        vis_init.capture_screen_image(path + 'Init-XYZ.png', do_render=True)
        vis_xyz.capture_screen_image(path + 'xyz.png', do_render=True)
        vis_parallax.capture_screen_image(path + 'parallax.png', do_render=True)

        if not vis_gt.poll_events() or not vis_init.poll_events() or not vis_xyz.poll_events() or not vis_parallax.poll_events():
            break

    vis_gt.destroy_window()
    vis_init.destroy_window()
    vis_xyz.destroy_window()
    vis_parallax.destroy_window()

str1 = '/'
str2 = '/Ground Truth/'
for path in paths:
    parent_dir = os.path.dirname(path)
    new_path1 = "%s%s" % (parent_dir, str1)
    base_dir = os.path.dirname(parent_dir)
    new_path2 = "%s%s" % (base_dir, str2)
    linked_display(new_path1,new_path2)