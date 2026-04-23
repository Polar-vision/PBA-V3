import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argv_path

file_path = argv_path.paths[0]
directory = os.path.dirname(file_path)

# 数据
pba_data = np.loadtxt(directory + '/pba_sigma_spectral_log.txt')
sba_data = np.loadtxt(directory + '/sba_sigma_spectral_log.txt')

pba_data1 = np.loadtxt(directory + '/pba_cond_log.txt')
sba_data1 = np.loadtxt(directory + '/sba_cond_log.txt')

pba_error  = pba_data1[:, 0]
pba_rho  = pba_data1[:, 1]
pba_mu   = pba_data1[:, 2]
pba_smax = pba_data1[:, 3]
pba_smin = pba_data1[:, 4]
pba_cond = pba_data1[:, 5]
pba_binf = pba_data1[:, 6]
pba_rsc = pba_data1[:, 7]
pba_rmc = pba_data1[:, 8]
pba_b_lips = pba_data1[:, 9]
pba_b_dir_valid = pba_data1[:, 10]
pba_nLmIt = pba_data1[:, 11]

sba_error = sba_data1[:, 0]
sba_rho  = sba_data1[:, 1]
sba_mu   = sba_data1[:, 2]
sba_smax = sba_data1[:, 3]
sba_smin = sba_data1[:, 4]
sba_cond = sba_data1[:, 5]
sba_binf = sba_data1[:, 6]
sba_rsc = sba_data1[:, 7]
sba_rmc = sba_data1[:, 8]
sba_b_lips = sba_data1[:, 9]
sba_b_dir_valid = sba_data1[:, 10]
sba_nLmIt = sba_data1[:, 11]

# 应用对数变换
pba_data_log = np.log10(pba_data + 1e-10)
sba_data_log = np.log10(sba_data + 1e-10)

# 获取迭代次数
n_iterations_pba = pba_data.shape[0]-1
n_iterations_sba = sba_data.shape[0]-1

print(f"PBA iterations: {n_iterations_pba}")
print(f"SBA iterations: {n_iterations_sba}")

# 统一最终迭代（两者取最小）
final_iter = min(n_iterations_pba, n_iterations_sba)
mid_iter = int((final_iter+1) / 2)

# 计算条件数（最大奇异值/最小奇异值）
def compute_condition_number(sigma_data):
    """计算条件数（最大奇异值/最小奇异值）"""
    condition_numbers = []
    for sigma in sigma_data:
        if len(sigma) > 0 and np.min(sigma) > 0:
            cond_num = sigma[0] / sigma[-1]  # 最大/最小
            condition_numbers.append(cond_num)
        else:
            condition_numbers.append(np.inf)
    return np.array(condition_numbers)

# 计算条件数
pba_condition = compute_condition_number(pba_data)
sba_condition = compute_condition_number(sba_data)

# 三个关键迭代阶段
stages = ['Initial', 'Middle', 'Final']
indices = [0, mid_iter, final_iter]

print(f"\n关键迭代:")
print(f"  Initial: iter {indices[0]}")
print(f"  Middle: iter {indices[1]}")
print(f"  Final: iter {indices[2]}")

# 风格设置
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica"],
    "font.size": 14,
    "axes.linewidth": 0.8,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "legend.frameon": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "figure.constrained_layout.use": True,
})

# 创建三个子图（1x3布局）
fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=300)

# 颜色定义
pba_color = '#2c7bb6'  # 蓝色
sba_color = '#d7191c'  # 红色

# 为每个子图设置不同的线型风格
stage_styles = {
    'Initial': {'linewidth': 1.5, 'alpha': 0.7, 'linestyle': '-', 'marker': 'o', 'markersize': 3, 'markevery': 0.1},
    'Middle': {'linewidth': 1.8, 'alpha': 0.85, 'linestyle': '-', 'marker': 's', 'markersize': 3, 'markevery': 0.1},
    'Final': {'linewidth': 2.0, 'alpha': 1.0, 'linestyle': '-', 'marker': '^', 'markersize': 4, 'markevery': 0.1}
}

# 为每个子图绘制对比
for idx, (stage, ax) in enumerate(zip(stages, axes)):
    pba_iter_idx = indices[idx] if indices[idx] < n_iterations_pba else n_iterations_pba
    sba_iter_idx = indices[idx] if indices[idx] < n_iterations_sba else n_iterations_sba
    
    # 获取对数变换后的数据
    pba_sigma_log = pba_data_log[pba_iter_idx, :]
    sba_sigma_log = sba_data_log[sba_iter_idx, :]
    pba_rank = np.arange(1, len(pba_sigma_log) + 1)
    sba_rank = np.arange(1, len(sba_sigma_log) + 1)
    
    # 获取原始数据用于条件数
    pba_sigma_orig = pba_data[pba_iter_idx, :]
    sba_sigma_orig = sba_data[sba_iter_idx, :]
    
    # 绘制PBA曲线
    style = stage_styles[stage]
    ax.plot(pba_rank, pba_sigma_log,
            color=pba_color,
            linewidth=style['linewidth'],
            linestyle=style['linestyle'],
            alpha=style['alpha'],
            marker=style['marker'],
            markersize=style['markersize'],
            markevery=max(1, len(pba_rank)//15),
            label='PBA')
    
    # 绘制SBA曲线
    ax.plot(sba_rank, sba_sigma_log,
            color=sba_color,
            linewidth=style['linewidth'],
            linestyle='--',
            alpha=style['alpha'],
            marker=style['marker'],
            markersize=style['markersize'],
            markevery=max(1, len(sba_rank)//15),
            label='SBA')
    
    # 创建右y轴用于条件数
    ax2 = ax.twinx()
    
    # 计算当前迭代的条件数
    pba_cond = pba_condition[pba_iter_idx]
    sba_cond = sba_condition[sba_iter_idx]
    
    # 在右y轴上绘制条件数（使用条形图更直观）
    x_pos = [len(pba_rank) * 0.2, len(pba_rank) * 0.8]
    cond_values = [pba_cond, sba_cond]
    colors = [pba_color, sba_color]
    labels = ['PBA', 'SBA']
    
    # 绘制条件数条形图
    bars = ax2.bar(x_pos, cond_values, width=len(pba_rank)*0.1, 
                   color=colors, alpha=0.6, edgecolor='black', linewidth=0.5)
    
    # 添加数值标签
    for bar, val, label in zip(bars, cond_values, labels):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{label}\n{val:.1e}',
                ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    # 设置右y轴
    ax2.set_ylabel('Condition Number', fontsize=14, color='black', fontweight='bold')
    # ax2.tick_params(axis='y', labelsize=0)
    # ax2.set_yscale('log')  # 条件数用对数坐标
    ax2.tick_params(axis='y', labelsize=0, length=0)  # 隐藏刻度标签和刻度线

    
    # 设置左y轴
    ax.set_ylabel('log10(Singular Value)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Rank', fontsize=14, fontweight='bold')
    
    # 设置子图标题
    ax.set_title(f'{stage} Stage\n(iter {pba_iter_idx}/{n_iterations_pba} vs {sba_iter_idx}/{n_iterations_sba})', 
                fontsize=14, fontweight='bold', pad=15)
    
    # 添加网格
    ax.grid(True, which='both', linestyle=':', linewidth=0.3, alpha=0.3)
    
    # 美化边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_linewidth(0.8)
    
    # 设置y轴范围（左y轴）
    all_data = np.concatenate([pba_sigma_log, sba_sigma_log])
    y_min = np.min(all_data) - 0.5
    y_max = np.max(all_data) + 0.5
    ax.set_ylim(bottom=y_min, top=y_max)
    
    # 设置x轴范围
    ax.set_xlim(left=-len(sba_rank)/10, right=max(len(pba_rank), len(sba_rank)) * 1.05)
    
    # 添加图例（左y轴）
    ax.legend(loc='best', fontsize=12, frameon=True, fancybox=True, 
             shadow=True, framealpha=0.8)
    
    # # 添加条件数信息文本框
    # textstr = f'PBA Cond#: {pba_cond:.2e}\nSBA Cond#: {sba_cond:.2e}'
    # ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=7,
    #         verticalalignment='bottom', horizontalalignment='right',
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8, pad=2))

# 调整布局
plt.tight_layout()

# 保存图形
output_pdf = directory + '/singular_spectrum_log_comparison.pdf'
output_png = directory + '/singular_spectrum_log_comparison.png'

plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')

print(f"\nFigures saved to:")
print(f"  {output_pdf}")
print(f"  {output_png}")

# 打印统计信息
print(f"\n条件数统计:")
print(f"  PBA - Initial: {pba_condition[0]:.2e}, Middle: {pba_condition[mid_iter]:.2e}, Final: {pba_condition[-1]:.2e}")
print(f"  SBA - Initial: {sba_condition[0]:.2e}, Middle: {sba_condition[mid_iter]:.2e}, Final: {sba_condition[-1]:.2e}")

# plt.show()

# 可选：创建第四个对比图（最终迭代详细对比）
def create_final_comparison():
    """创建最终迭代的详细对比图"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    # 获取最后一次迭代的数据
    pba_final = pba_data_log[-1, :]
    sba_final = sba_data_log[-1, :]
    pba_rank = np.arange(1, len(pba_final) + 1)
    sba_rank = np.arange(1, len(sba_final) + 1)
    
    # 绘制曲线
    ax.plot(pba_rank, pba_final, color=pba_color, linewidth=2.0, 
           marker='o', markersize=4, markevery=0.05, label=f'PBA (iter {n_iterations_pba})')
    ax.plot(sba_rank, sba_final, color=sba_color, linewidth=2.0, 
           linestyle='--', marker='s', markersize=4, markevery=0.05, 
           label=f'SBA (iter {n_iterations_sba})')
    
    # 创建右y轴
    ax2 = ax.twinx()
    
    # 计算最终条件数
    pba_cond_final = pba_condition[-1]
    sba_cond_final = sba_condition[-1]
    
    # 添加条件数标注
    ax2.axhline(y=pba_cond_final, color=pba_color, linestyle=':', linewidth=1.5, 
               alpha=0.7, label=f'PBA Cond#: {pba_cond_final:.2e}')
    ax2.axhline(y=sba_cond_final, color=sba_color, linestyle=':', linewidth=1.5, 
               alpha=0.7, label=f'SBA Cond#: {sba_cond_final:.2e}')
    
    # 设置坐标轴
    ax.set_xlabel('Rank', fontsize=14, fontweight='bold')
    ax.set_ylabel('log10(Singular Value)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Condition Number', fontsize=14, fontweight='bold')
    # ax2.set_yscale('log')
    
    # 设置标题
    ax.set_title(f'Final Iteration Comparison\n(PBA: {n_iterations_pba} iters, SBA: {n_iterations_sba} iters)', 
                fontsize=14, fontweight='bold')
    
    # 网格和美化
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    
    # 合并图例
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='best', fontsize=14, 
             frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    
    output_pdf = directory + '/final_iteration_comparison.pdf'
    output_png = directory + '/final_iteration_comparison.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"Final comparison saved to {output_pdf}")

# 创建最终对比图
create_final_comparison()

# 可选：创建条件数演化图
def plot_condition_evolution():
    """绘制条件数随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba+1)
    iterations_sba = np.arange(0, n_iterations_sba+1)
    
    ax.plot(iterations_pba, pba_condition, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_condition, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Condition Number", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Condition Number Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_condition[0], pba_condition[mid_iter], pba_condition[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_condition[0], sba_condition[mid_iter], sba_condition[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/condition_number_evolution.pdf'
    output_png = directory + '/condition_number_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"Condition number evolution saved to {output_pdf}")

# 可选：创建阻尼因子演化图
def plot_damping_evolution():
    """绘制阻尼因子随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_mu, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_mu, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Damping Factor", fontsize=14, fontweight='bold')
#     ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Damping Factor Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_mu[0], pba_mu[mid_iter], pba_mu[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_mu[0], sba_mu[mid_iter], sba_mu[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/damping_factor_evolution.pdf'
    output_png = directory + '/damping_factor_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"Damping factor evolution saved to {output_pdf}")

def plot_damping_evolution_log():
    """绘制阻尼因子随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    # 处理零值：将0替换为机器精度或极小值
    eps = 1e-16  # 或使用 np.finfo(float).eps
    pba_mu_safe = np.array(pba_mu, dtype=float)
    sba_mu_safe = np.array(sba_mu, dtype=float)
    pba_mu_safe[pba_mu_safe <= 0] = eps
    sba_mu_safe[sba_mu_safe <= 0] = eps
    
    ax.plot(iterations_pba, pba_mu_safe, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_mu_safe, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Damping Factor (log scale)", fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Damping Factor Evolution (log scale)', fontsize=14, fontweight='bold')
    
    # 标注关键点（同样使用安全版本）
    # 注意：如果原值为0，标注时显示极小值，位置会掉到图底部
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_mu_safe[0], pba_mu_safe[mid_iter], pba_mu_safe[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_mu_safe[0], sba_mu_safe[mid_iter], sba_mu_safe[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/damping_factor_evolution_log_scale.pdf'
    output_png = directory + '/damping_factor_evolution_log_scale.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    print(f"Damping factor evolution log scale saved to {output_pdf}")

def plot_rho_evolution():
    """绘制rho因子随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_rho, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_rho, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Gain Ratio", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Gain Ratio Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_rho[0], pba_rho[mid_iter], pba_rho[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_rho[0], sba_rho[mid_iter], sba_rho[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/rho_evolution.pdf'
    output_png = directory + '/rho_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"rho evolution saved to {output_pdf}")

# 可选：创建最小奇异值演化图
def plot_sigma_min_evolution():
    """绘制最小特征值随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_smin, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_smin, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Minimum Singular Value", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Minimum Singular Value Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_smin[0], pba_smin[mid_iter], pba_smin[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_smin[0], sba_smin[mid_iter], sba_smin[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/smin_evolution.pdf'
    output_png = directory + '/smin_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"smin evolution saved to {output_pdf}")

def plot_sigma_min_evolution_log_scale():
    """绘制最小特征值随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_smin, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_smin, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Minimum Singular Value (log scale)", fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Minimum Singular Value Evolution (log scale)', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_smin[0], pba_smin[mid_iter], pba_smin[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_smin[0], sba_smin[mid_iter], sba_smin[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/smin_evolution_log_scale.pdf'
    output_png = directory + '/smin_evolution_log_scale.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"smin evolution log scale saved to {output_pdf}")

# 可选：创建最大奇异值演化图
def plot_sigma_max_evolution():
    """绘制最大特征值随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_smax, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_smax, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Maximum Singular Value", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Maximum Singular Value Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_smax[0], pba_smax[mid_iter], pba_smax[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_smax[0], sba_smax[mid_iter], sba_smax[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/smax_evolution.pdf'
    output_png = directory + '/smax_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"smax evolution saved to {output_pdf}")


# 可选：创建收敛演化图
def plot_error_evolution():
    """绘制目标函数随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_error, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_error, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("log(MSE)", fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('MSE Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_error[0], pba_error[mid_iter], pba_error[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_error[0], sba_error[mid_iter], sba_error[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/MSE_evolution.pdf'
    output_png = directory + '/MSE_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"MSE evolution saved to {output_pdf}")

# 可选：创建梯度最大分量演化图
def plot_binf_evolution():
    """绘制梯度最大分量随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_binf, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_binf, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Maximum Gradient Component", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Maximum Gradient Component Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_binf[0], pba_binf[mid_iter], pba_binf[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_binf[0], sba_binf[mid_iter], sba_binf[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/binf_evolution.pdf'
    output_png = directory + '/binf_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"binf evolution saved to {output_pdf}")

def plot_binf_evolution_log_scale():
    """绘制梯度最大分量随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_binf, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_binf, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Maximum Gradient Component (log scale)", fontsize=14, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Maximum Gradient Component Evolution', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_binf[0], pba_binf[mid_iter], pba_binf[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_binf[0], sba_binf[mid_iter], sba_binf[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/binf_evolution_log_scale.pdf'
    output_png = directory + '/binf_evolution_log_scale.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"binf_log_scale evolution saved to {output_pdf}")

# 可选：创建状态变量的相对变化量演化图
def plot_rsc_evolution():
    """绘制状态变量的相对变化量随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba+1)
    iterations_sba = np.arange(0, n_iterations_sba+1)
    
    ax.plot(iterations_pba, pba_rsc, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_rsc, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Relative Change of State Vector", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Evolution of the Relative Change of State Vector', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_rsc[0], pba_rsc[mid_iter], pba_rsc[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_rsc[0], sba_rsc[mid_iter], sba_rsc[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/rsc_evolution.pdf'
    output_png = directory + '/rsc_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"rsc evolution saved to {output_pdf}")

# 可选：创建重投影误差相对变化量演化图
def plot_rmc_evolution():
    """绘制重投影误差相对变化量随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(1, n_iterations_pba + 1)
    iterations_sba = np.arange(1, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_rmc[1:], color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_rmc[1:], color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Relative Change of Projection Error", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Evolution of the Relative Change of Projection Error', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([mid_iter, n_iterations_pba], 
              [pba_rmc[mid_iter], pba_rmc[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([mid_iter, n_iterations_sba], 
              [sba_rmc[mid_iter], sba_rmc[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/rmc_evolution.pdf'
    output_png = directory + '/rmc_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"rmc evolution saved to {output_pdf}")


def plot_b_lips_evolution():
    """绘制梯度的Lipschitz连续性随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(1, n_iterations_pba + 1)
    iterations_sba = np.arange(1, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_b_lips[1:], color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_b_lips[1:], color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Gradient Lipschitz", fontsize=14, fontweight='bold')
#     ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Gradient Lipschitz', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([mid_iter, n_iterations_pba], 
              [pba_b_lips[mid_iter], pba_rmc[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([mid_iter, n_iterations_sba], 
              [sba_b_lips[mid_iter], sba_rmc[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/b_lips_evolution.pdf'
    output_png = directory + '/b_lips_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"b_lips evolution saved to {output_pdf}")

def plot_b_lips_evolution_log_scale():
    """绘制梯度的Lipschitz连续性随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(1, n_iterations_pba + 1)
    iterations_sba = np.arange(1, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_b_lips[1:], color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_b_lips[1:], color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    ax.set_yscale('log')
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Gradient Lipschitz (log scale)", fontsize=14, fontweight='bold')
    
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Gradient Lipschitz (log scale)', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([mid_iter, n_iterations_pba], 
              [pba_b_lips[mid_iter], pba_b_lips[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([mid_iter, n_iterations_sba], 
              [sba_b_lips[mid_iter], sba_b_lips[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/b_lips_evolution_log_scale.pdf'
    output_png = directory + '/b_lips_evolution_log_scale.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"b_lips evolution log scale saved to {output_pdf}")

def plot_b_dir_valid_evolution():
    """绘制梯度方向的有效性随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba+1)
    iterations_sba = np.arange(0, n_iterations_sba+1)
    
    ax.plot(iterations_pba, pba_b_dir_valid, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_b_dir_valid, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    ax.set_ylabel("Efficacy of the Gradient Direction", fontsize=14, fontweight='bold')
    # ax.set_yscale('log')
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Evolution of the Efficacy of the Gradient Direction', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_b_dir_valid[0], pba_b_dir_valid[mid_iter], pba_b_dir_valid[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_b_dir_valid[0], sba_b_dir_valid[mid_iter], sba_b_dir_valid[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/b_dir_valid_evolution.pdf'
    output_png = directory + '/b_dir_valid_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"b_dir_valid evolution saved to {output_pdf}")

def plot_nLimIte_evolution():
    """绘制当前阻尼参数（或信任区域半径）下，算法尝试更新参数的总次数随迭代的演化"""
    fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
    
    iterations_pba = np.arange(0, n_iterations_pba + 1)
    iterations_sba = np.arange(0, n_iterations_sba + 1)
    
    ax.plot(iterations_pba, pba_nLmIt, color=pba_color, linewidth=2.0, 
            marker='o', markersize=4, markevery=0.1, label='PBA')
    ax.plot(iterations_sba, sba_nLmIt, color=sba_color, linewidth=2.0, 
            linestyle='--', marker='s', markersize=4, markevery=0.1, label='SBA')
    
    ax.set_xlabel("Iteration", fontsize=14, fontweight='bold')
    # 改进1：使用更规范的Y轴名称
    ax.set_ylabel("Number of update attempts", fontsize=14, fontweight='bold')
    
    # 改进2：强制Y轴显示整数刻度
    from matplotlib.ticker import MaxNLocator
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
    # 改进3：设置Y轴范围从0开始（如果数据最小值大于0）
    y_min = min(np.min(pba_nLmIt), np.min(sba_nLmIt))
    y_max = max(np.max(pba_nLmIt), np.max(sba_nLmIt))
    ax.set_ylim(bottom=max(0, y_min - 1), top=y_max + 1)
    
    # ax.set_yscale('log')  # 整数数据通常不推荐log坐标
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    ax.legend(frameon=False, loc='best', fontsize=14)
    ax.set_title('Number of Parameter Update Attempts per Iteration', fontsize=14, fontweight='bold')
    
    # 标注关键点
    ax.scatter([0, mid_iter, n_iterations_pba], 
              [pba_nLmIt[0], pba_nLmIt[mid_iter], pba_nLmIt[-1]], 
              color=pba_color, s=50, zorder=5, edgecolors='black')
    ax.scatter([0, mid_iter, n_iterations_sba], 
              [sba_nLmIt[0], sba_nLmIt[mid_iter], sba_nLmIt[-1]], 
              color=sba_color, s=50, zorder=5, edgecolors='black')
    
    plt.tight_layout()
    
    output_pdf = directory + '/nLmIt_evolution.pdf'
    output_png = directory + '/nLmIt_evolution.png'
    plt.savefig(output_pdf, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(output_png, format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    print(f"nLmIt evolution saved to {output_pdf}")

# plot_condition_evolution()
# plot_damping_evolution()
# plot_damping_evolution_log()
# plot_rho_evolution()
# plot_sigma_min_evolution()
# plot_sigma_min_evolution_log_scale()
# plot_sigma_max_evolution()

# plot_error_evolution()
# plot_binf_evolution()
# plot_binf_evolution_log_scale()
# plot_rsc_evolution()
# plot_rmc_evolution()
# plot_b_lips_evolution()
plot_b_lips_evolution_log_scale()
# plot_b_dir_valid_evolution()
# plot_nLimIte_evolution()