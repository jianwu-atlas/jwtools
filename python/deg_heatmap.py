#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
deg_heatmap.py
────────────────────────────────────────────────────────────────
通用的「Up/Down 双面板 DEG presence heatmap」绘图工具，整合自
TE / ERV / gene 三套独立脚本的公共逻辑，抽象为单一可复用函数
plot_combined_deg_heatmap()。

核心特性（与三套原始脚本完全一致，逻辑未做任何简化）：
  1. 三层分类：upper(shared≥2 celltype) / middle(shared≥2 age-group
     comparison within same celltype) / lower(unique)
  2. Upper tier 两级排序：breadth(涉及celltype数)降序 + 组内
     hamming距离 + average linkage + optimal_ordering 层次聚类
  3. 混合坐标变换 (transAxes + transData)，避免热图外部标注
     (左侧DEG计数、顶部celltype分组名)意外撑大/错位 axes 边框
  4. Up/Down 面板高度按各自真实DEG总数比例分配
     (GridSpec height_ratios)，并施加 MIN_PANEL_FRAC 保底比例
     (仅显示层面保护，不影响任何实际数据/CSV导出)
  5. PPT参考图风格：Panel A(Up,红色) 在上，Panel B(Down,蓝色) 在下，
     各自独立表头(cell type + O/Y-M/Y-O/M)
  6. Illustrator 兼容性设置 (pdf.fonttype=42)，PNG+PDF双保存，
     所有中间数据(三层分类结果、breadth分组、面板高度比例)均导出CSV

用法示例见文件末尾 `if __name__ == '__main__':` 部分。
────────────────────────────────────────────────────────────────
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# ============================================================
# 默认配色主题（up=红色系, down=蓝色系）——可通过参数覆盖
# ============================================================
DEFAULT_THEME_UP = dict(
    color_upper='#b81111', color_middle='#7a0d0d',
    cmap_low=(0.88, 0.88, 0.88), cmap_high=(0.72, 0.07, 0.07),
    title='Upregulated', suffix='up',
)
DEFAULT_THEME_DOWN = dict(
    color_upper='#14407A', color_middle='#0d2a52',
    cmap_low=(0.88, 0.88, 0.88), cmap_high=(0.08, 0.30, 0.65),
    title='Downregulated', suffix='down',
)

DEFAULT_AGE_LABEL_MAP = {'OvsY': 'O/Y', 'MvsY': 'M/Y', 'OvsM': 'O/M'}


# ============================================================
# 内部辅助函数（均为纯函数，不依赖全局状态，便于单元测试）
# ============================================================
def _parse_col_factory(suffix):
    """根据列名后缀(up/down)生成解析函数: 'OvsY_O_up_Mono' -> ('Mono', 'OvsY')"""
    pattern = re.compile(rf'(OvsY|MvsY|OvsM)_\w+_{suffix}_(.+)')
    def parse_col(col):
        m = pattern.match(col)
        if m:
            return m.group(2), m.group(1)   # celltype, age
        return col, 'ZZZ'
    return parse_col


def _expected_conditions(suffix, celltype_priority):
    """按 celltype_priority 顺序 × 3个age比较，生成期望的完整列名列表"""
    age_map = {'OvsY': 'O', 'MvsY': 'M', 'OvsM': 'O'}
    expected = []
    for ct in celltype_priority:
        for age_str in ['OvsY', 'MvsY', 'OvsM']:
            expected.append(f'{age_str}_{age_map[age_str]}_{suffix}_{ct}')
    return expected


def _first_one_col_index(row, cols):
    """返回该行第一个值为1的列的索引位置，用于 lower/unique tier 的阶梯排序"""
    for i, c in enumerate(cols):
        if row[c] == 1:
            return i
    return len(cols)


def _cluster_upper_tier_by_breadth(df_upper, col_to_ct):
    """
    Upper tier 两级排序：
      Level 1(宏观): 按该DEG涉及的cell type数量降序排列 -> 漏斗形宏观结构
      Level 2(微观): 数量相同的行内部用 hamming距离+average linkage+
                      optimal_ordering 层次聚类，让共享模式相近的行聚拢
    返回: (排序后的DataFrame, [(start_row, end_row, n_celltype), ...])
    """
    if len(df_upper) < 2:
        return df_upper, [(0, len(df_upper), None)]

    cols = list(df_upper.columns)
    n_ct_per_row = []
    for _, row in df_upper.iterrows():
        cols_on = [c for c in cols if row[c] == 1]
        cts = set(col_to_ct[c] for c in cols_on)
        n_ct_per_row.append(len(cts))

    df_tagged = df_upper.copy()
    df_tagged['_n_ct'] = n_ct_per_row

    ordered_blocks, breadth_boundaries, running_idx = [], [], 0
    for n_ct in sorted(df_tagged['_n_ct'].unique(), reverse=True):
        block = df_tagged[df_tagged['_n_ct'] == n_ct].drop(columns=['_n_ct'])
        if len(block) >= 2:
            try:
                dist = pdist(block.values.astype(float), metric='hamming')
                Z = linkage(dist, method='average', optimal_ordering=True)
                block = block.iloc[leaves_list(Z)]
            except Exception as e:
                print(f"  [warn] breadth={n_ct} 组内聚类失败，保持原序: {e}")
        ordered_blocks.append(block)
        breadth_boundaries.append((running_idx, running_idx + len(block), n_ct))
        running_idx += len(block)

    return pd.concat(ordered_blocks), breadth_boundaries


def _process_direction(csv_path, suffix, order_celltype_priority, celltype_rank):
    """
    读取单一方向(up/down)的 presence matrix csv，执行完整的三层分类
    (upper/middle/lower) + upper两级breadth聚类 + middle/lower内部排序。
    返回包含画图所需全部中间结果的字典。
    """
    parse_col = _parse_col_factory(suffix)
    all_expected = _expected_conditions(suffix, order_celltype_priority)

    df = pd.read_csv(csv_path, index_col=0)
    if 'Freq' in df.columns:
        df = df.drop(columns=['Freq'])

    conditions_full = [c for c in all_expected if c in df.columns]
    if len(conditions_full) == 0:
        raise ValueError(
            f"[{suffix}] 未匹配到任何列，请检查 CSV 列名格式是否为 "
            f"'OvsY_O_{suffix}_<celltype>' 且 order_celltype_priority 与生成"
            f"presence matrix时使用的顺序一致。CSV路径: {csv_path}"
        )
    df = df[conditions_full].copy()

    col_to_ct  = {c: parse_col(c)[0] for c in conditions_full}
    col_to_grp = {c: parse_col(c)[1] for c in conditions_full}

    def classify_gene(row):
        cols_on = [c for c in conditions_full if row[c] == 1]
        cts  = set(col_to_ct[c]  for c in cols_on)
        grps = set(col_to_grp[c] for c in cols_on)
        if len(cts) >= 2:
            return 'shared_celltype'
        elif len(cts) == 1 and len(grps) >= 2:
            return 'shared_group'
        else:
            return 'unique'

    tier = df.apply(classify_gene, axis=1)
    freq = df.sum(axis=1)

    df_upper  = df[tier == 'shared_celltype'].copy()
    df_middle = df[tier == 'shared_group'].copy()
    df_lower  = df[tier == 'unique'].copy()

    assert len(df_upper) + len(df_middle) + len(df_lower) == len(df), \
        f"[{suffix}] 三层分类行数总和不等于总行数"

    df_upper_s, upper_breadth_boundaries = _cluster_upper_tier_by_breadth(df_upper, col_to_ct)

    def middle_sort_key(row):
        cols_on = [c for c in conditions_full if row[c] == 1]
        ct = col_to_ct[cols_on[0]]
        ct_rank = celltype_rank.get(ct, 999)
        n_grp_on = len(cols_on)
        first_idx = _first_one_col_index(row, conditions_full)
        return (ct_rank, -n_grp_on, first_idx)

    if len(df_middle) > 0:
        # 注意：df.apply(axis=1) 对返回tuple的函数会自动展开成多列DataFrame，
        # 这里用列表推导式+sorted手动排序以规避该坑
        middle_keys = [middle_sort_key(row) for _, row in df_middle.iterrows()]
        order_idx = sorted(range(len(middle_keys)), key=lambda i: middle_keys[i])
        df_middle_s = df_middle.iloc[order_idx]
    else:
        df_middle_s = df_middle

    if len(df_lower) > 0:
        df_lower = df_lower.copy()
        df_lower['_ord'] = df_lower.apply(
            lambda r: _first_one_col_index(r, conditions_full), axis=1)
        df_lower_s = df_lower.sort_values('_ord').drop(columns=['_ord'])
    else:
        df_lower_s = df_lower

    middle_counts_per_ct = {ct: 0 for ct in order_celltype_priority}
    if len(df_middle_s) > 0:
        for _, row in df_middle_s.iterrows():
            cols_on = [c for c in conditions_full if row[c] == 1]
            ct = col_to_ct[cols_on[0]]
            middle_counts_per_ct[ct] += 1

    lower_counts_per_col = (df_lower_s[conditions_full].sum(axis=0)
                            if len(df_lower_s) > 0
                            else pd.Series(0, index=conditions_full))

    n_upper, n_middle, n_lower = len(df_upper_s), len(df_middle_s), len(df_lower_s)
    n_total = n_upper + n_middle + n_lower

    mat = np.vstack([
        df_upper_s.values.astype(float)  if n_upper  > 0 else np.empty((0, len(conditions_full))),
        df_middle_s.values.astype(float) if n_middle > 0 else np.empty((0, len(conditions_full))),
        df_lower_s.values.astype(float)  if n_lower  > 0 else np.empty((0, len(conditions_full))),
    ])

    return dict(
        conditions_full=conditions_full, parse_col=parse_col,
        mat=mat, n_upper=n_upper, n_middle=n_middle, n_lower=n_lower, n_total=n_total,
        upper_breadth_boundaries=upper_breadth_boundaries,
        middle_counts_per_ct=middle_counts_per_ct,
        lower_counts_per_col=lower_counts_per_col,
        df_upper_s=df_upper_s, df_middle_s=df_middle_s, df_lower_s=df_lower_s,
        freq=freq,
    )


def _plot_panel(ax, data, theme, panel_label, title_text, order_celltype_priority,
                age_label_map):
    """在给定 ax 上绘制完整的三层热图面板（含独立表头、独立左侧标注）"""
    conditions_full = data['conditions_full']
    parse_col = data['parse_col']
    mat = data['mat']
    n_upper, n_middle = data['n_upper'], data['n_middle']
    n_cols = len(conditions_full)

    cmap_custom = LinearSegmentedColormap.from_list(
        'panel_cmap', [theme['cmap_low'], theme['cmap_high']], N=256)

    im = ax.imshow(mat, aspect='auto', cmap=cmap_custom, vmin=0, vmax=1,
                    interpolation='nearest', rasterized=True)

    # 锁定 xlim/ylim，防止外部标注(负坐标)撑大/错位 axes 边框
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(max(len(mat), 1) - 0.5, -0.5)
    ax.set_autoscale_on(False)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([age_label_map[parse_col(c)[1]] for c in conditions_full], fontsize=6)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.tick_params(axis='x', pad=2)
    ax.set_yticks([])

    ax.axhline(y=n_upper - 0.5, color='black', linewidth=1.1, linestyle='--')
    if n_middle > 0:
        ax.axhline(y=n_upper + n_middle - 0.5, color='black', linewidth=1.1, linestyle='--')

    for start_row, end_row, n_ct in data['upper_breadth_boundaries'][:-1]:
        ax.axhline(y=end_row - 0.5, color='gray', linewidth=0.4, linestyle=':', alpha=0.6)

    # 左侧标注：混合坐标(transAxes+transData) + 仅纯文字(无线条)
    trans_left = blended_transform_factory(ax.transAxes, ax.transData)
    TEXT_X_FRAC = -0.06

    ax.text(TEXT_X_FRAC, (n_upper - 1) / 2 if n_upper > 0 else -2,
            f'Num. of\nDEGs\n(n={n_upper})',
            transform=trans_left, va='center', ha='right', fontsize=7, clip_on=False)

    if n_middle > 0:
        middle_center_y = n_upper + (n_middle - 1) / 2
        ax.text(TEXT_X_FRAC, middle_center_y,
                f'Num. of\nDEGs\n(n={n_middle})',
                transform=trans_left, va='center', ha='right', fontsize=6.5,
                clip_on=False, color=theme['color_middle'])

    # 顶部 cell type 分组标签（混合坐标）
    trans_top = blended_transform_factory(ax.transData, ax.transAxes)
    CELLTYPE_Y_FRAC = 1.10
    group_size = 3
    n_groups_block = n_cols // group_size
    y_middle_label = n_upper + n_middle - 0.3 if n_middle > 0 else n_upper - 0.3

    for gi in range(n_groups_block):
        ct = order_celltype_priority[gi] if gi < len(order_celltype_priority) else None
        if ct is None:
            continue
        cnt = data['middle_counts_per_ct'].get(ct, 0)
        if cnt > 0:
            start = gi * group_size
            center = start + (group_size - 1) / 2
            ax.text(center, y_middle_label, str(cnt),
                    ha='center', va='bottom', fontsize=6,
                    color=theme['color_middle'], fontweight='bold', clip_on=False)

    for j, col in enumerate(conditions_full):
        cnt = int(data['lower_counts_per_col'].get(col, 0))
        if cnt > 0:
            ax.text(j, len(mat) + max(len(mat), 1) * 0.012, str(cnt),
                    ha='center', va='bottom', fontsize=6,
                    color=theme['color_upper'], fontweight='bold', rotation=90, clip_on=False)

    for gi in range(n_groups_block):
        start = gi * group_size
        end = start + group_size - 1
        center = (start + end) / 2
        ct_name = order_celltype_priority[gi] if gi < len(order_celltype_priority) else ''
        ax.text(center, CELLTYPE_Y_FRAC, ct_name, transform=trans_top,
                ha='center', va='bottom', fontsize=7.5, fontweight='bold', clip_on=False)
        if end < n_cols - 1:
            ax.axvline(x=end + 0.5, color='black', linewidth=0.6, alpha=0.5)

    ax.set_title(title_text, fontsize=10, pad=54, color=theme['color_upper'])
    ax.text(-0.02, 1.16, panel_label, transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='bottom', ha='right', clip_on=False)

    return im


# ============================================================
# 主函数：唯一需要对外暴露的公共 API
# ============================================================
def plot_combined_deg_heatmap(
    up_csv,
    down_csv,
    order_celltype_priority,
    data_type='TE',
    out_prefix=None,
    theme_up=None,
    theme_down=None,
    age_label_map=None,
    min_panel_frac=0.30,
    fig_width_per_col=0.42,
    fig_width_min=10.0,
    fig_height=15.0,
    dpi=300,
    save_png=True,
    save_pdf=True,
    save_csv=True,
    verbose=True,
):
    """
    生成 Up(红) / Down(蓝) 双面板 DEG presence heatmap，整合自 TE/ERV/gene
    三套原始脚本的通用逻辑。适用于任何"celltype × 3个age-group比较
    (OvsY/MvsY/OvsM) 的二值 presence matrix"数据。

    Parameters
    ----------
    up_csv, down_csv : str
        分别指向 up/down 方向的 full_gene_presence_matrix.csv 路径。
        要求列名格式为 '{OvsY|MvsY|OvsM}_{O|M}_{up|down}_{celltype}'。
    order_celltype_priority : list[str]
        cell type 的展示顺序（必须与生成 presence matrix 时使用的顺序
        完全一致，否则顶部分组标签会错位）。
    data_type : str, default 'TE'
        用于拼接标题/文件名，如 'TE' / 'ERV' / 'genes'。
        注意：TE/ERV建议用单数，gene集合建议用复数('genes')。
    out_prefix : str or None
        输出文件（png/pdf/csv）的路径前缀（不含扩展名）。
        默认为 up_csv 所在目录下的 f'full_fig_combined_{data_type}_v3'。
    theme_up, theme_down : dict or None
        配色主题字典，需包含 keys:
        color_upper, color_middle, cmap_low, cmap_high, title, suffix。
        默认分别为红色系(up)/蓝色系(down)。
    age_label_map : dict or None
        年龄组比较标签映射，默认 {'OvsY':'O/Y','MvsY':'M/Y','OvsM':'O/M'}。
    min_panel_frac : float, default 0.30
        up/down 面板各自的最小视觉高度占比（仅显示层面保底，不影响
        实际数据/CSV导出）。当两方向DEG总数差异悬殊时防止某面板
        被压缩到无法阅读。
    fig_width_per_col, fig_width_min, fig_height : float
        画布尺寸参数（英寸）。fig_width = max(fig_width_min,
        n_cols * fig_width_per_col)。
    dpi : int
        PNG/PDF 保存分辨率。
    save_png, save_pdf, save_csv : bool
        是否保存对应格式的输出文件。
    verbose : bool
        是否打印处理进度信息。

    Returns
    -------
    dict
        包含 'fig'(matplotlib Figure对象), 'data_up', 'data_down'
        (各自完整的中间处理结果字典), 'frac_up', 'frac_down'
        (最终应用的面板高度比例)。

    Examples
    --------
    >>> result = plot_combined_deg_heatmap(
    ...     up_csv='/path/to/2_TE/1_up/full_gene_presence_matrix.csv',
    ...     down_csv='/path/to/2_TE/2_down/full_gene_presence_matrix.csv',
    ...     order_celltype_priority=[
    ...         "Mono", "Mac", "DC", "NK", "T", "B",
    ...         "Epi", "EC", "Fib", "SMC"
    ...     ],
    ...     data_type='TE',
    ... )
    """
    # ---------- 参数初始化 ----------
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype']  = 42
    matplotlib.rcParams['font.family']  = 'Arial'

    theme_up   = theme_up   or DEFAULT_THEME_UP
    theme_down = theme_down or DEFAULT_THEME_DOWN
    age_label_map = age_label_map or DEFAULT_AGE_LABEL_MAP
    celltype_rank = {ct: i for i, ct in enumerate(order_celltype_priority)}

    if out_prefix is None:
        out_dir = os.path.dirname(os.path.dirname(up_csv))  # up_csv所在目录的上一级
        out_prefix = os.path.join(out_dir, f'full_fig_combined_{data_type}_v3')

    # ---------- 数据处理 ----------
    data_up   = _process_direction(up_csv,   theme_up['suffix'],   order_celltype_priority, celltype_rank)
    data_down = _process_direction(down_csv, theme_down['suffix'], order_celltype_priority, celltype_rank)

    if verbose:
        print(f"[up]   Upper={data_up['n_upper']}  Middle={data_up['n_middle']}  Lower={data_up['n_lower']}")
        print(f"[down] Upper={data_down['n_upper']}  Middle={data_down['n_middle']}  Lower={data_down['n_lower']}")

    n_total_up, n_total_down = data_up['n_total'], data_down['n_total']
    raw_total = n_total_up + n_total_down
    frac_up   = max(n_total_up   / raw_total, min_panel_frac)
    frac_down = max(n_total_down / raw_total, min_panel_frac)
    s = frac_up + frac_down
    frac_up, frac_down = frac_up / s, frac_down / s

    if verbose:
        print(f"📐 面板高度比例: Up={frac_up:.2%} (n_total={n_total_up})  "
              f"Down={frac_down:.2%} (n_total={n_total_down})")

    # ---------- 画图 ----------
    n_cols = len(data_up['conditions_full'])
    fig_w = max(fig_width_min, n_cols * fig_width_per_col)

    fig = plt.figure(figsize=(fig_w, fig_height))
    gs = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[frac_up, frac_down],
        top=0.90, bottom=0.06, left=0.14, right=0.93,
        hspace=0.55
    )
    ax_up   = fig.add_subplot(gs[0, 0])
    ax_down = fig.add_subplot(gs[1, 0])

    im_up = _plot_panel(
        ax_up, data_up, theme_up, 'A',
        f"{theme_up['title']} {data_type} (Aging)",
        order_celltype_priority, age_label_map)
    im_down = _plot_panel(
        ax_down, data_down, theme_down, 'B',
        f"{theme_down['title']} {data_type} (Aging)",
        order_celltype_priority, age_label_map)

    cbar_up = fig.colorbar(im_up, ax=ax_up, shrink=0.35, pad=0.015)
    cbar_up.set_label('DEG presence\n(1 = detected)', fontsize=6.5)
    cbar_up.ax.tick_params(labelsize=5.5)

    cbar_down = fig.colorbar(im_down, ax=ax_down, shrink=0.35, pad=0.015)
    cbar_down.set_label('DEG presence\n(1 = detected)', fontsize=6.5)
    cbar_down.ax.tick_params(labelsize=5.5)

    # ---------- 保存图片 ----------
    if save_png:
        fig.savefig(f'{out_prefix}.png', dpi=dpi, bbox_inches='tight')
    if save_pdf:
        fig.savefig(f'{out_prefix}.pdf', dpi=dpi, bbox_inches='tight')

    # ---------- 保存原始数据 CSV ----------
    if save_csv:
        for direction, data in [('up', data_up), ('down', data_down)]:
            labels = []
            for start, end, n_ct in data['upper_breadth_boundaries']:
                labels.extend([n_ct] * (end - start))

            data['df_upper_s'].assign(
                tier='shared_celltype',
                Freq=data['freq'].loc[data['df_upper_s'].index],
                n_celltype_sharing=labels
            ).to_csv(f'{out_prefix}_{direction}_upper_data.csv')

            if data['n_middle'] > 0:
                data['df_middle_s'].assign(
                    tier='shared_group', Freq=data['freq'].loc[data['df_middle_s'].index]
                ).to_csv(f'{out_prefix}_{direction}_middle_data.csv')

            data['df_lower_s'].assign(
                tier='unique', Freq=data['freq'].loc[data['df_lower_s'].index]
            ).to_csv(f'{out_prefix}_{direction}_lower_data.csv')

        pd.DataFrame({
            'direction': ['up', 'up', 'up', 'down', 'down', 'down'],
            'tier': ['upper', 'middle', 'lower'] * 2,
            'n_DEGs': [data_up['n_upper'], data_up['n_middle'], data_up['n_lower'],
                       data_down['n_upper'], data_down['n_middle'], data_down['n_lower']],
        }).to_csv(f'{out_prefix}_tier_summary_combined.csv', index=False)

        pd.DataFrame({
            'direction': ['up', 'down'],
            'n_total_real': [n_total_up, n_total_down],
            'raw_fraction': [n_total_up / raw_total, n_total_down / raw_total],
            'final_panel_fraction_with_floor': [frac_up, frac_down],
            'min_panel_frac_applied': [min_panel_frac, min_panel_frac],
        }).to_csv(f'{out_prefix}_panel_height_ratio.csv', index=False)

    if verbose:
        print(f"✅ 整合图已保存: {out_prefix}.png / .pdf (A=Up在上, B=Down在下)")
        print(f"✅ 面板高度已按真实DEG总数比例(带{min_panel_frac:.0%}保底)分配")

    plt.close(fig)
    return dict(fig=fig, data_up=data_up, data_down=data_down,
                frac_up=frac_up, frac_down=frac_down)


# ============================================================
# 使用示例（脱敏路径，实际使用时替换为你自己的数据路径）
# ============================================================
if __name__ == '__main__':

    CELLTYPE_ORDER = [
        "Mono", "Mac", "DC",   # Myeloid
        "NK", "T", "B",        # Lymphoid
        "Epi",                 # Epithelial
        "EC", "Fib", "SMC"     # Stromal
    ]

    # 示例1: TE
    plot_combined_deg_heatmap(
        up_csv='./data/TE/1_up/full_gene_presence_matrix.csv',
        down_csv='./data/TE/2_down/full_gene_presence_matrix.csv',
        order_celltype_priority=CELLTYPE_ORDER,
        data_type='TE',
    )

    # 示例2: ERV
    plot_combined_deg_heatmap(
        up_csv='./data/ERV/1_up/full_gene_presence_matrix.csv',
        down_csv='./data/ERV/2_down/full_gene_presence_matrix.csv',
        order_celltype_priority=CELLTYPE_ORDER,
        data_type='ERV',
    )

    # 示例3: genes（注意用复数）
    plot_combined_deg_heatmap(
        up_csv='./data/gene/1_up/full_gene_presence_matrix.csv',
        down_csv='./data/gene/2_down/full_gene_presence_matrix.csv',
        order_celltype_priority=CELLTYPE_ORDER,
        data_type='genes',
    )
