#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
example_usage.py
────────────────────────────────────────────────────────────────
plot_combined_deg_heatmap() 的最小可运行示例（脱敏路径）。
使用前请将 up_csv / down_csv 替换为你自己数据的实际路径。
────────────────────────────────────────────────────────────────
"""

import sys
import os

# 若未安装为包，直接从上级目录导入
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from deg_heatmap import plot_combined_deg_heatmap

CELLTYPE_ORDER = [
    "Mono", "Mac", "DC",   # Myeloid
    "NK", "T", "B",        # Lymphoid
    "Epi",                 # Epithelial
    "EC", "Fib", "SMC"     # Stromal
]

result = plot_combined_deg_heatmap(
    up_csv='./data/example_up_presence_matrix.csv',
    down_csv='./data/example_down_presence_matrix.csv',
    order_celltype_priority=CELLTYPE_ORDER,
    data_type='TE',
)

print("Panel height fractions:", result['frac_up'], result['frac_down'])

