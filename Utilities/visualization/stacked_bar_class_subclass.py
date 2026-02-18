# ═══════════════════════════════════════════════════════════════════════════
# FUNCTION          | plot_stacked_bar_class_subclass, plot_nested_pie_chart
# PURPOSE           | Stacked bar: Class(X) → Subclass(stacked Y) visualization
# INPUTS            | motifs(List[Dict]), title(str), figsize(tuple), exclude_classes(List[str])
# OUTPUTS           | plt.Figure object
# COLOR SCHEME      | UNIFIED_MOTIF_COLORS base + alpha gradient (0.4-1.0) for subclasses
# ORDERING          | Biological taxonomy via get_all_classes_taxonomy_order()
# ANNOTATIONS       | Count on segment (if ≥3), total above bar (n=X)
# COMPATIBILITY     | plot_nested_pie_chart() → wrapper redirect to stacked bar
# ═══════════════════════════════════════════════════════════════════════════

import matplotlib.pyplot as plt; import numpy as np; from collections import Counter; from matplotlib.colors import to_rgba
from Utilities.config.colors import UNIFIED_MOTIF_COLORS; from Utilities.config.motif_taxonomy import get_all_classes_taxonomy_order, get_subclasses_for_class

def plot_stacked_bar_class_subclass(motifs, title="Class → Subclass Composition", figsize=(14, 8), exclude_classes=None):
    if not motifs: fig, ax = plt.subplots(figsize=(8, 6)); ax.text(0.5, 0.5, 'No motifs', ha='center', va='center', fontsize=16); ax.axis('off'); return fig
    motifs = [m for m in motifs if m.get('Class') not in (exclude_classes or [])]
    if not motifs: fig, ax = plt.subplots(figsize=(8, 6)); ax.text(0.5, 0.5, 'All excluded', ha='center', va='center', fontsize=16); ax.axis('off'); return fig
    
    class_totals = Counter(); class_subclass = {}
    for m in motifs: cls, sub = m.get('Class', 'Unknown'), m.get('Subclass', 'Other'); class_totals[cls] += 1; class_subclass.setdefault(cls, Counter())[sub] += 1
    
    taxonomy = get_all_classes_taxonomy_order(); ordered = [c for c in taxonomy if c in class_totals] + [c for c in class_totals if c not in taxonomy]
    fig, ax = plt.subplots(figsize=figsize); x_pos = np.arange(len(ordered)); bottom = np.zeros(len(ordered))
    
    for cls_idx, cls in enumerate(ordered):
        base_rgba = to_rgba(UNIFIED_MOTIF_COLORS.get(cls, '#808080')); subs = class_subclass.get(cls, {})
        try: ordered_subs = [s for s in get_subclasses_for_class(cls) if s in subs]
        except: ordered_subs = sorted(subs.keys())
        
        for idx, sub in enumerate(ordered_subs):
            cnt = subs[sub]; alpha = 0.4 + 0.6 * (idx / max(1, len(ordered_subs) - 1)) if len(ordered_subs) > 1 else 1.0
            ax.bar(x_pos[cls_idx], cnt, 0.7, bottom=bottom[cls_idx], color=(*base_rgba[:3], alpha), edgecolor='white', linewidth=1.5)
            if cnt >= 3: ax.text(x_pos[cls_idx], bottom[cls_idx] + cnt/2, f"{cnt}", ha='center', va='center', fontsize=9, fontweight='bold', color='white')
            bottom[cls_idx] += cnt
        ax.text(x_pos[cls_idx], bottom[cls_idx] + 1, f"n={class_totals[cls]}", ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    ax.set_xticks(x_pos); ax.set_xticklabels(ordered, rotation=45, ha='right', fontsize=12, fontweight='bold')
    for lbl, cls in zip(ax.get_xticklabels(), ordered): lbl.set_color(UNIFIED_MOTIF_COLORS.get(cls, '#808080'))
    ax.set_ylabel('Motif Count', fontsize=14, fontweight='bold'); ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout(); return fig

def plot_nested_pie_chart(motifs, title="Class → Subclass", **kwargs): return plot_stacked_bar_class_subclass(motifs, title=title, **kwargs)
