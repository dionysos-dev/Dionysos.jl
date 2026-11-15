#!/usr/bin/env python3
"""Plot prompt-weighted SAT-and-implication rates from saved benchmark outputs.

No LLM, Spot, or controller calls. Recounts individual prompt outcomes, checks
the summary, and writes PDF/SVG/PNG plus an exact-data/provenance JSON file.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile


def verified_curve(source):
    source = Path(source).resolve()
    data = json.loads(source.read_text())
    if not data['complete']:
        raise ValueError('Input analysis is incomplete')
    n = data['prompt_count']
    prompts = data['prompt_results']
    if len(prompts) != n or len({p['prompt_index'] for p in prompts}) != n:
        raise ValueError('Prompt count or prompt identities do not match')
    rows = sorted(data['summary'], key=lambda row: row['k'])
    if [row['k'] for row in rows] != list(range(1, 11)):
        raise ValueError('Expected complete conjunction budgets k=1,...,10')
    curve = []
    for prompt in prompts:
        if [c['k'] for c in prompt['conjunctions']] != list(range(1,len(prompt['conjunctions'])+1)):
            raise ValueError('Saved cumulative conjunctions are not consecutive')
    for row in rows:
        k = row['k']
        good = 0
        for prompt in prompts:
            conjunctions = prompt['conjunctions']
            if not conjunctions:
                continue  # No candidate: abstention, retained in denominator.
            # The saved list contains one item per available distinct cluster.
            # After saturation, reuse the final cumulative conjunction.
            c = conjunctions[min(k, len(conjunctions)) - 1]
            logical_good = c['satisfiable'] is True and c['implies_reference'] is True
            if logical_good != c['logical_good']:
                raise ValueError('Per-prompt logical-good flag disagrees with SAT and implication')
            good += logical_good
        rate = good / n
        if row['total_prompts'] != n or row['counts']['logical_good'] != good:
            raise ValueError(f'Per-prompt recount disagrees with summary at k={k}')
        if not math.isclose(row['rates']['logical_good'], rate, abs_tol=1e-12):
            raise ValueError(f'Saved rate disagrees with exact counts at k={k}')
        curve.append(dict(k=k, logical_good=good, denominator=n, rate=rate, percentage=100*rate))
    return dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        model=data['model'], temperature=data['temperature'], repetitions=data['repetitions'],
        prompt_count=n, denominator_policy='all NL prompts, including no-candidate and unsatisfiable cases',
        event='SAT(Phi_k) and Phi_k implies reference',
        conjunction_policy='SAT-filtered candidates; frequency-ranked semantic classes; saturation at min(k,m)',
        curve=curve)


def plot_semantic_inclusion(source, output_directory):
    exact = verified_curve(source)
    # Keep rendering usable in a headless terminal with a writable font cache.
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'ltl-figure-matplotlib'))
    os.environ.setdefault('XDG_CACHE_HOME', str(Path(tempfile.gettempdir())/'ltl-figure-cache'))
    Path(os.environ['XDG_CACHE_HOME']).mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter

    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=True)
    k = [row['k'] for row in exact['curve']]
    percent = [row['percentage'] for row in exact['curve']]
    model_name = {'deepseek-v4-pro':'DeepSeek V4 Pro'}.get(exact['model'],exact['model'])
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.labelsize': 10, 'axes.titlesize': 10,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    })
    paths = []
    # Detail view has explicit 65–85% ticks; full-scale version supplies 0–100%.
    for name, ylim, ticks in [
            ('semantic_inclusion', (65,85), [65,70,75,80,85]),
            ('semantic_inclusion_full_scale', (0,100), [0,20,40,60,80,100])]:
        fig, ax = plt.subplots(figsize=(5.2,3.2), constrained_layout=True)
        ax.plot(k, percent, color='#1764AB', marker='o', markersize=4.5,
                linewidth=1.8, zorder=3)
        ax.set(xlabel=r'Semantic-cluster budget, $k$',
               ylabel='Non-trivial semantic inclusion',
               title=f"{model_name}: $T={exact['temperature']}$, $R={exact['repetitions']}$",
               xlim=(0.65,10.35), ylim=ylim, xticks=k, yticks=ticks)
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100,decimals=0))
        ax.grid(axis='y', color='#D6DCE2', linewidth=0.6, zorder=0)
        ax.spines[['top','right']].set_visible(False)
        ax.annotate(f'{percent[0]:.2f}%', (k[0],percent[0]), xytext=(8,-18),
                    textcoords='offset points', color='#1764AB', ha='left')
        ax.annotate(f'{percent[-1]:.2f}%', (k[-1],percent[-1]), xytext=(-4,10),
                    textcoords='offset points', color='#1764AB', ha='right')
        for suffix in ('pdf','svg','png'):
            path=output/f'{name}.{suffix}'
            fig.savefig(path,dpi=300,facecolor='white')
            paths.append(str(path))
        plt.close(fig)
    exact['matplotlib_version']=matplotlib.__version__
    exact['plot_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output/'semantic_inclusion_data.json').write_text(json.dumps(exact,indent=2)+'\n')
    caption=(f"On the full task-specific NL-to-LTL benchmark ({exact['prompt_count']:,} prompts), "
        'the prompt-weighted rate of non-trivial semantic inclusion, '
        r'$\operatorname{SAT}(\Phi_k(p))\land\Phi_k(p)\models\varphi^\star$, '
        f'increases from {percent[0]:.2f}% at k=1 to {percent[-1]:.2f}% at k=10. '
        f"DeepSeek V4 Pro is sampled R={exact['repetitions']} times per prompt at T={exact['temperature']}. "
        'Individually invalid or unsatisfiable candidates are discarded; semantic-equivalence '
        'classes are ranked by frequency and one representative per class is conjoined. '
        'The conjunction uses min(k,m(p)) available classes. All prompts remain in the denominator.\n')
    (output/'caption.txt').write_text(caption)
    for row in exact['curve']:
        print(f"k={row['k']:2d}: {row['logical_good']}/{row['denominator']} = {row['percentage']:.2f}%")
    for path in paths:
        print(path)
    return exact,paths


if __name__=='__main__':
    root=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=root/'Results/deepseek-v4-pro-task-specific-full-T0.7-K10-conjunctions/summary.json')
    parser.add_argument('--output',type=Path,default=root/'Figures/semantic-inclusion-r10')
    args=parser.parse_args()
    plot_semantic_inclusion(args.input,args.output)
