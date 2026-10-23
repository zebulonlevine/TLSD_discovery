"""Where every script in this repo reads and writes.

Two roots, both laid out like the original analysis tree
(GenomeWide_results/, Sublibrary_results/, fixed_NotI/, ...):

    data   inputs: deposited counts, library annotations, upstream stage outputs
    out    everything a script writes

For an end-to-end run set both to the same directory. To verify a single stage
against published outputs, point `data` at the published tree and `out` at a
scratch directory, so nothing published can be overwritten.

Resolution order for each root: environment variable (TLSD_DATA, TLSD_OUT),
then config/paths.yaml at the repo root. In-repo resources (category tables,
configs, library design tables) are addressed with repo().

Usage in a script two levels below the repo root:

    import os, sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
    from common import paths
    counts = paths.inp('Sublibrary_results/results/counts_for_mageck.txt')
    outdir = paths.outp('fixed_NotI/scores')
    cats   = paths.repo('screen_figures/categories/manyCats.csv')

Importing this module also exports TLSD_DATA / TLSD_OUT / TLSD_REPO to the
environment, so YAML configs can use ${TLSD_DATA} etc. (expanded on load).
"""
import os

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
_CONFIG = os.path.join(REPO, 'config', 'paths.yaml')


def _from_config(key):
    if not os.path.exists(_CONFIG):
        return None
    import yaml
    with open(_CONFIG) as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg.get(key)


def _root(env, key):
    v = os.environ.get(env) or _from_config(key)
    if not v:
        raise RuntimeError(
            f'{env} is not set and config/paths.yaml has no "{key}". '
            'Copy config/paths.example.yaml to config/paths.yaml and edit it, '
            f'or export {env}=/path/to/root.')
    return os.path.abspath(os.path.expanduser(os.path.expandvars(v)))


def data_root():
    return _root('TLSD_DATA', 'data')


def out_root():
    return _root('TLSD_OUT', 'out')


def inp(*parts):
    """Path to an input under the data root (read only)."""
    return os.path.join(data_root(), *parts)


def outp(*parts, mkdir=False):
    """Path under the output root; mkdir=True creates the parent directory."""
    p = os.path.join(out_root(), *parts)
    if mkdir:
        os.makedirs(os.path.dirname(p) if os.path.splitext(p)[1] else p, exist_ok=True)
    return p


def repo(*parts):
    """Path to a resource shipped inside this repository."""
    return os.path.join(REPO, *parts)


def expand(s):
    """Expand ${TLSD_DATA}/${TLSD_OUT}/${TLSD_REPO} (and ~) in a config string."""
    return os.path.expanduser(os.path.expandvars(s)) if isinstance(s, str) else s


def _export():
    os.environ.setdefault('TLSD_REPO', REPO)
    for env, fn in (('TLSD_DATA', data_root), ('TLSD_OUT', out_root)):
        try:
            os.environ.setdefault(env, fn())
        except RuntimeError:
            pass          # unresolved roots only fail when a script actually asks for them


_export()
