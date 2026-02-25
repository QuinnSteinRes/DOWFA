# dowfa

This is a simplified, **non-GUI** calculation library derived from the original Tkinter prototype.

It provides:
- Legacy stock disposal pathway calculator (`simulate_legacy_stock`)
- WVA pathway calculator (`simulate_wva`)

## Install (editable)
From the folder containing this repo:

```bash
python -m pip install -e .
```

## Quick start
```bash
python examples/run_legacy_stock.py
python examples/run_wva.py
```

## Notes
- All parameter names and equations follow the original scripts, but UI/API code has been removed.
- Monte Carlo sampling uses NumPy's Generator; set `seed=` for reproducibility.


## Table 2 printer

```bash
python examples/generate_table2.py
```
\n\n## Temperature model (Supplementary)\n\n- Reproduce Marin-style curves: `python examples/temperature_validate_marin.py`\n- Legacy stock pulse temperature: `python examples/temperature_legacy_pulse.py`\n
- WVA route temperature: `python examples/temperature_wva_routes.py`

- Legacy continuous (50y then stop): `python examples/temperature_legacy_pulse.py`
- Legacy continuous comparison plot: `python examples/temperature_legacy_continuous_compare.py`

- WVA route temperature plots (Marin scale): `python examples/temperature_wva_plots_marin_scale.py`
