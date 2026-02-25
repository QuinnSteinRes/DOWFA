"""Utility helpers for sampling and validation."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Tuple, Optional
import numpy as np

def rng(seed: Optional[int] = None) -> np.random.Generator:
    return np.random.default_rng(seed)

def uniform_samples(r: np.random.Generator, lo: float, hi: float, n: int) -> np.ndarray:
    if hi < lo:
        raise ValueError(f"Invalid range: {lo}..{hi}")
    return r.uniform(lo, hi, n)

def agent_key(agent: str) -> str:
    a = agent.strip().lower()
    if a in {"des", "desflurane"}:
        return "des"
    if a in {"sev", "sevo", "sevoflurane"}:
        return "sev"
    if a in {"iso", "isoflurane"}:
        return "iso"
    raise ValueError("agent must be one of: des, sev, iso")

def gwp_horizon_key(gwp: str | int) -> str:
    if isinstance(gwp, int):
        if gwp in (20, 100, 500):
            return f"GWP_{gwp}"
        raise ValueError("gwp int must be 20, 100, or 500")
    g = str(gwp).strip().upper().replace(" ", "")
    if g in {"20", "GWP20", "GWP_20"}:
        return "GWP_20"
    if g in {"100", "GWP100", "GWP_100"}:
        return "GWP_100"
    if g in {"500", "GWP500", "GWP_500"}:
        return "GWP_500"
    raise ValueError("gwp must be 20/100/500 or GWP_20/GWP_100/GWP_500")
