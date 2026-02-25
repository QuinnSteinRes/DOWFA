"""DOWFA calculation library.

Main entry points:
- simulate_legacy_stock (dowfa.legacy_stock)
- simulate_wva (dowfa.wva)
"""
from .legacy_stock import simulate_legacy_stock, LegacyStockResult
from .wva import simulate_wva, WVAResult

from .legacy_emissions import legacy_route_emissions_pulse, LegacyEmissions

from .wva_emissions import wva_route_emissions, WVAEmissions

from .legacy_emissions import legacy_route_emissions_continuous
