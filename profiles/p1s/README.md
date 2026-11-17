# P1S / PETG print profiles

| Variant | Bambu Studio project | Estimate | MakerWorld |
|---|---|---|---|
| Knob cup | [3MF](knob-cup-p1s-petg.3mf) | 4 h 18 min / 64.15 g | [Published](https://makerworld.com/en/models/3630615-scarlett-solo-4th-gen-10-inch-rack-mount#profileId-4163291) |
| Centered solid | [3MF](centered-solid-p1s-petg.3mf) | 3 h 57 min / 57.69 g | Draft |
| Centered front vents | [3MF](centered-front-vents-p1s-petg.3mf) | 3 h 53 min / 55.56 g | Draft |

All variants are **84.2 mm deep**. Estimates use Bambu Studio 02.09.01.70.

- P1S, 0.4 mm nozzle; PETG, 250 °C nozzle / 75 °C textured PEI.
- 0.2 mm layers, two walls, 20% gyroid, five top/bottom layers.
- Front face down, rotated 45°; normal supports from the build plate.
- No brim or purge tower. Remove supports and the Scarlett rubber feet before installation.

These editable projects include geometry, settings and thumbnails. Slice before printing.

Validation: [knob cup](validation.json), [centered solid](centered-solid-validation.json), [centered front vents](centered-front-vents-validation.json). Each checks the mesh against its STL, all 421 layers, bed bounds, temperatures, support toolpaths and volumetric flow. The knob profile also checks the cup membrane and magnet pocket.

To reproduce, run `python scripts/prepare_p1s.py <variant>`, open `build/p1s/<variant>/<variant>-project.3mf` in Bambu Studio, slice and save. Use the portable `settings.json` with its preset override list. Slice the saved project into `build/p1s/<variant>/verified`, then run `python scripts/verify_p1s.py <variant>`.
