# Bundled fonts

Self-hosted so the GUI keeps its typography offline — it runs as a local
server, and a CDN font would silently fall back to a system face on a
machine with no network.

| File | Family | Used for | Source |
|---|---|---|---|
| `Jost.woff2` | Jost (variable) | App title, card headers | [google/fonts](https://github.com/google/fonts/tree/main/ofl/jost) |
| `SpaceGrotesk.woff2` | Space Grotesk (variable) | UI text, network labels | [google/fonts](https://github.com/google/fonts/tree/main/ofl/spacegrotesk) |
| `JetBrainsMono.woff2` | JetBrains Mono (variable) | Alignment viewer, mutation ticks | [google/fonts](https://github.com/google/fonts/tree/main/ofl/jetbrainsmono) |

Jost is the display face because it descends from Futura, which Paul
Renner drew alongside the Bauhaus — the closest OFL type to that
lettering.

All three are licensed under the SIL Open Font License 1.1; see `OFL.txt`.
Converted from the upstream TTFs to woff2 with `fontTools`.
