# SPL Monitor: project notes

## Hardware / setups
- **Windows:** Focusrite Scarlett + reference mic, or a mic on the **Waves SoundGrid** network via the SoundGrid ASIO driver (pick the driver channel).
- **Mac:** USB measurement mic (UMIK / UMM-6 style). This is the preferred setup.
- **Display:** the Mac runs the meter and turns on **Web view**. The **Axis Proton** (SuperRack SoundGrid host) shows it in Edge (`--app=http://<mac-ip>:8080`) over the regular LAN, not the SoundGrid port. Nothing is installed on the Proton.
- Companion runs on a separate Mac, and the GainStage OSC bridge runs on the Proton.

## Features (current)
- Main meter: LAS by default (A/C/Z, Fast/Slow selectable)
- Tiles: LAeq 1 min, LAeq long (15 min default), LAeq since reset, LAFmax, LCpeak
- "Red at" level in the toolbar. Yellow starts N dB below red (default 5, 0 = off).
- Thumb: green up below the red level, red down at or above it (3 s hold). Can follow the live meter or LAeq 1 min.
- Mini mode: small frameless box, draggable, 📌 keep-on-top, double-click to return
- Logging: CSV (UTF-8 BOM for Excel) every 1 s / 10 s / 1 min. Each row has the interval LAeq, LCeq, LAFmax, LCpeak, end-of-interval LAS, LAeq long and the calibration info.
- Calibration: acoustic calibrator (94/114 dB) measured A-weighted over 3 s, or a manual offset. Saved in the config.
- Web view: built-in read-only HTTP server. `/` serves the page and `/data` serves JSON polled every 200 ms. Default port 8080.
- Device auto-pick: USB measurement mic, then Scarlett (WASAPI > MME/Core Audio). ASIO enabled. The ASIO path opens only the selected channel.

## Design decisions
- The user asked for a **simple meter + logging**. The RTA and history graph were removed on purpose.
- A/C weighting: bilinear transform with a pre-warped 12.2 kHz pole, normalized to 0 dB at 1 kHz. Within IEC 61672 Class 1 tolerance up to 12.5 kHz at 48 kHz.
- Time weighting: exponential, 125 ms / 1 s. Leq is computed from per-second energy sums (2 h kept).
- Config: `%APPDATA%\SPLMonitor\config.json` on Windows, `~/Library/Application Support/SPLMonitor/config.json` on Mac.
- Mac .app build adds NSMicrophoneUsageDescription and an ad-hoc codesign; without it the mic reads silence.

## Tested so far (in a Linux sandbox, simulated input)
- Weighting curves vs IEC values; 1 kHz calibration reads back 94.0 dB
- GUI, mini mode, pin toggle, logging intervals, web page rendered in Chrome at small and full-screen sizes
- Mocked SoundGrid device list: a 64-channel picker opens the correct ASIO channel
- **Not yet tested on real Windows, a real Mac, a Scarlett, SoundGrid or a USB mic.**

## Open ideas
- Load a USB mic calibration file (sensitivity factor) so calibration works without a calibrator
- Send the level to Companion (OSC/HTTP) to show on a Stream Deck button
