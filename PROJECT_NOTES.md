# SPL Monitor: project notes

## Hardware / setups
- **Windows:** Focusrite Scarlett + reference mic, or a mic on the **Waves SoundGrid** network via the SoundGrid ASIO driver (pick the driver channel).
- **Mac:** USB measurement mic (UMIK / UMM-6 style). This is the preferred setup.
- **Display:** the Mac runs the meter and turns on **Web view**. The **Axis Proton** (SuperRack SoundGrid host) shows it in Edge (`--app=http://<mac-ip>:8080`) over the regular LAN, not the SoundGrid port. Nothing is installed on the Proton.
- Companion runs on a separate Mac, and the GainStage OSC bridge runs on the Proton.

## Features (current, v2.0)
- See README.md for the full user-facing list. v2.0 (Nov 2026) added: calibration wizard + per-mic profiles (calibrator /
  UMIK cal file with FIR correction + Sens Factor / manual), LCeq, LA10/50/90, true-peak LCpeak, 1/3-octave RTA +
  spectrogram, OSHA PEL/HC + NIOSH dosimeter, theater calibration (pink noise per channel, presets, auto run), room EQ
  helper, session history/markers/alarms, session CSV + PDF report, scheduled logging, auto-save to OneDrive in dated
  folders, Planning Center Services Live (item stamps, per-item levels), Companion (OSC custom variables + button paint,
  HTTP /api/json), NDI output (ctypes to the installed NDI runtime), menu-bar meter (pyobjc NSStatusItem), dock-based
  dashboard, light/dark following macOS, in-app Guide.
- Web view: `/` page, `/data` JSON for the page, `/api/json` and `/api/level` for other software.

## Design decisions
- v1 was deliberately a simple meter + logging. In Nov 2026 Eric asked for the full pro feature set, so the RTA and
  history came back as optional dock panels; the big meter stays the centre of the window.
- UI choices (Eric, 2026-11-19): dashboard with dockable panels, follow macOS light/dark, big number + scrolling history,
  menu bar "● 87.3". AI chat assistant was dropped in favour of the built-in Guide.
- NDI uses ctypes against the user's installed NDI runtime because cyndilib has no Intel-Mac wheel (Universal build).
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
- Test NDI, the menu-bar item and Planning Center Live against the real thing on Eric's Mac
