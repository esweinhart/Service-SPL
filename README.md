# ServiceSPL 2.0

*(formerly ServiceSPL: your settings and mic profiles carry over automatically)*

A pro sound level meter for services, shows, theaters and mix rooms. It runs as a normal Mac app (no Python needed)
and on Windows.

**What it does**

- Big live meter (LAS by default) with a thumbs up/down, a 5-minute scrolling history and yellow/red levels
- LAeq 1 min / 15 min / since reset, LCeq, LAFmax, true-peak LCpeak, LA10 / LA50 / LA90
- Calibration wizard with a **profile per mic**: acoustic calibrator, UMIK-1/UMIK-2 cal files, or manual numbers
- **1/3-octave RTA** and a scrolling **spectrogram**
- **Theater calibration**: built-in pink noise generator, per-channel 85 dBC (cinema) targets, auto run
- **Room EQ helper**: measures several seats, compares with a target curve and tells you what to cut
- **OSHA / NIOSH noise dosimeter** with action-level and limit warnings
- **Session logging**: CSV every 1 s / 10 s / 1 min, markers, alarms, a **PDF report**, session CSV export
- **Scheduled logging** (e.g. every Sunday 08:30–12:30) and **auto-save to OneDrive** in dated folders
- **Planning Center Services Live**: every log row is stamped with the current item, with levels per item
- **Bitfocus Companion**: live level and color on a Stream Deck button
- **NDI** output of just the meter graphic, for your video switcher
- **Mac menu bar** readout (`● 87.3`, dot goes green/yellow/red), a mini always-on-top meter, and a **web view**
- Dashboard layout: every panel can be shown, hidden, dragged, stacked or floated. Follows macOS light/dark.
- A built-in **Guide** (toolbar **? Guide**) walks you through all of this inside the app

---

## 1. Install

### Mac (no Python needed)
1. Open this repo's **Releases** page. **latest** is the released version, **preview** is the newest test build of
   work in progress.
2. Download **ServiceSPL-Mac-AppleSilicon.dmg** (or the .zip). It's a Universal app, so it runs natively on M1/M2/M3/M4 *and* on
   Intel Macs. (**ServiceSPL-Mac-Intel.dmg** is an Intel-only build, if you prefer.)
3. Open the .dmg and drag **ServiceSPL** onto the Applications folder next to it.
4. First launch: **right-click → Open → Open**. If macOS still refuses, go to System Settings → Privacy & Security,
   scroll down and click **Open Anyway**. (The app is ad-hoc signed, not notarized, which is why macOS asks once.)
5. Click **Allow** when it asks for the microphone. Missed it? System Settings → Privacy & Security → Microphone →
   turn on ServiceSPL. Without this, macOS feeds the app silence.

A fresh build is made automatically every time the code changes (GitHub Actions, see `.github/workflows/build.yml`).

### Windows
Download **ServiceSPL-Windows.zip** from the same Releases page, unzip and run **ServiceSPL.exe**.

### From source (any OS)
```
pip install numpy scipy PySide6 sounddevice          # + pyobjc-framework-Cocoa on a Mac for the menu-bar text
python spl_monitor.py
```

---

## 2. First run in 60 seconds
1. Pick your mic in the input menu on the toolbar. A UMIK / UMM-6 / EMM-6 is picked automatically, then a Scarlett.
2. Press **▶ Start**.
3. Press **Calibrate…** and follow the wizard (next section).
4. Set **Red at** to your house limit. The reading turns yellow 5 dB before it and red at it, and the thumb flips down.
5. Hit **? Guide** any time for step-by-step help.

There are two test inputs at the bottom of the input menu: **Simulated pink noise** and **Simulated 94 dB calibrator**.
Handy for trying things without a mic.

---

## 3. Calibration and mic profiles
Every mic gets its own **profile** (toolbar menu next to Calibrate). Switch mics and the right calibration comes with
it. Tick *"pick this profile automatically when this input device is selected"* in the wizard and it switches by
itself.

**Calibrate…** offers three methods:

| Method | When | What happens |
|---|---|---|
| **Acoustic calibrator** | Best accuracy, any mic | Seat a 94 or 114 dB calibrator (1 kHz or 250 Hz), wait for *steady*, press **Measure (5 s)**. Measured unweighted, so it works with either frequency. |
| **Mic calibration file** | USB mics like miniDSP UMIK-1 / UMIK-2, or any REW-style `.txt` | Loads the frequency-response correction (applied as a linear-phase FIR to everything: meter, RTA, logs). If the file has a `Sens Factor`, the absolute level is set from it too. |
| **Manual** | You know the numbers | Either "0 dBFS sine = X dB SPL", or mic sensitivity (mV/Pa) + preamp gain + interface dBu at 0 dBFS. |

**UMIK tips**
- Download your mic's file from miniDSP with its serial number. There are two: **0°** (mic pointing at the source) and
  **90°** (mic pointing at the ceiling, typical for rooms and theaters).
- For the Sens Factor to be right, the Mac's input volume for the UMIK must be at **100 %** (Audio MIDI Setup →
  UMIK → Input → volume all the way up). Set the format to 48 kHz while you're there.
- A cal file plus a quick calibrator check afterwards is the gold standard.

After any gain change on an analog mic/interface, **recalibrate**. The status bar shows the active profile, how it was
calibrated, and turns yellow when the calibration is more than 30 days old.

---

## 4. The meter and statistics

| Tile | Meaning |
|---|---|
| Big number | Live level, LAS by default (weighting A/C/Z and Fast/Slow in Settings → Meter) |
| History strip | The last 5 minutes (1–30 min in Settings → Display); white line is LAeq 1 min, red dashed line is *Red at* |
| LAeq 1 min / 15 min | A-weighted average over the last minute / the long window (Settings → Meter) |
| LAeq since reset | Average since you pressed **Reset** |
| LCeq 15 min | C-weighted average: the bass energy. "LC − LA" above ~15 dB means a very bass-heavy mix |
| LAFmax | Loudest A-weighted Fast reading since reset |
| LCpeak | Highest C-weighted **true** peak since reset (4× oversampled). Turns red at 135 dB |
| LA10 / LA90 | Level exceeded 10 % / 90 % of the time over the long window. LA90 is the background (HVAC, room noise). Small text shows LA50 and the since-reset LA90 |

**Reset** clears Leq/max/peak/Ln. It doesn't touch the dose or the session history.

### Standards
- A and C weighting: IEC 61672-1 Class 1 tolerance (checked automatically at 44.1, 48 and 96 kHz, see `tests/`).
- Fast 125 ms / Slow 1 s exponential time weighting, checked with the IEC tone-burst and decay-rate tests.
- 1/3-octave bands: base-10 centre frequencies, 6th-order Butterworth bands (IEC 61260 style), 20 Hz – 20 kHz.
- This is a software meter. Its accuracy depends on your mic and calibration; it isn't a certified Class 1 instrument.

---

## 5. Dashboard panels
Use the **View** menu to show/hide panels. Drag a panel's title bar to move it, drop it on another to stack them as
tabs, or double-click to float it. The layout is saved when you quit.

### RTA · 1/3 octave + spectrogram
Weighting Z/A/C, speed Fast/Slow/**Leq** (average since *Reset peaks / Leq*). White ticks are peak hold. The
spectrogram scrolls the last ~20 seconds, 20 Hz at the bottom to 20 kHz at the top.

### Dosimeter · OSHA / NIOSH
| Meter | Criterion | Exchange | Threshold | Notes |
|---|---|---|---|---|
| OSHA PEL | 90 dBA / 8 h | 5 dB | 90 dB | The legal permissible exposure limit |
| OSHA HC | 90 dBA / 8 h | 5 dB | 80 dB | Hearing conservation; **50 % = action level** |
| NIOSH | 85 dBA / 8 h | 3 dB | 80 dB | Recommended limit (88 dBA = 4 h, 100 dBA = 15 min) |

Each shows dose %, TWA, an 8-hour projection and the time left to 100 %. Warnings fire at 50 % and 100 %, when the
pace projects over 100 %, and if any peak hits 140 dB. **Reset dose** at the start of a shift or call. The dose keeps
running across Stop/Start. Guidance only; use a certified dosimeter for compliance records.

### Theater calibration
1. **Output**: your interface (or *Simulated* to try it). **Format**: 2.0, 5.1, 7.1, 7.1.4. **Room** preset:

   | Preset | Screen | Surround / height | LFE | Noise |
   |---|---|---|---|---|
   | Cinema / dub stage | 85 dBC | 82 dBC | 91 dBC* | -20 dBFS |
   | Mid-size mix room | 82 | 79 | 88* | -20 dBFS |
   | Small room / nearfield | 79 | 76 | 85* | -20 dBFS |
   | Home theater | 75 | 75 | 85* | -30 dBFS |

   *LFE targets are a starting point: confirm on the RTA (about +10 dB in-band over the screen channels).
2. Check the **Out** column matches your patch (SMPTE order L R C LFE Ls Rs …). Out numbers and targets are editable.
3. Mic at the main listening position, pointing up (90° cal file).
4. **▶ Play**, adjust the channel until the big number says **on target ✓** (±0.5 dB), **Next ▶**. Or press
   **Auto run**: every channel plays 7 s (4 s settle + 3 s capture) and the table shows what to turn up or down.

Noise is pink, -20 dBFS RMS (AES17: a full-scale sine = 0 dBFS), band-limited 500 Hz – 2 kHz by default
(SMPTE/Dolby). Full-band and LFE-band noise are in the Noise menu. The LFE channel is always band-limited to 20–120 Hz.

### Session · markers · alarms · export
- Time history of the whole session (last 10 min / hour / 4 h / everything), markers in yellow, alarms as red ticks.
- **◆ Marker** or press **M**. **Shift+M** asks for a note ("Sermon", "Band in").
- **Alarms**: each rule (metric + threshold) fires once when crossed: red banner, beep (optional), log entry. It
  re-arms 1 dB below. Metrics: LAS, LAeq 1 min, LAeq long, LCeq long, LCpeak, LAFmax, OSHA HC dose %, NIOSH dose %.
  Defaults: LAeq 15 min ≥ 100, LCpeak ≥ 135.
- **Export CSV…** (one row per second for the whole session), **PDF report…**, **New session**, **Schedule…**.
- Venue and Operator fields go on the report.

### Service · Planning Center Live
1. At **https://api.planningcenteronline.com/oauth/applications** create a **Personal Access Token**.
2. Paste the **Application ID** and **Secret**, press **Connect** (on a Mac the secret is saved in your Keychain).
3. Pick the service type and plan (today's is selected) and tick **Follow Live**.
4. As Services Live advances, the current item shows in big type, stamps every log row (`Item` column), becomes an
   event, and gets its own LAeq / LCeq / LCpeak in the table. An item whose LAeq goes over **Item target** turns red
   and raises an alarm. The PDF report gets a per-item table.

No Planning Center? **Start item manually…** does the same with names you type.

---

## 6. Logging, schedules and OneDrive
- **● Log** writes a CSV (opens in Excel) every 1 s / 10 s / 1 min: LAeq, LCeq, LAFmax, LCpeak, LAS, LAeq long,
  LA10, LA90, OSHA HC and NIOSH dose %, the current service item, markers/events, and the calibration used.
- **Schedule…** (Session panel or Tools menu): add windows with days of the week (or a one-off date), start and stop
  time. Inside a window the meter starts if needed, a new session begins and logging runs; at the end logging stops
  (and the report is saved, see below). If you stop logging by hand during a window, it stays stopped until the next one.
- **Auto-save folder** (Settings → Logging & reports): press **Find OneDrive** (or Browse). Logs, the PDF report and the
  session CSV then go to

  `<OneDrive>/ServiceSPL/2026/2026-11-23 - Sunday AM/`

  (the plan name is added when Planning Center is followed) and OneDrive syncs them. With *"Save the PDF report +
  session CSV automatically when logging stops"* ticked, every stop of logging saves both.

---

## 7. Show it elsewhere

### Mac menu bar and mini meter
The menu bar shows `● 87.3` with the dot in green / yellow / red. Click it for Leq, peak, dose and quick controls
(start/stop, marker, logging, mini meter). Turn it off in Settings → Display.
**Mini ⤡** shrinks to a small box with the thumb and the number. 📌 keeps it on top, ⤢ or double-click goes back.

### Bitfocus Companion (Stream Deck)
1. In Companion, create custom variables: `spl_level`, `spl_leq1m`, `spl_leqlong`, `spl_peak`, `spl_dose`,
   `spl_state`, `spl_color`, `spl_item`. Make sure OSC is enabled in Companion (default port **12321**).
2. In ServiceSPL: Settings → Integrations → tick **Send levels to Bitfocus Companion**, set the Companion computer's IP
   (127.0.0.1 if it's the same Mac).
3. Put `$(custom:spl_level)` on a button. `spl_state` is `good` / `warn` / `bad` for feedback rules, `spl_color` is
   the hex color.
4. Or the zero-setup way: fill in **Paint this button** with a location like `1/0/0` (page/row/column). ServiceSPL then
   sets that button's text (level + weighting) and background color directly, 4 times a second.

Prefer polling? With **Web view** on, `http://<mac-ip>:8080/api/json` returns all values as JSON and
`/api/level` returns just the number.

### NDI
Press **NDI** on the toolbar. A source called **ServiceSPL** appears on the network with just the meter graphic (thumb,
big number, LAeq 1 min, current item, level bar) at 30 fps. Settings → Integrations: name, size (720p, 1080p, 540p or a
1280×360 lower-third strip) and **transparent background** for keying over video.
NDI needs the free **NDI Tools** (or NDI runtime) installed on the computer: https://ndi.video/tools. The app finds
`libndi` automatically; if not, point *NDI library* at it.

### Web view (e.g. the Axis Proton)
1. Click **Web view**. A box shows the address, e.g. `http://192.168.1.50:8080`. Allow incoming connections if macOS asks.
2. On the Proton, open Edge at that address: a big thumb and number that update 5×/second and fit any window.
3. Clean window: a desktop shortcut with
   `"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe" --app=http://192.168.1.50:8080`
- Use the regular LAN, not the SoundGrid port. Give the Mac a fixed IP (DHCP reservation).
- It's read-only. If the Mac app stops, the page shows "no connection" and reconnects by itself.

---

## 8. Room EQ helper (Help → Room EQ helper, or ? Guide → Room EQ)
1. Pick the output and channel of the speaker you're tuning, **▶ Play pink noise** (full-band, -20 dBFS) and bring it to
   a comfortable 80–85 dBC.
2. Mic where people sit. **Measure position (10 s)**. Do 3–5 seats; they're energy-averaged.
3. Pick a target: Flat, Live/worship house curve (+6 dB under 200 Hz, gentle HF roll-off), Cinema X-curve (SMPTE ST 202),
   small-room X-curve, or Hi-fi/studio tilt.
4. The chart shows each band above/below the target (green within ±2 dB). The list tells you what to do, e.g.
   *Cut 125 Hz by 5.2 dB (bell, Q ≈ 4.3)*. Boosts are capped at 3 dB, and anything more than 6 dB low is flagged as a
   room null to fix by moving the speaker or sub rather than by boosting.
5. Make the change, **Clear**, measure again.

---

## 9. Input setups
- **USB measurement mic (Mac)**: picked automatically. Audio MIDI Setup: 48 kHz, input volume 100 % for a cal file.
- **Focusrite Scarlett**: 48V on if the mic needs it, **Air off**, pick Ch 1 or Ch 2. Don't touch the gain after calibrating.
- **Waves SoundGrid (Windows)**: patch the mic to a SoundGrid Driver channel in SuperRack / SoundGrid Studio, choose
  *Waves SoundGrid [ASIO]* and that channel. The ASIO path opens only that channel. Recalibrate when you switch between
  Scarlett and SoundGrid.

---

## 10. Troubleshooting
| Problem | Fix |
|---|---|
| Reads 20–30 dB, doesn't move | macOS is giving the app silence: System Settings → Privacy & Security → Microphone → ServiceSPL on. Then quit and reopen. (The rename means macOS asks for the mic once more.) |
| "ServiceSPL can't be opened" | Right-click → Open → Open, or Privacy & Security → Open Anyway. |
| Levels way off with a UMIK cal file | Input volume must be 100 % in Audio MIDI Setup, and use the right file (0° vs 90°). Check with a calibrator. |
| NOT CALIBRATED in the status bar | Run Calibrate… for the active profile, or pick the right profile in the toolbar. |
| "dropouts" in the status bar | The computer couldn't keep up or the USB mic hiccupped. Close heavy apps; try another USB port. |
| No menu-bar number | Settings → Display → menu-bar box. (If the Mac build lacks the Cocoa bridge, a colored icon is shown instead.) |
| Companion button doesn't change | Create the custom variables first; check OSC is enabled in Companion and the IP/port match; for *Paint this button*, use page/row/column like `1/0/0`. |
| NDI button says runtime not found | Install NDI Tools (ndi.video/tools), or browse to `libndi.dylib` in Settings → Integrations. |
| Planning Center: "Couldn't connect" | Use a *Personal Access Token* (Application ID + Secret), not your login. Check the Mac has internet. |
| Planning Center doesn't follow | Someone has to be running **Services Live** for that plan; check you picked today's plan. |
| Nothing in OneDrive | Settings → Logging: auto-save folder set? Is OneDrive signed in and syncing? Files land under `ServiceSPL/<year>/<date>`. |
| Schedule didn't start | The app must be open (it can be minimized or in mini mode). Check the window's days/date and that it's ticked On. |
| Theater cal: no sound | Pick the right Output and check the Out numbers match your patch; for ASIO, only that one channel is opened. |
| Web view: page won't load | Same network? Allow incoming connections on the Mac. Try another port in Settings → Integrations. |

---

## 11. Files and data
- Settings and mic profiles: `~/Library/Application Support/ServiceSPL/config.json` (Mac),
  `%APPDATA%\ServiceSPL\config.json` (Windows). On first launch, settings from the old `SPLMonitor` folder are read in. The Planning Center secret is in the macOS Keychain.
- Logs: `Documents/SPL Logs` by default, or the auto-save folder.
- Tests: `python tests/test_iec61672.py` (weighting, time weighting, Leq, peak) and `python tests/test_engine.py`
  (true peak, RTA, cal files, dose, pink noise, Ln, EQ, OSC, schedule). Both run in CI before every Mac build.
