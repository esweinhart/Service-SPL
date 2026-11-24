# ServiceSPL 2.0

*(formerly ServiceSPL: your settings and mic profiles carry over automatically)*

A pro sound level meter for services, shows, theaters and mix rooms. It runs as a normal Apple Silicon Mac app (no Python
needed) and on Windows.

**What it does**

- Big live meter (LAS by default) with a thumbs up/down, a 5-minute scrolling history and yellow/red levels
- LAeq 1 min / 15 min / since reset, LCeq, LAFmax, true-peak LCpeak, LA10 / LA50 / LA90
- Calibration wizard with a **profile per mic**: acoustic calibrator, UMIK-1/UMIK-2 cal files, or manual numbers
- **1/3-octave RTA** and a scrolling **spectrogram**
- **Pink noise from your speakers** for calibration and room EQ (Left + Right, one side, or any output)
- **Room EQ helper**: measures several seats, compares with a target curve and tells you what to cut
- **OSHA / NIOSH noise dosimeter** with action-level and limit warnings
- **Session logging**: CSV every 1 s / 10 s / 1 min, markers, alarms, a **PDF report**, session CSV export
- **Scheduled logging** (e.g. every Sunday 8:30 AM – 12:30 PM, with a *Logging Complete* animation at the end) and **auto-save to OneDrive** in dated folders
- **Planning Center Services Live**: every log row is stamped with the current item, with levels per item
- **Bitfocus Companion**: live level and color on a Stream Deck button
- **NDI** output of just the meter graphic, for your video switcher
- **Mac menu bar** readout (`● 87.3`, dot goes green/yellow/red), a mini always-on-top meter, and a **web view**
- Dashboard layout: every panel can be shown, hidden, dragged or stacked, always inside the main window. Every
  popup (settings, calibration, schedule, guide, messages) opens inside the main window too; only the Mini meter,
  the Web view (a browser) and NDI live outside it. Follows macOS light/dark.
- A second toolbar row has one button per docked panel (RTA | Spectrogram, Session), plus **Show all**, **¼ screen** and **Fit**
- **Dosimeter** and **Planning Center** open as pop-up cards over the dashboard (drag the title to move, corner to
  resize, ✕ to close). They never leave the main window and reopen where you left them. Their buttons sit on the
  right of the bottom bar
- Toolbar buttons can be **dragged into your own order**: click the **🔒** (end of the top row) to unlock, drag,
  click again to lock. Dividers drag just like buttons while unlocked. Right-click a toolbar any time to **add a
  divider**, or right-click a divider to **remove** it (dragging one off the bar removes it too)
- **Mini**, **Web view** and **NDI** sit on the left of the bottom bar
- Drag the thin lines between the meter and the panels (they light up on hover) to resize them
- **Help → Keyboard shortcuts…** (⌘/): ⌘R start/stop, ⌘L log, ⌘⇧R reset, M / ⇧M marker, ⌘E CSV, ⌘P PDF,
  ⌘1-⌘5 panels, ⌘⇧M Mini, ⌘, Settings
- Little **!** marks next to readings and settings: hover over one for a plain-English explanation of LAeq, LAFmax, LCpeak, L10/L90, A/C/Z weighting, Fast/Slow, dose, TWA, dBFS and more
- A built-in **Guide** (**Help → Guide…**) walks you through all of this inside the app

---

## 1. Install

### Mac (no Python needed)
1. Open this repo's **Releases** page. **latest** is the released version, **preview** is the newest test build of
   work in progress.
2. Download **ServiceSPL-Mac-AppleSilicon.dmg** (or the .zip) for M1/M2/M3/M4 and newer Macs. Only the Apple Silicon
   Mac app is built now. Intel Macs can still run it from source (below).
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
3. Open **Settings → Meter → Calibration wizard…** (or Tools → Calibration wizard…) and follow the wizard (next section).
4. Set **Red at** to your house limit. The reading turns yellow 5 dB before it and red at it, and the thumb flips down.
5. **Help → Guide…** any time for step-by-step help.

There are two test inputs at the bottom of the input menu: **Simulated pink noise** and **Simulated 94 dB calibrator**.
Handy for trying things without a mic.

---

## 3. Calibration and mic profiles
Every mic gets its own **profile** (toolbar menu next to Calibrate). Switch mics and the right calibration comes with
it. Tick *"pick this profile automatically when this input device is selected"* in the wizard and it switches by
itself.

The **Calibration wizard** (Settings → Meter) offers these methods:

| Method | When | What happens |
|---|---|---|
| **Acoustic calibrator** | Best accuracy, any mic | Seat a 94 or 114 dB calibrator (1 kHz or 250 Hz), wait for *steady*, press **Measure (5 s)**. Measured unweighted, so it works with either frequency. |
| **Mic calibration file** | USB mics like miniDSP UMIK-1 / UMIK-2, analog measurement mics like the **Audix TM1**, Dayton EMM-6, or any REW-style `.txt` / `.cal` / `.frd` / `.csv` | Loads the frequency-response correction (applied as a linear-phase FIR to everything: meter, RTA, logs). If the file has a `Sens Factor` (UMIK), the absolute level is set from it too. |
| **Pink noise from my speakers + a reference SPL meter** | No calibrator, but you have a trusted handheld meter | ServiceSPL plays pink noise out of your speakers (**▶ Play pink noise** on that page), hold the meter next to the mic, set the **same weighting (A/C/Z) and speed** as the handheld, type its reading and press **Measure (10 s)**. |
| **Manual** | You know the numbers | Either "0 dBFS sine = X dB SPL", or mic sensitivity (mV/Pa) + preamp gain + interface dBu at 0 dBFS. |

**Calibrate with pink noise from your speakers** (also **Tools → Calibrate with pink noise from my speakers…**)
1. Set the mic's preamp gain and phantom power where you'll leave them for the service. If you change the gain
   later, calibrate again.
2. Put the mic (e.g. the TM1) on a stand at the main listening position, pointing at the ceiling. Hold or mount
   your handheld meter right next to it, at the same height and pointing the same way.
3. On the wizard page, pick **Speakers out of** (the output that feeds your console or speakers; the Mac's current
   sound output is picked for you), **Channels** (Left + Right, Left only, Right only or any output number) and
   **Level** (start at **-20 dBFS** with the speakers turned down), then press **▶ Play pink noise**.
   **Using your soundboard's generator instead?** Set **Pink noise comes from** to *My soundboard / console generator
   plays it*: the playback controls hide and ServiceSPL just listens.
4. Set **Reference SPL (target)** (85 dB by default, anything from 60 to 110 dB; 75-80 dB is kinder to small rooms)
   and bring the speakers up until the handheld reads about that. Keep the room quiet: no HVAC blasts, talking or
   walking around.
5. Set the handheld and ServiceSPL to the **same weighting and speed**. **Z (flat) + Slow** is the easiest match,
   or C on both if the handheld has no Z.
6. Press **Start meter** if it isn't running, and wait until both readings settle.
7. Type in what the handheld reads, then press **Measure (10 s)**. Don't move anything while it measures.
8. Press **Next**. The reading in the same weighting should match the handheld. Then **Save**.

**UMIK tips**
- Download your mic's file from miniDSP with its serial number. There are two: **0°** (mic pointing at the source) and
  **90°** (mic pointing at the ceiling, typical for rooms and theaters).
- For the Sens Factor to be right, the Mac's input volume for the UMIK must be at **100 %** (Audio MIDI Setup →
  UMIK → Input → volume all the way up). Set the format to 48 kHz while you're there.
- A cal file plus a quick calibrator check afterwards is the gold standard.

**Audix TM1 (and other analog mics)**
- In the Calibration wizard, choose **Mic calibration file** and pick the TM1's file. The serial number, model and the mV/Pa
  sensitivity are read from it, and the profile is named after the mic (e.g. "TM1 #210345").
- Analog mics have no Sens Factor, because the level depends on your preamp gain. So the file fixes the *frequency
  response*, and then you set the *absolute level* with the **acoustic calibrator** or the pink noise
  **Match a reference SPL meter** step on the same profile. The wizard goes straight there after you save the file.
  Do it again whenever you change the preamp gain.
- Phantom power on, and keep the preamp gain where you'll leave it for the service.

**Compare like with like.** A-weighting reads lower than C or Z on pink noise and bass-heavy sound, often by
3 to 10 dB. The big number's title (e.g. **LAS · A-weighted · Slow**) tells you what you're looking at, and the
wizard's last page shows the reading in the dashboard's weighting with LAS / LCS / LZF underneath.

After any gain change on an analog mic/interface, **recalibrate**. The bottom bar shows the active profile, how it was
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
tabs. Panels stay inside the main window (no floating windows). The layout is saved when you quit.

**Any window size.** Drag the window as small as about a quarter of a laptop screen (or press **¼ screen**) or as
big as a 4K/5K display. The toolbars wrap onto extra rows when the window is narrow, the meter drops the small tiles
and then the history graph as it gets short, and panels shrink down to their title bar and scroll. Hide panels you
don't need for the tiniest window, or use **Mini ⤡** for a floating always-on-top box.

### RTA · 1/3 octave
Weighting Z/A/C, speed Fast/Slow/**Leq** (average since **Reset**). White ticks are peak hold.

### Spectrogram
Its own panel, so you can place it anywhere or close it. It scrolls the last ~20 seconds, 20 Hz at the bottom to
20 kHz at the top, and brighter means louder. It's good for spotting feedback rings and hum.

### Dosimeter · OSHA / NIOSH
| Meter | Criterion | Exchange | Threshold | Notes |
|---|---|---|---|---|
| OSHA PEL | 90 dBA / 8 h | 5 dB | 90 dB | The legal permissible exposure limit |
| OSHA HC | 90 dBA / 8 h | 5 dB | 80 dB | Hearing conservation; **50 % = action level** |
| NIOSH | 85 dBA / 8 h | 3 dB | 80 dB | Recommended limit (88 dBA = 4 h, 100 dBA = 15 min) |

Each shows dose %, TWA, an 8-hour projection and the time left to 100 %. Warnings fire at 50 % and 100 %, when the
pace projects over 100 %, and if any peak hits 140 dB. **Reset dose** at the start of a shift or call. The dose keeps
running across Stop/Start. Guidance only; use a certified dosimeter for compliance records.

### Session · markers · alarms · export
- Time history of the whole session (last 10 min / hour / 4 h / everything), markers in yellow, alarms as red ticks.
- **◆ Marker** or press **M**. **Shift+M** asks for a note ("Sermon", "Band in").
- **Alarms**: each rule (metric + threshold) fires once when crossed: red banner, beep (optional), log entry. It
  re-arms 1 dB below. Metrics: LAS, LAeq 1 min, LAeq long, LCeq long, LCpeak, LAFmax, OSHA HC dose %, NIOSH dose %.
  Defaults: LAeq 15 min ≥ 100, LCpeak ≥ 135.
- **Export CSV…** (one row per second for the whole session), **PDF report…**, **New session**, **Schedule…**.
- Venue and Operator fields go on the report.

### Service · Planning Center Live
ServiceSPL follows **Services Live**: when whoever runs Live advances to the next item, a new item starts, the log rows
get stamped with it, and the levels are kept per item. Setup takes about 3 minutes, once.

**1. Make a Personal Access Token (in your browser)**
1. Sign in to Planning Center with an account that can **see the plans** you want to follow (Services *Viewer* or higher
   on that service type). The token gets the same access as that account.
2. Go to **https://api.planningcenteronline.com/oauth/applications** (the panel's **Open token page** button goes there).
3. Under **Personal Access Tokens**, click **New Personal Access Token**.
4. Description: `ServiceSPL`. Keep the newest **Services** API version offered. Click **Submit**.
5. Copy the **Application ID** and the **Secret**. Keep the page open until step 2 works.

**2. Connect ServiceSPL**
1. **View → Service · Planning Center Live** to show the panel.
2. Paste the Application ID into the first box and the Secret into the second, then press **Connect**. The status line
   says *Connected · N service types*. On a Mac the secret is saved in your Keychain, so this is a one-time step.
3. Pick the **service type** (left drop-down). The right one lists recent and upcoming plans; today's is selected.
4. Set **Item target LAeq**. An item whose average goes over it turns red and raises an alarm (95 dB is a common
   worship target).

**3. On the day**
1. Start the meter and logging (or let the **Schedule** do it).
2. Tick **Follow Live**. The status says *Following Live: …*.
3. Someone opens **Services Live** for that plan, **takes control** and advances items. ServiceSPL checks every 3 s and
   shows **Now: <item> LAeq …**.
4. Every log row has an `Item` column. Afterwards the panel table and the PDF report list LAeq / LCeq / LAFmax / LCpeak
   per item.

**Troubleshooting**
| Symptom | Fix |
|---|---|
| `CERTIFICATE_VERIFY_FAILED` (ServiceSPL 2.0 on a Mac) | Fixed in 2.1: the app now carries its own certificate bundle. Update the app; your keys were fine. If it still happens, check the Mac's date/time and any web filter. |
| `401` / Unauthorized | Wrong Application ID or Secret, or the token was revoked. Make a new token. |
| Connected, but no service types or plans | The account behind the token can't see them. Give it Viewer access in Services. |
| Following, but the item never changes | Nobody has taken control in Services Live, or a different plan is picked. Check the plan's date. |
| Couldn't connect / timed out | The Mac needs internet access to `api.planningcenteronline.com`. |

The **? Setup help** button in the panel opens these same steps inside the app.
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
Also in the app: **Help → Guide → Stream Deck (Companion)**.

1. Install **Bitfocus Companion** and press **Launch GUI**.
2. Companion **Settings** → OSC: enable the **OSC Listener** on port **12321**.
3. Companion **Variables → Custom Variables**: add `spl_level`, `spl_leq1m`, `spl_leqlong`, `spl_peak`, `spl_dose`,
   `spl_state`, `spl_color`, `spl_item` (plus `spl_d1` … `spl_d4` if you build your own big-number layout).
   `spl_state` is `good` / `warn` / `bad` for feedback rules, `spl_color` is the hex color.
4. ServiceSPL: Settings → Integrations → tick **Send levels to Bitfocus Companion (OSC)**, set the Companion computer
   to `127.0.0.1` (or the Companion machine's IP) and the port to `12321`.
5. Put `$(custom:spl_level)` on a button with text size **Auto**. Test with Companion's **Emulator**, no Stream Deck
   needed: press Start in ServiceSPL and watch the button update.

6. Or the zero-setup way: Settings → Integrations → **Button location**, a page/row/column like `1/0/0`, and pick
   **Paint button(s)**:
   - **One button**: that button shows the level + weighting (e.g. `87.3` / `LAS`) and turns green / yellow / red.
   - **Four buttons in a row**: that button and the 3 to its right show the reading big, one piece per key:
     `8` · `7` · `.3` · `dBA` (or `10` · `2` · `.4` · `dBA` above 100). All 4 turn green / yellow / red together.
     In Companion, set those 4 buttons' **text size to Auto** (or the largest size) so each piece fills its key.

   ServiceSPL updates them 4 times a second. Building your own layout instead? The 4 pieces are also in the custom
   variables `spl_d1` … `spl_d4` (create them first).

Good to know: the ports must match; many apps can send to Companion's port but only one app per computer can listen
on it; across machines, allow UDP 12321 through the firewall.

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
2. On the Proton, open Edge at that address: a big thumb and number that update 5×/second and fit any window
   (side by side when wide, thumb on top when tall). Browsers won't make a window narrower than about 400-500 px;
   an OBS / ProPresenter browser source can be any size.
3. Clean window: a desktop shortcut with
   `"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe" --app=http://192.168.1.50:8080`
- Use the regular LAN, not the SoundGrid port. Give the Mac a fixed IP (DHCP reservation).
- It's read-only. If the Mac app stops, the page shows "no connection" and reconnects by itself.

---

## 8. Room EQ helper (Help → Room EQ helper)
1. Pick the output and channels of the speakers you're tuning, **▶ Play pink noise** (full-band, -20 dBFS) and bring it to
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
| NOT CALIBRATED in the bottom bar | Run the Calibration wizard (Settings → Meter) for the active profile, or pick the right profile in the toolbar. |
| "dropouts" in the bottom bar | The computer couldn't keep up or the USB mic hiccupped. Close heavy apps; try another USB port. |
| No menu-bar number | Settings → Display → menu-bar box. (If the Mac build lacks the Cocoa bridge, a colored icon is shown instead.) |
| Companion button doesn't change | Create the custom variables first (including `spl_d1`…`spl_d4` if you use them); check OSC is enabled in Companion and the IP/port match; for *Paint this button*, use page/row/column like `1/0/0`. |
| NDI button says runtime not found | Install NDI Tools (ndi.video/tools), or browse to `libndi.dylib` in Settings → Integrations. |
| Planning Center problems | See the troubleshooting table in the Planning Center section above, or **? Setup help** in the panel. |
| Nothing in OneDrive | Settings → Logging: auto-save folder set? Is OneDrive signed in and syncing? Files land under `ServiceSPL/<year>/<date>`. |
| Schedule didn't start | The app must be open (it can be minimized or in mini mode). Check the window's days/date and that it's ticked On. |
| Pink noise: no sound | Pick the right **Speakers out of** device and channels; check the Mac's output volume and that the speakers are up. For ASIO, only the chosen channels are opened. |
| Web view: page won't load | Same network? Allow incoming connections on the Mac. Try another port in Settings → Integrations. |

---

## 11. Files and data

### Back up / move your settings (File → Export settings… / Import settings…)
Everything goes in one `.servicespl.json` file: every mic profile (with its calibration offset **and** the cal-file
curve, so you don't need the original file again), meter and alarm settings, theme, panel layout and window size,
Companion, NDI, web view, the logging schedule, the OneDrive / auto-save folder and Planning Center.

- **Updating the app:** you don't need to do anything. Settings live outside the app (below), so dragging a new
  ServiceSPL into Applications keeps them. Exporting first is a cheap safety net.
- **New Mac, or a second booth computer:** Export on the old one, Import on the new one. The app restarts with
  everything loaded.
- **Planning Center secret:** Export asks whether to include it. If you say yes, the file contains a password-like
  secret, so keep it private (don't email it around or leave it on a shared drive). If you say no, paste the secret
  again after importing.
- After updating, macOS may ask for **microphone** permission again, because each build is signed fresh. Click Allow.

- Settings and mic profiles: `~/Library/Application Support/ServiceSPL/config.json` (Mac),
  `%APPDATA%\ServiceSPL\config.json` (Windows). On first launch, settings from the old `SPLMonitor` folder are read in. The Planning Center secret is in the macOS Keychain.
- Logs: `Documents/SPL Logs` by default, or the auto-save folder.
- Tests: `python tests/test_iec61672.py` (weighting, time weighting, Leq, peak) and `python tests/test_engine.py`
  (true peak, RTA, cal files incl. TM1/CSV/decimal-comma formats, dose, pink noise, Ln, EQ, OSC, schedule, settings backup). Both run in CI before every Mac build.
