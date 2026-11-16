# SPL Monitor

A simple SPL meter with logging, a thumbs-up/down indicator, a mini always-on-top mode,
and a **Web view** for showing the meter on another computer's screen (like the Axis Proton).
It runs on **Windows** (Scarlett or Waves SoundGrid) and **Mac** (a USB measurement mic such as a miniDSP UMIK or Dayton UMM-6).

## Mac app, no Python needed (easiest)
1. Go to this repo's **Releases → latest** and download **SPL-Monitor-Mac-AppleSilicon.zip** (M1/M2/M3/M4) or **SPL-Monitor-Mac-Intel.zip** (older Macs). Not sure? Apple menu → About This Mac: "Chip: Apple…" means Apple Silicon.
2. Double-click the zip, then drag **SPL Monitor** into Applications.
3. First launch: **right-click → Open → Open**. (If macOS still says it can't be opened, go to System Settings → Privacy & Security, scroll down and click **Open Anyway**.)
4. Click **Allow** when it asks for the microphone. If you missed it: System Settings → Privacy & Security → Microphone → turn on SPL Monitor.
A fresh build is made automatically every time the code changes. Windows users can grab **SPL-Monitor-Windows.zip** from the same place.

## Mac setup with Python (the old way)
1. Install Python 3.12 or 3.13 from python.org. Then open the Python folder in Applications and run **Install Certificates.command** once.
2. Unzip, then **right-click `run.command` → Open** (a plain double-click gets blocked the first time). The first run takes a few minutes.
3. When macOS asks whether **Terminal** can use the microphone, click **Allow**. Leave the Terminal window open while the meter runs.
4. Plug in the USB mic. The app picks it automatically. Click **Start**, then calibrate (below).
5. In **Audio MIDI Setup**, set the USB mic to 48 kHz. Don't change its input volume after you calibrate.

Want a normal Mac app you can double-click? Run **`build_mac_app.command`** once (right-click → Open). You'll get `dist/SPL Monitor.app`; drag it to Applications. The first time, right-click → Open, then allow microphone access.

## Showing it on the Proton (all the processing stays on the Mac)
1. On the Mac, click **Web view**. A box shows the address, for example `http://192.168.1.50:8080`.
   If macOS asks to *"accept incoming network connections"*, click **Allow**.
2. On the Proton, open **Edge**, which is already installed, so nothing new goes on the Proton, and go to that address.
   You'll get a big thumb and number that update 5 times a second and resize to fit any window.
3. For a clean window without browser bars, make a desktop shortcut on the Proton with this target:
   `"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe" --app=http://192.168.1.50:8080`
   Size it and park it in a corner of the SuperRack screen.
- The Proton reaches the Mac over your **regular LAN**, the same network your GainStage OSC bridge uses to talk to Companion. Don't use the SoundGrid port.
- Give the Mac a fixed IP address (a DHCP reservation) so the shortcut keeps working.
- If the Mac app stops, the page shows "no connection" and a gray thumb. It reconnects automatically.
- The web view is read-only. Anyone on the network who opens the page can see the reading, but nobody can change anything.

## Windows setup (one time)
1. Install Python 3.10+ from python.org and tick **"Add Python to PATH"**.
2. Double-click **run.bat**. The first launch takes about a minute; after that it opens right away.
   (Or run **build_exe.bat** once to get a standalone `dist\SPLMonitor.exe` for the booth PC.)

## Choosing the input
- **Scarlett:** the app picks it automatically. Turn on 48V if your mic needs it, turn **Air** off, and choose Ch 1 or Ch 2.
- **Waves SoundGrid:** in your SoundGrid host app (SuperRack / SoundGrid Studio), patch the mic's I/O channel to a **SoundGrid Driver** channel. In SPL Monitor, choose **Waves SoundGrid [ASIO]** and pick that driver channel number.
- Click **Start**.

## Calibrate (first time, and after any gain change)
Set the preamp gain so a loud service won't clip. Put the calibrator on the mic, click **Calibrate…**, enter 94 or 114 and wait 3 seconds.
⚠️ After that, leave the gain alone. On SoundGrid that includes the preamp gain in the host app.
The calibration is saved per PC. If you switch between the Scarlett and SoundGrid, recalibrate after switching.

## Red level and thumbs
- **Red at** (in the toolbar) sets your limit. At or above it, the reading turns red and the thumb flips down 👎. Below it, you get a green thumbs up 👍.
- Yellow starts 5 dB below red. You can change that gap, or turn yellow off, in Settings.
- The thumbs-down holds for 3 seconds so it doesn't flicker. In Settings you can have the thumb and color follow **LAeq 1 min** instead of the live meter, which gives a steadier "how are we doing" read.

## Mini mode
Click **Mini ⤡** to shrink to a small box with the thumb and the live number. Drag it anywhere.
📌 keeps it on top of other windows. ⤢ (or a double-click) goes back to the full view.

## Logging
Press **Reset** at the start of the service, then **● Log**. Press it again to stop.
Logs go to `Documents\SPL Logs\SPL_<date>_<time>.csv` and open directly in Excel. The interval is set in Settings: 1 s, 10 s or 1 min.

| Column | Meaning |
|---|---|
| LAeq / LCeq | Average A / C level over the interval |
| LAFmax | Loudest moment (A, Fast) in the interval |
| LCpeak | Highest C-weighted peak in the interval |
| LAS | Meter reading at the end of the interval |
| LAeq 15 min | Rolling 15-minute average |
| Calibration | Calibration level and date |

Settings are saved in `%APPDATA%\SPLMonitor\config.json`. A/C weighting is within IEC 61672 Class 1 tolerance up to 12.5 kHz at 48 kHz.
