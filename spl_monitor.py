"""
ServiceSPL - pro sound level meter for Mac and Windows.

Main meter + Leq / Ln statistics, calibration wizard with per-mic profiles (UMIK-style cal files),
1/3-octave RTA + spectrogram, OSHA/NIOSH dosimeter, theater calibration (pink noise generator),
session logging with markers and alarms (CSV + PDF report), mini meter, menu-bar meter, web view.

Requires: numpy, scipy, PySide6, sounddevice
Optional on Mac: pyobjc-framework-Cocoa (live number in the menu bar)
"""

import sys, os, re, json, time, csv, queue, math, datetime, socket, threading, base64, collections, array, copy
import struct, subprocess, ctypes, urllib.request, urllib.parse
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
import numpy as np
from scipy import signal
from PySide6 import QtCore, QtGui, QtWidgets, QtSvg

os.environ.setdefault("SD_ENABLE_ASIO", "1")  # show ASIO drivers (Waves SoundGrid, Focusrite ASIO)
try:
    import sounddevice as sd
except Exception:  # PortAudio missing -> simulator only
    sd = None

APP_NAME = "ServiceSPL"
APP_VERSION = "2.0"
if sys.platform == "darwin":
    _CFG_BASE = os.path.expanduser("~/Library/Application Support")
else:
    _CFG_BASE = os.environ.get("APPDATA", os.path.expanduser("~"))
CFG_DIR = os.path.join(_CFG_BASE, "ServiceSPL")
CFG_PATH = os.path.join(CFG_DIR, "config.json")
OLD_CFG_PATH = os.path.join(_CFG_BASE, "SPLMonitor", "config.json")   # settings from the "SPL Monitor" days
SIM_NAME = "Simulated pink noise (test)"
SIM_CAL_NAME = "Simulated 94 dB calibrator (test)"
SIM_OUT = "Simulated (no audio out)"

ALARM_METRICS = ["LAS (live)", "LAeq 1 min", "LAeq long", "LCeq long", "LCpeak", "LAFmax",
                 "OSHA HC dose %", "NIOSH dose %"]

DEFAULTS = dict(
    device=None, channel=1, warn_gap=5.0, alarm=100.0, main_weight="A", main_time="S",
    long_leq_min=15, log_interval=1, thumb_src="live", web_port=8080, mini_on_top=True, mini_pos=None,
    log_dir=os.path.join(os.path.expanduser("~"), "Documents", "SPL Logs"),
    profiles={"Default": dict(offset=120.0, calibrated=False, cal_info="", method="", date="",
                              cal_points=None, cal_file="", use_curve=True, sens_factor=None,
                              serial="", device_hint="", auto_select=False)},
    profile="Default",
    theme="auto", meter_style="history", history_secs=300,
    menubar=True, menubar_style="dot",
    alarm_beep=True,
    alarm_rules=[dict(on=True, metric="LAeq long", thr=100.0),
                 dict(on=True, metric="LCpeak", thr=135.0),
                 dict(on=False, metric="LCeq long", thr=115.0)],
    rta_weight="Z", rta_speed="F", rta_spectro=True, spec_top=100.0, spec_range=70.0,
    theater=dict(layout="5.1", preset="Cinema / dub stage", noise="SMPTE band-limited (500 Hz - 2 kHz)",
                 level=-20.0, out_device=None, outputs={}, targets={}),
    session_view=600, venue="", operator="",
    companion_on=False, companion_host="127.0.0.1", companion_port=12321, companion_button="",
    ndi_name="ServiceSPL", ndi_size="1280x720", ndi_transparent=False, ndi_lib="",
    autosave_dir="", autosave_reports=True, schedule=[],
    pco_app_id="", pco_target=95.0, pco_service_type=None,
    win_geom=None, dock_state=None,
)

# ---------------------------------------------------------------- DSP ----

def weighting_sos(kind, fs):
    """IEC 61672 A / C weighting as a digital SOS filter (normalised to 0 dB @ 1 kHz)."""
    if kind == "Z":
        return None
    f1, f2, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217
    w1, w2, w3 = (2 * np.pi * f for f in (f1, f2, f3))
    w4 = 2 * fs * np.tan(np.pi * f4 / fs)  # pre-warped HF pole -> within Class 1 tolerance to 12.5 kHz
    den = np.polymul([1, 2 * w4, w4 ** 2], [1, 2 * w1, w1 ** 2])
    if kind == "A":
        num = [w4 ** 2, 0, 0, 0, 0]
        den = np.polymul(np.polymul(den, [1, w3]), [1, w2])
    else:  # C
        num = [w4 ** 2, 0, 0]
    b, a = signal.bilinear(num, den, fs)
    sos = signal.tf2sos(b, a)
    _, h = signal.sosfreqz(sos, worN=[1000.0], fs=fs)
    sos[0, :3] /= abs(h[0])
    return sos


def weighting_db(kind, f):
    """Analog IEC 61672 weighting in dB (used to weight RTA bands for display)."""
    f = np.asarray(f, dtype=float)
    if kind == "Z":
        return np.zeros_like(f)
    f1, f2, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217

    def c(f):
        return 20 * np.log10(f4 ** 2 * f ** 2 / ((f ** 2 + f1 ** 2) * (f ** 2 + f4 ** 2)))

    def a(f):
        return 20 * np.log10(f4 ** 2 * f ** 4 / ((f ** 2 + f1 ** 2) * np.sqrt(f ** 2 + f2 ** 2)
                                                * np.sqrt(f ** 2 + f3 ** 2) * (f ** 2 + f4 ** 2)))
    g = a if kind == "A" else c
    return g(f) - g(np.array(1000.0))


class TruePeak:
    """4x oversampled peak detector so LCpeak catches inter-sample peaks."""
    def __init__(self, up=4, half_taps=12):
        self.up = up
        self.h = signal.firwin(2 * up * half_taps + 1, 0.9 / up) * up
        self.keep = len(self.h) // up + 1
        self.tail = np.zeros(self.keep)

    def peak_sq(self, x):
        xx = np.concatenate([self.tail, x])
        self.tail = xx[-self.keep:]
        y = signal.upfirdn(self.h, xx, up=self.up)[len(self.h):len(xx) * self.up]
        return float(np.max(y * y)) if len(y) else 0.0


def parse_cal_file(path):
    """Read a UMIK / REW-style mic calibration file.
    Header like:  "Sens Factor =-0.7716dB, AGain =18dB, SERNO: 7023270"
    Data lines:   <freq Hz> <gain dB> [phase]  (space, tab or comma separated)."""
    pts, sens, serial, again = [], None, "", None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip().strip('"').strip()
            m = re.search(r"Sens\s*Factor\s*=\s*([-+]?\d+(?:\.\d+)?)", s, re.I)
            if m:
                sens = float(m.group(1))
            m = re.search(r"SERNO\s*:\s*([\w-]+)", s, re.I)
            if m:
                serial = m.group(1)
            m = re.search(r"AGain\s*=\s*([-+]?\d+(?:\.\d+)?)", s, re.I)
            if m:
                again = float(m.group(1))
            if s and (s[0].isdigit() or (s[0] in "+-." and len(s) > 1 and (s[1].isdigit() or s[1] == "."))):
                tok = [t for t in re.split(r"[\s,;]+", s) if t]
                try:
                    fr, g = float(tok[0]), float(tok[1])
                except (ValueError, IndexError):
                    continue
                if fr > 0 and (not pts or fr > pts[-1][0]):
                    pts.append([fr, g])
    if len(pts) < 2:
        raise ValueError("No frequency / dB pairs found in that file.")
    return dict(points=pts, sens_factor=sens, serial=serial, again=again, name=os.path.basename(path))


class CalFilter:
    """Linear-phase FIR that applies the inverse of a mic calibration curve."""
    def __init__(self, points, fs):
        f = np.array([p[0] for p in points], float)
        g = np.array([p[1] for p in points], float)
        n = int(fs / 12) | 1
        nf = 2 ** int(math.ceil(math.log2(n))) + 1
        grid = np.linspace(0, fs / 2, nf)
        corr = -np.interp(np.log10(np.clip(grid, f[0], f[-1])), np.log10(f), g)
        corr = np.clip(corr, -20, 20)
        self.h = signal.firwin2(n, grid, 10 ** (corr / 20), fs=fs, nfreqs=nf)
        self.hist = np.zeros(n - 1)

    def process(self, x):
        xx = np.concatenate([self.hist, x])
        self.hist = xx[-(len(self.h) - 1):]
        return signal.fftconvolve(xx, self.h, mode="valid")


RTA_LABELS = ["20", "25", "31.5", "40", "50", "63", "80", "100", "125", "160", "200", "250", "315", "400",
              "500", "630", "800", "1k", "1.25k", "1.6k", "2k", "2.5k", "3.15k", "4k", "5k", "6.3k", "8k",
              "10k", "12.5k", "16k", "20k"]


class ThirdOctaveBank:
    """IEC 61260 style 1/3-octave bank (base-10, 6th-order Butterworth bands), 20 Hz - 20 kHz."""
    def __init__(self, fs):
        self.fs = fs
        self.fc, self.labels, self.sos = [], [], []
        for i, n in enumerate(range(-17, 14)):
            fc = 1000.0 * 10 ** (n / 10)
            lo, hi = fc * 10 ** (-1 / 20), fc * 10 ** (1 / 20)
            if hi >= 0.49 * fs:
                break
            self.fc.append(fc)
            self.labels.append(RTA_LABELS[i])
            self.sos.append(signal.butter(3, [lo, hi], btype="bandpass", fs=fs, output="sos"))
        self.zi = [np.zeros((s.shape[0], 2)) for s in self.sos]
        nb = len(self.fc)
        self.fast = np.full(nb, 1e-20)
        self.slow = np.full(nb, 1e-20)
        self.reset_leq()

    def reset_leq(self):
        self.leq_sum = np.zeros(len(self.fc))
        self.leq_n = 0

    def process(self, x):
        n = len(x)
        dt = n / self.fs
        af, a_s = 1 - math.exp(-dt / 0.125), 1 - math.exp(-dt / 1.0)
        for i, sos in enumerate(self.sos):
            y, self.zi[i] = signal.sosfilt(sos, x, zi=self.zi[i])
            e = float(np.dot(y, y))
            self.leq_sum[i] += e
            e /= n
            self.fast[i] += af * (e - self.fast[i])
            self.slow[i] += a_s * (e - self.slow[i])
        self.leq_n += n

    def levels(self, speed):
        if speed == "Leq":
            return self.leq_sum / max(self.leq_n, 1)
        return self.fast if speed == "F" else self.slow


class Dosimeter:
    """Noise dose per OSHA / NIOSH style criteria (criterion level, exchange rate, threshold)."""
    def __init__(self, name, criterion, exchange, threshold, desc):
        self.name, self.criterion, self.exchange, self.threshold, self.desc = name, criterion, exchange, threshold, desc
        self.reset()

    def reset(self):
        self.dose = 0.0          # fraction, 1.0 = 100 %
        self.elapsed = 0.0
        self.peak_over = False

    def add(self, level, dt):
        self.elapsed += dt
        if np.isfinite(level) and level >= self.threshold:
            t_allowed = 8 * 3600.0 / 2 ** ((level - self.criterion) / self.exchange)
            self.dose += dt / t_allowed

    def twa(self):
        if self.dose <= 0:
            return float("nan")
        return self.criterion + self.exchange / math.log10(2) * math.log10(self.dose)

    def projected(self):
        """Projected 8-hour dose in percent at the average rate so far."""
        return self.dose * 100 * 8 * 3600 / self.elapsed if self.elapsed > 60 else float("nan")

    def time_to_full(self):
        if self.dose >= 1:
            return 0.0
        rate = self.dose / self.elapsed if self.elapsed else 0
        return (1 - self.dose) / rate if rate > 0 else float("inf")


def make_dosimeters():
    return [Dosimeter("OSHA PEL", 90, 5, 90, "90 dBA criterion · 5 dB exchange · 90 dB threshold · Slow"),
            Dosimeter("OSHA HC", 90, 5, 80, "Hearing conservation / action level · 80 dB threshold · Slow"),
            Dosimeter("NIOSH", 85, 3, 80, "85 dBA REL · 3 dB exchange · 80 dB threshold")]


class StatLevels:
    """0.1 dB histogram of LAF samples -> L10 / L50 / L90 since reset."""
    LO, HI, STEP = -60.0, 200.0, 0.1

    def __init__(self):
        self.hist = np.zeros(int((self.HI - self.LO) / self.STEP) + 1, np.int64)
        self.n = 0

    def add(self, levels):
        if len(levels):
            idx = np.clip(((np.asarray(levels) - self.LO) / self.STEP).astype(int), 0, len(self.hist) - 1)
            np.add.at(self.hist, idx, 1)
            self.n += len(levels)

    def ln(self, pct):
        if not self.n:
            return float("nan")
        c = np.cumsum(self.hist[::-1])
        i = int(np.searchsorted(c, pct / 100.0 * self.n))
        return self.HI - i * self.STEP


def ln_of(samples, pct):
    """Level exceeded pct % of the time."""
    return float(np.percentile(samples, 100 - pct)) if len(samples) else float("nan")


class Engine:
    """All metering maths. Feed it mono float blocks via process()."""
    HIST_MAX = 7200          # 2 h of per-second data
    LN_RATE = 10             # LAF samples per second for L10 / L90

    def __init__(self, fs, offset, long_min, cal_points=None):
        self.fs = fs
        self.offset = offset
        self.long_min = long_min
        self.calfilt = CalFilter(cal_points, fs) if cal_points else None
        self.sos = {k: weighting_sos(k, fs) for k in "ACZ"}
        self.zi = {k: (np.zeros((s.shape[0], 2)) if s is not None else None)
                   for k, s in self.sos.items()}
        self.alpha = {"F": 1 - math.exp(-1 / (0.125 * fs)), "S": 1 - math.exp(-1 / (1.0 * fs))}
        self.tz = {(k, t): np.zeros(1) for k in "ACZ" for t in "FS"}
        self.cur = {(k, t): 1e-20 for k in "ACZ" for t in "FS"}
        self.hist_sum = {k: [] for k in "ACZ"}
        self.hist_n = []
        self.laf_hist = collections.deque(maxlen=self.HIST_MAX)
        self.records = []
        self.tp = TruePeak()
        self.ln_step = max(1, fs // self.LN_RATE)
        self.samples = 0
        self.rta = None          # ThirdOctaveBank when the RTA panel is visible
        self.extra = []          # more banks (room EQ measurements)
        self.tap = None          # list of corrected blocks for the spectrogram
        self.dose = make_dosimeters()
        self.cal_target = 0
        self.cal_sum = 0.0
        self.cal_n = 0
        self.cal_weight = "Z"
        self.reset_stats()
        self._new_second()

    # --- helpers
    def L(self, ms):
        return 10 * math.log10(max(ms, 1e-20)) + self.offset

    def reset_stats(self):
        self.tot_sum = {k: 0.0 for k in "ACZ"}
        self.tot_n = 0
        self.max_af = 1e-20
        self.max_pk = 1e-20
        self.stat = StatLevels()
        self.reset_time = time.time()

    def reset_dose(self):
        for d in self.dose:
            d.reset()

    def _new_second(self):
        self.sec_sum = {k: 0.0 for k in "ACZ"}
        self.sec_n = 0
        self.sec_pk = 1e-20
        self.sec_maxaf = 1e-20
        self.sec_laf = []

    # --- main processing
    def process(self, x):
        x = np.asarray(x, dtype=np.float64)
        n = len(x)
        if n == 0:
            return
        if self.calfilt is not None:
            x = self.calfilt.process(x)
        if self.tap is not None:
            self.tap.append(x)
            del self.tap[:-400]
        if self.rta is not None:
            self.rta.process(x)
        for b in self.extra:
            b.process(x)
        blk_a = 0.0
        for k in "ACZ":
            if self.sos[k] is None:
                xw = x
            else:
                xw, self.zi[k] = signal.sosfilt(self.sos[k], x, zi=self.zi[k])
            sq = xw * xw
            for t in "FS":
                a = self.alpha[t]
                y, self.tz[(k, t)] = signal.lfilter([a], [1.0, a - 1.0], sq, zi=self.tz[(k, t)])
                self.cur[(k, t)] = float(y[-1])
                if k == "A" and t == "F":
                    self.sec_maxaf = max(self.sec_maxaf, float(y.max()))
                    self.max_af = max(self.max_af, self.sec_maxaf)
                    i0 = (-self.samples) % self.ln_step
                    picks = y[i0::self.ln_step]
                    if len(picks):
                        self.sec_laf.extend((10 * np.log10(np.maximum(picks, 1e-20)) + self.offset).tolist())
            s = float(sq.sum())
            self.sec_sum[k] += s
            self.tot_sum[k] += s
            if k == "A":
                blk_a = s / n
            if k == "C":
                self.sec_pk = max(self.sec_pk, float(sq.max()), self.tp.peak_sq(xw))
            if k == self.cal_weight and self.cal_target:
                self.cal_sum += s
                self.cal_n += n
        dt = n / self.fs
        las, laeq = self.L(self.cur[("A", "S")]), self.L(blk_a)
        self.dose[0].add(las, dt)
        self.dose[1].add(las, dt)
        self.dose[2].add(laeq, dt)
        self.samples += n
        self.sec_n += n
        self.tot_n += n
        self.max_pk = max(self.max_pk, self.sec_pk)
        if self.L(self.sec_pk) >= 140:
            for d in self.dose:
                d.peak_over = True
        if self.sec_n >= self.fs:
            self._close_second()

    def _close_second(self):
        for k in "ACZ":
            self.hist_sum[k].append(self.sec_sum[k])
            del self.hist_sum[k][:-self.HIST_MAX]
        self.hist_n.append(self.sec_n)
        del self.hist_n[:-self.HIST_MAX]
        laf = np.array(self.sec_laf, np.float32)
        self.laf_hist.append(laf)
        self.stat.add(laf)
        L, c, sn = self.L, self.cur, self.sec_n
        lw = int(self.long_min * 60)
        win = self.ln_window(lw)
        self.records.append({
            "time": time.time(),
            "LAS": L(c["A", "S"]), "LAF": L(c["A", "F"]),
            "LCS": L(c["C", "S"]), "LZS": L(c["Z", "S"]),
            "LAeq1s": L(self.sec_sum["A"] / sn), "LCeq1s": L(self.sec_sum["C"] / sn),
            "LZeq1s": L(self.sec_sum["Z"] / sn),
            "LAeq1m": self.leq("A", 60), "LCeq1m": self.leq("C", 60),
            "LAeqLong": self.leq("A", lw), "LCeqLong": self.leq("C", lw),
            "LCpeak1s": L(self.sec_pk), "LAFmax1s": L(self.sec_maxaf),
            "LAF10": laf,
            "L10": ln_of(win, 10), "L50": ln_of(win, 50), "L90": ln_of(win, 90),
            "doseOSHA": self.dose[1].dose * 100, "doseNIOSH": self.dose[2].dose * 100,
            "dosePEL": self.dose[0].dose * 100,
        })
        self._new_second()

    def ln_window(self, secs):
        n = min(secs, len(self.laf_hist))
        if not n:
            return np.zeros(0)
        return np.concatenate([self.laf_hist[-i] for i in range(n, 0, -1)])

    def leq(self, k, secs):
        if not self.hist_n:
            return float("nan")
        return self.L(sum(self.hist_sum[k][-secs:]) / sum(self.hist_n[-secs:]))

    def level(self, k, t):
        return self.L(self.cur[(k, t)])

    def leq_total(self, k):
        return self.L(self.tot_sum[k] / self.tot_n) if self.tot_n else float("nan")

    # --- calibration
    def begin_cal(self, seconds=5.0, weight="Z"):
        self.cal_sum, self.cal_n = 0.0, 0
        self.cal_weight = weight
        self.cal_target = int(seconds * self.fs)

    @property
    def cal_done(self):
        return bool(self.cal_target) and self.cal_n >= self.cal_target

    @property
    def cal_progress(self):
        return min(1.0, self.cal_n / self.cal_target) if self.cal_target else 0.0

    def finish_cal(self):
        """Returns the raw (un-offset) level of the calibration period."""
        ms = self.cal_sum / max(self.cal_n, 1)
        self.cal_target = 0
        return 10 * math.log10(max(ms, 1e-20))


def make_pink(fs, lo, hi, rms_dbfs, seconds=10.0, seed=1):
    """Seamless loop of band-limited pink noise. Level is AES17 dBFS RMS (full-scale sine = 0 dBFS)."""
    n = 2 ** int(round(math.log2(fs * seconds)))
    rng = np.random.default_rng(seed)
    f = np.fft.rfftfreq(n, 1 / fs)
    spec = rng.standard_normal(len(f)) + 1j * rng.standard_normal(len(f))
    mag = np.zeros(len(f))
    band = (f >= lo) & (f <= hi)
    mag[band] = 1 / np.sqrt(f[band])
    x = np.fft.irfft(spec * mag, n)
    x *= (10 ** (rms_dbfs / 20) / math.sqrt(2)) / np.sqrt(np.mean(x * x))
    return x.astype(np.float32)


NOISE_TYPES = {"SMPTE band-limited (500 Hz - 2 kHz)": (500, 2000),
               "Full-band pink (20 Hz - 20 kHz)": (20, 20000),
               "LFE band (20 - 120 Hz)": (20, 120)}


class NoisePlayer:
    """Plays a looping noise buffer into one output channel."""
    def __init__(self):
        self.stream = None
        self.buf = None
        self.pos = 0
        self.ch = 0
        self.gain = 1.0
        self.lock = threading.Lock()

    @property
    def playing(self):
        return self.buf is not None

    def set_buffer(self, buf):
        with self.lock:
            self.buf, self.pos = buf, 0

    def read(self, n):
        with self.lock:
            if self.buf is None:
                return np.zeros(n, np.float32)
            idx = (self.pos + np.arange(n)) % len(self.buf)
            self.pos = (self.pos + n) % len(self.buf)
            return self.buf[idx] * self.gain

    def open(self, device, ch):
        """device None = no audio out (simulator reads the noise directly)."""
        self.close()
        self.ch = ch
        if device is None or sd is None:
            return 48000
        info = sd.query_devices(device)
        fs = int(info["default_samplerate"])
        api = sd.query_hostapis(info["hostapi"])["name"]

        def cb(out, frames, t, status):
            out.fill(0)
            col = 0 if "ASIO" in api else self.ch
            out[:, col] = self.read(frames)

        if "ASIO" in api:
            self.stream = sd.OutputStream(device=device, channels=1, samplerate=fs, dtype="float32",
                                          callback=cb, extra_settings=sd.AsioSettings(channel_selectors=[ch]))
        else:
            if ch + 1 > info["max_output_channels"]:
                raise RuntimeError(f"That output only has {info['max_output_channels']} channel(s).")
            self.stream = sd.OutputStream(device=device, channels=info["max_output_channels"], samplerate=fs,
                                          dtype="float32", callback=cb)
        self.stream.start()
        return fs

    def close(self):
        if self.stream is not None:
            try:
                self.stream.stop(); self.stream.close()
            except Exception:
                pass
        self.stream = None

    def stop(self):
        self.set_buffer(None)
        self.close()


class SimSource:
    """Test input paced by the wall clock: pink noise with a slow swell, a 94 dB calibrator tone,
    or the theater noise generator routed through a fake room (per-channel speaker gains)."""
    B = [0.049922035, -0.095993537, 0.050612699, -0.004408786]
    A = [1, -2.494956002, 2.017265875, -0.522189400]
    ROOM = [0.0, -1.2, 0.8, 4.0, -2.5, -1.8, 1.5, -3.0, 0.4, -0.9, 1.1, -1.6]   # dB per output channel

    def __init__(self, fs, mode="pink"):
        self.fs = fs
        self.mode = mode
        self.t0 = time.monotonic()
        self.done = 0
        self.zi = np.zeros(3)
        self.inject = None           # NoisePlayer when the theater generator is simulated

    def read(self):
        target = int((time.monotonic() - self.t0) * self.fs)
        n = min(target - self.done, self.fs)
        if n <= 0:
            return []
        t = (self.done + np.arange(n)) / self.fs
        self.done = target
        if self.inject is not None and self.inject.playing:
            g = 10 ** ((self.ROOM[self.inject.ch % len(self.ROOM)] - 12) / 20)
            return [self.inject.read(n).astype(np.float64) * g + 1e-4 * np.random.randn(n)]
        if self.mode == "cal":   # 1 kHz tone that reads 94.0 dB with the default 120 dB offset
            return [math.sqrt(2) * 10 ** ((94 - 120) / 20) * np.sin(2 * np.pi * 1000 * t)
                    + 1e-5 * np.random.randn(n)]
        p, self.zi = signal.lfilter(self.B, self.A, np.random.randn(n), zi=self.zi)
        env = 0.05 * (1 + 0.8 * np.sin(2 * np.pi * t / 20))
        return [p * env]


# --------------------------------------------------------------- theme ----

PALETTES = {
    "dark": dict(bg="#121417", panel="#1b1e23", border="#2a2e35", btn="#2a2e35", hover="#353a42",
                 btn_border="#3a3f47", field="#22262c", fg="#e8e8e8", dim="#8a9099", track="#0c0d0f",
                 grid="#262a31", green="#3ddc84", yellow="#ffcc33", red="#ff4d4d", accent="#4da3ff",
                 checked="#1f6f43", rec="#7a1f1f", banner="#5c1a1a", off="#3a3f47"),
    "light": dict(bg="#eef0f3", panel="#ffffff", border="#d3d7dd", btn="#e6e9ed", hover="#d9dde3",
                  btn_border="#c2c8cf", field="#ffffff", fg="#1b1e23", dim="#5f6670", track="#e1e4e8",
                  grid="#e3e6ea", green="#14a457", yellow="#c99300", red="#dc3030", accent="#1f6fd1",
                  checked="#c9ecd8", rec="#f5cccc", banner="#fbe0e0", off="#b9bfc7"),
}
T = dict(PALETTES["dark"])


def style_sheet():
    return f"""
QWidget {{ background:{T['bg']}; color:{T['fg']}; font-family:'Segoe UI','Helvetica Neue',Arial; font-size:10pt; }}
QFrame#tile {{ background:{T['panel']}; border:1px solid {T['border']}; border-radius:6px; }}
QFrame#tile QLabel {{ background:transparent; }}
QPushButton, QToolButton {{ background:{T['btn']}; border:1px solid {T['btn_border']}; border-radius:4px; padding:5px 12px; }}
QPushButton:hover, QToolButton:hover {{ background:{T['hover']}; }}
QPushButton:checked {{ background:{T['checked']}; border-color:{T['green']}; }}
QPushButton#log:checked {{ background:{T['rec']}; border-color:{T['red']}; }}
QPushButton:disabled {{ color:{T['dim']}; }}
QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit, QPlainTextEdit {{ background:{T['field']}; border:1px solid {T['btn_border']};
    border-radius:4px; padding:3px 6px; }}
QTableWidget, QListWidget {{ background:{T['panel']}; border:1px solid {T['border']}; gridline-color:{T['border']}; }}
QHeaderView::section {{ background:{T['btn']}; color:{T['dim']}; border:none; padding:3px; }}
QDockWidget {{ color:{T['dim']}; font-weight:600; }}
QDockWidget::title {{ background:{T['panel']}; padding:5px 8px; border-bottom:1px solid {T['border']}; }}
QProgressBar {{ background:{T['track']}; border:1px solid {T['border']}; border-radius:4px; text-align:center; }}
QStatusBar {{ color:{T['dim']}; }}
QMenu {{ background:{T['panel']}; border:1px solid {T['border']}; }}
QMenu::item:selected {{ background:{T['hover']}; }}
QToolTip {{ background:{T['panel']}; color:{T['fg']}; border:1px solid {T['border']}; }}
"""


def system_is_dark():
    try:
        return QtGui.QGuiApplication.styleHints().colorScheme() != QtCore.Qt.ColorScheme.Light
    except Exception:
        return True


LOG_INTERVALS = {"1 second": 1, "10 seconds": 10, "1 minute": 60}
MONO = "font-family:Menlo,Consolas,'Courier New',monospace;"


def fmt(v, nd=1):
    return "--.-" if v is None or not np.isfinite(v) else f"{v:5.{nd}f}"


def hms(secs):
    secs = int(max(secs, 0))
    return f"{secs // 3600}:{secs // 60 % 60:02d}:{secs % 60:02d}"


def energy_mean(vals):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    return 10 * math.log10(np.mean(10 ** (v / 10))) if len(v) else float("nan")


# ------------------------------------------------------------- widgets ----

class Tile(QtWidgets.QFrame):
    def __init__(self, title, size=22, tip=""):
        super().__init__()
        self.setObjectName("tile")
        if tip:
            self.setToolTip(tip)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 6, 10, 6)
        lay.setSpacing(0)
        self.t = QtWidgets.QLabel(title)
        self.v = QtWidgets.QLabel()
        self.v.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.v.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Preferred)
        self.t.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Preferred)
        self.s = QtWidgets.QLabel("")
        self.s.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Preferred)
        self.s.setAlignment(QtCore.Qt.AlignRight)
        self.size = size
        self.color = None
        for w in (self.t, self.v, self.s):
            lay.addWidget(w)
        self.refresh_theme()
        self.set(None)

    def refresh_theme(self):
        self.t.setStyleSheet(f"color:{T['dim']}; font-size:9pt; background:transparent;")
        self.s.setStyleSheet(f"color:{T['dim']}; font-size:8pt; background:transparent;")
        self._paint_value()

    def _paint_value(self):
        self.v.setStyleSheet(f"color:{self.color or T['fg']}; font-size:{self.size}pt; font-weight:600;"
                             f"{MONO} background:transparent;")

    def set_size(self, pt):
        pt = int(pt)
        if pt != self.size:
            self.size = pt
            self._paint_value()

    def set(self, val, color=None, sub=None, title=None, text=None):
        if text is None and (val is None or not np.isfinite(val)):
            color = T["dim"]          # stopped / no data yet: a quiet grey "--.-"
        self.v.setText(text if text is not None else fmt(val))
        if color != self.color:
            self.color = color
            self._paint_value()
        if sub is not None:
            self.s.setText(sub)
        if title:
            self.t.setText(title)


class LevelBar(QtWidgets.QWidget):
    """Horizontal bar meter with 2-second peak hold."""
    def __init__(self, lo=30, hi=130):
        super().__init__()
        self.lo, self.hi = lo, hi
        self.level = self.peak = lo
        self.peak_t = 0
        self.warn, self.alarm = 95, 100
        self.setMinimumHeight(40)

    def setLevel(self, v):
        v = v if np.isfinite(v) else self.lo
        self.level = v
        now = time.time()
        if v >= self.peak or now - self.peak_t > 2:
            self.peak, self.peak_t = v, now
        self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        r = self.rect().adjusted(8, 4, -8, -16)
        p.fillRect(r, QtGui.QColor(T["track"]))
        frac = lambda v: (min(max(v, self.lo), self.hi) - self.lo) / (self.hi - self.lo)
        x = lambda v: r.left() + frac(v) * r.width()
        g = QtGui.QLinearGradient(r.left(), 0, r.right(), 0)
        fw, fa = frac(self.warn), frac(self.alarm)
        for pos, col in ((0, T["green"]), (max(fw - 1e-3, 0), T["green"]), (fw, T["yellow"]),
                         (max(fa - 1e-3, fw), T["yellow"]), (fa, T["red"]), (1, T["red"])):
            g.setColorAt(pos, QtGui.QColor(col))
        p.fillRect(QtCore.QRectF(r.left(), r.top(), x(self.level) - r.left(), r.height()), g)
        p.setPen(QtGui.QPen(QtGui.QColor(T["fg"]), 2))
        p.drawLine(QtCore.QPointF(x(self.peak), r.top()), QtCore.QPointF(x(self.peak), r.bottom()))
        p.setPen(QtGui.QColor(T["dim"]))
        f = p.font(); f.setPointSize(8); p.setFont(f)
        for v in range(self.lo, self.hi + 1, 10):
            xx = x(v)
            p.drawLine(QtCore.QPointF(xx, r.bottom()), QtCore.QPointF(xx, r.bottom() + 3))
            p.drawText(QtCore.QRectF(xx - 20, r.bottom() + 2, 40, 14), QtCore.Qt.AlignCenter, str(v))


class HistoryStrip(QtWidgets.QWidget):
    """Scrolling level history (default 5 minutes) under the big number."""
    def __init__(self, secs=300):
        super().__init__()
        self.secs = secs
        self.pts = collections.deque()
        self.warn, self.alarm = 95, 100
        self.lo, self.hi = 40, 120
        self.leq_pts = collections.deque()
        self._last = 0
        self.setMinimumHeight(50)

    def add(self, v, leq1m=None):
        now = time.time()
        if now - self._last < 0.1:
            return
        self._last = now
        self.pts.append((now, v if np.isfinite(v) else float("nan")))
        if leq1m is not None:
            self.leq_pts.append((now, leq1m))
        while self.pts and self.pts[0][0] < now - self.secs:
            self.pts.popleft()
        while self.leq_pts and self.leq_pts[0][0] < now - self.secs:
            self.leq_pts.popleft()
        self.update()

    def clear(self):
        self.pts.clear(); self.leq_pts.clear(); self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        r = QtCore.QRectF(self.rect().adjusted(34, 6, -8, -16))
        p.fillRect(r, QtGui.QColor(T["track"]))
        vals = [v for _, v in self.pts if np.isfinite(v)]
        if vals:   # auto-range in 10 dB steps, always showing the red line
            lo = min(min(vals), self.warn) - 5
            hi = max(max(vals), self.alarm) + 3
            self.lo, self.hi = 10 * math.floor(lo / 10), 10 * math.ceil(hi / 10)
        span = max(self.hi - self.lo, 10)
        y = lambda v: r.bottom() - (min(max(v, self.lo), self.hi) - self.lo) / span * r.height()
        now = time.time()
        x = lambda t: r.right() - (now - t) / self.secs * r.width()
        f = p.font(); f.setPointSize(8); p.setFont(f)
        for v in range(int(self.lo), int(self.hi) + 1, 10):
            p.setPen(QtGui.QPen(QtGui.QColor(T["grid"]), 1))
            p.drawLine(QtCore.QPointF(r.left(), y(v)), QtCore.QPointF(r.right(), y(v)))
            p.setPen(QtGui.QColor(T["dim"]))
            p.drawText(QtCore.QRectF(0, y(v) - 7, 30, 14), QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, str(v))
        for m in range(0, int(self.secs) + 1, 60):
            xx = r.right() - m / self.secs * r.width()
            p.drawText(QtCore.QRectF(xx - 25, r.bottom() + 1, 50, 14), QtCore.Qt.AlignCenter,
                       "now" if m == 0 else f"-{m // 60}m")
        if len(self.pts) > 1:
            path = QtGui.QPainterPath()
            path.moveTo(x(self.pts[0][0]), r.bottom())
            for t, v in self.pts:
                path.lineTo(x(t), y(v if np.isfinite(v) else self.lo))
            path.lineTo(x(self.pts[-1][0]), r.bottom())
            g = QtGui.QLinearGradient(0, r.bottom(), 0, r.top())
            fw = (self.warn - self.lo) / span
            fa = (self.alarm - self.lo) / span
            for pos, col in ((0, T["green"]), (fw - 1e-3, T["green"]), (fw, T["yellow"]),
                             (fa - 1e-3, T["yellow"]), (fa, T["red"]), (1, T["red"])):
                c = QtGui.QColor(col); c.setAlpha(150)
                g.setColorAt(min(max(pos, 0), 1), c)
            p.fillPath(path, g)
        if len(self.leq_pts) > 1:
            p.setPen(QtGui.QPen(QtGui.QColor(T["fg"]), 1.5))
            pl = QtGui.QPolygonF([QtCore.QPointF(x(t), y(v)) for t, v in self.leq_pts if np.isfinite(v)])
            p.drawPolyline(pl)
        p.setPen(QtGui.QPen(QtGui.QColor(T["red"]), 1, QtCore.Qt.DashLine))
        p.drawLine(QtCore.QPointF(r.left(), y(self.alarm)), QtCore.QPointF(r.right(), y(self.alarm)))
        p.setPen(QtGui.QColor(T["dim"]))
        p.drawText(r.adjusted(6, 2, 0, 0), QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop, "— LAeq 1 min")


THUMB_SVG = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100">
<circle cx="50" cy="50" r="48" fill="{bg}"/>
<g transform="{xf}" fill="#ffffff">
  <rect x="20" y="47" width="13" height="33" rx="3"/>
  <path d="M37 50 L37 76 Q37 80 41 80 L68 80 Q74 80 75 74 L79 56 Q80 49 73 49 L58 49
           L60 33 Q61 22 53 21 Q49 21 48 26 L44 42 Q42 47 37 50 Z"/>
</g>
<g transform="{xf}" stroke="{bg}" stroke-width="2.5" stroke-linecap="round">
  <line x1="62" y1="58" x2="77" y2="58"/>
  <line x1="62" y1="66" x2="76" y2="66"/>
  <line x1="62" y1="73" x2="75" y2="73"/>
</g>
</svg>"""


class Thumb(QtWidgets.QWidget):
    """Green thumbs-up below the red level, red thumbs-down at/above it, grey when stopped."""
    def __init__(self, size=150):
        super().__init__()
        self.setFixedSize(size, size)
        self.r = {k: QtSvg.QSvgRenderer(QtCore.QByteArray(THUMB_SVG.format(bg=bg, xf=xf).encode()))
                  for k, bg, xf in (("good", "#2ea866", ""), ("bad", "#e03b3b", "rotate(180 50 50)"),
                                    ("off", "#8a9099", ""))}
        self.state = "off"

    def set_size(self, px):
        px = int(px)
        if px != self.width():
            self.setFixedSize(px, px)

    def setState(self, st):
        if st != self.state:
            self.state = st
            self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        self.r[self.state].render(p, QtCore.QRectF(self.rect()))


SPEC_LUT = None


def spec_lut():
    """256-entry dark-blue -> purple -> orange -> yellow colour map (ARGB32)."""
    global SPEC_LUT
    if SPEC_LUT is None:
        anchors = np.array([[0, 0, 4], [40, 11, 84], [101, 21, 110], [159, 42, 99], [212, 72, 66],
                            [245, 125, 21], [250, 193, 39], [252, 255, 164]], float)
        xs = np.linspace(0, 1, len(anchors))
        t = np.linspace(0, 1, 256)
        rgb = np.stack([np.interp(t, xs, anchors[:, i]) for i in range(3)], 1).astype(np.uint32)
        SPEC_LUT = (0xFF000000 | (rgb[:, 0] << 16) | (rgb[:, 1] << 8) | rgb[:, 2]).astype(np.uint32)
    return SPEC_LUT


class RtaWidget(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.labels, self.levels = [], np.zeros(0)
        self.peaks = np.zeros(0)
        self.lo, self.hi = 20, 120
        self.setMinimumHeight(90)

    def set_data(self, labels, levels):
        if len(labels) != len(self.labels):
            self.peaks = np.full(len(labels), -999.0)
        self.labels, self.levels = labels, levels
        self.peaks = np.maximum(self.peaks, levels)
        self.update()

    def reset_peaks(self):
        self.peaks = np.full(len(self.labels), -999.0)

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        r = QtCore.QRectF(self.rect().adjusted(32, 6, -6, -18))
        p.fillRect(r, QtGui.QColor(T["track"]))
        f = p.font(); f.setPointSize(8); p.setFont(f)
        if len(self.levels):
            hi = max(np.nanmax(self.levels) + 10, 60)
            self.hi = 10 * math.ceil(hi / 10)
            self.lo = self.hi - 90
        span = self.hi - self.lo
        y = lambda v: r.bottom() - (min(max(v, self.lo), self.hi) - self.lo) / span * r.height()
        for v in range(int(self.lo), int(self.hi) + 1, 10):
            p.setPen(QtGui.QColor(T["grid"]))
            p.drawLine(QtCore.QPointF(r.left(), y(v)), QtCore.QPointF(r.right(), y(v)))
            p.setPen(QtGui.QColor(T["dim"]))
            p.drawText(QtCore.QRectF(0, y(v) - 7, 28, 14), QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, str(v))
        n = len(self.levels)
        if not n:
            p.setPen(QtGui.QColor(T["dim"]))
            p.drawText(r, QtCore.Qt.AlignCenter, "Press Start to see the spectrum")
            return
        bw = r.width() / n
        g = QtGui.QLinearGradient(0, r.bottom(), 0, r.top())
        g.setColorAt(0, QtGui.QColor(T["accent"]).darker(130))
        g.setColorAt(1, QtGui.QColor(T["accent"]).lighter(130))
        for i, (lab, v) in enumerate(zip(self.labels, self.levels)):
            x0 = r.left() + i * bw
            p.fillRect(QtCore.QRectF(x0 + 1, y(v), bw - 2, r.bottom() - y(v)), g)
            p.setPen(QtGui.QPen(QtGui.QColor(T["fg"]), 2))
            p.drawLine(QtCore.QPointF(x0 + 1, y(self.peaks[i])), QtCore.QPointF(x0 + bw - 1, y(self.peaks[i])))
            if i % 3 == 2 or n <= 12:
                p.setPen(QtGui.QColor(T["dim"]))
                p.drawText(QtCore.QRectF(x0 - 12, r.bottom() + 2, bw + 24, 14), QtCore.Qt.AlignCenter, lab)


class SpectrogramWidget(QtWidgets.QWidget):
    """Scrolling log-frequency spectrogram (20 Hz - 20 kHz)."""
    W, H = 480, 160

    def __init__(self):
        super().__init__()
        self.img = np.zeros((self.H, self.W), np.uint8)
        self.top, self.range = 100.0, 70.0
        self.fs = None
        self.buf = np.zeros(0)
        self.setMinimumHeight(60)

    def setup(self, fs):
        self.fs = fs
        self.nfft = 8192 if fs <= 50000 else 16384
        self.hop = fs // 25
        self.win = np.hanning(self.nfft)
        self.norm = 2.0 / self.win.sum() ** 2
        f = np.fft.rfftfreq(self.nfft, 1 / fs)
        self.freqs = f
        top = min(20000, fs / 2 * 0.98)
        edges = np.geomspace(20, top, self.H + 1)
        self.lo_i = np.searchsorted(f, edges[:-1])
        self.hi_i = np.searchsorted(f, edges[1:])
        self.centre = np.sqrt(edges[:-1] * edges[1:])
        self.buf = np.zeros(0)

    def feed(self, blocks, offset):
        if not blocks or self.fs is None:
            return
        self.buf = np.concatenate([self.buf] + blocks)[-self.fs * 2:]
        cols = []
        while len(self.buf) >= self.nfft:
            seg = self.buf[:self.nfft]
            self.buf = self.buf[self.hop:]
            pw = np.abs(np.fft.rfft(seg * self.win)) ** 2 * self.norm
            cum = np.concatenate([[0], np.cumsum(pw)])
            cnt = self.hi_i - self.lo_i
            band = np.where(cnt > 0, (cum[self.hi_i] - cum[self.lo_i]) / np.maximum(cnt, 1),
                            np.interp(self.centre, self.freqs, pw))
            db = 10 * np.log10(np.maximum(band, 1e-20)) + offset
            cols.append(np.clip((db - (self.top - self.range)) / self.range * 255, 0, 255)[::-1])
        if cols:
            k = min(len(cols), self.W)
            self.img = np.roll(self.img, -k, axis=1)
            self.img[:, -k:] = np.array(cols[-k:]).T.astype(np.uint8)
            self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        r = self.rect().adjusted(32, 4, -6, -4)
        rgb = np.ascontiguousarray(spec_lut()[self.img])
        self._q = QtGui.QImage(rgb.data, self.W, self.H, self.W * 4, QtGui.QImage.Format_RGB32)
        self._keep = rgb
        p.setRenderHint(QtGui.QPainter.SmoothPixmapTransform)
        p.drawImage(QtCore.QRectF(r), self._q)
        p.setPen(QtGui.QColor(T["dim"]))
        f = p.font(); f.setPointSize(8); p.setFont(f)
        top = min(20000, (self.fs or 48000) / 2 * 0.98)
        for fr, lab in ((50, "50"), (100, "100"), (500, "500"), (1000, "1k"), (5000, "5k"), (10000, "10k")):
            yy = r.bottom() - math.log(fr / 20) / math.log(top / 20) * r.height()
            p.drawText(QtCore.QRectF(0, yy - 7, 28, 14), QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, lab)


# ------------------------------------------------------------- session ----

class Session:
    """Per-second history for the whole session (compact columns), plus markers and events."""
    COLS = ["time", "LAeq1s", "LCeq1s", "LZeq1s", "LAFmax1s", "LCpeak1s", "LAS", "LAeq1m", "LAeqLong",
            "LCeqLong", "doseOSHA", "doseNIOSH", "dosePEL"]
    MAX_ROWS = 7 * 24 * 3600

    def __init__(self):
        self.start = time.time()
        self.cols = {k: array.array("d") for k in self.COLS}
        self.laf = array.array("f")
        self.markers = []     # (time, text)
        self.events = []      # (time, text)
        self.items = []       # (start time, title): Planning Center / manual service items
        self.cal_info = ""

    def __len__(self):
        return len(self.cols["time"])

    def add(self, rec):
        if len(self) >= self.MAX_ROWS:
            return
        for k in self.COLS:
            self.cols[k].append(float(rec[k]))
        self.laf.extend(rec["LAF10"].tolist())

    def col(self, k):
        return np.frombuffer(self.cols[k], dtype=np.float64) if len(self) else np.zeros(0)

    def summary(self):
        if not len(self):
            return {}
        la, lc = self.col("LAeq1s"), self.col("LCeq1s")
        laf = np.frombuffer(self.laf, dtype=np.float32) if len(self.laf) else np.zeros(0)
        lw = self.col("LAeqLong")
        return dict(dur=self.col("time")[-1] - self.col("time")[0] + 1,
                    LAeq=energy_mean(la), LCeq=energy_mean(lc),
                    LAFmax=float(np.nanmax(self.col("LAFmax1s"))), LCpeak=float(np.nanmax(self.col("LCpeak1s"))),
                    L10=ln_of(laf, 10), L50=ln_of(laf, 50), L90=ln_of(laf, 90),
                    maxLAeq1m=float(np.nanmax(self.col("LAeq1m"))), maxLong=float(np.nanmax(lw)),
                    doseOSHA=self.col("doseOSHA")[-1], doseNIOSH=self.col("doseNIOSH")[-1],
                    dosePEL=self.col("dosePEL")[-1])


def paint_session_graph(p, r, sess, window, alarm, warn, small_font=8):
    """Draw the session time history into rect r (used on screen and in the PDF)."""
    r = QtCore.QRectF(r)
    p.save()
    p.setRenderHint(QtGui.QPainter.Antialiasing)
    p.fillRect(r, QtGui.QColor(T["track"]))
    f = p.font(); f.setPointSizeF(small_font); p.setFont(f)
    t = sess.col("time")
    if len(t) < 2:
        p.setPen(QtGui.QColor(T["dim"]))
        p.drawText(r, QtCore.Qt.AlignCenter, "Session history appears here once the meter is running")
        p.restore()
        return
    t1 = t[-1]
    t0 = t[0] if not window else max(t[0], t1 - window)
    sel = t >= t0
    tt, la, l1 = t[sel], sess.col("LAeq1s")[sel], sess.col("LAeq1m")[sel]
    fin = la[np.isfinite(la)]
    lo = 10 * math.floor((min(fin.min() if len(fin) else 40, warn) - 5) / 10)
    hi = 10 * math.ceil((max(fin.max() if len(fin) else 100, alarm) + 3) / 10)
    span = max(t1 - t0, 1)
    fm = QtGui.QFontMetricsF(p.font())
    lm = fm.horizontalAdvance("120 ") + 4
    g = r.adjusted(lm, 4, -4, -fm.height() - 4)
    x = lambda v: g.left() + (v - t0) / span * g.width()
    y = lambda v: g.bottom() - (min(max(v, lo), hi) - lo) / (hi - lo) * g.height()
    for v in range(int(lo), int(hi) + 1, 10):
        p.setPen(QtGui.QPen(QtGui.QColor(T["grid"]), 1))
        p.drawLine(QtCore.QPointF(g.left(), y(v)), QtCore.QPointF(g.right(), y(v)))
        p.setPen(QtGui.QColor(T["dim"]))
        p.drawText(QtCore.QRectF(r.left(), y(v) - fm.height() / 2, lm - 4, fm.height()),
                   QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, str(v))
    step = next(s for s in (60, 300, 600, 900, 1800, 3600, 7200, 14400, 86400) if span / s <= 8)
    tk = math.ceil(t0 / step) * step
    while tk <= t1:
        p.drawText(QtCore.QRectF(x(tk) - 40, g.bottom() + 2, 80, fm.height()), QtCore.Qt.AlignCenter,
                   datetime.datetime.fromtimestamp(tk).strftime("%H:%M"))
        tk += step
    # bin to pixels (energy mean) for the area, keep 1-min line as is
    npx = max(int(g.width()), 2)
    bins = np.clip(((tt - t0) / span * (npx - 1)).astype(int), 0, npx - 1)
    e = np.where(np.isfinite(la), 10 ** (la / 10), 0)
    s = np.bincount(bins, e, npx)
    c = np.bincount(bins, None, npx)
    have = c > 0
    lv = np.full(npx, np.nan)
    lv[have] = 10 * np.log10(np.maximum(s[have] / c[have], 1e-20))
    path = QtGui.QPainterPath()
    started = False
    for i in np.nonzero(have)[0]:
        px = g.left() + i
        if not started:
            path.moveTo(px, g.bottom()); started = True
        path.lineTo(px, y(lv[i]))
    if started:
        path.lineTo(g.left() + np.nonzero(have)[0][-1], g.bottom())
        col = QtGui.QColor(T["accent"]); col.setAlpha(140)
        p.fillPath(path, col)
    pts = [QtCore.QPointF(x(a), y(b)) for a, b in zip(tt[::max(1, len(tt) // npx)], l1[::max(1, len(tt) // npx)])
           if np.isfinite(b)]
    if len(pts) > 1:
        p.setPen(QtGui.QPen(QtGui.QColor(T["fg"]), 1.5))
        p.drawPolyline(QtGui.QPolygonF(pts))
    p.setPen(QtGui.QPen(QtGui.QColor(T["red"]), 1, QtCore.Qt.DashLine))
    p.drawLine(QtCore.QPointF(g.left(), y(alarm)), QtCore.QPointF(g.right(), y(alarm)))
    last_x, row = -1e9, 0
    for mt, txt in sess.markers:
        if t0 <= mt <= t1:
            mx = x(mt)
            row = row + 1 if mx - last_x < 120 else 0     # stagger labels that would overlap
            last_x = mx
            p.setPen(QtGui.QPen(QtGui.QColor(T["yellow"]), 1))
            p.drawLine(QtCore.QPointF(mx, g.top()), QtCore.QPointF(mx, g.bottom()))
            p.drawText(QtCore.QPointF(mx + 3, g.top() + fm.ascent() + (row % 4) * fm.height()), txt[:24])
    for et, _ in sess.events:
        if t0 <= et <= t1:
            p.setPen(QtGui.QPen(QtGui.QColor(T["red"]), 2))
            p.drawLine(QtCore.QPointF(x(et), g.bottom() - 6), QtCore.QPointF(x(et), g.bottom()))
    p.setPen(QtGui.QColor(T["dim"]))
    p.drawText(g.adjusted(6, 2, 0, 0), QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop,
               "▮ LAeq 1 s   — LAeq 1 min   ┆ markers")
    p.restore()


class SessionGraph(QtWidgets.QWidget):
    def __init__(self, main):
        super().__init__()
        self.main = main
        self.setMinimumHeight(150)

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        c = self.main.cfg
        paint_session_graph(p, self.rect().adjusted(2, 2, -2, -2), self.main.session, c["session_view"],
                            c["alarm"], c["alarm"] - c["warn_gap"])


# -------------------------------------------------------------- panels ----

class MeterPanel(QtWidgets.QWidget):
    """Central widget: big number + thumb, history strip / bar, and the statistics tiles."""
    def __init__(self, main):
        super().__init__()
        self.main = main
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(8, 6, 8, 6)
        self.banner = QtWidgets.QPushButton("")
        self.banner.setVisible(False)
        self.banner.setToolTip("Click to dismiss")
        self.banner.clicked.connect(lambda: self.banner.setVisible(False))
        v.addWidget(self.banner)
        self.main_tile = Tile("LAS", size=96)
        self.thumb = Thumb(130)
        ml = self.main_tile.layout()
        ml.removeWidget(self.main_tile.v)
        hr = QtWidgets.QHBoxLayout()
        hr.addWidget(self.thumb); hr.addWidget(self.main_tile.v, 1)
        ml.insertLayout(1, hr)
        self.history = HistoryStrip(main.cfg["history_secs"])
        self.bar = LevelBar()
        ml.addWidget(self.history, 1)
        ml.addWidget(self.bar)
        v.addWidget(self.main_tile, 3)
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(6)
        self.t_1m = Tile("LAeq 1 min", tip="A-weighted equivalent level over the last minute")
        self.t_long = Tile("LAeq 15 min", tip="A-weighted equivalent level over the long window (Settings)")
        self.t_tot = Tile("LAeq since reset")
        self.t_clong = Tile("LCeq 15 min", tip="C-weighted equivalent level over the long window (bass energy)")
        self.t_max = Tile("LAFmax", tip="Highest A-weighted Fast level since reset")
        self.t_pk = Tile("LCpeak", tip="Highest C-weighted true peak since reset (4x oversampled)")
        self.t_l10 = Tile("LA10", tip="Level exceeded 10 % of the time over the long window (LAF)")
        self.t_l90 = Tile("LA90", tip="Level exceeded 90 % of the time over the long window: the background")
        self.grid = grid
        self.cols = 0
        v.addLayout(grid)
        self.banner_until = 0
        self.apply_style()
        self.setMinimumSize(260, 220)
        self.reflow()

    def reflow(self):
        """Scale the big number, thumb and tiles to the window, and re-grid the tiles when narrow / short."""
        w, h = max(self.width(), 260), max(self.height(), 220)
        cols = 4 if w >= 720 else 2
        small = [self.t_1m, self.t_long, self.t_tot, self.t_clong, self.t_max, self.t_pk, self.t_l10, self.t_l90]
        show = small if h >= 520 or cols == 4 else [self.t_1m, self.t_long, self.t_max, self.t_pk]
        if h < 330:
            show = [self.t_1m, self.t_long] if cols == 2 else small[:4]
        key = (cols, len(show))
        if key != self.cols:
            self.cols = key
            for t in small:
                self.grid.removeWidget(t)
                t.setVisible(t in show)
            for i, t in enumerate(show):
                self.grid.addWidget(t, i // cols, i % cols)
        rows = math.ceil(len(show) / cols)
        tiles_h = rows * (70 if h >= 520 else 52)
        top_h = max(h - tiles_h - 20, 120)
        self.main_tile.set_size(max(20, min(96, top_h * 0.13, w * 0.075)))
        self.thumb.set_size(max(36, min(130, top_h * 0.22, w * 0.12)))
        hist = self.main.cfg["meter_style"] == "history"
        self.history.setVisible(hist and top_h > 200)
        self.bar.setVisible(not hist or top_h <= 200)
        for t in small:
            t.set_size(max(11, min(22, w / cols / 14, h / 30)))
            t.s.setVisible(h >= 420)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self.reflow()

    def tiles(self):
        return [self.main_tile, self.t_1m, self.t_long, self.t_tot, self.t_clong,
                self.t_max, self.t_pk, self.t_l10, self.t_l90]

    def apply_style(self):
        hist = self.main.cfg["meter_style"] == "history"
        self.history.setVisible(hist)
        self.bar.setVisible(not hist)
        if getattr(self, "grid", None) is not None:
            self.cols = 0
            self.reflow()
        self.banner.setStyleSheet(f"QPushButton {{ background:{T['banner']}; color:{T['red']}; border:1px solid "
                                  f"{T['red']}; font-weight:700; font-size:12pt; padding:8px; }}")

    def flash(self, text):
        self.banner.setText("⚠  " + text)
        self.banner.setVisible(True)
        self.banner_until = time.time() + 20

    def tick(self):
        if self.banner.isVisible() and time.time() > self.banner_until:
            self.banner.setVisible(False)


class RtaPanel(QtWidgets.QWidget):
    def __init__(self, main):
        super().__init__()
        self.main = main
        c = main.cfg
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(6, 6, 6, 6)
        row = QtWidgets.QHBoxLayout()
        self.w = QtWidgets.QComboBox(); self.w.addItems(["Z", "A", "C"]); self.w.setCurrentText(c["rta_weight"])
        self.sp = QtWidgets.QComboBox(); self.sp.addItems(["Fast", "Slow", "Leq"])
        self.sp.setCurrentText({"F": "Fast", "S": "Slow", "Leq": "Leq"}[c["rta_speed"]])
        rst = QtWidgets.QPushButton("Reset"); rst.setToolTip("Reset peak hold and the Leq average")
        rst.clicked.connect(self.reset)
        self.spec = QtWidgets.QCheckBox("Spectro"); self.spec.setChecked(c["rta_spectro"])
        self.spec.setToolTip("Show the scrolling spectrogram")
        self.spec.toggled.connect(self._spec_toggled)
        for lab, wdg in (("Weighting", self.w), ("Speed", self.sp)):
            row.addWidget(QtWidgets.QLabel(lab)); row.addWidget(wdg)
        row.addWidget(rst); row.addStretch(); row.addWidget(self.spec)
        v.addLayout(row)
        self.rta = RtaWidget()
        v.addWidget(self.rta, 3)
        self.spectro = SpectrogramWidget()
        self.spectro.top, self.spectro.range = c["spec_top"], c["spec_range"]
        v.addWidget(self.spectro, 2)
        self.spectro.setVisible(c["rta_spectro"])
        self.w.currentTextChanged.connect(self._save)
        self.sp.currentTextChanged.connect(self._save)
        self.fs = None

    def _save(self):
        self.main.cfg["rta_weight"] = self.w.currentText()
        self.main.cfg["rta_speed"] = {"Fast": "F", "Slow": "S", "Leq": "Leq"}[self.sp.currentText()]
        self.rta.reset_peaks()

    def _spec_toggled(self, on):
        self.main.cfg["rta_spectro"] = on
        self.spectro.setVisible(on)

    def reset(self):
        self.rta.reset_peaks()
        e = self.main.engine
        if e and e.rta:
            e.rta.reset_leq()

    def active(self):
        return self.isVisible() and self.visibleRegion().isEmpty() is False

    def tick(self, e):
        """Called every GUI tick with the running engine (or None)."""
        if e is None or not self.active():
            if e is not None:
                e.rta, e.tap = None, None
            return
        if e.rta is None:
            e.rta = ThirdOctaveBank(e.fs)
            self.rta.reset_peaks()
        want_tap = self.spectro.isVisible()
        if want_tap and e.tap is None:
            e.tap = []
            self.spectro.setup(e.fs)
        elif not want_tap:
            e.tap = None
        if self.spectro.fs != e.fs and want_tap:
            self.spectro.setup(e.fs)
        c = self.main.cfg
        ms = e.rta.levels(c["rta_speed"])
        lv = 10 * np.log10(np.maximum(ms, 1e-20)) + e.offset + weighting_db(c["rta_weight"], e.rta.fc)
        self.rta.set_data(e.rta.labels, lv)
        if e.tap:
            blocks, e.tap[:] = list(e.tap), []
            top = float(np.nanmax(lv)) if len(lv) else 100
            self.spectro.top = 10 * math.ceil((top + 5) / 10)
            self.spectro.feed(blocks, e.offset)


class DosePanel(QtWidgets.QWidget):
    def __init__(self, main):
        super().__init__()
        self.main = main
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(6, 6, 6, 6)
        top = QtWidgets.QHBoxLayout()
        self.elapsed = QtWidgets.QLabel("Not running")
        rst = QtWidgets.QPushButton("Reset dose")
        rst.setToolTip("Start a new exposure period (e.g. a new shift or show)")
        rst.clicked.connect(self.reset)
        top.addWidget(self.elapsed); top.addStretch(); top.addWidget(rst)
        v.addLayout(top)
        self.rows = []
        for d in make_dosimeters():
            box = QtWidgets.QFrame(); box.setObjectName("tile")
            g = QtWidgets.QGridLayout(box)
            g.setContentsMargins(10, 6, 10, 6)
            name = QtWidgets.QLabel(f"<b>{d.name}</b>")
            desc = QtWidgets.QLabel(d.desc)
            desc.setWordWrap(True)
            bar = QtWidgets.QProgressBar(); bar.setRange(0, 1000); bar.setFormat("%v")
            bar.setTextVisible(True)
            info = QtWidgets.QLabel("")
            info.setWordWrap(True)
            status = QtWidgets.QLabel("")
            g.addWidget(name, 0, 0); g.addWidget(status, 0, 1, QtCore.Qt.AlignRight)
            g.addWidget(desc, 1, 0, 1, 2)
            g.addWidget(bar, 2, 0, 1, 2)
            g.addWidget(info, 3, 0, 1, 2)
            v.addWidget(box)
            self.rows.append((desc, bar, info, status))
        note = QtWidgets.QLabel("For guidance only. Use a certified Type 2 dosimeter for compliance records.")
        note.setWordWrap(True)
        self.note = note
        v.addWidget(note)
        v.addStretch()
        self.refresh_theme()

    def refresh_theme(self):
        for desc, bar, info, status in self.rows:
            desc.setStyleSheet(f"color:{T['dim']}; font-size:8pt; background:transparent;")
            info.setStyleSheet(f"color:{T['fg']}; {MONO} font-size:9pt; background:transparent;")
        self.note.setStyleSheet(f"color:{T['dim']}; font-size:8pt;")

    def reset(self):
        if self.main.engine:
            self.main.engine.reset_dose()
        self.main.dose_flags = set()
        self.main.add_event("Dose reset")

    def tick(self, e):
        if e is None:
            return
        self.elapsed.setText(f"Exposure time  {hms(e.dose[0].elapsed)}")
        for d, (desc, bar, info, status) in zip(e.dose, self.rows):
            pct = d.dose * 100
            bar.setValue(int(min(pct, 100) * 10))
            bar.setFormat(f"{pct:.1f} %")
            col = T["red"] if pct >= 100 else T["yellow"] if pct >= 50 else T["green"]
            bar.setStyleSheet(f"QProgressBar::chunk {{ background:{col}; border-radius:3px; }}")
            proj = d.projected()
            ttf = d.time_to_full()
            info.setText(f"TWA {fmt(d.twa())} dBA   ·   8 h projection {fmt(proj, 0)} %   ·   "
                         f"100 % in {'—' if not np.isfinite(ttf) else 'reached' if ttf == 0 else hms(ttf)}")
            if d.peak_over:
                st, sc = "PEAK ≥ 140 dB", T["red"]
            elif pct >= 100:
                st, sc = "OVER LIMIT", T["red"]
            elif pct >= 50:
                st, sc = "ACTION LEVEL", T["yellow"]
            elif np.isfinite(proj) and proj >= 100:
                st, sc = "ON PACE TO EXCEED", T["yellow"]
            else:
                st, sc = "OK", T["green"]
            status.setText(st)
            status.setStyleSheet(f"color:{sc}; font-weight:700; background:transparent;")


THEATER_LAYOUTS = {
    "2.0": ["L", "R"],
    "5.1": ["L", "R", "C", "LFE", "Ls", "Rs"],
    "7.1": ["L", "R", "C", "LFE", "Lss", "Rss", "Lrs", "Rrs"],
    "7.1.4": ["L", "R", "C", "LFE", "Lss", "Rss", "Lrs", "Rrs", "Ltf", "Rtf", "Ltr", "Rtr"],
}
THEATER_PRESETS = {   # target dBC per channel kind, noise level dBFS
    "Cinema / dub stage": dict(screen=85, surround=82, lfe=91, height=82, level=-20),
    "Mid-size mix room": dict(screen=82, surround=79, lfe=88, height=79, level=-20),
    "Small room / nearfield": dict(screen=79, surround=76, lfe=85, height=76, level=-20),
    "Home theater (75 dBC)": dict(screen=75, surround=75, lfe=85, height=75, level=-30),
}


def chan_kind(name):
    if name in ("L", "R", "C"):
        return "screen"
    if name == "LFE":
        return "lfe"
    if name.endswith(("tf", "tr")):
        return "height"
    return "surround"


class TheaterPanel(QtWidgets.QWidget):
    """Per-channel pink-noise level calibration (C-weighted, Slow)."""
    TOL = 0.5

    def __init__(self, main):
        super().__init__()
        self.main = main
        self.tc = main.cfg["theater"]
        self.player = NoisePlayer()
        self.cur = -1
        self.auto = None          # (channel index, t_start) while auto-sequencing
        self.results = {}
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(6, 6, 6, 6)
        f = QtWidgets.QGridLayout()
        self.out = QtWidgets.QComboBox()
        self.layout_cb = QtWidgets.QComboBox(); self.layout_cb.addItems(list(THEATER_LAYOUTS))
        self.preset = QtWidgets.QComboBox(); self.preset.addItems(list(THEATER_PRESETS))
        self.noise = QtWidgets.QComboBox(); self.noise.addItems(list(NOISE_TYPES))
        self.level = QtWidgets.QDoubleSpinBox(); self.level.setRange(-40, -6); self.level.setSuffix(" dBFS RMS")
        self.level.setToolTip("AES17 level: a full-scale sine reads 0 dBFS. Cinema standard is -20 dBFS.")
        f.addWidget(QtWidgets.QLabel("Output"), 0, 0); f.addWidget(self.out, 0, 1, 1, 3)
        f.addWidget(QtWidgets.QLabel("Format"), 1, 0); f.addWidget(self.layout_cb, 1, 1)
        f.addWidget(QtWidgets.QLabel("Room"), 1, 2); f.addWidget(self.preset, 1, 3)
        f.addWidget(QtWidgets.QLabel("Noise"), 2, 0); f.addWidget(self.noise, 2, 1)
        f.addWidget(QtWidgets.QLabel("Level"), 2, 2); f.addWidget(self.level, 2, 3)
        v.addLayout(f)
        self.big = QtWidgets.QLabel("--.- dBC")
        self.big.setAlignment(QtCore.Qt.AlignCenter)
        self.hint = QtWidgets.QLabel("Pick a channel and press Play")
        self.hint.setAlignment(QtCore.Qt.AlignCenter)
        v.addWidget(self.big); v.addWidget(self.hint)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Ch", "Out", "Target", "Measured", "Adjust", ""])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.itemChanged.connect(self._item_changed)
        v.addWidget(self.table, 1)
        b = QtWidgets.QHBoxLayout()
        self.play = QtWidgets.QPushButton("▶ Play"); self.play.setCheckable(True)
        self.play.clicked.connect(self.toggle_play)
        prev = QtWidgets.QPushButton("◀ Prev"); prev.clicked.connect(lambda: self.step(-1))
        nxt = QtWidgets.QPushButton("Next ▶"); nxt.clicked.connect(lambda: self.step(1))
        cap = QtWidgets.QPushButton("Capture"); cap.setToolTip("Store the last 3 s LCeq for this channel")
        cap.clicked.connect(self.capture)
        self.auto_btn = QtWidgets.QPushButton("Auto run"); self.auto_btn.setCheckable(True)
        self.auto_btn.setToolTip("Play each channel, wait 4 s to settle, capture 3 s, move on")
        self.auto_btn.clicked.connect(self.toggle_auto)
        for w in (self.play, prev, nxt, cap, self.auto_btn):
            b.addWidget(w)
        v.addLayout(b)
        self.populate_outputs()
        self.layout_cb.setCurrentText(self.tc["layout"])
        self.preset.setCurrentText(self.tc["preset"])
        self.noise.setCurrentText(self.tc["noise"])
        self.level.setValue(self.tc["level"])
        self.layout_cb.currentTextChanged.connect(self.build_table)
        self.preset.currentTextChanged.connect(self.apply_preset)
        self.noise.currentTextChanged.connect(self._regen)
        self.level.valueChanged.connect(self._regen)
        self.out.currentIndexChanged.connect(self._out_changed)
        self.table.currentCellChanged.connect(lambda r, *_: self.select(r))
        self.build_table()
        self.refresh_theme()

    def refresh_theme(self):
        self.big.setStyleSheet(f"font-size:40pt; font-weight:700; {MONO}")
        self.hint.setStyleSheet(f"color:{T['dim']}; font-size:12pt;")

    def populate_outputs(self):
        self.out.blockSignals(True)
        self.out.clear()
        if sd is not None:
            try:
                apis = sd.query_hostapis()
                for i, d in enumerate(sd.query_devices()):
                    if d["max_output_channels"] > 0:
                        self.out.addItem(f"{d['name']}  [{apis[d['hostapi']]['name']}]", i)
            except Exception:
                pass
        self.out.addItem(SIM_OUT, None)
        i = self.out.findText(self.tc.get("out_device") or "")
        self.out.setCurrentIndex(i if i >= 0 else 0)
        self.out.blockSignals(False)

    def channels(self):
        return THEATER_LAYOUTS[self.layout_cb.currentText()]

    def build_table(self):
        self.tc["layout"] = self.layout_cb.currentText()
        self.table.blockSignals(True)
        names = self.channels()
        self.table.setRowCount(len(names))
        pr = THEATER_PRESETS[self.preset.currentText()]
        for i, n in enumerate(names):
            out = self.tc["outputs"].get(n, i + 1)
            tgt = self.tc["targets"].get(n, pr[chan_kind(n)])
            items = [n, str(out), f"{tgt:.1f}", "", "", ""]
            for c, txt in enumerate(items):
                it = QtWidgets.QTableWidgetItem(txt)
                if c not in (1, 2):
                    it.setFlags(it.flags() & ~QtCore.Qt.ItemIsEditable)
                it.setTextAlignment(QtCore.Qt.AlignCenter)
                self.table.setItem(i, c, it)
            if n in self.results:
                self._show_result(i, self.results[n])
        self.table.blockSignals(False)

    def apply_preset(self, name):
        self.tc["preset"] = name
        self.tc["targets"] = {}
        pr = THEATER_PRESETS[name]
        self.level.setValue(pr["level"])
        self.build_table()

    def _item_changed(self, it):
        n = self.table.item(it.row(), 0).text()
        try:
            val = float(it.text())
        except ValueError:
            return
        if it.column() == 1:
            self.tc["outputs"][n] = max(1, int(val))
            if it.row() == self.cur and self.player.playing:
                self.start_noise()
        elif it.column() == 2:
            self.tc["targets"][n] = val
        self.main.save_cfg()

    def target(self, row):
        try:
            return float(self.table.item(row, 2).text())
        except (ValueError, AttributeError):
            return 85.0

    def out_ch(self, row):
        try:
            return max(1, int(float(self.table.item(row, 1).text()))) - 1
        except (ValueError, AttributeError):
            return row

    def _out_changed(self):
        self.tc["out_device"] = self.out.currentText()
        if self.player.playing:
            self.start_noise()

    def _regen(self):
        self.tc["noise"] = self.noise.currentText()
        self.tc["level"] = self.level.value()
        if self.player.playing:
            self.start_noise()

    def select(self, row):
        if row < 0:
            return
        self.cur = row
        self.table.selectRow(row)
        if self.player.playing:
            self.start_noise()

    def step(self, d):
        n = self.table.rowCount()
        self.select((max(self.cur, 0) + d) % n)

    def start_noise(self):
        if self.cur < 0:
            self.select(0)
        ch = self.out_ch(self.cur)
        name = self.channels()[self.cur]
        lo, hi = NOISE_TYPES[self.noise.currentText()]
        if name == "LFE" and hi > 200:
            lo, hi = 20, 120          # always band-limit the LFE channel
        try:
            fs = self.player.open(self.out.currentData(), ch)
        except Exception as ex:
            self.play.setChecked(False)
            QtWidgets.QMessageBox.critical(self, APP_NAME, f"Could not open the output:\n{ex}")
            return
        self.player.set_buffer(make_pink(fs, lo, min(hi, fs * 0.45), self.level.value()))
        self.player.ch = ch
        if self.main.sim is not None:
            self.main.sim.inject = self.player if self.out.currentData() is None else None
        self.play.setChecked(True)
        self.play.setText("■ Stop")
        self.t_start = time.time()

    def stop_noise(self):
        self.player.stop()
        self.play.setChecked(False)
        self.play.setText("▶ Play")

    def toggle_play(self, on):
        if on:
            if not self.main.running():
                self.main.run_btn.click()
            self.start_noise()
        else:
            self.auto = None
            self.auto_btn.setChecked(False)
            self.stop_noise()

    def toggle_auto(self, on):
        if on:
            if not self.main.running():
                self.main.run_btn.click()
            self.select(0)
            self.start_noise()
            self.auto = time.time()
        else:
            self.auto = None
            self.stop_noise()

    def capture(self):
        e = self.main.engine
        if e is None or self.cur < 0:
            return
        val = e.leq("C", 3)
        name = self.channels()[self.cur]
        self.results[name] = val
        self._show_result(self.cur, val)
        self.main.add_event(f"Theater cal {name}: {val:.1f} dBC (target {self.target(self.cur):.1f})")

    def _show_result(self, row, val):
        d = self.target(row) - val
        ok = abs(d) <= self.TOL
        self.table.blockSignals(True)
        self.table.item(row, 3).setText(f"{val:.1f}")
        self.table.item(row, 4).setText("OK" if ok else f"{'▲ up' if d > 0 else '▼ down'} {abs(d):.1f} dB")
        self.table.item(row, 5).setText("✓" if ok else "✗")
        col = QtGui.QColor(T["green"] if ok else T["yellow"] if abs(d) <= 2 else T["red"])
        for c in (3, 4, 5):
            self.table.item(row, c).setForeground(col)
        self.table.blockSignals(False)

    def tick(self, e):
        if e is None or not self.player.playing or self.cur < 0:
            if not self.player.playing and self.big.text() != "--.- dBC":
                self.big.setText("--.- dBC")
                self.refresh_theme()
            return
        lv = e.level("C", "S")
        d = self.target(self.cur) - lv
        col = T["green"] if abs(d) <= self.TOL else T["yellow"] if abs(d) <= 2 else T["red"]
        self.big.setText(f"{fmt(lv).strip()} dBC")
        self.big.setStyleSheet(f"font-size:40pt; font-weight:700; {MONO} color:{col};")
        name = self.channels()[self.cur]
        if abs(d) <= self.TOL:
            msg = f"{name}: on target ✓"
        else:
            msg = f"{name}: turn {'UP' if d > 0 else 'DOWN'} {abs(d):.1f} dB"
        if self.auto is not None:
            el = time.time() - self.auto
            msg += f"   ·   auto: {'settling' if el < 4 else 'capturing'} {el:.0f}/7 s"
            if el >= 7:
                self.capture()
                if self.cur + 1 >= self.table.rowCount():
                    self.auto_btn.setChecked(False)
                    self.toggle_auto(False)
                    self.hint.setText("Auto run finished. Adjust the red/yellow channels and re-run.")
                    return
                self.step(1)
                self.auto = time.time()
        self.hint.setText(msg)


class SessionPanel(QtWidgets.QWidget):
    def __init__(self, main):
        super().__init__()
        self.main = main
        c = main.cfg
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(6, 6, 6, 6)
        top = QtWidgets.QHBoxLayout()
        self.info = QtWidgets.QLabel("")
        self.view = QtWidgets.QComboBox()
        self.views = {"Last 10 min": 600, "Last hour": 3600, "Last 4 hours": 14400, "Whole session": 0}
        self.view.addItems(list(self.views))
        self.view.setCurrentText(next((k for k, s in self.views.items() if s == c["session_view"]), "Last 10 min"))
        self.view.currentTextChanged.connect(self._view)
        mk = QtWidgets.QPushButton("◆ Marker"); mk.setToolTip("Drop a marker (shortcut M)")
        mk.clicked.connect(lambda: main.add_marker())
        new = QtWidgets.QPushButton("New session"); new.clicked.connect(lambda: main.new_session())
        ecsv = QtWidgets.QPushButton("Export CSV…"); ecsv.clicked.connect(lambda: main.export_csv())
        epdf = QtWidgets.QPushButton("PDF report…"); epdf.clicked.connect(lambda: main.export_pdf())
        sch = QtWidgets.QPushButton("Schedule…"); sch.setToolTip("Start and stop logging by time (e.g. every Sunday)")
        sch.clicked.connect(lambda: main.edit_schedule())
        top.addWidget(self.info, 1)
        for w in (self.view, mk, new, ecsv, epdf, sch):
            top.addWidget(w)
        v.addLayout(top)
        meta = QtWidgets.QHBoxLayout()
        self.venue = QtWidgets.QLineEdit(c["venue"]); self.venue.setPlaceholderText("Venue / show / room")
        self.oper = QtWidgets.QLineEdit(c["operator"]); self.oper.setPlaceholderText("Operator")
        self.venue.textChanged.connect(lambda s: c.__setitem__("venue", s))
        self.oper.textChanged.connect(lambda s: c.__setitem__("operator", s))
        meta.addWidget(self.venue, 2); meta.addWidget(self.oper, 1)
        v.addLayout(meta)
        split = QtWidgets.QHBoxLayout()
        self.graph = SessionGraph(main)
        split.addWidget(self.graph, 3)
        side = QtWidgets.QVBoxLayout()
        side.addWidget(QtWidgets.QLabel("Markers & events"))
        self.events = QtWidgets.QListWidget()
        side.addWidget(self.events, 1)
        side.addWidget(QtWidgets.QLabel("Alarms"))
        self.rules = QtWidgets.QTableWidget(0, 3)
        self.rules.setHorizontalHeaderLabels(["On", "Metric", "≥ dB / %"])
        self.rules.verticalHeader().setVisible(False)
        self.rules.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        side.addWidget(self.rules, 1)
        rb = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("+"); add.clicked.connect(self.add_rule)
        rem = QtWidgets.QPushButton("−"); rem.clicked.connect(self.remove_rule)
        self.beep = QtWidgets.QCheckBox("Beep"); self.beep.setChecked(c["alarm_beep"])
        self.beep.toggled.connect(lambda on: c.__setitem__("alarm_beep", on))
        rb.addWidget(add); rb.addWidget(rem); rb.addStretch(); rb.addWidget(self.beep)
        side.addLayout(rb)
        split.addLayout(side, 2)
        v.addLayout(split, 1)
        self.load_rules()

    def _view(self, t):
        self.main.cfg["session_view"] = self.views[t]
        self.graph.update()

    def load_rules(self):
        self.rules.blockSignals(True)
        rs = self.main.cfg["alarm_rules"]
        self.rules.setRowCount(len(rs))
        for i, r in enumerate(rs):
            cb = QtWidgets.QCheckBox(); cb.setChecked(r["on"])
            cb.toggled.connect(self.save_rules)
            cell = QtWidgets.QWidget(); h = QtWidgets.QHBoxLayout(cell); h.setContentsMargins(6, 0, 0, 0)
            h.addWidget(cb)
            self.rules.setCellWidget(i, 0, cell)
            m = QtWidgets.QComboBox(); m.addItems(ALARM_METRICS); m.setCurrentText(r["metric"])
            m.currentTextChanged.connect(self.save_rules)
            self.rules.setCellWidget(i, 1, m)
            s = QtWidgets.QDoubleSpinBox(); s.setRange(0, 200); s.setDecimals(1); s.setValue(r["thr"])
            s.valueChanged.connect(self.save_rules)
            self.rules.setCellWidget(i, 2, s)
        self.rules.resizeRowsToContents()
        self.rules.blockSignals(False)

    def save_rules(self, *_):
        rs = []
        for i in range(self.rules.rowCount()):
            cb = self.rules.cellWidget(i, 0).findChild(QtWidgets.QCheckBox)
            rs.append(dict(on=cb.isChecked(), metric=self.rules.cellWidget(i, 1).currentText(),
                           thr=self.rules.cellWidget(i, 2).value()))
        self.main.cfg["alarm_rules"] = rs
        self.main.alarm_armed = {}
        self.main.save_cfg()

    def add_rule(self):
        self.main.cfg["alarm_rules"].append(dict(on=True, metric="LAeq 1 min", thr=self.main.cfg["alarm"]))
        self.load_rules(); self.save_rules()

    def remove_rule(self):
        r = self.rules.currentRow()
        if r < 0:
            r = self.rules.rowCount() - 1
        if r >= 0:
            del self.main.cfg["alarm_rules"][r]
            self.load_rules(); self.save_rules()

    def add_event_row(self, t, text, kind):
        it = QtWidgets.QListWidgetItem(f"{datetime.datetime.fromtimestamp(t):%H:%M:%S}  "
                                       f"{'◆' if kind == 'marker' else '⚠' if kind == 'alarm' else '·'} {text}")
        if kind == "alarm":
            it.setForeground(QtGui.QColor(T["red"]))
        elif kind == "marker":
            it.setForeground(QtGui.QColor(T["yellow"]))
        self.events.addItem(it)
        self.events.scrollToBottom()

    def tick(self):
        s = self.main.session
        n = len(s)
        self.info.setText(f"Session {datetime.datetime.fromtimestamp(s.start):%b %d %H:%M}  ·  "
                          f"{hms(n)} recorded  ·  {len(s.markers)} markers")
        if self.graph.isVisible():
            self.graph.update()


# ---------------------------------------------------- calibration wizard ----

class CalWizard(QtWidgets.QDialog):
    """Step-by-step calibration into a named mic profile."""
    def __init__(self, main):
        super().__init__(main)
        self.main = main
        self.setWindowTitle("Calibration wizard")
        self.setMinimumSize(560, 430)
        self.result = None        # dict(offset, info, method, points, ...)
        self.readings = collections.deque(maxlen=20)
        v = QtWidgets.QVBoxLayout(self)
        self.title = QtWidgets.QLabel()
        self.title.setStyleSheet("font-size:15pt; font-weight:700;")
        v.addWidget(self.title)
        self.stack = QtWidgets.QStackedWidget()
        v.addWidget(self.stack, 1)
        b = QtWidgets.QHBoxLayout()
        self.back = QtWidgets.QPushButton("◀ Back"); self.back.clicked.connect(self.go_back)
        self.next = QtWidgets.QPushButton("Next ▶"); self.next.clicked.connect(self.go_next)
        self.next.setDefault(True)
        cancel = QtWidgets.QPushButton("Cancel"); cancel.clicked.connect(self.reject)
        b.addWidget(cancel); b.addStretch(); b.addWidget(self.back); b.addWidget(self.next)
        v.addLayout(b)
        self._build_profile_page()
        self._build_method_page()
        self._build_calibrator_page()
        self._build_file_page()
        self._build_manual_page()
        self._build_done_page()
        self.history = [0]
        self.show_page(0)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.poll)
        self.timer.start(100)

    # pages ---------------------------------------------------------------
    def _page(self, title):
        w = QtWidgets.QWidget()
        w.setProperty("title", title)
        lay = QtWidgets.QVBoxLayout(w)
        self.stack.addWidget(w)
        return lay

    @staticmethod
    def _note(text):
        l = QtWidgets.QLabel(text)
        l.setWordWrap(True)
        l.setStyleSheet(f"color:{T['dim']};")
        return l

    def _build_profile_page(self):
        lay = self._page("1 · Which microphone?")
        lay.addWidget(self._note("Each mic gets its own profile, so swapping mics never mixes up calibrations."))
        self.p_existing = QtWidgets.QRadioButton("Update an existing profile")
        self.p_combo = QtWidgets.QComboBox(); self.p_combo.addItems(list(self.main.cfg["profiles"]))
        self.p_combo.setCurrentText(self.main.cfg["profile"])
        self.p_new = QtWidgets.QRadioButton("Create a new profile")
        self.p_name = QtWidgets.QLineEdit(); self.p_name.setPlaceholderText("e.g. UMIK-1 #7023270, ECM8000 on Scarlett")
        self.p_existing.setChecked(True)
        self.p_auto = QtWidgets.QCheckBox("Pick this profile automatically when this input device is selected")
        dev = self.main.dev.currentText()
        self.p_auto.setToolTip(dev)
        for w in (self.p_existing, self.p_combo, self.p_new, self.p_name, self.p_auto):
            lay.addWidget(w)
        lay.addWidget(self._note(f"Current input: {dev}"))
        lay.addStretch()

    def _build_method_page(self):
        lay = self._page("2 · How do you want to calibrate?")
        self.m_cal = QtWidgets.QRadioButton("Acoustic calibrator (most accurate)")
        self.m_file = QtWidgets.QRadioButton("Load the mic's calibration file (UMIK-1 / UMIK-2 / REW .txt)")
        self.m_man = QtWidgets.QRadioButton("Enter numbers by hand (mic sensitivity or known offset)")
        self.m_cal.setChecked(True)
        lay.addWidget(self.m_cal)
        lay.addWidget(self._note("Put a 94 or 114 dB calibrator on the mic. Takes 5 seconds."))
        lay.addWidget(self.m_file)
        lay.addWidget(self._note("Applies the frequency correction, and the Sens Factor for absolute level "
                                 "when your input volume is at 100 %. You can fine-tune with a calibrator afterwards."))
        lay.addWidget(self.m_man)
        lay.addWidget(self._note("Mic mV/Pa plus interface gain, or type an offset you already know."))
        lay.addStretch()

    def _build_calibrator_page(self):
        lay = self._page("3 · Acoustic calibrator")
        f = QtWidgets.QFormLayout()
        self.c_level = QtWidgets.QComboBox(); self.c_level.setEditable(True)
        self.c_level.addItems(["94.0", "114.0", "104.0", "124.0"])
        self.c_freq = QtWidgets.QComboBox(); self.c_freq.addItems(["1 kHz", "250 Hz"])
        f.addRow("Calibrator level (dB SPL)", self.c_level)
        f.addRow("Calibrator frequency", self.c_freq)
        lay.addLayout(f)
        lay.addWidget(self._note("1. Set your interface gain where you'll leave it. Don't touch it afterwards.\n"
                                 "2. Seat the calibrator fully on the mic and switch it on.\n"
                                 "3. Wait for 'steady', then press Measure."))
        self.c_live = QtWidgets.QLabel("--.- dB")
        self.c_live.setAlignment(QtCore.Qt.AlignCenter)
        self.c_live.setStyleSheet(f"font-size:30pt; font-weight:700; {MONO}")
        self.c_stab = QtWidgets.QLabel("")
        self.c_stab.setAlignment(QtCore.Qt.AlignCenter)
        lay.addWidget(self.c_live); lay.addWidget(self.c_stab)
        row = QtWidgets.QHBoxLayout()
        self.c_start = QtWidgets.QPushButton("Start meter"); self.c_start.clicked.connect(self.main.run_btn.click)
        self.c_meas = QtWidgets.QPushButton("Measure (5 s)"); self.c_meas.clicked.connect(self.measure)
        self.c_prog = QtWidgets.QProgressBar(); self.c_prog.setRange(0, 100); self.c_prog.setValue(0)
        row.addWidget(self.c_start); row.addWidget(self.c_meas); row.addWidget(self.c_prog, 1)
        lay.addLayout(row)
        self.c_res = QtWidgets.QLabel("")
        lay.addWidget(self.c_res)
        lay.addStretch()

    def _build_file_page(self):
        lay = self._page("3 · Mic calibration file")
        row = QtWidgets.QHBoxLayout()
        self.f_path = QtWidgets.QLineEdit(); self.f_path.setReadOnly(True)
        br = QtWidgets.QPushButton("Choose file…"); br.clicked.connect(self.choose_file)
        row.addWidget(self.f_path, 1); row.addWidget(br)
        lay.addLayout(row)
        lay.addWidget(self._note("UMIK files come in a 0° and a 90° version. Use 0° when the mic points at the "
                                 "source, 90° when it points at the ceiling (typical for a room or theater)."))
        self.f_info = QtWidgets.QLabel("")
        self.f_info.setWordWrap(True)
        lay.addWidget(self.f_info)
        self.f_curve = QtWidgets.QCheckBox("Apply the frequency-response correction"); self.f_curve.setChecked(True)
        self.f_sens = QtWidgets.QCheckBox("Set absolute level from the Sens Factor "
                                          "(macOS input volume for this mic must be at 100 %)")
        self.f_sens.setChecked(True)
        lay.addWidget(self.f_curve); lay.addWidget(self.f_sens)
        self.f_plot = QtWidgets.QLabel()
        self.f_plot.setMinimumHeight(120)
        lay.addWidget(self.f_plot, 1)
        self.file_data = None

    def _build_manual_page(self):
        lay = self._page("3 · Enter numbers")
        self.n_offset = QtWidgets.QRadioButton("I know the offset: a 0 dBFS full-scale sine equals")
        self.n_sens = QtWidgets.QRadioButton("Calculate from mic sensitivity and interface specs")
        self.n_offset.setChecked(True)
        self.n_off = QtWidgets.QDoubleSpinBox(); self.n_off.setRange(60, 200); self.n_off.setDecimals(2)
        self.n_off.setSuffix(" dB SPL")
        prof = self.main.active_profile()
        self.n_off.setValue(prof["offset"] - 3.01)
        lay.addWidget(self.n_offset); lay.addWidget(self.n_off)
        lay.addWidget(self.n_sens)
        f = QtWidgets.QFormLayout()
        self.n_mv = QtWidgets.QDoubleSpinBox(); self.n_mv.setRange(0.01, 1000); self.n_mv.setDecimals(2)
        self.n_mv.setValue(10.0); self.n_mv.setSuffix(" mV/Pa")
        self.n_gain = QtWidgets.QDoubleSpinBox(); self.n_gain.setRange(-20, 80); self.n_gain.setSuffix(" dB")
        self.n_gain.setValue(40)
        self.n_fs = QtWidgets.QDoubleSpinBox(); self.n_fs.setRange(-20, 40); self.n_fs.setValue(9.0)
        self.n_fs.setSuffix(" dBu"); self.n_fs.setDecimals(1)
        f.addRow("Mic sensitivity", self.n_mv)
        f.addRow("Preamp gain", self.n_gain)
        f.addRow("Input level for 0 dBFS", self.n_fs)
        lay.addLayout(f)
        self.n_res = QtWidgets.QLabel("")
        lay.addWidget(self.n_res)
        for w in (self.n_off, self.n_mv, self.n_gain, self.n_fs):
            w.valueChanged.connect(self.update_manual)
        self.n_offset.toggled.connect(self.update_manual)
        lay.addWidget(self._note("Preamp gain is often a guess on analog knobs, so a calibrator check afterwards "
                                 "is a good idea."))
        lay.addStretch()
        self.update_manual()

    def _build_done_page(self):
        lay = self._page("4 · Check and save")
        self.d_sum = QtWidgets.QLabel("")
        self.d_sum.setWordWrap(True)
        lay.addWidget(self.d_sum)
        self.d_live = QtWidgets.QLabel("")
        self.d_live.setAlignment(QtCore.Qt.AlignCenter)
        self.d_live.setStyleSheet(f"font-size:30pt; font-weight:700; {MONO}")
        lay.addWidget(self.d_live)
        lay.addWidget(self._note("This is the live reading (Z, Fast) with the new calibration. "
                                 "With the calibrator still on, it should match its level."))
        lay.addStretch()

    # navigation ----------------------------------------------------------
    def show_page(self, i):
        self.stack.setCurrentIndex(i)
        self.title.setText(self.stack.widget(i).property("title"))
        self.back.setEnabled(len(self.history) > 1)
        self.next.setText("Save ✓" if i == 5 else "Next ▶")
        if i == 5:
            self.fill_summary()

    def go_back(self):
        if len(self.history) > 1:
            self.history.pop()
            self.show_page(self.history[-1])

    def go_next(self):
        i = self.stack.currentIndex()
        if i == 0:
            if self.p_new.isChecked() and not self.p_name.text().strip():
                self.p_name.setFocus(); return
            nxt = 1
        elif i == 1:
            nxt = 2 if self.m_cal.isChecked() else 3 if self.m_file.isChecked() else 4
        elif i in (2, 3, 4):
            if not self.make_result(i):
                return
            nxt = 5
        else:
            self.save()
            return
        self.history.append(nxt)
        self.show_page(nxt)

    # actions -------------------------------------------------------------
    def base_profile(self):
        name = self.p_name.text().strip() if self.p_new.isChecked() else self.p_combo.currentText()
        prof = copy.deepcopy(self.main.cfg["profiles"].get(name) or DEFAULTS["profiles"]["Default"])
        return name, prof

    def measure(self):
        e = self.main.engine
        if e is None:
            return
        e.begin_cal(5.0, "Z")
        self.c_meas.setEnabled(False)
        self.c_res.setText("Measuring… keep it steady")

    def poll(self):
        e = self.main.engine
        running = e is not None and self.main.running()
        self.c_start.setVisible(not running)
        self.c_meas.setEnabled(running and not (e and e.cal_target))
        if not running:
            self.c_live.setText("--.- dB")
            self.c_stab.setText("Meter is stopped: press Start meter")
            return
        lz = e.level("Z", "F")
        self.readings.append(lz)
        self.c_live.setText(f"{fmt(lz).strip()} dB")
        spread = max(self.readings) - min(self.readings) if len(self.readings) > 5 else 99
        self.c_stab.setText("● steady" if spread < 0.3 else "○ waiting for a steady tone…")
        self.c_stab.setStyleSheet(f"color:{T['green'] if spread < 0.3 else T['yellow']}; font-weight:600;")
        if e.cal_target:
            self.c_prog.setValue(int(e.cal_progress * 100))
        if e.cal_done:
            raw = e.finish_cal()
            try:
                ref = float(self.c_level.currentText())
            except ValueError:
                ref = 94.0
            self.cal_raw = (raw, ref)
            off = ref - raw
            self.c_res.setText(f"Done. Offset {off:.2f} dB (was {e.offset:.2f}, change {off - e.offset:+.2f} dB). "
                               "Press Next.")
            self.c_prog.setValue(100)
        if self.stack.currentIndex() == 5 and self.result:
            raw_now = lz - e.offset
            self.d_live.setText(f"{raw_now + self.result['offset']:.1f} dB")

    def choose_file(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Mic calibration file", os.path.expanduser("~"),
                                                        "Calibration files (*.txt *.cal *.frd *.csv);;All files (*)")
        if not path:
            return
        try:
            self.file_data = parse_cal_file(path)
        except Exception as ex:
            self.f_info.setText(f"Couldn't read that file: {ex}")
            self.file_data = None
            return
        d = self.file_data
        self.f_path.setText(path)
        pts = np.array(d["points"])
        self.f_info.setText(f"{len(pts)} points, {pts[0, 0]:.0f} Hz – {pts[-1, 0]:.0f} Hz, "
                            f"correction {pts[:, 1].min():+.1f} to {pts[:, 1].max():+.1f} dB"
                            + (f"<br>Sens Factor {d['sens_factor']:+.3f} dB" if d["sens_factor"] is not None else
                               "<br>No Sens Factor in this file: absolute level stays as is (use a calibrator).")
                            + (f" · serial {d['serial']}" if d["serial"] else "")
                            + (f" · AGain {d['again']:g} dB" if d["again"] is not None else ""))
        self.f_sens.setEnabled(d["sens_factor"] is not None)
        self.f_sens.setChecked(d["sens_factor"] is not None)
        if self.p_new.isChecked() and not self.p_name.text().strip() and d["serial"]:
            self.p_name.setText(f"Mic #{d['serial']}")
        self.draw_curve(pts)

    def draw_curve(self, pts):
        w, h = 500, 120
        pm = QtGui.QPixmap(w, h); pm.fill(QtGui.QColor(T["track"]))
        p = QtGui.QPainter(pm); p.setRenderHint(QtGui.QPainter.Antialiasing)
        lo, hi = min(-6, pts[:, 1].min() - 1), max(6, pts[:, 1].max() + 1)
        x = lambda f: math.log10(f / 10) / math.log10(2400) * w
        y = lambda d: h - (d - lo) / (hi - lo) * h
        p.setPen(QtGui.QPen(QtGui.QColor(T["grid"]), 1))
        p.drawLine(QtCore.QPointF(0, y(0)), QtCore.QPointF(w, y(0)))
        p.setPen(QtGui.QColor(T["dim"]))
        for fr, lab in ((100, "100"), (1000, "1k"), (10000, "10k")):
            p.drawText(QtCore.QPointF(x(fr) + 2, h - 3), lab)
        p.setPen(QtGui.QPen(QtGui.QColor(T["accent"]), 2))
        p.drawPolyline(QtGui.QPolygonF([QtCore.QPointF(x(max(f, 10)), y(d)) for f, d in pts]))
        p.end()
        self.f_plot.setPixmap(pm)

    def update_manual(self, *_):
        off = self.manual_offset()
        self.n_res.setText(f"→ offset {off:.2f} dB  (1 Pa / 94 dB reads {94:.0f} at "
                           f"{94 - off:.1f} dBFS RMS)")

    def manual_offset(self):
        if self.n_offset.isChecked():
            return self.n_off.value() + 3.01      # full-scale sine RMS is -3.01 in engine units
        dbu = 20 * math.log10(self.n_mv.value() / 1000 / 0.7746)
        raw = dbu + self.n_gain.value() - self.n_fs.value() - 3.01
        return 94 - raw

    def make_result(self, page):
        _, prof = self.base_profile()
        now = datetime.datetime.now()
        if page == 2:
            if not getattr(self, "cal_raw", None):
                self.c_res.setText("Press Measure first.")
                return False
            raw, ref = self.cal_raw
            self.result = dict(offset=ref - raw, method="calibrator",
                               info=f"{ref:.1f} dB calibrator · {now:%b %d %H:%M}",
                               cal_points=prof.get("cal_points"), cal_file=prof.get("cal_file", ""),
                               sens_factor=prof.get("sens_factor"), serial=prof.get("serial", ""))
        elif page == 3:
            d = self.file_data
            if not d:
                self.f_info.setText("Choose a calibration file first.")
                return False
            off = prof["offset"]
            what = "level unchanged"
            if self.f_sens.isChecked() and d["sens_factor"] is not None:
                off = 94.0 - d["sens_factor"]
                what = f"Sens {d['sens_factor']:+.2f} dB"
            self.result = dict(offset=off, method="cal file",
                               info=f"{d['name']} ({what}) · {now:%b %d %H:%M}",
                               cal_points=d["points"] if self.f_curve.isChecked() else None,
                               cal_file=d["name"], sens_factor=d["sens_factor"], serial=d["serial"])
            if not self.f_sens.isChecked() and prof.get("calibrated"):
                self.result["info"] = prof["cal_info"] + f" + {d['name']} curve"
        else:
            off = self.manual_offset()
            self.result = dict(offset=off, method="manual",
                               info=f"manual offset {off:.2f} dB · {now:%b %d %H:%M}",
                               cal_points=prof.get("cal_points"), cal_file=prof.get("cal_file", ""),
                               sens_factor=prof.get("sens_factor"), serial=prof.get("serial", ""))
        return True

    def fill_summary(self):
        name, prof = self.base_profile()
        r = self.result
        self.d_sum.setText(f"<b>Profile:</b> {name}<br><b>Method:</b> {r['method']}<br>"
                           f"<b>Offset:</b> {r['offset']:.2f} dB (was {prof['offset']:.2f}, "
                           f"change {r['offset'] - prof['offset']:+.2f} dB)<br>"
                           f"<b>Frequency correction:</b> {r['cal_file'] if r['cal_points'] else 'none'}")

    def save(self):
        name, prof = self.base_profile()
        prof.update(self.result)
        prof.update(calibrated=True, cal_info=self.result["info"], date=datetime.datetime.now().isoformat(),
                    use_curve=bool(self.result["cal_points"]))
        if self.p_auto.isChecked():
            prof.update(auto_select=True, device_hint=self.main.dev.currentText())
        prof.pop("info", None)
        self.main.cfg["profiles"][name] = prof
        self.main.set_profile(name)
        self.main.add_event(f"Calibrated '{name}': offset {prof['offset']:.2f} dB ({prof['method']})")
        self.accept()

    def done(self, r):
        self.timer.stop()
        e = self.main.engine
        if e is not None:
            e.cal_target = 0
        super().done(r)


# ---------------------------------------------------- mini + menu bar ----

class MiniMeter(QtWidgets.QWidget):
    """Small frameless box: thumb + live reading. Drag to move, pin to keep on top, ⤢ to go back."""
    def __init__(self, main):
        super().__init__(None, QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)
        self.main = main
        self.setWindowTitle(APP_NAME)
        self.setObjectName("mini")
        self.setAttribute(QtCore.Qt.WA_StyledBackground)
        h = QtWidgets.QHBoxLayout(self)
        h.setContentsMargins(8, 6, 4, 6)
        self.thumb = Thumb(64)
        h.addWidget(self.thumb)
        col = QtWidgets.QVBoxLayout(); col.setSpacing(0)
        self.val = QtWidgets.QLabel("--.-")
        self.val.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.sub = QtWidgets.QLabel("")
        self.sub.setAlignment(QtCore.Qt.AlignRight)
        col.addWidget(self.val); col.addWidget(self.sub)
        h.addLayout(col, 1)
        btns = QtWidgets.QVBoxLayout(); btns.setSpacing(0)
        self.pin = QtWidgets.QToolButton(); self.pin.setText("📌"); self.pin.setCheckable(True)
        self.pin.setToolTip("Keep on top"); self.pin.setChecked(main.cfg["mini_on_top"])
        self.pin.toggled.connect(self.set_on_top)
        big = QtWidgets.QToolButton(); big.setText("⤢"); big.setToolTip("Back to full view")
        big.clicked.connect(main.exit_mini)
        btns.addWidget(self.pin); btns.addWidget(big); btns.addStretch()
        h.addLayout(btns)
        self.setFixedSize(290, 84)
        self._drag = None
        self.refresh_theme()
        self.set_on_top(self.pin.isChecked(), show=False)

    def refresh_theme(self):
        self.setStyleSheet(f"QWidget#mini {{ background:{T['panel']}; border:1px solid {T['btn_border']}; border-radius:8px; }}"
                           f"QToolButton {{ background:transparent; border:none; color:{T['dim']}; font-size:11pt; padding:0; }}"
                           f"QToolButton:checked, QToolButton:hover {{ color:{T['fg']}; }}")
        self.sub.setStyleSheet(f"color:{T['dim']}; font-size:8pt; background:transparent;")

    def set_on_top(self, on, show=True):
        self.main.cfg["mini_on_top"] = on
        self.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint, on)
        if show:
            self.show()
        self.main.save_cfg()

    def update_view(self, value, color, state, sub):
        self.val.setText(fmt(value))
        self.val.setStyleSheet(f"color:{color}; font-size:34pt; font-weight:600; {MONO} background:transparent;")
        self.sub.setText(sub)
        self.thumb.setState(state)

    def mousePressEvent(self, e):
        if e.button() == QtCore.Qt.LeftButton:
            self._drag = e.globalPosition().toPoint() - self.frameGeometry().topLeft()

    def mouseMoveEvent(self, e):
        if self._drag is not None:
            self.move(e.globalPosition().toPoint() - self._drag)

    def mouseReleaseEvent(self, e):
        self._drag = None
        self.main.cfg["mini_pos"] = [self.x(), self.y()]
        self.main.save_cfg()

    def mouseDoubleClickEvent(self, e):
        self.main.exit_mini()

    def closeEvent(self, e):
        e.ignore()
        self.main.exit_mini()


_status_click = [None]
_NSStatusTarget = None
if sys.platform == "darwin":
    try:
        import objc  # noqa: F401
        from AppKit import (NSStatusBar, NSVariableStatusItemLength, NSColor, NSFont, NSFontAttributeName,
                            NSForegroundColorAttributeName, NSMutableAttributedString, NSAttributedString)
        from Foundation import NSObject

        class _NSStatusTarget(NSObject):
            def clicked_(self, sender):
                if _status_click[0]:
                    _status_click[0]()
    except Exception:
        _NSStatusTarget = None


class MenuBarMeter:
    """Mac: live '● 87.3' text in the menu bar (green / yellow / red dot). Elsewhere: a tray icon."""
    def __init__(self, main):
        self.main = main
        self.item = None
        self.tray = None
        self.menu = QtWidgets.QMenu()
        self.lines = [self.menu.addAction("") for _ in range(4)]
        for a in self.lines:
            a.setEnabled(False)
        self.menu.addSeparator()
        self.menu.addAction("Show ServiceSPL", main.exit_mini)
        self.menu.addAction("Mini meter", main.enter_mini)
        self.a_run = self.menu.addAction("Start", main.run_btn.click)
        self.menu.addAction("Drop marker", lambda: main.add_marker())
        self.a_log = self.menu.addAction("Start logging", main.log_btn.click)
        self.menu.addSeparator()
        self.menu.addAction("Quit", QtWidgets.QApplication.quit)
        self.menu.aboutToShow.connect(self.fill_menu)
        self._last = ("", "")
        if _NSStatusTarget is not None:
            try:
                self.item = NSStatusBar.systemStatusBar().statusItemWithLength_(NSVariableStatusItemLength)
                self.target = _NSStatusTarget.alloc().init()
                _status_click[0] = self.popup
                btn = self.item.button()
                btn.setTarget_(self.target)
                btn.setAction_("clicked:")
                self.font = NSFont.monospacedDigitSystemFontOfSize_weight_(13.0, 0.23)
            except Exception:
                self.item = None
        if self.item is None and QtWidgets.QSystemTrayIcon.isSystemTrayAvailable():
            self.tray = QtWidgets.QSystemTrayIcon()
            self.tray.setContextMenu(self.menu)
            self.tray.activated.connect(lambda reason: self.popup()
                                        if reason == QtWidgets.QSystemTrayIcon.Trigger else None)
            self.tray.show()
        self.set(None, None)

    def popup(self):
        self.menu.popup(QtGui.QCursor.pos())

    def fill_menu(self):
        m = self.main
        e = m.engine
        self.a_run.setText("Stop" if m.running() else "Start")
        self.a_log.setText("Stop logging" if m.log_writer else "Start logging")
        vals = ["Not running", "", "", ""]
        if e is not None and m.running():
            w, t = m.cfg["main_weight"], m.cfg["main_time"]
            vals = [f"L{w}{t}  {fmt(e.level(w, t)).strip()} dB",
                    f"LAeq 1 min {fmt(m.last_rec.get('LAeq1m')).strip()}  ·  {m.cfg['long_leq_min']} min "
                    f"{fmt(m.last_rec.get('LAeqLong')).strip()}",
                    f"LCpeak {fmt(e.L(e.max_pk)).strip()}  ·  LAFmax {fmt(e.L(e.max_af)).strip()}",
                    f"Dose OSHA {e.dose[1].dose * 100:.1f} %  ·  NIOSH {e.dose[2].dose * 100:.1f} %"]
        for a, s in zip(self.lines, vals):
            a.setText(s)
            a.setVisible(bool(s))

    def set(self, value, color):
        text = "--.-" if value is None or not np.isfinite(value) else f"{value:.1f}"
        col = color or T["off"]
        if (text, col) == self._last:
            return
        self._last = (text, col)
        if self.item is not None:
            try:
                q = QtGui.QColor(col)
                s = NSMutableAttributedString.alloc().init()
                s.appendAttributedString_(NSAttributedString.alloc().initWithString_attributes_(
                    "● ", {NSFontAttributeName: self.font,
                           NSForegroundColorAttributeName: NSColor.colorWithSRGBRed_green_blue_alpha_(
                               q.redF(), q.greenF(), q.blueF(), 1.0)}))
                s.appendAttributedString_(NSAttributedString.alloc().initWithString_attributes_(
                    text, {NSFontAttributeName: self.font,
                           NSForegroundColorAttributeName: NSColor.labelColor()}))
                self.item.button().setAttributedTitle_(s)
            except Exception:
                pass
        elif self.tray is not None:
            pm = QtGui.QPixmap(64, 64); pm.fill(QtCore.Qt.transparent)
            p = QtGui.QPainter(pm); p.setRenderHint(QtGui.QPainter.Antialiasing)
            p.setBrush(QtGui.QColor(col)); p.setPen(QtCore.Qt.NoPen)
            p.drawEllipse(2, 2, 60, 60)
            p.setPen(QtGui.QColor("white"))
            f = p.font(); f.setPixelSize(30); f.setBold(True); p.setFont(f)
            p.drawText(pm.rect(), QtCore.Qt.AlignCenter, text.split(".")[0])
            p.end()
            self.tray.setIcon(QtGui.QIcon(pm))
            self.tray.setToolTip(f"{APP_NAME}: {text} dB")

    def remove(self):
        if self.item is not None:
            try:
                NSStatusBar.systemStatusBar().removeStatusItem_(self.item)
            except Exception:
                pass
            self.item = None
        if self.tray is not None:
            self.tray.hide()
            self.tray = None


# -------------------------------------------------------- integrations ----

def osc_message(addr, *args):
    """Encode one OSC 1.0 message (int / float / string arguments)."""
    def pad(b):
        return b + b"\0" * (4 - len(b) % 4)
    tags, payload = ",", b""
    for a in args:
        if isinstance(a, bool):
            a = int(a)
        if isinstance(a, int):
            tags += "i"; payload += struct.pack(">i", a)
        elif isinstance(a, float):
            tags += "f"; payload += struct.pack(">f", a)
        else:
            tags += "s"; payload += pad(str(a).encode("utf-8"))
    return pad(addr.encode()) + pad(tags.encode()) + payload


class CompanionLink:
    """Pushes the live level to Bitfocus Companion over OSC (UDP, default port 12321):
    custom variables spl_level, spl_leq1m, spl_leqlong, spl_peak, spl_dose, spl_state, spl_color, spl_item
    and, optionally, the text + background colour of one button."""
    def __init__(self):
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.last = 0.0
        self.error = ""

    def send(self, cfg, d):
        now = time.time()
        if now - self.last < 0.25:
            return
        self.last = now
        host, port = cfg["companion_host"] or "127.0.0.1", int(cfg["companion_port"])
        msgs = [osc_message(f"/custom-variable/spl_{k}/value", str(v)) for k, v in (
            ("level", d["level"]), ("leq1m", d["leq1m"]), ("leqlong", d["leqlong"]), ("peak", d["peak"]),
            ("dose", d["dose"]), ("state", d["state"]), ("color", d["color"]), ("item", d.get("item", "")))]
        loc = (cfg.get("companion_button") or "").strip()
        if re.fullmatch(r"\d+/\d+/\d+", loc):
            q = QtGui.QColor(d["color"])
            msgs.append(osc_message(f"/location/{loc}/style/text", f"{d['level']}\n{d['label']}"))
            msgs.append(osc_message(f"/location/{loc}/style/bgcolor", q.red(), q.green(), q.blue()))
            msgs.append(osc_message(f"/location/{loc}/style/color", 255, 255, 255))
        try:
            for m in msgs:
                self.sock.sendto(m, (host, port))
            self.error = ""
        except OSError as e:
            self.error = str(e)


def secret_get(account):
    """Read a secret from the macOS Keychain (or the config fallback elsewhere)."""
    if sys.platform == "darwin":
        try:
            for service in ("ServiceSPL", "SPL Monitor"):
                r = subprocess.run(["security", "find-generic-password", "-s", service, "-a", account, "-w"],
                                   capture_output=True, text=True, timeout=5)
                if r.returncode == 0:
                    return r.stdout.strip()
            return ""
        except Exception:
            return ""
    return ""


def secret_set(account, value):
    if sys.platform == "darwin":
        try:
            subprocess.run(["security", "add-generic-password", "-U", "-s", "ServiceSPL", "-a", account, "-w", value],
                           capture_output=True, timeout=5)
            return True
        except Exception:
            return False
    return False


class PcoClient:
    """Planning Center Services API (personal access token: application id + secret)."""
    BASE = "https://api.planningcenteronline.com/services/v2"

    def __init__(self, app_id, secret):
        self.auth = "Basic " + base64.b64encode(f"{app_id}:{secret}".encode()).decode()

    def get(self, path, **params):
        url = self.BASE + path + ("?" + urllib.parse.urlencode(params) if params else "")
        req = urllib.request.Request(url, headers={"Authorization": self.auth, "Accept": "application/json",
                                                   "User-Agent": f"ServiceSPL/{APP_VERSION}"})
        with urllib.request.urlopen(req, timeout=10) as r:
            return json.load(r)

    def service_types(self):
        d = self.get("/service_types", per_page=100)
        return [(x["id"], x["attributes"].get("name") or x["id"]) for x in d.get("data", [])]

    def plans(self, st):
        out = []
        for flt, order, n in (("past", "-sort_date", 2), ("future", "sort_date", 6)):
            d = self.get(f"/service_types/{st}/plans", filter=flt, order=order, per_page=n)
            for x in d.get("data", []):
                a = x["attributes"]
                label = " · ".join(s for s in (a.get("dates") or a.get("short_dates") or "", a.get("title") or "",
                                               a.get("series_title") or "") if s)
                out.append((x["id"], label or x["id"], a.get("sort_date") or ""))
        out.sort(key=lambda p: p[2])
        return out

    def live_item(self, st, plan):
        """Current Live item: (item_id, title) or (None, None) when Live isn't running."""
        d = self.get(f"/service_types/{st}/plans/{plan}/live", include="current_item_time,items")
        rel = (d.get("data") or {}).get("relationships", {}).get("current_item_time", {}).get("data")
        if not rel:
            return None, None
        inc = {(x["type"], x["id"]): x for x in d.get("included", [])}
        it = inc.get(("ItemTime", rel["id"]))
        item_ref = (it or {}).get("relationships", {}).get("item", {}).get("data") if it else None
        if not item_ref:
            return rel["id"], "Live item"
        item = inc.get(("Item", item_ref["id"]))
        title = item["attributes"].get("title") if item else None
        if not title:
            try:
                title = self.get(f"/service_types/{st}/plans/{plan}/items/{item_ref['id']}")["data"]["attributes"]["title"]
            except Exception:
                title = "Live item"
        return item_ref["id"], title


class PcoFollower:
    """Background poller for Services Live. The GUI reads .current and .error."""
    def __init__(self):
        self.client = None
        self.st = self.plan = None
        self.current = (None, None)
        self.error = ""
        self.on = False
        self._stop = threading.Event()
        threading.Thread(target=self._run, daemon=True).start()

    def follow(self, client, st, plan):
        self.client, self.st, self.plan = client, st, plan
        self.current = (None, None)
        self.on = True

    def stop(self):
        self.on = False
        self.current = (None, None)

    def _run(self):
        while not self._stop.wait(3.0):
            if not (self.on and self.client):
                continue
            try:
                self.current = self.client.live_item(self.st, self.plan)
                self.error = ""
            except Exception as e:
                self.error = str(e)


def item_stats(sess, t0, t1):
    t = sess.col("time")
    if not len(t):
        return None
    i0, i1 = np.searchsorted(t, t0), np.searchsorted(t, t1)
    if i1 <= i0:
        return None
    return dict(LAeq=energy_mean(sess.col("LAeq1s")[i0:i1]), LCeq=energy_mean(sess.col("LCeq1s")[i0:i1]),
                LAFmax=float(np.nanmax(sess.col("LAFmax1s")[i0:i1])),
                LCpeak=float(np.nanmax(sess.col("LCpeak1s")[i0:i1])), dur=i1 - i0)


def find_onedrive():
    """Best guess at the user's synced OneDrive folder."""
    cands = []
    if sys.platform == "darwin":
        import glob
        cands += sorted(glob.glob(os.path.expanduser("~/Library/CloudStorage/OneDrive*")))
    for k in ("OneDriveCommercial", "OneDriveConsumer", "OneDrive"):
        if os.environ.get(k):
            cands.append(os.environ[k])
    cands.append(os.path.expanduser("~/OneDrive"))
    return next((c for c in cands if os.path.isdir(c)), "")


# ----------------------------------------------------------------- NDI ----

class NDIlib_send_create_t(ctypes.Structure):
    _fields_ = [("p_ndi_name", ctypes.c_char_p), ("p_groups", ctypes.c_char_p),
                ("clock_video", ctypes.c_bool), ("clock_audio", ctypes.c_bool)]


class NDIlib_video_frame_v2_t(ctypes.Structure):
    _fields_ = [("xres", ctypes.c_int), ("yres", ctypes.c_int), ("FourCC", ctypes.c_int),
                ("frame_rate_N", ctypes.c_int), ("frame_rate_D", ctypes.c_int),
                ("picture_aspect_ratio", ctypes.c_float), ("frame_format_type", ctypes.c_int),
                ("timecode", ctypes.c_int64), ("p_data", ctypes.c_void_p),
                ("line_stride_in_bytes", ctypes.c_int), ("p_metadata", ctypes.c_char_p),
                ("timestamp", ctypes.c_int64)]


def ndi_library_candidates(custom=""):
    c = [custom] if custom else []
    for v in ("NDI_RUNTIME_DIR_V6", "NDI_RUNTIME_DIR_V5", "NDI_RUNTIME_DIR_V4"):
        d = os.environ.get(v)
        if d:
            c += [os.path.join(d, n) for n in ("libndi.dylib", "Processing.NDI.Lib.x64.dll", "libndi.so")]
    if sys.platform == "darwin":
        c += ["/usr/local/lib/libndi.dylib", "/opt/homebrew/lib/libndi.dylib",
              "/Library/NDI SDK for Apple/lib/macOS/libndi.dylib",
              "/Library/NDI Advanced SDK for Apple/lib/macOS/libndi_advanced.dylib",
              "/Applications/NDI Video Monitor.app/Contents/Frameworks/libndi.dylib",
              "/Applications/NDI Tools/NDI Video Monitor.app/Contents/Frameworks/libndi.dylib",
              "/Applications/NDI Scan Converter.app/Contents/Frameworks/libndi.dylib"]
    elif sys.platform.startswith("win"):
        for v in ("6", "5"):
            c.append(rf"C:\Program Files\NDI\NDI {v} Runtime\v{v}\Processing.NDI.Lib.x64.dll")
            c.append(rf"C:\Program Files\NDI\NDI {v} Tools\Runtime\Processing.NDI.Lib.x64.dll")
        c.append("Processing.NDI.Lib.x64.dll")
    else:
        c += ["libndi.so.6", "libndi.so.5", "libndi.so"]
    return c


class NdiSender:
    """Sends rendered meter frames as an NDI source, using the installed NDI runtime (no Python package)."""
    FOURCC_BGRA = ord("B") | ord("G") << 8 | ord("R") << 16 | ord("A") << 24

    def __init__(self, name, custom_lib=""):
        self.lib = None
        err = []
        for path in ndi_library_candidates(custom_lib):
            if os.path.isabs(path) and not os.path.exists(path):
                continue
            try:
                self.lib = ctypes.CDLL(path)
                self.path = path
                break
            except OSError as e:
                err.append(f"{path}: {e}")
        if self.lib is None:
            raise RuntimeError("NDI runtime not found. Install the free NDI Tools (ndi.video/tools) "
                               "or point Settings → Integrations → NDI library at libndi.")
        L = self.lib
        L.NDIlib_initialize.restype = ctypes.c_bool
        L.NDIlib_send_create.restype = ctypes.c_void_p
        L.NDIlib_send_create.argtypes = [ctypes.POINTER(NDIlib_send_create_t)]
        L.NDIlib_send_send_video_v2.argtypes = [ctypes.c_void_p, ctypes.POINTER(NDIlib_video_frame_v2_t)]
        L.NDIlib_send_destroy.argtypes = [ctypes.c_void_p]
        if not L.NDIlib_initialize():
            raise RuntimeError("NDI could not start on this computer (CPU not supported?).")
        self._name = name.encode()
        cs = NDIlib_send_create_t(self._name, None, False, False)
        self.inst = L.NDIlib_send_create(ctypes.byref(cs))
        if not self.inst:
            raise RuntimeError("Could not create the NDI sender.")

    def send(self, img, fps=30):
        img = img.convertToFormat(QtGui.QImage.Format_ARGB32)   # BGRA in memory on little-endian CPUs
        w, h = img.width(), img.height()
        ptr = img.constBits()
        buf = (ctypes.c_char * img.sizeInBytes()).from_buffer_copy(bytes(ptr))
        fr = NDIlib_video_frame_v2_t(w, h, self.FOURCC_BGRA, fps * 1000, 1000, w / h, 1, 2 ** 63 - 1,
                                     ctypes.cast(buf, ctypes.c_void_p), img.bytesPerLine(), None, 0)
        self.lib.NDIlib_send_send_video_v2(self.inst, ctypes.byref(fr))

    def close(self):
        if self.lib and getattr(self, "inst", None):
            self.lib.NDIlib_send_destroy(self.inst)
            self.inst = None


_THUMB_RENDER = {}


def paint_meter_card(p, r, d, transparent=False):
    """The 'just the SPL meter' graphic used for NDI: thumb, big number, label line and a level bar."""
    r = QtCore.QRectF(r)
    p.save()
    p.setRenderHint(QtGui.QPainter.Antialiasing)
    p.setCompositionMode(QtGui.QPainter.CompositionMode_Source)
    p.fillRect(r, QtGui.QColor(0, 0, 0, 0) if transparent else QtGui.QColor("#121417"))
    p.setCompositionMode(QtGui.QPainter.CompositionMode_SourceOver)
    if transparent:
        p.setBrush(QtGui.QColor(18, 20, 23, 200)); p.setPen(QtCore.Qt.NoPen)
        p.drawRoundedRect(r.adjusted(r.width() * 0.03, r.height() * 0.08, -r.width() * 0.03, -r.height() * 0.08),
                          r.height() * 0.06, r.height() * 0.06)
    H = r.height()
    pad = H * 0.12
    ts = H * 0.5
    st = d.get("state", "off")
    if st not in _THUMB_RENDER:
        bg, xf = {"good": ("#2ea866", ""), "bad": ("#e03b3b", "rotate(180 50 50)")}.get(st, ("#8a9099", ""))
        _THUMB_RENDER[st] = QtSvg.QSvgRenderer(QtCore.QByteArray(THUMB_SVG.format(bg=bg, xf=xf).encode()))
    _THUMB_RENDER[st].render(p, QtCore.QRectF(r.left() + pad, r.top() + (H - ts) / 2 - H * 0.05, ts, ts))
    f = QtGui.QFont("Menlo"); f.setStyleHint(QtGui.QFont.Monospace); f.setBold(True); f.setPixelSize(int(H * 0.42))
    p.setFont(f); p.setPen(QtGui.QColor(d.get("color", "#e8e8e8")))
    num = QtCore.QRectF(r.left() + pad + ts, r.top() + H * 0.12, r.width() - 2 * pad - ts, H * 0.5)
    p.drawText(num, QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, d.get("level", "--.-"))
    f2 = QtGui.QFont("Helvetica"); f2.setPixelSize(int(H * 0.06)); p.setFont(f2)
    p.setPen(QtGui.QColor("#c9ced6"))
    sub_r = QtCore.QRectF(num.left(), num.bottom(), num.width(), H * 0.1)
    sub = QtGui.QFontMetrics(f2).elidedText(d.get("sub", ""), QtCore.Qt.ElideRight, int(sub_r.width()))
    p.drawText(sub_r, QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, sub)
    br = QtCore.QRectF(num.left(), r.top() + H * 0.76, num.width(), H * 0.06)
    p.setPen(QtCore.Qt.NoPen); p.setBrush(QtGui.QColor("#0c0d0f")); p.drawRoundedRect(br, 4, 4)
    try:
        lv = float(d.get("level"))
    except (TypeError, ValueError):
        lv = 30.0
    frac = min(max((lv - 40) / 80, 0), 1)
    p.setBrush(QtGui.QColor(d.get("color", "#3ddc84")))
    p.drawRoundedRect(QtCore.QRectF(br.left(), br.top(), br.width() * frac, br.height()), 4, 4)
    p.restore()


# ---------------------------------------------------- schedule + guide ----

DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]


def schedule_window(entry, now):
    """True when 'now' falls inside this schedule entry."""
    if not entry.get("on", True):
        return False
    try:
        h0, m0 = map(int, entry["start"].split(":"))
        h1, m1 = map(int, entry["stop"].split(":"))
    except (KeyError, ValueError):
        return False
    if entry.get("date"):
        if now.strftime("%Y-%m-%d") != entry["date"]:
            return False
    elif now.weekday() not in entry.get("days", []):
        return False
    a, b, x = h0 * 60 + m0, h1 * 60 + m1, now.hour * 60 + now.minute
    return a <= x < b if a < b else False


class ScheduleDialog(QtWidgets.QDialog):
    """Edit recurring (weekly) or one-off logging windows."""
    def __init__(self, cfg, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Logging schedule")
        self.setMinimumWidth(720)
        self.cfg = cfg
        v = QtWidgets.QVBoxLayout(self)
        v.addWidget(QtWidgets.QLabel("Logging starts and stops by itself inside these windows. The meter starts "
                                     "automatically if it isn't running. Leave Date empty for every week."))
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["On", "Days", "Date (one-off)", "Start", "Stop", "New session + report"])
        self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        self.table.verticalHeader().setVisible(False)
        v.addWidget(self.table)
        b = QtWidgets.QHBoxLayout()
        add = QtWidgets.QPushButton("+ Add window"); add.clicked.connect(lambda: self.add_row({}))
        rem = QtWidgets.QPushButton("− Remove"); rem.clicked.connect(lambda: self.table.removeRow(self.table.currentRow()))
        b.addWidget(add); b.addWidget(rem); b.addStretch()
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        b.addWidget(bb)
        v.addLayout(b)
        for e in cfg["schedule"]:
            self.add_row(e)

    def add_row(self, e):
        r = self.table.rowCount()
        self.table.insertRow(r)
        on = QtWidgets.QCheckBox(); on.setChecked(e.get("on", True))
        self.table.setCellWidget(r, 0, on)
        days = QtWidgets.QWidget(); h = QtWidgets.QHBoxLayout(days); h.setContentsMargins(2, 0, 2, 0); h.setSpacing(2)
        sel = e.get("days", [6])
        for i, dname in enumerate(DAYS):
            cb = QtWidgets.QCheckBox(dname[:2]); cb.setChecked(i in sel); h.addWidget(cb)
        self.table.setCellWidget(r, 1, days)
        date = QtWidgets.QLineEdit(e.get("date", "")); date.setPlaceholderText("YYYY-MM-DD")
        self.table.setCellWidget(r, 2, date)
        for c, key, dflt in ((3, "start", "09:00"), (4, "stop", "12:30")):
            te = QtWidgets.QTimeEdit(QtCore.QTime.fromString(e.get(key, dflt), "HH:mm")); te.setDisplayFormat("HH:mm")
            self.table.setCellWidget(r, c, te)
        rep = QtWidgets.QCheckBox(); rep.setChecked(e.get("report", True))
        self.table.setCellWidget(r, 5, rep)
        self.table.resizeRowsToContents()

    def apply(self):
        out = []
        for r in range(self.table.rowCount()):
            days = [i for i, cb in enumerate(self.table.cellWidget(r, 1).findChildren(QtWidgets.QCheckBox)) if cb.isChecked()]
            out.append(dict(on=self.table.cellWidget(r, 0).isChecked(), days=days,
                            date=self.table.cellWidget(r, 2).text().strip(),
                            start=self.table.cellWidget(r, 3).time().toString("HH:mm"),
                            stop=self.table.cellWidget(r, 4).time().toString("HH:mm"),
                            report=self.table.cellWidget(r, 5).isChecked()))
        self.cfg["schedule"] = out


EQ_TARGETS = {
    "Flat": lambda f: np.zeros_like(f),
    "Live / worship house curve": lambda f: 6 * np.clip((np.log10(200) - np.log10(f)) / (np.log10(200) - np.log10(40)), 0, 1)
                                        - np.where(f > 1000, np.log2(np.maximum(f, 1000) / 1000), 0),
    "Cinema X-curve (SMPTE ST 202)": lambda f: -3 * np.log2(np.maximum(f, 2000) / 2000),
    "X-curve, small room (from 4 kHz)": lambda f: -1.5 * np.log2(np.maximum(f, 4000) / 4000),
    "Hi-fi / studio (gentle tilt)": lambda f: 4 * np.clip((np.log10(120) - np.log10(f)) / (np.log10(120) - np.log10(30)), 0, 1)
                                          - 0.5 * np.log2(f / 1000),
}


def eq_suggestions(fc, labels, measured, target, lo_hz=40, hi_hz=16000, tol=2.0, max_boost=3.0):
    """Compare a measured 1/3-octave response with a target curve (aligned on 250 Hz - 4 kHz)
    and return (normalized deviation, [suggestion strings])."""
    fc = np.asarray(fc)
    ref = (fc >= 250) & (fc <= 4000)
    dev = measured - target
    dev = dev - np.mean(dev[ref]) if ref.any() else dev
    tips = []
    for f, lab, dv in zip(fc, labels, dev):
        if f < lo_hz or f > hi_hz or abs(dv) < tol:
            continue
        if dv > 0:
            tips.append(f"Cut {lab} Hz by {dv:.1f} dB (bell, Q ≈ 4.3)")
        elif dv < -6:
            tips.append(f"{lab} Hz is {-dv:.1f} dB low: likely a room null. Move the speaker or sub, don't boost it")
        else:
            tips.append(f"Boost {lab} Hz by {min(-dv, max_boost):.1f} dB (bell, Q ≈ 4.3)"
                        + (" (capped, boosts eat headroom)" if -dv > max_boost else ""))
    return dev, tips


class EqPlot(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.fc, self.dev, self.target = [], None, None
        self.setMinimumHeight(170)

    def set_data(self, fc, labels, dev):
        self.fc, self.labels, self.dev = list(fc), labels, dev
        self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        r = QtCore.QRectF(self.rect().adjusted(34, 6, -6, -18))
        p.fillRect(r, QtGui.QColor(T["track"]))
        f = p.font(); f.setPointSize(8); p.setFont(f)
        y = lambda v: r.center().y() - max(min(v, 12), -12) / 12 * r.height() / 2
        for v in (-12, -6, -2, 0, 2, 6, 12):
            p.setPen(QtGui.QPen(QtGui.QColor(T["dim"] if v == 0 else T["grid"]), 1,
                                QtCore.Qt.SolidLine if v in (0, 6, -6, 12, -12) else QtCore.Qt.DotLine))
            p.drawLine(QtCore.QPointF(r.left(), y(v)), QtCore.QPointF(r.right(), y(v)))
            p.setPen(QtGui.QColor(T["dim"]))
            p.drawText(QtCore.QRectF(0, y(v) - 7, 30, 14), QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter, f"{v:+d}")
        if self.dev is None:
            p.drawText(r, QtCore.Qt.AlignCenter, "Measure to see how the room compares with the target")
            return
        n = len(self.dev)
        bw = r.width() / n
        for i, dv in enumerate(self.dev):
            col = T["green"] if abs(dv) < 2 else T["yellow"] if abs(dv) < 4 else T["red"]
            x0 = r.left() + i * bw
            p.fillRect(QtCore.QRectF(x0 + 1, min(y(0), y(dv)), bw - 2, abs(y(dv) - y(0))), QtGui.QColor(col))
            if i % 3 == 2:
                p.setPen(QtGui.QColor(T["dim"]))
                p.drawText(QtCore.QRectF(x0 - 12, r.bottom() + 2, bw + 24, 14), QtCore.Qt.AlignCenter, self.labels[i])
        p.setPen(QtGui.QColor(T["dim"]))
        p.drawText(r.adjusted(6, 2, 0, 0), QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop,
                   "dB above (+) / below (−) the target curve · green within ±2 dB")


GUIDE_PAGES = [
    ("Welcome", """<h2>Welcome to ServiceSPL</h2>
<p>This guide walks you through setting up and using the meter. Pick a topic on the left, or press <b>Next</b>.</p>
<ol><li><b>Pick your mic</b> in the toolbar input menu (a UMIK or other USB measurement mic is picked automatically).</li>
<li>Press <b>▶ Start</b>. The big number is the live level; the strip under it is the last few minutes.</li>
<li><b>Calibrate</b> once per mic (next page), so the numbers are real dB SPL.</li>
<li>Set <b>Red at</b> to your house limit. Everything turns yellow 5 dB before that and red at it.</li></ol>
<p>Panels around the meter (RTA, dosimeter, theater cal, session, Planning Center) can be dragged, stacked, floated or
hidden from the <b>View</b> menu. The layout is remembered.</p>"""),
    ("Calibrate your mic", """<h2>Calibrate your mic</h2>
<p>Each mic gets a <b>profile</b>. Run <b>Calibrate…</b> and pick one of three ways:</p>
<ul><li><b>Acoustic calibrator</b> (best): put a 94 or 114 dB calibrator on the mic, wait for “steady”, press Measure.</li>
<li><b>Cal file</b>: load the UMIK-1/UMIK-2 file from miniDSP (enter your serial number on their site). This fixes the
frequency response and, with the Mac's input volume for the mic at 100 %, sets the level from the Sens Factor.
Use the <b>90°</b> file when the mic points at the ceiling.</li>
<li><b>Manual</b>: type the mic's mV/Pa and your interface gain, or an offset you already know.</li></ul>
<p>Tick “pick automatically” and the profile follows that mic whenever you select it. The status bar shows which profile
is active and warns when a calibration is over 30 days old.</p>"""),
    ("Monitoring a show", """<h2>Monitoring a service or show</h2>
<ul><li><b>LAeq 1 min / 15 min</b> are the averages most venue limits use (e.g. 100 dB LAeq 15 min).</li>
<li><b>LCeq</b> shows the bass energy; LC − LA above ~15 dB means the mix is very bass heavy.</li>
<li><b>LA10</b> is the loud moments, <b>LA90</b> the background (room noise, HVAC).</li>
<li><b>LCpeak</b> is the highest instantaneous peak. Keep it under 135 dB for audiences, never above 140.</li></ul>
<p>Alarms (Session panel) fire once when a value crosses its threshold: a red banner, a beep and an entry in the log.
Press <b>M</b> (or ◆ Marker) to drop markers like “sermon” or “band in”.</p>"""),
    ("Hearing dose", """<h2>Hearing dose (OSHA / NIOSH)</h2>
<p>The dosimeter adds up exposure since its last reset. 100 % is the daily limit:</p>
<ul><li><b>NIOSH</b>: 85 dBA for 8 h, every +3 dB halves the time (88 dBA = 4 h, 100 dBA = 15 min).</li>
<li><b>OSHA HC</b>: 90 dBA for 8 h with a 5 dB exchange; 50 % is the action level for hearing protection programs.</li>
<li><b>OSHA PEL</b>: the legal limit (ignores anything under 90 dBA).</li></ul>
<p>Reset the dose at the start of a shift or call time. It survives Stop/Start. Guidance only, not a certified dosimeter.</p>"""),
    ("Theater calibration", """<h2>Theater / room level calibration</h2>
<ol><li>Open the <b>Theater calibration</b> panel and pick your interface as the Output, the format (5.1, 7.1, 7.1.4…) and
the room preset (cinema: 85 dBC screen, 82 dBC surrounds).</li>
<li>Check the <b>Out</b> column matches your patch (L R C LFE Ls Rs is the SMPTE order).</li>
<li>Put the mic at the main listening position, pointing at the ceiling (90° cal file).</li>
<li>Press <b>▶ Play</b>, adjust the amp/processor until it says <b>on target ✓</b>, then <b>Next ▶</b>.
Or press <b>Auto run</b> to measure every channel in one go and fix the red ones afterwards.</li></ol>
<p>Noise is -20 dBFS RMS (AES17) band-limited 500 Hz – 2 kHz, the SMPTE/Dolby standard. LFE is always band-limited to
20 – 120 Hz; its target is a starting point, confirm it on the RTA (about +10 dB in-band over the screen channels).</p>"""),
    ("Room EQ", "__EQ__"),
    ("Logging & reports", """<h2>Logging, schedules and reports</h2>
<ul><li><b>● Log</b> writes a CSV every 1 s / 10 s / 1 min with Leq, max, peak, LA10/90, dose, markers and the current
Planning Center item.</li>
<li><b>Schedule…</b> (Session panel) starts and stops logging by itself, e.g. every Sunday 08:30 – 12:30.</li>
<li><b>Auto-save folder</b> (Settings → Logging): point it at your OneDrive folder. Logs, the session CSV and the PDF
report are saved there in dated folders (YYYY / YYYY-MM-DD - plan) and OneDrive syncs them.</li>
<li><b>PDF report</b>: summary, time history, markers, alarms and per-item levels for the service.</li></ul>"""),
    ("Companion, NDI & web", """<h2>Stream Deck (Companion), NDI and web view</h2>
<ul><li><b>Companion</b>: enable it in Settings → Integrations. ServiceSPL sends OSC to Companion (port 12321) and sets
custom variables <code>spl_level</code>, <code>spl_leq1m</code>, <code>spl_leqlong</code>, <code>spl_peak</code>,
<code>spl_dose</code>, <code>spl_state</code>, <code>spl_color</code>, <code>spl_item</code>. Create those custom variables
in Companion, then use <code>$(custom:spl_level)</code> on a button. Or give a button location (page/row/column) and
ServiceSPL paints that button's text and color directly.</li>
<li><b>NDI</b>: the <b>NDI</b> toolbar button sends just the meter graphic as an NDI source named “ServiceSPL” for your
switcher. Needs the free NDI Tools/runtime installed. Tick “transparent background” for a key-able overlay.</li>
<li><b>Web view</b>: shows the big number in any browser on the network; <code>/api/json</code> and <code>/api/level</code>
give the numbers to other software.</li></ul>"""),
    ("Planning Center", """<h2>Planning Center Services Live: setup</h2>
<p>ServiceSPL follows <b>Services Live</b>: whenever whoever runs Live advances to the next item, ServiceSPL starts a new
item, stamps the log with it and keeps levels per item. One-time setup, about 3 minutes:</p>
<h3>1 · Make a Personal Access Token (in your web browser)</h3>
<ol><li>Sign in to Planning Center with an account that can <b>see the plans</b> you want to follow (Services
<i>Viewer</i> or higher on that service type). The token gets the same access as that account.</li>
<li>Open <b>https://api.planningcenteronline.com/oauth/applications</b> (or press <b>Open token page</b> in the panel).</li>
<li>Scroll to <b>Personal Access Tokens</b> and click <b>New Personal Access Token</b>.</li>
<li>Description: <code>ServiceSPL</code>. For the API version, keep the newest <b>Services</b> version offered. Click
<b>Submit</b>.</li>
<li>You'll see an <b>Application ID</b> and a <b>Secret</b>. Copy both now (keep the page open until step 2 works).</li></ol>
<h3>2 · Connect ServiceSPL</h3>
<ol><li>Show the panel: <b>View → Service · Planning Center Live</b>.</li>
<li>Paste the <b>Application ID</b> into the first box and the <b>Secret</b> into the second, then press <b>Connect</b>.
The status line should say <i>Connected · N service types</i>. On a Mac the secret goes into your Keychain, so you only do
this once.</li>
<li>Pick the <b>service type</b> (e.g. "Sunday Services") in the left drop-down. The right drop-down lists recent and
upcoming plans; today's plan is picked for you.</li>
<li>Set <b>Item target LAeq</b>: an item whose average goes over this turns red and raises an alarm (95 dB is a common
worship target).</li></ol>
<h3>3 · On the day</h3>
<ol><li>Start the meter (and logging, or let the <b>Schedule</b> do it).</li>
<li>Tick <b>Follow Live</b>. The status line says <i>Following Live: …</i>.</li>
<li>Someone opens <b>Services Live</b> for that plan, <b>takes control</b> and advances items as the service goes.
ServiceSPL checks every 3 seconds; <b>Now: …</b> shows the current item and its running LAeq.</li>
<li>Afterwards the panel table and the PDF report list LAeq / LCeq / LAFmax / LCpeak for every item.</li></ol>
<h3>If it doesn't work</h3>
<ul><li><b>401 / Unauthorized</b>: the ID or secret is wrong or the token was revoked. Make a new token.</li>
<li><b>Connected, but no service types or plans</b>: the account behind the token can't see them. Give it Viewer access.</li>
<li><b>Following, but the item never changes</b>: nobody has taken control in Services Live, or a different plan is
selected. Check the date in the plan drop-down.</li>
<li><b>Couldn't connect / timed out</b>: the Mac needs internet access to api.planningcenteronline.com.</li>
<li>No Planning Center at all? <b>Start item manually…</b> does the same with names you type.</li></ul>"""),
]


class GuideWizard(QtWidgets.QDialog):
    """Built-in step-by-step help, with an interactive room EQ page driven by the live RTA."""
    def __init__(self, main, start="Welcome"):
        super().__init__(main)
        self.main = main
        self.setWindowTitle("ServiceSPL guide")
        self.resize(940, 640)
        h = QtWidgets.QHBoxLayout(self)
        self.nav = QtWidgets.QListWidget()
        self.nav.setFixedWidth(190)
        self.nav.addItems([t for t, _ in GUIDE_PAGES])
        h.addWidget(self.nav)
        right = QtWidgets.QVBoxLayout()
        self.stack = QtWidgets.QStackedWidget()
        for title, html in GUIDE_PAGES:
            if html == "__EQ__":
                self.stack.addWidget(self._eq_page())
            else:
                tb = QtWidgets.QTextBrowser(); tb.setOpenExternalLinks(True); tb.setHtml(html)
                tb.setStyleSheet("font-size:11pt;")
                self.stack.addWidget(tb)
        right.addWidget(self.stack, 1)
        b = QtWidgets.QHBoxLayout()
        self.action = QtWidgets.QPushButton(""); self.action.clicked.connect(self.do_action)
        prev = QtWidgets.QPushButton("◀ Back"); prev.clicked.connect(lambda: self.nav.setCurrentRow(max(0, self.nav.currentRow() - 1)))
        nxt = QtWidgets.QPushButton("Next ▶"); nxt.clicked.connect(lambda: self.nav.setCurrentRow(min(self.nav.count() - 1, self.nav.currentRow() + 1)))
        close = QtWidgets.QPushButton("Close"); close.clicked.connect(self.accept)
        b.addWidget(self.action); b.addStretch(); b.addWidget(prev); b.addWidget(nxt); b.addWidget(close)
        right.addLayout(b)
        h.addLayout(right, 1)
        self.nav.currentRowChanged.connect(self._page)
        self.nav.setCurrentRow(max(0, [t for t, _ in GUIDE_PAGES].index(start) if start in [t for t, _ in GUIDE_PAGES] else 0))
        self.timer = QtCore.QTimer(self); self.timer.timeout.connect(self.poll); self.timer.start(200)
        self.capture = None
        self.positions = []

    ACTIONS = {"Welcome": "▶ Start the meter", "Calibrate your mic": "Open the calibration wizard…",
               "Monitoring a show": "Show the Session panel", "Hearing dose": "Show the dosimeter",
               "Theater calibration": "Show the theater panel", "Logging & reports": "Edit the schedule…",
               "Companion, NDI & web": "Open Settings…", "Planning Center": "Show the Planning Center panel"}

    def _page(self, i):
        self.stack.setCurrentIndex(i)
        title = GUIDE_PAGES[i][0]
        self.action.setText(self.ACTIONS.get(title, ""))
        self.action.setVisible(title in self.ACTIONS)

    def do_action(self):
        t = GUIDE_PAGES[self.nav.currentRow()][0]
        m = self.main
        if t == "Welcome" and not m.running():
            m.run_btn.click()
        elif t == "Calibrate your mic":
            m.calibrate()
        elif t in ("Monitoring a show", "Hearing dose", "Theater calibration", "Planning Center"):
            key = {"Monitoring a show": "session", "Hearing dose": "dose", "Theater calibration": "theater",
                   "Planning Center": "pco"}[t]
            m.docks[key].show(); m.docks[key].raise_()
        elif t == "Logging & reports":
            m.edit_schedule()
        elif t == "Companion, NDI & web":
            m.settings()

    def _eq_page(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        intro = QtWidgets.QLabel(
            "<h2>Room EQ with the RTA</h2>"
            "<ol><li>Pick the output/channel for the speaker you're tuning and press <b>Play pink noise</b> "
            "(full-band, -20 dBFS). Bring it up to a comfortable 80–85 dBC.</li>"
            "<li>Put the mic where people listen. Press <b>Measure position</b> (10 s average). Measure 3–5 seats "
            "and the app averages them, which stops you EQing one seat's quirks.</li>"
            "<li>Make the suggested cuts on your processor/console EQ, then <b>Clear</b> and measure again. "
            "Prefer cuts; big dips are usually room nulls that EQ can't fix.</li></ol>")
        intro.setWordWrap(True)
        v.addWidget(intro)
        row = QtWidgets.QHBoxLayout()
        self.eq_out = QtWidgets.QComboBox()
        th = self.main.theater
        for i in range(th.out.count()):
            self.eq_out.addItem(th.out.itemText(i), th.out.itemData(i))
        self.eq_out.setCurrentIndex(th.out.currentIndex())
        self.eq_ch = QtWidgets.QSpinBox(); self.eq_ch.setRange(1, 64); self.eq_ch.setPrefix("Out ")
        self.eq_target = QtWidgets.QComboBox(); self.eq_target.addItems(list(EQ_TARGETS))
        self.eq_target.currentTextChanged.connect(self.recompute)
        self.eq_play = QtWidgets.QPushButton("▶ Play pink noise"); self.eq_play.setCheckable(True)
        self.eq_play.clicked.connect(self.toggle_noise)
        row.addWidget(self.eq_out, 2); row.addWidget(self.eq_ch); row.addWidget(QtWidgets.QLabel("Target"))
        row.addWidget(self.eq_target, 2); row.addWidget(self.eq_play)
        v.addLayout(row)
        row2 = QtWidgets.QHBoxLayout()
        self.eq_meas = QtWidgets.QPushButton("Measure position (10 s)"); self.eq_meas.clicked.connect(self.measure)
        clr = QtWidgets.QPushButton("Clear"); clr.clicked.connect(self.clear)
        self.eq_status = QtWidgets.QLabel("No measurements yet")
        self.eq_prog = QtWidgets.QProgressBar(); self.eq_prog.setRange(0, 100); self.eq_prog.setValue(0)
        row2.addWidget(self.eq_meas); row2.addWidget(clr); row2.addWidget(self.eq_prog, 1); row2.addWidget(self.eq_status)
        v.addLayout(row2)
        self.eq_plot = EqPlot()
        v.addWidget(self.eq_plot, 1)
        self.eq_tips = QtWidgets.QListWidget()
        self.eq_tips.setMaximumHeight(130)
        v.addWidget(self.eq_tips)
        self.eq_player = NoisePlayer()
        self.eq_bank = None
        return w

    def toggle_noise(self, on):
        m = self.main
        if on:
            if not m.running():
                m.run_btn.click()
            try:
                fs = self.eq_player.open(self.eq_out.currentData(), self.eq_ch.value() - 1)
            except Exception as ex:
                self.eq_play.setChecked(False)
                QtWidgets.QMessageBox.critical(self, APP_NAME, f"Could not open the output:\n{ex}")
                return
            self.eq_player.set_buffer(make_pink(fs, 20, min(20000, fs * 0.45), -20.0, seed=7))
            self.eq_player.ch = self.eq_ch.value() - 1
            if m.sim is not None and self.eq_out.currentData() is None:
                m.sim.inject = self.eq_player
            self.eq_play.setText("■ Stop noise")
        else:
            self.eq_player.stop()
            if m.sim is not None and m.sim.inject is self.eq_player:
                m.sim.inject = None
            self.eq_play.setText("▶ Play pink noise")

    def measure(self):
        e = self.main.engine
        if e is None or not self.main.running():
            self.eq_status.setText("Press Start first")
            return
        self.eq_bank = ThirdOctaveBank(e.fs)
        e.extra.append(self.eq_bank)
        self.capture = time.time()
        self.eq_meas.setEnabled(False)

    def clear(self):
        self.positions = []
        self.eq_plot.set_data([], [], None)
        self.eq_tips.clear()
        self.eq_status.setText("Cleared")

    def poll(self):
        e = self.main.engine
        if self.capture is None or e is None:
            return
        el = time.time() - self.capture
        self.eq_prog.setValue(int(min(el / 10, 1) * 100))
        if el >= 10:
            self.capture = None
            if self.eq_bank in e.extra:
                e.extra.remove(self.eq_bank)
            self.eq_meas.setEnabled(True)
            lv = 10 * np.log10(np.maximum(self.eq_bank.levels("Leq"), 1e-20)) + e.offset
            self.positions.append(lv)
            self.eq_status.setText(f"{len(self.positions)} position(s) averaged")
            self.recompute()

    def recompute(self, *_):
        if not self.positions or self.eq_bank is None:
            return
        avg = 10 * np.log10(np.mean([10 ** (p / 10) for p in self.positions], axis=0))
        fc = np.array(self.eq_bank.fc)
        tgt = EQ_TARGETS[self.eq_target.currentText()](fc)
        dev, tips = eq_suggestions(fc, self.eq_bank.labels, avg, tgt)
        self.eq_plot.set_data(fc, self.eq_bank.labels, dev)
        self.eq_tips.clear()
        self.eq_tips.addItems(tips or ["Within ±2 dB of the target everywhere. Nice room! 🎉"])

    def done(self, r):
        self.timer.stop()
        self.toggle_noise(False)
        super().done(r)


class PcoPanel(QtWidgets.QWidget):
    """Planning Center Services Live: follow the current item, per-item levels."""
    def __init__(self, main):
        super().__init__()
        self.main = main
        c = main.cfg
        self.client = None
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(6, 6, 6, 6)
        g = QtWidgets.QGridLayout()
        self.app_id = QtWidgets.QLineEdit(c["pco_app_id"]); self.app_id.setPlaceholderText("Application ID")
        self.secret = QtWidgets.QLineEdit(); self.secret.setEchoMode(QtWidgets.QLineEdit.Password)
        self.secret.setPlaceholderText("Secret")
        self.secret.setText(secret_get("pco") or c.get("pco_secret", ""))
        con = QtWidgets.QPushButton("Connect"); con.clicked.connect(self.connect_pco)
        tok = QtWidgets.QPushButton("Open token page")
        tok.setToolTip("Create a Personal Access Token in your browser")
        tok.clicked.connect(lambda: QtGui.QDesktopServices.openUrl(
            QtCore.QUrl("https://api.planningcenteronline.com/oauth/applications")))
        hlp = QtWidgets.QPushButton("? Setup help"); hlp.clicked.connect(lambda: main.guide("Planning Center"))
        self.st = QtWidgets.QComboBox(); self.plan = QtWidgets.QComboBox()
        self.st.currentIndexChanged.connect(self.load_plans)
        self.follow = QtWidgets.QCheckBox("Follow Live"); self.follow.toggled.connect(self.toggle_follow)
        self.target = QtWidgets.QDoubleSpinBox(); self.target.setRange(60, 130); self.target.setValue(c["pco_target"])
        self.target.setPrefix("Item target LAeq "); self.target.setSuffix(" dB")
        self.target.valueChanged.connect(lambda x: c.__setitem__("pco_target", x))
        g.addWidget(self.app_id, 0, 0); g.addWidget(self.secret, 0, 1); g.addWidget(con, 0, 2)
        g.addWidget(self.st, 1, 0); g.addWidget(self.plan, 1, 1); g.addWidget(self.follow, 1, 2)
        g.addWidget(self.target, 2, 0)
        hb = QtWidgets.QHBoxLayout(); hb.addWidget(tok); hb.addWidget(hlp)
        g.addLayout(hb, 2, 1, 1, 2)
        v.addLayout(g)
        self.status = QtWidgets.QLabel("Not connected")
        self.now = QtWidgets.QLabel("—")
        self.now.setWordWrap(True)
        v.addWidget(self.status); v.addWidget(self.now)
        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Item", "Start", "Time", "LAeq", "LCeq", "LCpeak"])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        v.addWidget(self.table, 1)
        self.manual = QtWidgets.QPushButton("Start item manually…")
        self.manual.setToolTip("No Planning Center? Name the current segment yourself")
        self.manual.clicked.connect(self.manual_item)
        v.addWidget(self.manual)
        self.follower = PcoFollower()
        self.last_id = None
        self.flagged = set()
        self.refresh_theme()

    def refresh_theme(self):
        self.status.setStyleSheet(f"color:{T['dim']};")
        self.now.setStyleSheet("font-size:15pt; font-weight:700;")

    def connect_pco(self):
        a, s = self.app_id.text().strip(), self.secret.text().strip()
        if not (a and s):
            self.status.setText("Enter the Application ID and Secret of a Personal Access Token")
            return
        self.main.cfg["pco_app_id"] = a
        if not secret_set("pco", s):
            self.main.cfg["pco_secret"] = s
        self.main.save_cfg()
        self.client = PcoClient(a, s)
        try:
            sts = self.client.service_types()
        except Exception as e:
            hint = " (check the Application ID and Secret)" if "401" in str(e) else ""
            self.status.setText(f"Couldn't connect: {e}{hint}")
            return
        self.st.blockSignals(True); self.st.clear()
        for i, n in sts:
            self.st.addItem(n, i)
        want = self.st.findData(self.main.cfg.get("pco_service_type"))
        self.st.setCurrentIndex(want if want >= 0 else 0)
        self.st.blockSignals(False)
        self.status.setText(f"Connected · {len(sts)} service types")
        self.load_plans()

    def load_plans(self):
        if not self.client or self.st.currentData() is None:
            return
        self.main.cfg["pco_service_type"] = self.st.currentData()
        try:
            plans = self.client.plans(self.st.currentData())
        except Exception as e:
            self.status.setText(f"Couldn't load plans: {e}")
            return
        self.plan.clear()
        today = datetime.date.today().isoformat()
        pick = 0
        for i, (pid, label, sd_) in enumerate(plans):
            self.plan.addItem(label, pid)
            if sd_[:10] <= today:
                pick = i
            if sd_[:10] == today:
                pick = i
        self.plan.setCurrentIndex(pick)

    def toggle_follow(self, on):
        if on and self.client and self.plan.currentData():
            self.follower.follow(self.client, self.st.currentData(), self.plan.currentData())
            self.main.pco_plan_label = self.plan.currentText()
            self.status.setText(f"Following Live: {self.plan.currentText()}")
        else:
            self.follower.stop()
            if on:
                self.follow.setChecked(False)
                self.status.setText("Connect and pick a plan first")

    def manual_item(self):
        txt, ok = QtWidgets.QInputDialog.getText(self, "Current item", "Item name (e.g. Worship set, Sermon):")
        if ok and txt.strip():
            self.main.set_item(txt.strip())

    def tick(self):
        f = self.follower
        if f.on:
            iid, title = f.current
            if f.error:
                self.status.setText(f"Planning Center: {f.error}")
            if iid and iid != self.last_id:
                self.last_id = iid
                self.main.set_item(title)
        self.refresh_table()

    def refresh_table(self):
        s = self.main.session
        items = s.items
        if not items:
            self.now.setText("No item yet")
            return
        self.table.setRowCount(len(items))
        end = time.time()
        tgt = self.main.cfg["pco_target"]
        for i, (t0, title) in enumerate(items):
            t1 = items[i + 1][0] if i + 1 < len(items) else end
            st = item_stats(s, t0, t1)
            vals = [title, f"{datetime.datetime.fromtimestamp(t0):%H:%M}", hms(t1 - t0)]
            vals += [fmt(st[k]).strip() if st else "--.-" for k in ("LAeq", "LCeq", "LCpeak")]
            over = st is not None and np.isfinite(st["LAeq"]) and st["LAeq"] >= tgt
            for c, txt in enumerate(vals):
                it = self.table.item(i, c)
                if it is None:
                    it = QtWidgets.QTableWidgetItem(); self.table.setItem(i, c, it)
                it.setText(txt)
                it.setForeground(QtGui.QColor(T["red"] if over and c == 3 else T["fg"]))
            if over and i not in self.flagged and st["dur"] > 30:
                self.flagged.add(i)
                self.main.raise_alarm(f"'{title}' running LAeq {st['LAeq']:.1f} ≥ item target {tgt:.0f}")
        t0, title = items[-1]
        st = item_stats(s, t0, end)
        lv = st["LAeq"] if st else float("nan")
        col = T["red"] if np.isfinite(lv) and lv >= tgt else T["fg"]
        self.now.setText(f"Now: {title}   <span style='color:{col}; font-family:Menlo,monospace'>"
                         f"LAeq {fmt(lv).strip()}</span>")


# ---------------------------------------------------------- web view -----

WEB_PAGE = """<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>SPL</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
html,body{margin:0;height:100%;background:#121417;color:#e8e8e8;
  font-family:'Segoe UI','Helvetica Neue',Arial,sans-serif;overflow:hidden}
#wrap{display:flex;align-items:center;height:100%;padding:0 3vw;gap:3vw;box-sizing:border-box}
#thumb{height:min(70vh,28vw);flex:none}
#right{flex:1;text-align:right}
#val{font-family:Consolas,Menlo,monospace;font-weight:700;font-size:min(55vh,16vw);line-height:1;white-space:nowrap}
#sub{color:#8a9099;font-size:min(9vh,2.4vw);white-space:nowrap}
</style></head><body>
<div id="wrap"><div id="thumb"></div><div id="right"><div id="val">--.-</div><div id="sub">connecting…</div></div></div>
<script>
const SVG=__SVG__;
const BG={good:"#2ea866",bad:"#e03b3b",off:"#3a3f47"};
let last="";
function thumb(st){if(st===last)return;last=st;
 document.getElementById("thumb").innerHTML=SVG.replaceAll("{bg}",BG[st])
  .replaceAll("{xf}",st==="bad"?"rotate(180 50 50)":"").replace("<svg ","<svg height='100%' ");}
async function poll(){
 try{const r=await fetch("/data",{cache:"no-store"});const d=await r.json();
  const v=document.getElementById("val");
  v.textContent=d.value==null?"--.-":d.value.toFixed(1);v.style.color=d.color;
  document.getElementById("sub").textContent=d.sub;thumb(d.state);}
 catch(e){document.getElementById("sub").textContent="no connection to ServiceSPL";
  document.getElementById("val").style.color="#8a9099";thumb("off");}
 setTimeout(poll,200);}
poll();
</script></body></html>"""


class WebView:
    """Tiny read-only HTTP server: '/' = big-number page, '/data' = live JSON."""
    def __init__(self, port):
        self.state = {"value": None, "color": "#e8e8e8", "state": "off", "sub": "stopped"}
        self.api = {}
        page = WEB_PAGE.replace("__SVG__", json.dumps(THUMB_SVG)).encode()
        view = self

        class H(BaseHTTPRequestHandler):
            def do_GET(self):
                if self.path.startswith("/api/json"):
                    body, ctype = json.dumps(view.api).encode(), "application/json"
                elif self.path.startswith("/api/level"):
                    body, ctype = str(view.api.get("level", "--.-")).encode(), "text/plain; charset=utf-8"
                elif self.path.startswith("/data"):
                    body, ctype = json.dumps(view.state).encode(), "application/json"
                else:
                    body, ctype = page, "text/html; charset=utf-8"
                self.send_response(200)
                self.send_header("Content-Type", ctype)
                self.send_header("Cache-Control", "no-store")
                self.send_header("Access-Control-Allow-Origin", "*")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *a):
                pass

        self.httpd = ThreadingHTTPServer(("0.0.0.0", port), H)
        self.httpd.daemon_threads = True
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.port = port

    def urls(self):
        ips = set()
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            s.connect(("10.255.255.255", 1))
            ips.add(s.getsockname()[0])
            s.close()
        except Exception:
            pass
        try:
            for info in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET):
                ips.add(info[4][0])
        except Exception:
            pass
        ips = sorted(i for i in ips if not i.startswith("127."))
        return [f"http://{i}:{self.port}" for i in ips] or [f"http://localhost:{self.port}"]

    def stop(self):
        self.httpd.shutdown()
        self.httpd.server_close()


# ------------------------------------------------------------ settings ----

class SettingsDialog(QtWidgets.QDialog):
    def __init__(self, cfg, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumWidth(560)
        outer = QtWidgets.QVBoxLayout(self)
        tabs = QtWidgets.QTabWidget()
        outer.addWidget(tabs)

        def page(name):
            w = QtWidgets.QWidget(); f = QtWidgets.QFormLayout(w); tabs.addTab(w, name); return f

        # --- meter
        f = page("Meter")
        self.w = QtWidgets.QComboBox(); self.w.addItems(["A", "C", "Z"]); self.w.setCurrentText(cfg["main_weight"])
        self.t = QtWidgets.QComboBox(); self.t.addItems(["Fast", "Slow"])
        self.t.setCurrentText("Fast" if cfg["main_time"] == "F" else "Slow")
        self.warn = QtWidgets.QDoubleSpinBox(); self.warn.setRange(0, 30); self.warn.setValue(cfg["warn_gap"])
        self.warn.setSuffix(" dB below red"); self.warn.setSpecialValueText("off")
        self.alarm = QtWidgets.QDoubleSpinBox(); self.alarm.setRange(40, 140); self.alarm.setValue(cfg["alarm"])
        self.long = QtWidgets.QSpinBox(); self.long.setRange(1, 120); self.long.setValue(cfg["long_leq_min"]); self.long.setSuffix(" min")
        self.src = QtWidgets.QComboBox(); self.src.addItems(["Live meter", "LAeq 1 min"])
        self.src.setCurrentIndex(0 if cfg["thumb_src"] == "live" else 1)
        f.addRow("Meter weighting", self.w)
        f.addRow("Meter speed", self.t)
        f.addRow("Yellow warning starts", self.warn)
        f.addRow("Turn red / thumbs down at (dB)", self.alarm)
        f.addRow("Thumb / color follows", self.src)
        f.addRow("Long Leq window", self.long)

        # --- display
        f = page("Display")
        self.theme = QtWidgets.QComboBox(); self.theme.addItems(["Follow system", "Dark", "Light"])
        self.theme.setCurrentIndex({"auto": 0, "dark": 1, "light": 2}[cfg["theme"]])
        self.style = QtWidgets.QComboBox(); self.style.addItems(["Big number + history", "Big number + bar"])
        self.style.setCurrentIndex(0 if cfg["meter_style"] == "history" else 1)
        self.hist = QtWidgets.QSpinBox(); self.hist.setRange(1, 30); self.hist.setSuffix(" min")
        self.hist.setValue(max(1, cfg["history_secs"] // 60))
        self.menubar = QtWidgets.QCheckBox("Show the live level in the menu bar (● 87.3)")
        self.menubar.setChecked(cfg["menubar"])
        f.addRow("Theme", self.theme)
        f.addRow("Meter style", self.style)
        f.addRow("History length", self.hist)
        f.addRow("", self.menubar)

        # --- logging
        f = page("Logging && reports")
        self.li = QtWidgets.QComboBox(); self.li.addItems(list(LOG_INTERVALS))
        self.li.setCurrentText(next(k for k, v in LOG_INTERVALS.items() if v == cfg["log_interval"]))
        self.logdir = QtWidgets.QLineEdit(cfg["log_dir"])
        row = QtWidgets.QHBoxLayout(); row.addWidget(self.logdir)
        b = QtWidgets.QPushButton("Browse…"); b.clicked.connect(lambda: self._browse(self.logdir)); row.addWidget(b)
        self.auto = QtWidgets.QLineEdit(cfg["autosave_dir"])
        self.auto.setPlaceholderText("off (e.g. your OneDrive folder)")
        row2 = QtWidgets.QHBoxLayout(); row2.addWidget(self.auto)
        b2 = QtWidgets.QPushButton("Browse…"); b2.clicked.connect(lambda: self._browse(self.auto)); row2.addWidget(b2)
        od = QtWidgets.QPushButton("Find OneDrive"); od.clicked.connect(self._onedrive); row2.addWidget(od)
        self.auto_rep = QtWidgets.QCheckBox("Save the PDF report + session CSV automatically when logging stops")
        self.auto_rep.setChecked(cfg["autosave_reports"])
        f.addRow("Log interval", self.li)
        f.addRow("Log folder", row)
        f.addRow("Auto-save folder", row2)
        f.addRow("", self.auto_rep)
        f.addRow("", QtWidgets.QLabel("With an auto-save folder, logs and reports go to\n"
                                      "<folder>/ServiceSPL/YYYY/YYYY-MM-DD[ - plan]/ instead of the log folder."))

        # --- integrations
        f = page("Integrations")
        self.comp = QtWidgets.QCheckBox("Send levels to Bitfocus Companion (OSC)")
        self.comp.setChecked(cfg["companion_on"])
        self.ch = QtWidgets.QLineEdit(cfg["companion_host"])
        self.cp = QtWidgets.QSpinBox(); self.cp.setRange(1, 65535); self.cp.setValue(int(cfg["companion_port"]))
        self.cb = QtWidgets.QLineEdit(cfg["companion_button"]); self.cb.setPlaceholderText("optional page/row/column, e.g. 1/0/0")
        f.addRow(self.comp)
        f.addRow("Companion computer", self.ch)
        f.addRow("Companion OSC port", self.cp)
        f.addRow("Paint this button", self.cb)
        self.port = QtWidgets.QSpinBox(); self.port.setRange(1024, 65535); self.port.setValue(cfg["web_port"])
        f.addRow("Web view / API port", self.port)
        self.ndi_name = QtWidgets.QLineEdit(cfg["ndi_name"])
        self.ndi_size = QtWidgets.QComboBox(); self.ndi_size.addItems(["1280x720", "1920x1080", "960x540", "1280x360"])
        self.ndi_size.setCurrentText(cfg["ndi_size"])
        self.ndi_tr = QtWidgets.QCheckBox("Transparent background (for keying over video)")
        self.ndi_tr.setChecked(cfg["ndi_transparent"])
        self.ndi_lib = QtWidgets.QLineEdit(cfg["ndi_lib"]); self.ndi_lib.setPlaceholderText("auto-detect")
        row3 = QtWidgets.QHBoxLayout(); row3.addWidget(self.ndi_lib)
        b3 = QtWidgets.QPushButton("Browse…"); b3.clicked.connect(self._ndi_lib); row3.addWidget(b3)
        f.addRow("NDI source name", self.ndi_name)
        f.addRow("NDI size", self.ndi_size)
        f.addRow("", self.ndi_tr)
        f.addRow("NDI library", row3)

        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        outer.addWidget(bb)

    def _browse(self, edit):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose folder", edit.text() or os.path.expanduser("~"))
        if d:
            edit.setText(d)

    def _onedrive(self):
        d = find_onedrive()
        if d:
            self.auto.setText(d)
        else:
            QtWidgets.QMessageBox.information(self, APP_NAME, "Couldn't find a OneDrive folder. Install OneDrive, "
                                              "sign in, then press Browse and pick it (Finder → OneDrive).")

    def _ndi_lib(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "NDI library", "/", "NDI runtime (libndi*.dylib *.dll *.so*)")
        if path:
            self.ndi_lib.setText(path)

    def apply(self, cfg):
        cfg.update(main_weight=self.w.currentText(), main_time=self.t.currentText()[0],
                   warn_gap=self.warn.value(), alarm=self.alarm.value(), long_leq_min=self.long.value(),
                   log_interval=LOG_INTERVALS[self.li.currentText()], log_dir=self.logdir.text(),
                   thumb_src="live" if self.src.currentIndex() == 0 else "leq1m",
                   web_port=self.port.value(), theme=["auto", "dark", "light"][self.theme.currentIndex()],
                   meter_style="history" if self.style.currentIndex() == 0 else "bar",
                   history_secs=self.hist.value() * 60, menubar=self.menubar.isChecked(),
                   autosave_dir=self.auto.text().strip(), autosave_reports=self.auto_rep.isChecked(),
                   companion_on=self.comp.isChecked(), companion_host=self.ch.text().strip() or "127.0.0.1",
                   companion_port=self.cp.value(), companion_button=self.cb.text().strip(),
                   ndi_name=self.ndi_name.text().strip() or "ServiceSPL", ndi_size=self.ndi_size.currentText(),
                   ndi_transparent=self.ndi_tr.isChecked(), ndi_lib=self.ndi_lib.text().strip())


# ---------------------------------------------------------- PDF report ----

def write_pdf_report(path, sess, cfg, profile_name, prof, plan_label=""):
    """One-page (plus overflow) session report: summary, graph, markers and events."""
    saved = dict(T)
    T.update(PALETTES["light"])        # reports are always printed light
    try:
        w = QtGui.QPdfWriter(path)
        w.setPageSize(QtGui.QPageSize(QtGui.QPageSize.Letter))
        w.setResolution(144)
        w.setPageMargins(QtCore.QMarginsF(14, 14, 14, 14), QtGui.QPageLayout.Millimeter)
        w.setTitle("SPL session report")
        w.setCreator(f"{APP_NAME} {APP_VERSION}")
        p = QtGui.QPainter(w)
        W, H = w.width(), w.height()
        y = 0

        def text(s, size=10, bold=False, color="#1b1e23", x=0, width=None, align=QtCore.Qt.AlignLeft):
            nonlocal y
            f = QtGui.QFont("Helvetica", size); f.setBold(bold)
            p.setFont(f); p.setPen(QtGui.QColor(color))
            fm = QtGui.QFontMetricsF(f)
            p.drawText(QtCore.QRectF(x, y, width or W - x, fm.height() * 1.2), align | QtCore.Qt.AlignVCenter, s)
            return fm.height() * 1.25

        y += text("SPL Session Report", 20, True)
        y += text(" · ".join(x for x in (cfg.get("venue") or ("" if plan_label else "Untitled session"), plan_label,
                                          cfg.get("operator") or "") if x), 12, color="#5f6670")
        t = sess.col("time")
        start = datetime.datetime.fromtimestamp(t[0] if len(t) else sess.start)
        end = datetime.datetime.fromtimestamp(t[-1] if len(t) else time.time())
        y += text(f"{start:%A %B %d, %Y  %H:%M} – {end:%H:%M}   ·   Mic profile: {profile_name}   ·   "
                  f"{prof.get('cal_info') if prof.get('calibrated') else 'NOT CALIBRATED'}", 9, color="#5f6670")
        y += 10
        s = sess.summary()
        if s:
            items = [("Duration", hms(s["dur"])), ("LAeq", fmt(s["LAeq"])), ("LCeq", fmt(s["LCeq"])),
                     ("LC − LA", fmt(s["LCeq"] - s["LAeq"])), ("LAFmax", fmt(s["LAFmax"])),
                     ("LCpeak", fmt(s["LCpeak"])), ("LA10", fmt(s["L10"])), ("LA50", fmt(s["L50"])),
                     ("LA90", fmt(s["L90"])), ("Max LAeq 1 min", fmt(s["maxLAeq1m"])),
                     (f"Max LAeq {cfg['long_leq_min']} min", fmt(s["maxLong"])),
                     ("Dose OSHA HC / NIOSH", f"{s['doseOSHA']:.1f} % / {s['doseNIOSH']:.1f} %")]
            cols = 4
            cw = W / cols
            box_h = 70
            for i, (k, v) in enumerate(items):
                cx, cy = (i % cols) * cw, y + (i // cols) * (box_h + 8)
                p.setPen(QtGui.QPen(QtGui.QColor("#d3d7dd"), 1)); p.setBrush(QtGui.QColor("#f6f7f9"))
                p.drawRoundedRect(QtCore.QRectF(cx + 2, cy, cw - 8, box_h), 6, 6)
                f = QtGui.QFont("Helvetica", 8); p.setFont(f); p.setPen(QtGui.QColor("#5f6670"))
                p.drawText(QtCore.QRectF(cx + 10, cy + 4, cw - 20, 22), QtCore.Qt.AlignLeft, k)
                f = QtGui.QFont("Menlo", 15); f.setBold(True); p.setFont(f); p.setPen(QtGui.QColor("#1b1e23"))
                p.drawText(QtCore.QRectF(cx + 10, cy + 24, cw - 20, 40), QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter,
                           v.strip())
            y += math.ceil(len(items) / cols) * (box_h + 8) + 10
        gh = H * 0.30
        paint_session_graph(p, QtCore.QRectF(0, y, W, gh), sess, 0, cfg["alarm"], cfg["alarm"] - cfg["warn_gap"],
                            small_font=7)
        y += gh + 14
        if sess.items:
            y += text("Service items", 12, True)
            cx = [0, W * 0.42, W * 0.52, W * 0.62, W * 0.72, W * 0.82, W * 0.91]
            hdr = ["Item", "Start", "Time", "LAeq", "LCeq", "LAFmax", "LCpeak"]
            f = QtGui.QFont("Helvetica", 8); f.setBold(True); p.setFont(f); p.setPen(QtGui.QColor("#5f6670"))
            for x0, h in zip(cx, hdr):
                p.drawText(QtCore.QRectF(x0, y, W * 0.4, 24), QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter, h)
            y += 26
            tend = t[-1] + 1 if len(t) else time.time()
            for i, (t0, title) in enumerate(sess.items):
                if y > H - 30:
                    w.newPage(); y = 0
                t1 = sess.items[i + 1][0] if i + 1 < len(sess.items) else tend
                st = item_stats(sess, t0, t1)
                vals = [title[:48], f"{datetime.datetime.fromtimestamp(t0):%H:%M}", hms(t1 - t0)]
                vals += [fmt(st[k]).strip() if st else "--.-" for k in ("LAeq", "LCeq", "LAFmax", "LCpeak")]
                over = st is not None and np.isfinite(st["LAeq"]) and st["LAeq"] >= cfg.get("pco_target", 999)
                f = QtGui.QFont("Helvetica", 9); p.setFont(f)
                for j, (x0, v) in enumerate(zip(cx, vals)):
                    p.setPen(QtGui.QColor("#dc3030" if over and j == 3 else "#1b1e23"))
                    p.drawText(QtCore.QRectF(x0, y, W * 0.4, 24), QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter, v)
                y += 24
            y += 10
        rows = [(mt, "◆ " + txt) for mt, txt in sess.markers] + [(et, "· " + txt) for et, txt in sess.events]
        rows.sort()
        y += text("Markers and events", 12, True)
        if not rows:
            y += text("None", 9, color="#5f6670")
        for rt, txt in rows:
            if y > H - 30:
                w.newPage(); y = 0
            y += text(f"{datetime.datetime.fromtimestamp(rt):%H:%M:%S}   {txt}", 9)
        y += 8
        if y > H - 60:
            w.newPage(); y = 0
        y += text("Levels in dB re 20 µPa. A/C weighting and Fast/Slow per IEC 61672-1. LCpeak is 4x oversampled. "
                  "Dose figures are for guidance only.", 7, color="#5f6670")
        y += text(f"Generated by {APP_NAME} {APP_VERSION} on {datetime.datetime.now():%Y-%m-%d %H:%M}", 7, color="#5f6670")
        p.end()
    finally:
        T.clear(); T.update(saved)


# ---------------------------------------------------------- main window ---

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.setDockOptions(QtWidgets.QMainWindow.AnimatedDocks | QtWidgets.QMainWindow.AllowTabbedDocks)
        self.cfg = self.load_cfg()
        self.apply_theme()
        self.engine = self.stream = self.sim = None
        self.q = queue.Queue()
        self.xruns = 0
        self.log_file = self.log_writer = None
        self.log_rows = 0
        self.pending = []      # per-second records waiting to be logged
        self.pending_notes = []
        self.session = Session()
        self.alarm_armed = {}
        self.dose_flags = set()
        self.last_rec = {}
        self.bad_until = 0.0
        self.last_1m = float("nan")
        self.mini = None
        self.web = None
        self.menubar_meter = None
        self.companion = CompanionLink()
        self.ndi = None
        self.card = {"level": "--.-", "color": "#8a9099", "state": "off", "sub": "stopped"}
        self.current_item = ""
        self.pco_plan_label = ""
        self.sched_key = None
        self.sched_skip = set()
        self._sched_stopping = False
        self._last_sec = 0
        self._build_ui()
        self.populate_devices()
        self.restore_layout()
        self.update_menubar_meter()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.tick)
        self.timer.start(50)
        self.ndi_timer = QtCore.QTimer(self)
        self.ndi_timer.timeout.connect(self.send_ndi)
        try:
            QtGui.QGuiApplication.styleHints().colorSchemeChanged.connect(lambda *_: self.on_theme_change())
        except Exception:
            pass

    # ---------- config
    def load_cfg(self):
        cfg = copy.deepcopy(DEFAULTS)
        saved = {}
        for path in (CFG_PATH, OLD_CFG_PATH):      # carry settings + mic profiles over from SPL Monitor
            try:
                with open(path) as f:
                    saved = json.load(f)
                break
            except Exception:
                continue
        if saved and "profiles" not in saved:      # v1 config: turn the single calibration into a profile
            saved["profiles"] = {"Default": dict(DEFAULTS["profiles"]["Default"],
                                                 offset=saved.get("cal_offset", 120.0),
                                                 calibrated=saved.get("calibrated", False),
                                                 cal_info=saved.get("cal_info", ""),
                                                 method="calibrator" if saved.get("calibrated") else "")}
        for k in ("cal_offset", "calibrated", "cal_info"):
            saved.pop(k, None)
        cfg.update(saved)
        for k in ("theater",):
            merged = copy.deepcopy(DEFAULTS[k]); merged.update(cfg.get(k) or {}); cfg[k] = merged
        if cfg["profile"] not in cfg["profiles"]:
            cfg["profile"] = next(iter(cfg["profiles"]))
        return cfg

    def save_cfg(self):
        try:
            os.makedirs(CFG_DIR, exist_ok=True)
            tmp = CFG_PATH + ".tmp"
            with open(tmp, "w") as f:
                json.dump(self.cfg, f, indent=2)
            os.replace(tmp, CFG_PATH)
        except Exception:
            pass

    def active_profile(self):
        return self.cfg["profiles"][self.cfg["profile"]]

    # ---------- theme
    def apply_theme(self):
        mode = self.cfg["theme"]
        dark = system_is_dark() if mode == "auto" else mode == "dark"
        T.clear(); T.update(PALETTES["dark" if dark else "light"])
        QtWidgets.QApplication.instance().setStyleSheet(style_sheet())

    def on_theme_change(self):
        self.apply_theme()
        for t in self.meter.tiles():
            t.refresh_theme()
        self.meter.apply_style()
        for w in (self.dose_panel, self.theater, self.mini, self.pco_panel):
            if w is not None:
                w.refresh_theme()
        self.apply_cfg_to_ui()
        for w in self.findChildren(QtWidgets.QWidget):
            w.update()
        if self.menubar_meter:
            self.menubar_meter._last = ("", "")

    # ---------- UI
    def _build_ui(self):
        tb = QtWidgets.QToolBar("Main"); tb.setObjectName("toolbar"); tb.setMovable(False); self.addToolBar(tb)
        self.dev = QtWidgets.QComboBox(); self.dev.setMinimumWidth(150)
        self.dev.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.dev.setMinimumContentsLength(16)
        self.ch = QtWidgets.QComboBox(); self.ch.setMinimumWidth(80)
        self.dev.currentIndexChanged.connect(self.fill_channels)
        self.dev.currentIndexChanged.connect(self.auto_profile)
        self.run_btn = QtWidgets.QPushButton("▶  Start"); self.run_btn.setCheckable(True)
        self.run_btn.clicked.connect(self.toggle_run)
        self.prof_cb = QtWidgets.QComboBox(); self.prof_cb.setMinimumWidth(110)
        self.prof_cb.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.prof_cb.setMinimumContentsLength(10)
        self.prof_cb.setToolTip("Mic profile (calibration)")
        self.prof_cb.activated.connect(lambda i: self.set_profile(self.prof_cb.itemText(i)))
        cal = QtWidgets.QPushButton("Calibrate…"); cal.clicked.connect(self.calibrate)
        rst = QtWidgets.QPushButton("Reset"); rst.setToolTip("Reset Leq / Max / Peak / Ln (not the dose or session)")
        rst.clicked.connect(self.reset_stats)
        self.log_btn = QtWidgets.QPushButton("●  Log"); self.log_btn.setObjectName("log")
        self.log_btn.setCheckable(True); self.log_btn.clicked.connect(self.toggle_log)
        mk = QtWidgets.QPushButton("◆ Marker"); mk.setToolTip("Drop a marker in the session and log (M)")
        mk.clicked.connect(lambda: self.add_marker())
        for w in (self.dev, self.ch, self.run_btn, self.prof_cb, cal, rst, self.log_btn, mk):
            tb.addWidget(w)
        tb.addSeparator()
        lbl = QtWidgets.QLabel("  Red at ")
        self.red_lbl = lbl
        self.red_at = QtWidgets.QDoubleSpinBox()
        self.red_at.setRange(50, 130); self.red_at.setDecimals(1); self.red_at.setSingleStep(1)
        self.red_at.setSuffix(" dB"); self.red_at.setValue(self.cfg["alarm"])
        self.red_at.setToolTip("Readings at or above this level turn red and show thumbs down")
        self.red_at.valueChanged.connect(self.set_red_at)
        tb.addWidget(lbl); tb.addWidget(self.red_at)
        tb.addSeparator()
        mini_btn = QtWidgets.QPushButton("Mini ⤡")
        mini_btn.setToolTip("Shrink to a small always-on-top box")
        mini_btn.clicked.connect(self.enter_mini)
        tb.addWidget(mini_btn)
        self.web_btn = QtWidgets.QPushButton("Web view")
        self.web_btn.setCheckable(True)
        self.web_btn.setToolTip("Show the meter in a browser on another computer (e.g. the Proton)")
        self.web_btn.clicked.connect(self.toggle_web)
        tb.addWidget(self.web_btn)
        self.ndi_btn = QtWidgets.QPushButton("NDI")
        self.ndi_btn.setCheckable(True)
        self.ndi_btn.setToolTip("Send just the meter graphic as an NDI source for your video switcher")
        self.ndi_btn.clicked.connect(self.toggle_ndi)
        tb.addWidget(self.ndi_btn)
        st = QtWidgets.QPushButton("Settings…"); st.clicked.connect(self.settings)
        tb.addWidget(st)
        gd = QtWidgets.QPushButton("? Guide"); gd.setToolTip("Step-by-step guide and room EQ helper")
        gd.clicked.connect(lambda: self.guide())
        tb.addWidget(gd)
        self.cal_lbl = QtWidgets.QLabel()
        self.status = self.statusBar()
        self.status.addPermanentWidget(self.cal_lbl)

        self.meter = MeterPanel(self)
        self.setCentralWidget(self.meter)
        self.thumb = self.meter.thumb
        self.rta_panel = RtaPanel(self)
        self.dose_panel = DosePanel(self)
        self.theater = TheaterPanel(self)
        self.session_panel = SessionPanel(self)
        self.pco_panel = PcoPanel(self)
        self.docks = {}
        for name, title, w, area in (("rta", "RTA · 1/3 octave + spectrogram", self.rta_panel, QtCore.Qt.RightDockWidgetArea),
                                     ("dose", "Dosimeter · OSHA / NIOSH", self.dose_panel, QtCore.Qt.RightDockWidgetArea),
                                     ("theater", "Theater calibration", self.theater, QtCore.Qt.RightDockWidgetArea),
                                     ("session", "Session · markers · alarms · export", self.session_panel,
                                      QtCore.Qt.BottomDockWidgetArea),
                                     ("pco", "Service · Planning Center Live", self.pco_panel,
                                      QtCore.Qt.BottomDockWidgetArea)):
            d = QtWidgets.QDockWidget(title, self)
            d.setObjectName("dock_" + name)
            sc = QtWidgets.QScrollArea()          # panels scroll instead of forcing the window bigger
            sc.setWidgetResizable(True)
            sc.setFrameShape(QtWidgets.QFrame.NoFrame)
            sc.setWidget(w)
            sc.setMinimumSize(200, 80)
            d.setWidget(sc)
            d.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable | QtWidgets.QDockWidget.DockWidgetClosable
                          | QtWidgets.QDockWidget.DockWidgetFloatable)
            self.addDockWidget(area, d)
            self.docks[name] = d
        self.splitDockWidget(self.docks["rta"], self.docks["dose"], QtCore.Qt.Vertical)
        self.tabifyDockWidget(self.docks["dose"], self.docks["theater"])
        self.docks["dose"].raise_()
        self.splitDockWidget(self.docks["session"], self.docks["pco"], QtCore.Qt.Horizontal)
        if not self.cfg.get("dock_state"):
            QtCore.QTimer.singleShot(0, self.default_dock_sizes)

        mb = self.menuBar()
        fm = mb.addMenu("&File")
        fm.addAction("Export session CSV…", lambda: self.export_csv())
        fm.addAction("Export PDF report…", lambda: self.export_pdf())
        fm.addAction("New session", lambda: self.new_session())
        fm.addSeparator()
        fm.addAction("Settings…", self.settings)
        vm = mb.addMenu("&View")
        for d in self.docks.values():
            vm.addAction(d.toggleViewAction())
        vm.addSeparator()
        vm.addAction("Show all panels", lambda: [d.show() for d in self.docks.values()])
        vm.addAction("Hide all panels", lambda: [d.hide() for d in self.docks.values()])
        vm.addAction("Mini meter", self.enter_mini)
        vm.addAction("Quarter-screen size", self.quarter_size)
        vm.addAction("Fit to screen", lambda: (self.resize(10_000, 10_000), self.fit_to_screen()))
        tm = mb.addMenu("&Tools")
        tm.addAction("Calibration wizard…", self.calibrate)
        tm.addAction("Drop marker", lambda: self.add_marker(), QtGui.QKeySequence("M"))
        tm.addAction("Marker with note…", self.add_marker_note, QtGui.QKeySequence("Shift+M"))
        tm.addAction("Reset Leq / Max / Peak", self.reset_stats)
        tm.addAction("Reset dose", self.dose_panel.reset)
        tm.addAction("Logging schedule…", self.edit_schedule)
        tm.addAction("Start service item…", self.pco_panel.manual_item)
        hm = mb.addMenu("&Help")
        hm.addAction("Guide…", lambda: self.guide())
        hm.addAction("Room EQ helper…", lambda: self.guide("Room EQ"))
        hm.addAction("Theater calibration help", lambda: self.guide("Theater calibration"))
        hm.addSeparator()
        about = hm.addAction(f"About {APP_NAME}", self.about)
        about.setMenuRole(QtGui.QAction.AboutRole)
        self.refresh_profiles()
        self.apply_cfg_to_ui()

    def restore_layout(self):
        scr = (self.screen() or QtGui.QGuiApplication.primaryScreen()).availableGeometry()
        try:
            if self.cfg.get("win_geom"):
                self.restoreGeometry(QtCore.QByteArray(base64.b64decode(self.cfg["win_geom"])))
            else:
                self.resize(int(scr.width() * 0.85), int(scr.height() * 0.85))
                self.move(scr.center() - self.rect().center())
            if self.cfg.get("dock_state"):
                self.restoreState(QtCore.QByteArray(base64.b64decode(self.cfg["dock_state"])))
        except Exception:
            self.resize(int(scr.width() * 0.85), int(scr.height() * 0.85))
        self.fit_to_screen()

    def default_dock_sizes(self):
        w, h = self.width(), self.height()
        self.resizeDocks([self.docks["rta"]], [int(w * 0.32)], QtCore.Qt.Horizontal)
        self.resizeDocks([self.docks["session"], self.docks["pco"]], [int(w * 0.62), int(w * 0.38)], QtCore.Qt.Horizontal)
        self.resizeDocks([self.docks["rta"], self.docks["dose"]], [int(h * 0.5), int(h * 0.4)], QtCore.Qt.Vertical)
        self.resizeDocks([self.docks["session"]], [int(h * 0.28)], QtCore.Qt.Vertical)

    def quarter_size(self):
        """Half the screen's width and height (1/4 of its area), parked top-right."""
        scr = (self.screen() or QtGui.QGuiApplication.primaryScreen()).availableGeometry()
        self.showNormal()
        self.resize(scr.width() // 2, scr.height() // 2)
        self.move(scr.right() - self.frameGeometry().width() + 1, scr.top())

    def fit_to_screen(self):
        """Never open taller or wider than the screen (e.g. after moving from a bigger display)."""
        scr = (self.screen() or QtGui.QGuiApplication.primaryScreen()).availableGeometry()
        g = self.frameGeometry()
        w, h = min(g.width(), scr.width()), min(g.height(), scr.height())
        if (w, h) != (g.width(), g.height()):
            self.resize(w - (g.width() - self.width()), h - (g.height() - self.height()))
        g = self.frameGeometry()
        x = min(max(g.left(), scr.left()), scr.right() - g.width() + 1)
        y = min(max(g.top(), scr.top()), scr.bottom() - g.height() + 1)
        self.move(x, y)

    def refresh_profiles(self):
        self.prof_cb.blockSignals(True)
        self.prof_cb.clear()
        self.prof_cb.addItems(list(self.cfg["profiles"]))
        self.prof_cb.setCurrentText(self.cfg["profile"])
        self.prof_cb.blockSignals(False)

    def set_profile(self, name):
        if name not in self.cfg["profiles"]:
            return
        self.cfg["profile"] = name
        self.refresh_profiles()
        prof = self.active_profile()
        if self.engine is not None:
            self.engine.offset = prof["offset"]
            pts = prof.get("cal_points") if prof.get("use_curve") else None
            self.engine.calfilt = CalFilter(pts, self.engine.fs) if pts else None
        self.session.cal_info = prof.get("cal_info", "")
        self.apply_cfg_to_ui()
        self.save_cfg()

    def auto_profile(self):
        dev = self.dev.currentText()
        for name, p in self.cfg["profiles"].items():
            if p.get("auto_select") and p.get("device_hint") and p["device_hint"] == dev:
                if name != self.cfg["profile"]:
                    self.set_profile(name)
                    self.status.showMessage(f"Mic profile '{name}' selected for this input")
                return

    def apply_cfg_to_ui(self):
        c = self.cfg
        m = self.meter
        m.main_tile.set(None, title=f"L{c['main_weight']}{c['main_time']}   "
                                    f"{c['main_weight']}-weighted · {'Fast' if c['main_time'] == 'F' else 'Slow'}")
        m.t_long.set(None, title=f"LAeq {c['long_leq_min']} min")
        m.t_clong.set(None, title=f"LCeq {c['long_leq_min']} min")
        m.bar.warn, m.bar.alarm = c["alarm"] - c["warn_gap"], c["alarm"]
        m.history.warn, m.history.alarm = c["alarm"] - c["warn_gap"], c["alarm"]
        m.history.secs = c["history_secs"]
        m.apply_style()
        self.red_lbl.setStyleSheet(f"color:{T['red']}; font-weight:600;")
        self.red_at.blockSignals(True); self.red_at.setValue(c["alarm"]); self.red_at.blockSignals(False)
        prof = self.active_profile()
        if prof["calibrated"]:
            age = ""
            try:
                days = (datetime.datetime.now() - datetime.datetime.fromisoformat(prof["date"])).days
                age = f" · {days} days old" if days >= 30 else ""
            except Exception:
                pass
            curve = " · freq-corrected" if prof.get("use_curve") and prof.get("cal_points") else ""
            self.cal_lbl.setText(f"  CAL ✓ {self.cfg['profile']}: {prof['cal_info']}{curve}{age}  ")
            self.cal_lbl.setStyleSheet(f"color:{T['yellow'] if age else T['green']};")
        else:
            self.cal_lbl.setText(f"  {self.cfg['profile']}: NOT CALIBRATED  ")
            self.cal_lbl.setStyleSheet(f"color:{T['yellow']}; font-weight:600;")

    def fill_channels(self):
        """Channel list follows the selected device (2 for a Scarlett, up to 128+ for SoundGrid)."""
        dev = self.dev.currentData()
        n = 1
        if dev is not None and sd is not None:
            try:
                n = sd.query_devices(dev)["max_input_channels"]
            except Exception:
                pass
        want = self.ch.currentIndex() + 1 if self.ch.count() else self.cfg["channel"]
        self.ch.blockSignals(True)
        self.ch.clear()
        self.ch.addItems([f"Ch {i}" for i in range(1, n + 1)])
        self.ch.setCurrentIndex(min(max(want, 1), n) - 1)
        self.ch.blockSignals(False)

    def populate_devices(self):
        """List inputs; auto-pick a USB measurement mic, then the Scarlett, if nothing was saved."""
        self.dev.blockSignals(True)
        self.dev.clear()
        best, best_rank = -1, 99
        if sd is not None:
            try:
                apis = sd.query_hostapis()
                for i, d in enumerate(sd.query_devices()):
                    if d["max_input_channels"] < 1:
                        continue
                    api = apis[d["hostapi"]]["name"]
                    self.dev.addItem(f"{d['name']}  [{api}]", i)
                    nm = d["name"].lower()
                    if any(k in nm for k in ("umik", "umm-6", "umm6", "dayton", "emm-6", "xref", "measurement")):
                        rank = 0                     # USB measurement mic
                    elif any(k in nm for k in ("scarlett", "focusrite")):
                        rank = 1 if "WASAPI" in api else 2 if ("MME" in api or "Core Audio" in api) else 3
                    else:
                        rank = 99
                    if rank < best_rank and rank < 99:
                        best, best_rank = self.dev.count() - 1, rank
            except Exception as e:
                self.status.showMessage(f"Audio device scan failed: {e}")
        self.dev.addItem(SIM_NAME, "sim")
        self.dev.addItem(SIM_CAL_NAME, "simcal")
        saved = self.dev.findText(self.cfg["device"]) if self.cfg["device"] else -1
        self.dev.setCurrentIndex(saved if saved >= 0 else max(best, 0))
        self.dev.blockSignals(False)
        self.fill_channels()
        self.auto_profile()

    # ---------- run / stop
    def running(self):
        return self.engine is not None and (self.stream is not None or self.sim is not None)

    def toggle_run(self, on):
        if on:
            try:
                self.start()
                self.run_btn.setText("■  Stop")
            except Exception as e:
                self.run_btn.setChecked(False)
                QtWidgets.QMessageBox.critical(self, APP_NAME, f"Could not open input:\n{e}")
        else:
            self.stop()
            self.run_btn.setText("▶  Start")
            self.thumb.setState("off")
            if self.mini:
                self.mini.update_view(None, T["fg"], "off", "stopped")
            if self.web:
                self.web.state = {"value": None, "color": T["fg"], "state": "off", "sub": "stopped"}
            if self.menubar_meter:
                self.menubar_meter.set(None, None)
            self.status.showMessage("Stopped")

    def start(self):
        self.stop()
        dev = self.dev.currentData()
        ch = self.ch.currentIndex() + 1
        self.cfg.update(device=self.dev.currentText(), channel=ch)
        self.save_cfg()
        if dev in ("sim", "simcal"):
            fs = 48000
            self.sim = SimSource(fs, "cal" if dev == "simcal" else "pink")
            if self.theater.player.playing and self.theater.out.currentData() is None:
                self.sim.inject = self.theater.player
        else:
            info = sd.query_devices(dev)
            if ch > info["max_input_channels"]:
                raise RuntimeError(f"That device only has {info['max_input_channels']} input channel(s).")
            fs = int(info["default_samplerate"])
            idx = ch - 1
            self.q = queue.Queue()

            def cb(indata, frames, t, status):
                if status:
                    self.xruns += 1
                self.q.put(indata[:, idx].copy())

            api = sd.query_hostapis(info["hostapi"])["name"]
            if "ASIO" in api:
                # open just the one channel we need from a big ASIO device (e.g. SoundGrid 64 ch)
                idx = 0
                self.stream = sd.InputStream(device=dev, channels=1, samplerate=fs, blocksize=1024,
                                             dtype="float32", callback=cb,
                                             extra_settings=sd.AsioSettings(channel_selectors=[ch - 1]))
            else:
                self.stream = sd.InputStream(device=dev, channels=ch, samplerate=fs, blocksize=1024,
                                             dtype="float32", callback=cb)
            self.stream.start()
        prof = self.active_profile()
        pts = prof.get("cal_points") if prof.get("use_curve") else None
        old = self.engine
        self.engine = Engine(fs, prof["offset"], self.cfg["long_leq_min"], pts)
        if old is not None and old.fs == fs:      # keep the exposure going across stop / start
            self.engine.dose = old.dose
        self.session.cal_info = prof.get("cal_info", "")
        self.pending = []
        self.status.showMessage(f"Running @ {fs} Hz")

    def stop(self):
        if self.stream is not None:
            try:
                self.stream.stop(); self.stream.close()
            except Exception:
                pass
        self.stream = self.sim = None

    # ---------- actions
    def reset_stats(self):
        if self.engine:
            self.engine.reset_stats()
        self.meter.history.clear()
        self.add_event("Meter stats reset")

    def calibrate(self):
        CalWizard(self).exec()

    def settings(self):
        d = SettingsDialog(self.cfg, self)
        if d.exec():
            d.apply(self.cfg)
            if self.engine:
                self.engine.long_min = self.cfg["long_leq_min"]
            self.apply_theme()
            self.on_theme_change()
            self.update_menubar_meter()
            self.save_cfg()

    def update_menubar_meter(self):
        if self.cfg["menubar"] and self.menubar_meter is None:
            self.menubar_meter = MenuBarMeter(self)
        elif not self.cfg["menubar"] and self.menubar_meter is not None:
            self.menubar_meter.remove()
            self.menubar_meter = None

    def add_event(self, text, kind="event"):
        t = time.time()
        if kind == "marker":
            self.session.markers.append((t, text))
        else:
            self.session.events.append((t, text))
        self.pending_notes.append(("◆ " if kind == "marker" else "⚠ " if kind == "alarm" else "") + text)
        self.session_panel.add_event_row(t, text, kind)

    def add_marker(self, note=None):
        n = len(self.session.markers) + 1
        self.add_event(note or f"Marker {n}", "marker")
        self.status.showMessage(f"Marker {n} dropped", 3000)

    def add_marker_note(self):
        txt, ok = QtWidgets.QInputDialog.getText(self, "Marker", "Note for this marker:")
        if ok:
            self.add_marker(txt.strip() or None)

    def new_session(self, confirm=True):
        if confirm and len(self.session) > 60:
            r = QtWidgets.QMessageBox.question(self, APP_NAME, "Start a new session? The current history is cleared "
                                               "(export it first if you need it).")
            if r != QtWidgets.QMessageBox.Yes:
                return
        self.session = Session()
        self.session.cal_info = self.active_profile().get("cal_info", "")
        self.session_panel.events.clear()
        self.alarm_armed = {}
        self.pco_panel.flagged = set()
        self.pco_panel.table.setRowCount(0)
        if self.current_item:
            self.session.items.append((time.time(), self.current_item))
        self.session_panel.graph.update()

    def about(self):
        QtWidgets.QMessageBox.about(self, f"About {APP_NAME}",
            f"<h2>{APP_NAME} {APP_VERSION}</h2>"
            "<p>Sound level meter for services, shows and theaters: Leq, RTA, dosimetry, theater calibration, "
            "logging and reports, Planning Center, Companion and NDI.</p>"
            "<p>A/C weighting and Fast/Slow per IEC 61672-1 · 1/3 octaves per IEC 61260.</p>"
            "<p>Built with Python, Qt (PySide6), NumPy, SciPy and PortAudio. NDI® is a registered trademark of Vizrt NDI AB.</p>")

    def guide(self, page="Welcome"):
        GuideWizard(self, page).exec()

    def edit_schedule(self):
        d = ScheduleDialog(self.cfg, self)
        if d.exec():
            d.apply()
            self.save_cfg()

    def set_item(self, title):
        """Start a new service item (from Planning Center Live or by hand)."""
        self.current_item = title
        self.session.items.append((time.time(), title))
        self.add_event(f"Item: {title}")

    def out_dir(self):
        """Where logs and reports go: the dated auto-save (OneDrive) folder, or the plain log folder."""
        base = self.cfg["autosave_dir"]
        if not base:
            return self.cfg["log_dir"]
        day = datetime.datetime.fromtimestamp(self.session.start)
        name = day.strftime("%Y-%m-%d")
        if self.pco_plan_label:
            name += " - " + re.sub(r'[\\/:*?"<>|]+', "-", self.pco_plan_label)[:60].strip()
        return os.path.join(base, "ServiceSPL", day.strftime("%Y"), name)

    def export_csv(self, path=None):
        if not len(self.session):
            if path is None:
                QtWidgets.QMessageBox.information(self, APP_NAME, "Nothing recorded yet: press Start first.")
            return
        os.makedirs(self.out_dir(), exist_ok=True)
        if path is None:
            default = os.path.join(self.out_dir(), datetime.datetime.fromtimestamp(self.session.start)
                                   .strftime("Session_%Y-%m-%d_%H%M.csv"))
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export session CSV", default, "CSV (*.csv)")
            if not path:
                return
        s = self.session
        notes = collections.defaultdict(list)
        for t, txt in s.markers:
            notes[int(t)].append("◆ " + txt)
        for t, txt in s.events:
            notes[int(t)].append(txt)
        keys = sorted(notes)
        it_t = [t for t, _ in s.items]
        with open(path, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f)
            w.writerow(["Timestamp", "LAeq 1s", "LCeq 1s", "LZeq 1s", "LAFmax", "LCpeak", "LAS", "LAeq 1 min",
                        f"LAeq {self.cfg['long_leq_min']} min", f"LCeq {self.cfg['long_leq_min']} min",
                        "LA10 1s", "LA90 1s", "Dose OSHA HC %", "Dose NIOSH %", "Item", "Markers / events",
                        "Calibration"])
            cal = s.cal_info or "UNCALIBRATED"
            laf = np.frombuffer(s.laf, dtype=np.float32)
            pos, ki = 0, 0
            for i in range(len(s)):
                row = {k: s.cols[k][i] for k in s.COLS}
                seg = laf[pos:pos + 10]; pos += 10
                tt = int(row["time"])
                txt = []
                while ki < len(keys) and keys[ki] <= tt:
                    txt += notes[keys[ki]]
                    ki += 1
                j = int(np.searchsorted(it_t, row["time"], side="right")) - 1
                item = s.items[j][1] if j >= 0 else ""
                w.writerow([datetime.datetime.fromtimestamp(row["time"]).isoformat(sep=" ", timespec="seconds")]
                           + [f"{row[k]:.1f}" for k in ("LAeq1s", "LCeq1s", "LZeq1s", "LAFmax1s", "LCpeak1s",
                                                        "LAS", "LAeq1m", "LAeqLong", "LCeqLong")]
                           + [f"{ln_of(seg, 10):.1f}", f"{ln_of(seg, 90):.1f}",
                              f"{row['doseOSHA']:.2f}", f"{row['doseNIOSH']:.2f}", item, " | ".join(txt), cal])
        self.status.showMessage(f"Exported {len(s)} rows to {path}")
        return path

    def export_pdf(self, path=None):
        if not len(self.session):
            if path is None:
                QtWidgets.QMessageBox.information(self, APP_NAME, "Nothing recorded yet: press Start first.")
            return
        os.makedirs(self.out_dir(), exist_ok=True)
        interactive = path is None
        if interactive:
            default = os.path.join(self.out_dir(), datetime.datetime.fromtimestamp(self.session.start)
                                   .strftime("SPL_Report_%Y-%m-%d_%H%M.pdf"))
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export PDF report", default, "PDF (*.pdf)")
            if not path:
                return
        try:
            write_pdf_report(path, self.session, self.cfg, self.cfg["profile"], self.active_profile(),
                             self.pco_plan_label)
            self.status.showMessage(f"Report saved to {path}")
            if interactive:
                QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(path))
        except Exception as e:
            if interactive:
                QtWidgets.QMessageBox.critical(self, APP_NAME, f"Could not write the report:\n{e}")
            else:
                self.add_event(f"Auto-save of the report failed: {e}")
        return path

    def autosave(self):
        """Report + session CSV into the auto-save (OneDrive) folder."""
        stamp = datetime.datetime.now().strftime("%H%M")
        d = self.out_dir()
        try:
            os.makedirs(d, exist_ok=True)
            self.export_pdf(os.path.join(d, f"SPL_Report_{stamp}.pdf"))
            self.export_csv(os.path.join(d, f"Session_{stamp}.csv"))
            self.status.showMessage(f"Report and session CSV saved to {d}")
        except Exception as e:
            self.add_event(f"Auto-save failed: {e}")

    def toggle_log(self, on):
        if on:
            try:
                d = self.out_dir()
                os.makedirs(d, exist_ok=True)
                path = os.path.join(d, datetime.datetime.now().strftime("SPL_%Y-%m-%d_%H%M%S.csv"))
                self.log_file = open(path, "w", newline="", encoding="utf-8-sig")
                self.log_writer = csv.writer(self.log_file)
                self.log_writer.writerow(["Timestamp", "LAeq", "LCeq", "LAFmax", "LCpeak", "LAS",
                                          f"LAeq {self.cfg['long_leq_min']} min", "LA10", "LA90",
                                          "Dose OSHA HC %", "Dose NIOSH %", "Item", "Markers / events",
                                          "Calibration"])
                self.log_rows = 0
                self.pending = []
                self.pending_notes = []
                self.log_path = path
                self.log_btn.setText("●  Logging")
            except Exception as e:
                self.log_btn.setChecked(False)
                QtWidgets.QMessageBox.critical(self, APP_NAME, f"Could not start log:\n{e}")
        else:
            if self.log_file:
                self.log_file.close()
                self.status.showMessage(f"Saved {self.log_rows} rows to {self.log_path}")
            was = self.log_writer is not None
            self.log_file = self.log_writer = None
            self.log_btn.setText("●  Log")
            if self.sched_key and not self._sched_stopping:      # stopped by hand inside a schedule window
                self.sched_skip.add(self.sched_key)
                self.sched_key = None
            if was and self.cfg["autosave_dir"] and self.cfg["autosave_reports"]:
                self.autosave()

    def check_schedule(self):
        now = datetime.datetime.now()
        active, entry = None, None
        for i, e in enumerate(self.cfg["schedule"]):
            if schedule_window(e, now):
                active, entry = (i, now.strftime("%Y-%m-%d")), e
                break
        if active and active != self.sched_key and active not in self.sched_skip and not self.log_writer:
            if not self.running():
                self.run_btn.click()
                if not self.running():
                    return
            if entry.get("report", True):
                self.new_session(confirm=False)
            self.sched_key = active
            self.add_event(f"Scheduled logging started ({entry['start']}–{entry['stop']})")
            self.log_btn.click()
        elif self.sched_key and active != self.sched_key:
            self.add_event("Scheduled logging stopped")
            self._sched_stopping = True
            if self.log_writer:
                self.log_btn.click()
            self._sched_stopping = False
            self.sched_key = None

    def write_log(self):
        """Every N seconds write one row: Leq / max / peak / Ln over that interval, plus markers."""
        n = self.cfg["log_interval"]
        while len(self.pending) >= n:
            chunk, self.pending = self.pending[:n], self.pending[n:]
            en = lambda key: 10 * math.log10(sum(10 ** (r[key] / 10) for r in chunk) / n)
            last = chunk[-1]
            laf = np.concatenate([r["LAF10"] for r in chunk])
            ts = datetime.datetime.fromtimestamp(last["time"]).isoformat(sep=" ", timespec="seconds")
            prof = self.active_profile()
            cal = f"{self.cfg['profile']}: {prof['cal_info']}" if prof["calibrated"] else "UNCALIBRATED"
            notes, self.pending_notes = " | ".join(self.pending_notes), []
            self.log_writer.writerow([ts, f"{en('LAeq1s'):.1f}", f"{en('LCeq1s'):.1f}",
                                      f"{max(r['LAFmax1s'] for r in chunk):.1f}",
                                      f"{max(r['LCpeak1s'] for r in chunk):.1f}",
                                      f"{last['LAS']:.1f}", f"{last['LAeqLong']:.1f}",
                                      f"{ln_of(laf, 10):.1f}", f"{ln_of(laf, 90):.1f}",
                                      f"{last['doseOSHA']:.2f}", f"{last['doseNIOSH']:.2f}",
                                      self.current_item, notes, cal])
            self.log_file.flush()
            self.log_rows += 1

    # ---------- NDI
    def toggle_ndi(self, on):
        if on:
            try:
                self.ndi = NdiSender(self.cfg["ndi_name"], self.cfg["ndi_lib"])
            except Exception as e:
                self.ndi = None
                self.ndi_btn.setChecked(False)
                QtWidgets.QMessageBox.warning(self, APP_NAME, str(e))
                return
            self.ndi_btn.setText("NDI ✓")
            self.ndi_timer.start(33)
            self.status.showMessage(f"NDI source '{self.cfg['ndi_name']}' is live")
        else:
            self.ndi_timer.stop()
            if self.ndi:
                self.ndi.close()
            self.ndi = None
            self.ndi_btn.setText("NDI")

    def render_card(self):
        w, h = (int(x) for x in self.cfg["ndi_size"].split("x"))
        img = QtGui.QImage(w, h, QtGui.QImage.Format_ARGB32)
        img.fill(QtCore.Qt.transparent)
        p = QtGui.QPainter(img)
        paint_meter_card(p, QtCore.QRectF(0, 0, w, h), self.card, self.cfg["ndi_transparent"])
        p.end()
        return img

    def send_ndi(self):
        if self.ndi is None:
            return
        try:
            self.ndi.send(self.render_card(), 30)
        except Exception as e:
            self.ndi_btn.setChecked(False)
            self.toggle_ndi(False)
            self.status.showMessage(f"NDI stopped: {e}")

    # ---------- alarms
    def check_alarms(self, rec):
        e = self.engine
        vals = {"LAS (live)": rec["LAS"], "LAeq 1 min": rec["LAeq1m"], "LAeq long": rec["LAeqLong"],
                "LCeq long": rec["LCeqLong"], "LCpeak": rec["LCpeak1s"], "LAFmax": rec["LAFmax1s"],
                "OSHA HC dose %": rec["doseOSHA"], "NIOSH dose %": rec["doseNIOSH"]}
        for i, r in enumerate(self.cfg["alarm_rules"]):
            if not r["on"]:
                continue
            v = vals.get(r["metric"], float("nan"))
            armed = self.alarm_armed.get(i, True)
            if armed and np.isfinite(v) and v >= r["thr"]:
                self.alarm_armed[i] = False
                name = r["metric"].replace("long", f"{self.cfg['long_leq_min']} min")
                self.raise_alarm(f"{name} {v:.1f} ≥ {r['thr']:.1f}")
            elif not armed and np.isfinite(v) and v < r["thr"] - 1.0:
                self.alarm_armed[i] = True
        for d in e.dose[1:]:
            for lim in (50, 100):
                key = (d.name, lim)
                if d.dose * 100 >= lim and key not in self.dose_flags:
                    self.dose_flags.add(key)
                    self.raise_alarm(f"{d.name} noise dose reached {lim} %"
                                     + (" (action level)" if lim == 50 else " (limit)"))
            proj = d.projected()
            key = (d.name, "proj")
            if np.isfinite(proj) and proj >= 100 and d.elapsed > 600 and key not in self.dose_flags:
                self.dose_flags.add(key)
                self.add_event(f"{d.name}: on pace for {proj:.0f} % over 8 h")

    def raise_alarm(self, text):
        self.add_event(text, "alarm")
        self.meter.flash(text)
        if self.cfg["alarm_beep"]:
            QtWidgets.QApplication.beep()

    # ---------- web / mini
    def toggle_web(self, on):
        if on:
            try:
                self.web = WebView(self.cfg["web_port"])
            except OSError as e:
                self.web_btn.setChecked(False)
                QtWidgets.QMessageBox.critical(self, APP_NAME,
                    f"Could not start web view on port {self.cfg['web_port']}:\n{e}\n\nTry another port in Settings.")
                return
            urls = self.web.urls()
            self.web_btn.setText("Web view ✓")
            self.web_btn.setToolTip("Open in a browser: " + "  or  ".join(urls))
            box = QtWidgets.QMessageBox(self)
            box.setWindowTitle("Web view is on")
            box.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
            box.setText("On the other computer (e.g. the Proton), open a browser to:\n\n" + "\n".join(urls) +
                        "\n\nIf there's more than one address, use the one on the same network as that computer.")
            box.exec()
        else:
            if self.web:
                self.web.stop()
            self.web = None
            self.web_btn.setText("Web view")

    def enter_mini(self):
        if self.mini is None:
            self.mini = MiniMeter(self)
        pos = self.cfg.get("mini_pos")
        if pos:
            self.mini.move(*pos)
        else:
            g = self.geometry()
            self.mini.move(g.right() - 300, g.top() + 10)
        self.mini.thumb.setState(self.thumb.state)
        self.mini.show()
        self.hide()

    def exit_mini(self):
        if self.mini:
            self.mini.hide()
        self.show(); self.raise_(); self.activateWindow()

    def set_red_at(self, v):
        self.cfg["alarm"] = v
        self.apply_cfg_to_ui()
        self.meter.bar.update(); self.meter.history.update()
        self.save_cfg()

    def color_for(self, v):
        if v is None or not np.isfinite(v):
            return T["fg"]
        return (T["red"] if v >= self.cfg["alarm"] else T["yellow"]
                if v >= self.cfg["alarm"] - self.cfg["warn_gap"] else T["green"])

    # ---------- main loop
    def tick(self):
        e = self.engine
        self.meter.tick()
        now_s = int(time.time())
        if now_s != self._last_sec:
            self._last_sec = now_s
            self.check_schedule()
            self.pco_panel.tick()
        e = self.engine
        if not self.running():
            return
        blocks = self.sim.read() if self.sim else []
        if self.stream is not None:
            for _ in range(200):
                try:
                    blocks.append(self.q.get_nowait())
                except queue.Empty:
                    break
        for b in blocks:
            e.process(b)

        new, e.records = e.records, []
        for r in new:
            self.session.add(r)
            self.check_alarms(r)
        if new:
            self.last_rec = new[-1]
        if self.log_writer and new:
            self.pending.extend(new)
            self.write_log()

        m = self.meter
        main = e.level(self.cfg["main_weight"], self.cfg["main_time"])
        if new:
            self.last_1m = new[-1]["LAeq1m"]
        judged = main if self.cfg["thumb_src"] == "live" else self.last_1m
        now = time.time()
        if np.isfinite(judged) and judged >= self.cfg["alarm"]:
            self.bad_until = now + 3.0          # hold thumbs-down 3 s so it doesn't flicker
        self.thumb.setState("bad" if now < self.bad_until else "good")
        col = T["red"] if now < self.bad_until else self.color_for(main)
        m.main_tile.set(main, col)
        w, t = self.cfg["main_weight"], self.cfg["main_time"]
        sub = (f"L{w}{t}  ·  1 min {fmt(self.last_1m).strip()}  ·  red at {self.cfg['alarm']:.0f}"
               + ("  ·  ● REC" if self.log_writer else ""))
        if self.mini is not None and self.mini.isVisible():
            self.mini.update_view(main, col, self.thumb.state, sub)
        if self.web:
            self.web.state = {"value": round(main, 1) if np.isfinite(main) else None,
                              "color": col, "state": self.thumb.state, "sub": sub}
        if self.menubar_meter:
            self.menubar_meter.set(main, col)
        lr = self.last_rec
        lvl = fmt(main).strip()
        state = "bad" if now < self.bad_until else ("warn" if col == T["yellow"] else "good")
        self.card = {"level": lvl, "color": col if col != T["fg"] else "#e8e8e8", "state": self.thumb.state,
                     "sub": f"L{w}{t}  ·  LAeq 1 min {fmt(self.last_1m).strip()}"
                            + (f"  ·  {self.current_item}" if self.current_item else "")}
        if self.cfg["companion_on"] or self.web:
            api = {"level": lvl, "leq1m": fmt(lr.get("LAeq1m")).strip(), "leqlong": fmt(lr.get("LAeqLong")).strip(),
                   "lceqlong": fmt(lr.get("LCeqLong")).strip(), "peak": fmt(e.L(e.max_pk)).strip(),
                   "dose": f"{e.dose[2].dose * 100:.0f}%", "state": state,
                   "color": QtGui.QColor(PALETTES["dark"][{"bad": "red", "warn": "yellow", "good": "green"}[state]]).name(),
                   "label": f"L{w}{t}", "item": self.current_item, "red_at": self.cfg["alarm"],
                   "recording": bool(self.log_writer), "profile": self.cfg["profile"]}
            if self.web:
                self.web.api = api
            if self.cfg["companion_on"]:
                self.companion.send(self.cfg, api)
        m.bar.setLevel(main)
        m.history.add(main, self.last_1m)
        self.rta_panel.tick(e)
        self.theater.tick(e)
        if new:
            r = new[-1]
            n_s, lw = len(e.hist_n), self.cfg["long_leq_min"]
            filling = "" if n_s >= lw * 60 else f"filling {n_s / 60:.1f}/{lw} min"
            m.t_1m.set(r["LAeq1m"], self.color_for(r["LAeq1m"]))
            m.t_long.set(r["LAeqLong"], self.color_for(r["LAeqLong"]), sub=filling)
            m.t_clong.set(r["LCeqLong"], sub=f"LC − LA {fmt(r['LCeqLong'] - r['LAeqLong']).strip()}")
            m.t_tot.set(e.leq_total("A"), sub=hms(time.time() - e.reset_time))
            m.t_max.set(e.L(e.max_af))
            pk = e.L(e.max_pk)
            m.t_pk.set(pk, T["red"] if pk >= 135 else None)
            m.t_l10.set(r["L10"], sub=f"LA50 {fmt(r['L50']).strip()}")
            m.t_l90.set(r["L90"], sub=f"since reset {fmt(e.stat.ln(90)).strip()}")
            if self.docks["dose"].isVisible():
                self.dose_panel.tick(e)
            self.session_panel.tick()
            msg = f"Running @ {e.fs} Hz"
            if self.log_writer:
                msg += f"   ·   Logging every {self.cfg['log_interval']} s → {self.log_rows} rows"
            if self.xruns:
                msg += f"   ·   dropouts: {self.xruns}"
            if not self.status.currentMessage().startswith(("Marker", "Exported", "Report", "Mic profile")):
                self.status.showMessage(msg)

    def closeEvent(self, ev):
        self.stop()
        self.theater.player.stop()
        self.toggle_ndi(False)
        if self.web:
            self.web.stop()
        if self.mini:
            self.mini.deleteLater()
        if self.menubar_meter:
            self.menubar_meter.remove()
        if self.log_file:
            self.log_file.close()
        self.cfg["win_geom"] = base64.b64encode(bytes(self.saveGeometry())).decode()
        self.cfg["dock_state"] = base64.b64encode(bytes(self.saveState())).decode()
        self.save_cfg()
        super().closeEvent(ev)
        QtWidgets.QApplication.quit()


def main():
    QtWidgets.QApplication.setHighDpiScaleFactorRoundingPolicy(
        QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    app.setApplicationVersion(APP_VERSION)
    app.setQuitOnLastWindowClosed(False)    # mini meter / menu bar keep running while the window is hidden
    w = MainWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
