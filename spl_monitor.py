"""
SPL Monitor - simple SPL meter with logging (Windows).
Built for a Focusrite Scarlett + calibrated measurement mic.

Requires: numpy, scipy, PySide6, sounddevice
"""

import sys, os, json, time, csv, queue, math, datetime, socket, threading
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
import numpy as np
from scipy import signal
from PySide6 import QtCore, QtGui, QtWidgets, QtSvg

os.environ.setdefault("SD_ENABLE_ASIO", "1")  # show ASIO drivers (Waves SoundGrid, Focusrite ASIO)
try:
    import sounddevice as sd
except Exception:  # PortAudio missing -> simulator only
    sd = None

APP_NAME = "SPL Monitor"
if sys.platform == "darwin":
    CFG_DIR = os.path.expanduser("~/Library/Application Support/SPLMonitor")
else:
    CFG_DIR = os.path.join(os.environ.get("APPDATA", os.path.expanduser("~")), "SPLMonitor")
CFG_PATH = os.path.join(CFG_DIR, "config.json")
SIM_NAME = "Simulated pink noise (test)"

DEFAULTS = dict(
    device=None, channel=1, cal_offset=120.0, calibrated=False, cal_info="",
    warn_gap=5.0, alarm=100.0, main_weight="A", main_time="S",
    long_leq_min=15, log_interval=1, thumb_src="live", web_port=8080, mini_on_top=True, mini_pos=None,
    log_dir=os.path.join(os.path.expanduser("~"), "Documents", "SPL Logs"),
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



class Engine:
    """All metering maths. Feed it mono float blocks via process()."""
    HIST_MAX = 7200          # 2 h of per-second data

    def __init__(self, fs, offset, long_min):
        self.fs = fs
        self.offset = offset
        self.long_min = long_min
        self.sos = {k: weighting_sos(k, fs) for k in "ACZ"}
        self.zi = {k: (np.zeros((s.shape[0], 2)) if s is not None else None)
                   for k, s in self.sos.items()}
        self.alpha = {"F": 1 - math.exp(-1 / (0.125 * fs)), "S": 1 - math.exp(-1 / (1.0 * fs))}
        self.tz = {(k, t): np.zeros(1) for k in "ACZ" for t in "FS"}
        self.cur = {(k, t): 1e-20 for k in "ACZ" for t in "FS"}
        self.hist_sum = {k: [] for k in "ACZ"}
        self.hist_n = []
        self.records = []
        self.cal_target = 0
        self.cal_sum = 0.0
        self.cal_n = 0
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
        self.reset_time = time.time()

    def _new_second(self):
        self.sec_sum = {k: 0.0 for k in "ACZ"}
        self.sec_n = 0
        self.sec_pk = 1e-20
        self.sec_maxaf = 1e-20


    # --- main processing
    def process(self, x):
        x = np.asarray(x, dtype=np.float64)
        n = len(x)
        if n == 0:
            return
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
            s = float(sq.sum())
            self.sec_sum[k] += s
            self.tot_sum[k] += s
            if k == "C":
                self.sec_pk = max(self.sec_pk, float(sq.max()))
            if k == "A" and self.cal_target:
                self.cal_sum += s
                self.cal_n += n
        self.sec_n += n
        self.tot_n += n
        self.max_pk = max(self.max_pk, self.sec_pk)
        if self.sec_n >= self.fs:
            self._close_second()

    def _close_second(self):
        for k in "ACZ":
            self.hist_sum[k].append(self.sec_sum[k])
            del self.hist_sum[k][:-self.HIST_MAX]
        self.hist_n.append(self.sec_n)
        del self.hist_n[:-self.HIST_MAX]
        L, c, sn = self.L, self.cur, self.sec_n
        self.records.append({
            "time": time.time(),
            "LAS": L(c["A", "S"]), "LAF": L(c["A", "F"]),
            "LCS": L(c["C", "S"]), "LZS": L(c["Z", "S"]),
            "LAeq1s": L(self.sec_sum["A"] / sn), "LCeq1s": L(self.sec_sum["C"] / sn),
            "LZeq1s": L(self.sec_sum["Z"] / sn),
            "LAeq1m": self.leq("A", 60), "LCeq1m": self.leq("C", 60),
            "LAeqLong": self.leq("A", int(self.long_min * 60)),
            "LCpeak1s": L(self.sec_pk), "LAFmax1s": L(self.sec_maxaf),
        })
        self._new_second()

    def leq(self, k, secs):
        if not self.hist_n:
            return float("nan")
        return self.L(sum(self.hist_sum[k][-secs:]) / sum(self.hist_n[-secs:]))

    def level(self, k, t):
        return self.L(self.cur[(k, t)])

    def leq_total(self, k):
        return self.L(self.tot_sum[k] / self.tot_n) if self.tot_n else float("nan")


    # --- calibration
    def begin_cal(self, seconds=3.0):
        self.cal_sum, self.cal_n = 0.0, 0
        self.cal_target = int(seconds * self.fs)

    @property
    def cal_done(self):
        return self.cal_target and self.cal_n >= self.cal_target

    def finish_cal(self, ref_db):
        ms = self.cal_sum / max(self.cal_n, 1)
        self.cal_target = 0
        self.offset = ref_db - 10 * math.log10(max(ms, 1e-20))
        return self.offset


class SimSource:
    """Pink noise with a slow 20-second swell, paced by the wall clock."""
    B = [0.049922035, -0.095993537, 0.050612699, -0.004408786]
    A = [1, -2.494956002, 2.017265875, -0.522189400]

    def __init__(self, fs):
        self.fs = fs
        self.t0 = time.monotonic()
        self.done = 0
        self.zi = np.zeros(3)

    def read(self):
        target = int((time.monotonic() - self.t0) * self.fs)
        n = min(target - self.done, self.fs)
        if n <= 0:
            return []
        t = (self.done + np.arange(n)) / self.fs
        self.done = target
        p, self.zi = signal.lfilter(self.B, self.A, np.random.randn(n), zi=self.zi)
        env = 0.05 * (1 + 0.8 * np.sin(2 * np.pi * t / 20))
        return [p * env]


# ---------------------------------------------------------------- UI -----

BG, PANEL, FG, DIM = "#121417", "#1b1e23", "#e8e8e8", "#8a9099"
GREEN, YELLOW, RED = "#3ddc84", "#ffcc33", "#ff4d4d"

STYLE = f"""
QWidget {{ background:{BG}; color:{FG}; font-family:'Segoe UI','Helvetica Neue',Arial; font-size:10pt; }}
QFrame#tile {{ background:{PANEL}; border:1px solid #2a2e35; border-radius:6px; }}
QPushButton {{ background:#2a2e35; border:1px solid #3a3f47; border-radius:4px; padding:5px 12px; }}
QPushButton:hover {{ background:#353a42; }}
QPushButton:checked {{ background:#1f6f43; border-color:{GREEN}; }}
QPushButton#log:checked {{ background:#7a1f1f; border-color:{RED}; }}
QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit {{ background:#22262c; border:1px solid #3a3f47;
    border-radius:4px; padding:3px 6px; }}
QStatusBar {{ color:{DIM}; }}
"""
LOG_INTERVALS = {"1 second": 1, "10 seconds": 10, "1 minute": 60}


def fmt(v):
    return "--.-" if v is None or not np.isfinite(v) else f"{v:5.1f}"


class Tile(QtWidgets.QFrame):
    def __init__(self, title, size=24):
        super().__init__()
        self.setObjectName("tile")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 6, 10, 6)
        lay.setSpacing(0)
        self.t = QtWidgets.QLabel(title)
        self.t.setStyleSheet(f"color:{DIM}; font-size:9pt; background:transparent;")
        self.v = QtWidgets.QLabel()
        self.v.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.s = QtWidgets.QLabel("")
        self.s.setStyleSheet(f"color:{DIM}; font-size:8pt; background:transparent;")
        self.s.setAlignment(QtCore.Qt.AlignRight)
        self.size = size
        for w in (self.t, self.v, self.s):
            lay.addWidget(w)
        self.set(None)

    def set(self, val, color=FG, sub=None, title=None):
        self.v.setText(fmt(val))
        self.v.setStyleSheet(f"color:{color}; font-size:{self.size}pt; font-weight:600;"
                             f"font-family:Consolas,Menlo,'Courier New',monospace; background:transparent;")
        if sub is not None:
            self.s.setText(sub)
        if title:
            self.t.setText(title)


class LevelBar(QtWidgets.QWidget):
    """Horizontal bar meter, 40-120 dB, with 2-second peak hold."""
    def __init__(self, lo=40, hi=120):
        super().__init__()
        self.lo, self.hi = lo, hi
        self.level = self.peak = lo
        self.peak_t = 0
        self.warn, self.alarm = 95, 100
        self.setMinimumHeight(46)

    def setLevel(self, v):
        v = v if np.isfinite(v) else self.lo
        self.level = v
        now = time.time()
        if v >= self.peak or now - self.peak_t > 2:
            self.peak, self.peak_t = v, now
        self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        r = self.rect().adjusted(8, 4, -8, -18)
        p.fillRect(r, QtGui.QColor("#0c0d0f"))
        frac = lambda v: (min(max(v, self.lo), self.hi) - self.lo) / (self.hi - self.lo)
        x = lambda v: r.left() + frac(v) * r.width()
        g = QtGui.QLinearGradient(r.left(), 0, r.right(), 0)
        fw, fa = frac(self.warn), frac(self.alarm)
        for pos, col in ((0, GREEN), (max(fw - 1e-3, 0), GREEN), (fw, YELLOW),
                         (max(fa - 1e-3, fw), YELLOW), (fa, RED), (1, RED)):
            g.setColorAt(pos, QtGui.QColor(col))
        p.fillRect(QtCore.QRectF(r.left(), r.top(), x(self.level) - r.left(), r.height()), g)
        p.setPen(QtGui.QPen(QtGui.QColor("white"), 2))
        p.drawLine(QtCore.QPointF(x(self.peak), r.top()), QtCore.QPointF(x(self.peak), r.bottom()))
        p.setPen(QtGui.QColor(DIM))
        f = p.font(); f.setPointSize(8); p.setFont(f)
        for v in range(self.lo, self.hi + 1, 10):
            xx = x(v)
            p.drawLine(QtCore.QPointF(xx, r.bottom()), QtCore.QPointF(xx, r.bottom() + 4))
            p.drawText(QtCore.QRectF(xx - 20, r.bottom() + 3, 40, 14), QtCore.Qt.AlignCenter, str(v))


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
                                    ("off", "#3a3f47", ""))}
        self.state = "off"

    def setState(self, st):
        if st != self.state:
            self.state = st
            self.update()

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        self.r[self.state].render(p, QtCore.QRectF(self.rect()))


class MiniMeter(QtWidgets.QWidget):
    """Small frameless box: thumb + live reading. Drag to move, pin to keep on top, ⤢ to go back."""
    def __init__(self, main):
        super().__init__(None, QtCore.Qt.Window | QtCore.Qt.FramelessWindowHint)
        self.main = main
        self.setWindowTitle(APP_NAME)
        self.setObjectName("mini")
        self.setStyleSheet(f"QWidget#mini {{ background:{PANEL}; border:1px solid #3a3f47; border-radius:8px; }}"
                           "QToolButton { background:transparent; border:none; color:#8a9099; font-size:11pt; }"
                           "QToolButton:checked { color:#ffffff; }"
                           "QToolButton:hover { color:#ffffff; }")
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
        self.sub.setStyleSheet(f"color:{DIM}; font-size:8pt; background:transparent;")
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
        self.set_on_top(self.pin.isChecked(), show=False)

    def set_on_top(self, on, show=True):
        self.main.cfg["mini_on_top"] = on
        self.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint, on)
        if show:
            self.show()
        self.main.save_cfg()

    def update_view(self, value, color, state, sub):
        self.val.setText(fmt(value))
        self.val.setStyleSheet(f"color:{color}; font-size:34pt; font-weight:600;"
                               f"font-family:Consolas,Menlo,'Courier New',monospace; background:transparent;")
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
 catch(e){document.getElementById("sub").textContent="no connection to SPL Monitor";
  document.getElementById("val").style.color="#8a9099";thumb("off");}
 setTimeout(poll,200);}
poll();
</script></body></html>"""


class WebView:
    """Tiny read-only HTTP server: '/' = big-number page, '/data' = live JSON."""
    def __init__(self, port):
        self.state = {"value": None, "color": FG, "state": "off", "sub": "stopped"}
        page = WEB_PAGE.replace("__SVG__", json.dumps(THUMB_SVG)).encode()
        view = self

        class H(BaseHTTPRequestHandler):
            def do_GET(self):
                if self.path.startswith("/data"):
                    body, ctype = json.dumps(view.state).encode(), "application/json"
                else:
                    body, ctype = page, "text/html; charset=utf-8"
                self.send_response(200)
                self.send_header("Content-Type", ctype)
                self.send_header("Cache-Control", "no-store")
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


class SettingsDialog(QtWidgets.QDialog):
    def __init__(self, cfg, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        f = QtWidgets.QFormLayout(self)
        self.w = QtWidgets.QComboBox(); self.w.addItems(["A", "C", "Z"]); self.w.setCurrentText(cfg["main_weight"])
        self.t = QtWidgets.QComboBox(); self.t.addItems(["Fast", "Slow"])
        self.t.setCurrentText("Fast" if cfg["main_time"] == "F" else "Slow")
        self.warn = QtWidgets.QDoubleSpinBox(); self.warn.setRange(0, 30); self.warn.setValue(cfg["warn_gap"])
        self.warn.setSuffix(" dB below red"); self.warn.setSpecialValueText("off")
        self.alarm = QtWidgets.QDoubleSpinBox(); self.alarm.setRange(40, 140); self.alarm.setValue(cfg["alarm"])
        self.long = QtWidgets.QSpinBox(); self.long.setRange(1, 120); self.long.setValue(cfg["long_leq_min"]); self.long.setSuffix(" min")
        self.li = QtWidgets.QComboBox(); self.li.addItems(list(LOG_INTERVALS))
        self.li.setCurrentText(next(k for k, v in LOG_INTERVALS.items() if v == cfg["log_interval"]))
        self.off = QtWidgets.QDoubleSpinBox(); self.off.setRange(0, 200); self.off.setDecimals(2); self.off.setValue(cfg["cal_offset"])
        self.logdir = QtWidgets.QLineEdit(cfg["log_dir"])
        browse = QtWidgets.QPushButton("Browse…"); browse.clicked.connect(self._browse)
        row = QtWidgets.QHBoxLayout(); row.addWidget(self.logdir); row.addWidget(browse)
        f.addRow("Meter weighting", self.w)
        f.addRow("Meter speed", self.t)
        f.addRow("Yellow warning starts", self.warn)
        f.addRow("Turn red / thumbs down at (dB)", self.alarm)
        self.src = QtWidgets.QComboBox(); self.src.addItems(["Live meter", "LAeq 1 min"])
        self.src.setCurrentIndex(0 if cfg["thumb_src"] == "live" else 1)
        f.addRow("Thumb / color follows", self.src)
        f.addRow("Long Leq window", self.long)
        f.addRow("Log interval", self.li)
        f.addRow("Log folder", row)
        self.port = QtWidgets.QSpinBox(); self.port.setRange(1024, 65535); self.port.setValue(cfg["web_port"])
        f.addRow("Web view port", self.port)
        f.addRow("Calibration offset (dB)", self.off)
        bb = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        f.addRow(bb)
        self._off0 = cfg["cal_offset"]

    def _browse(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Log folder", self.logdir.text())
        if d:
            self.logdir.setText(d)

    def apply(self, cfg):
        cfg.update(main_weight=self.w.currentText(), main_time=self.t.currentText()[0],
                   warn_gap=self.warn.value(), alarm=self.alarm.value(), long_leq_min=self.long.value(),
                   log_interval=LOG_INTERVALS[self.li.currentText()], log_dir=self.logdir.text(),
                   thumb_src="live" if self.src.currentIndex() == 0 else "leq1m",
                   web_port=self.port.value())
        if abs(self.off.value() - self._off0) > 1e-6:
            cfg.update(cal_offset=self.off.value(), calibrated=True,
                       cal_info=f"manual offset {self.off.value():.2f} dB")


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1180, 480)
        self.cfg = self.load_cfg()
        self.engine = self.stream = self.sim = None
        self.q = queue.Queue()
        self.xruns = 0
        self.log_file = self.log_writer = None
        self.log_rows = 0
        self.pending = []      # per-second records waiting to be logged
        self.cal_ref = None
        self._build_ui()
        self.populate_devices()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.tick)
        self.timer.start(50)

    # ---------- config
    def load_cfg(self):
        cfg = dict(DEFAULTS)
        try:
            with open(CFG_PATH) as f:
                cfg.update(json.load(f))
        except Exception:
            pass
        return cfg

    def save_cfg(self):
        try:
            os.makedirs(CFG_DIR, exist_ok=True)
            with open(CFG_PATH, "w") as f:
                json.dump(self.cfg, f, indent=2)
        except Exception:
            pass

    # ---------- UI
    def _build_ui(self):
        tb = QtWidgets.QToolBar(); tb.setMovable(False); self.addToolBar(tb)
        self.dev = QtWidgets.QComboBox(); self.dev.setMinimumWidth(240)
        self.ch = QtWidgets.QComboBox(); self.ch.setMinimumWidth(90)
        self.dev.currentIndexChanged.connect(self.fill_channels)
        self.run_btn = QtWidgets.QPushButton("▶  Start"); self.run_btn.setCheckable(True)
        self.run_btn.clicked.connect(self.toggle_run)
        cal = QtWidgets.QPushButton("Calibrate…"); cal.clicked.connect(self.calibrate)
        rst = QtWidgets.QPushButton("Reset"); rst.setToolTip("Reset Leq / Max / Peak")
        rst.clicked.connect(self.reset_stats)
        self.log_btn = QtWidgets.QPushButton("●  Log"); self.log_btn.setObjectName("log")
        self.log_btn.setCheckable(True); self.log_btn.clicked.connect(self.toggle_log)
        st = QtWidgets.QPushButton("Settings…"); st.clicked.connect(self.settings)
        for w in (self.dev, self.ch, self.run_btn, cal, rst, self.log_btn, st):
            tb.addWidget(w)
        tb.addSeparator()
        lbl = QtWidgets.QLabel("  Red at ")
        lbl.setStyleSheet(f"color:{RED}; font-weight:600;")
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
        self.mini = None
        self.web = None
        self.cal_lbl = QtWidgets.QLabel()
        self.status = self.statusBar()
        self.status.addPermanentWidget(self.cal_lbl)

        central = QtWidgets.QWidget(); self.setCentralWidget(central)
        v = QtWidgets.QVBoxLayout(central)
        self.main_tile = Tile("LAS", size=96)
        self.thumb = Thumb()
        ml = self.main_tile.layout()
        ml.removeWidget(self.main_tile.v)
        hr = QtWidgets.QHBoxLayout()
        hr.addWidget(self.thumb); hr.addStretch(); hr.addWidget(self.main_tile.v)
        ml.insertLayout(1, hr)
        self.bad_until = 0.0
        self.last_1m = float("nan")
        v.addWidget(self.main_tile, 1)
        self.bar = LevelBar(); v.addWidget(self.bar)
        row = QtWidgets.QHBoxLayout()
        self.t_1m = Tile("LAeq 1 min")
        self.t_long = Tile("LAeq 15 min")
        self.t_tot = Tile("LAeq since reset")
        self.t_max = Tile("LAFmax")
        self.t_pk = Tile("LCpeak")
        for t in (self.t_1m, self.t_long, self.t_tot, self.t_max, self.t_pk):
            row.addWidget(t)
        v.addLayout(row)
        self.apply_cfg_to_ui()

    def apply_cfg_to_ui(self):
        c = self.cfg
        self.main_tile.set(None, title=f"L{c['main_weight']}{c['main_time']}   "
                                       f"{c['main_weight']}-weighted · {'Fast' if c['main_time'] == 'F' else 'Slow'}")
        self.t_long.set(None, title=f"LAeq {c['long_leq_min']} min")
        self.bar.warn, self.bar.alarm = c["alarm"] - c["warn_gap"], c["alarm"]
        if hasattr(self, "red_at"):
            self.red_at.blockSignals(True); self.red_at.setValue(c["alarm"]); self.red_at.blockSignals(False)
        if c["calibrated"]:
            self.cal_lbl.setText(f"  CAL ✓ {c['cal_info']}  ")
            self.cal_lbl.setStyleSheet(f"color:{GREEN};")
        else:
            self.cal_lbl.setText("  NOT CALIBRATED  ")
            self.cal_lbl.setStyleSheet(f"color:{YELLOW}; font-weight:600;")

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
        """List inputs; auto-pick the Scarlett (WASAPI, then MME) if nothing was saved."""
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
        self.dev.addItem(SIM_NAME, None)
        saved = self.dev.findText(self.cfg["device"]) if self.cfg["device"] else -1
        self.dev.setCurrentIndex(saved if saved >= 0 else max(best, 0))
        self.fill_channels()

    # ---------- run / stop
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
                self.mini.update_view(None, FG, "off", "stopped")
            if self.web:
                self.web.state = {"value": None, "color": FG, "state": "off", "sub": "stopped"}
            self.status.showMessage("Stopped")

    def start(self):
        self.stop()
        dev = self.dev.currentData()
        ch = self.ch.currentIndex() + 1
        self.cfg.update(device=self.dev.currentText(), channel=ch)
        self.save_cfg()
        if dev is None:
            fs = 48000
            self.sim = SimSource(fs)
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
        self.engine = Engine(fs, self.cfg["cal_offset"], self.cfg["long_leq_min"])
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

    def calibrate(self):
        if not self.engine or (self.stream is None and self.sim is None):
            QtWidgets.QMessageBox.information(self, APP_NAME, "Press Start first, then calibrate.")
            return
        ref, ok = QtWidgets.QInputDialog.getDouble(
            self, "Calibrate",
            "Put the calibrator on the mic and switch it on.\n"
            "Set the Scarlett gain where you'll leave it — don't touch it after calibrating.\n\n"
            "Calibrator level (dB SPL):", 94.0, 60.0, 140.0, 1)
        if ok:
            self.cal_ref = ref
            self.engine.begin_cal(3.0)
            self.status.showMessage("Calibrating… hold steady for 3 seconds")

    def settings(self):
        d = SettingsDialog(self.cfg, self)
        if d.exec():
            d.apply(self.cfg)
            if self.engine:
                self.engine.offset = self.cfg["cal_offset"]
                self.engine.long_min = self.cfg["long_leq_min"]
            self.apply_cfg_to_ui()
            self.save_cfg()

    def toggle_log(self, on):
        if on:
            try:
                os.makedirs(self.cfg["log_dir"], exist_ok=True)
                path = os.path.join(self.cfg["log_dir"],
                                    datetime.datetime.now().strftime("SPL_%Y-%m-%d_%H%M%S.csv"))
                self.log_file = open(path, "w", newline="", encoding="utf-8-sig")
                self.log_writer = csv.writer(self.log_file)
                self.log_writer.writerow(["Timestamp", "LAeq", "LCeq", "LAFmax", "LCpeak", "LAS",
                                          f"LAeq {self.cfg['long_leq_min']} min", "Calibration"])
                self.log_rows = 0
                self.pending = []
                self.log_path = path
                self.log_btn.setText("●  Logging")
            except Exception as e:
                self.log_btn.setChecked(False)
                QtWidgets.QMessageBox.critical(self, APP_NAME, f"Could not start log:\n{e}")
        else:
            if self.log_file:
                self.log_file.close()
                self.status.showMessage(f"Saved {self.log_rows} rows to {self.log_path}")
            self.log_file = self.log_writer = None
            self.log_btn.setText("●  Log")

    def write_log(self):
        """Every N seconds write one row: Leq / max / peak over that interval."""
        n = self.cfg["log_interval"]
        e = self.engine
        while len(self.pending) >= n:
            chunk, self.pending = self.pending[:n], self.pending[n:]
            en = lambda key: 10 * math.log10(sum(10 ** (r[key] / 10) for r in chunk) / n)
            last = chunk[-1]
            ts = datetime.datetime.fromtimestamp(last["time"]).isoformat(sep=" ", timespec="seconds")
            cal = self.cfg["cal_info"] if self.cfg["calibrated"] else "UNCALIBRATED"
            self.log_writer.writerow([ts, f"{en('LAeq1s'):.1f}", f"{en('LCeq1s'):.1f}",
                                      f"{max(r['LAFmax1s'] for r in chunk):.1f}",
                                      f"{max(r['LCpeak1s'] for r in chunk):.1f}",
                                      f"{last['LAS']:.1f}", f"{last['LAeqLong']:.1f}", cal])
            self.log_file.flush()
            self.log_rows += 1

    # ---------- main loop
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
        self.bar.warn, self.bar.alarm = v - self.cfg["warn_gap"], v
        self.bar.update()
        self.save_cfg()

    def color_for(self, v):
        if not np.isfinite(v):
            return FG
        return RED if v >= self.cfg["alarm"] else YELLOW if v >= self.cfg["alarm"] - self.cfg["warn_gap"] else GREEN

    def tick(self):
        e = self.engine
        if e is None or (self.stream is None and self.sim is None):
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

        if e.cal_done:
            off = e.finish_cal(self.cal_ref)
            self.cfg.update(cal_offset=off, calibrated=True,
                            cal_info=f"{self.cal_ref:.1f} dB · {datetime.datetime.now():%b %d %H:%M}")
            self.save_cfg(); self.apply_cfg_to_ui()
            self.status.showMessage(f"Calibrated (offset {off:.2f} dB)")

        new, e.records = e.records, []
        if self.log_writer and new:
            self.pending.extend(new)
            self.write_log()

        main = e.level(self.cfg["main_weight"], self.cfg["main_time"])
        if new:
            self.last_1m = new[-1]["LAeq1m"]
        judged = main if self.cfg["thumb_src"] == "live" else self.last_1m
        now = time.time()
        if np.isfinite(judged) and judged >= self.cfg["alarm"]:
            self.bad_until = now + 3.0          # hold thumbs-down 3 s so it doesn't flicker
        self.thumb.setState("bad" if now < self.bad_until else "good")
        col = RED if now < self.bad_until else self.color_for(main)
        self.main_tile.set(main, col)
        w, t = self.cfg["main_weight"], self.cfg["main_time"]
        sub = (f"L{w}{t}  ·  1 min {fmt(self.last_1m).strip()}  ·  red at {self.cfg['alarm']:.0f}"
               + ("  ·  ● REC" if self.log_writer else ""))
        if self.mini is not None and self.mini.isVisible():
            self.mini.update_view(main, col, self.thumb.state, sub)
        if self.web:
            self.web.state = {"value": round(main, 1) if np.isfinite(main) else None,
                              "color": col, "state": self.thumb.state, "sub": sub}
        self.bar.setLevel(main)
        if new:
            r = new[-1]
            n_s, lw = len(e.hist_n), self.cfg["long_leq_min"]
            self.t_1m.set(r["LAeq1m"], self.color_for(r["LAeq1m"]))
            self.t_long.set(r["LAeqLong"], self.color_for(r["LAeqLong"]),
                            sub="" if n_s >= lw * 60 else f"filling {n_s / 60:.1f}/{lw} min")
            el = int(time.time() - e.reset_time)
            self.t_tot.set(e.leq_total("A"), sub=f"{el // 3600}:{el // 60 % 60:02d}:{el % 60:02d}")
            self.t_max.set(e.L(e.max_af))
            self.t_pk.set(e.L(e.max_pk), RED if e.L(e.max_pk) >= 135 else FG)
            msg = f"Running @ {e.fs} Hz"
            if self.log_writer:
                msg += f"   ·   Logging every {self.cfg['log_interval']} s → {self.log_rows} rows"
            if self.xruns:
                msg += f"   ·   dropouts: {self.xruns}"
            self.status.showMessage(msg)

    def closeEvent(self, ev):
        self.stop()
        if self.web:
            self.web.stop()
        if self.mini:
            self.mini.deleteLater()
        if self.log_file:
            self.log_file.close()
        self.save_cfg()
        super().closeEvent(ev)


def main():
    QtWidgets.QApplication.setHighDpiScaleFactorRoundingPolicy(
        QtCore.Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    app.setStyleSheet(STYLE)
    w = MainWindow()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
