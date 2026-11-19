"""IEC 61672-1:2013 conformance checks for the SPL Monitor engine (run: python tests/test_iec61672.py)."""
import math, os, sys
import numpy as np
from scipy import signal
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import spl_monitor as sm

# Frequency, exact analog A and C (IEC 61672-1 Annex E formulas), Class 1 tolerance (+, -)
FREQS = [10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800,
         1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000]
TOL = {10: (3.5, 99), 12.5: (3.0, 99), 16: (2.5, 4.5), 20: (2.5, 2.5), 25: (2.5, 2.0), 31.5: (2.0, 1.5),
       40: (1.5, 1.5), 50: (1.5, 1.5), 63: (1.5, 1.5), 80: (1.5, 1.5), 100: (1.0, 1.0), 125: (1.0, 1.0),
       160: (1.0, 1.0), 200: (1.0, 1.0), 250: (1.0, 1.0), 315: (1.0, 1.0), 400: (1.0, 1.0), 500: (1.0, 1.0),
       630: (1.0, 1.0), 800: (1.0, 1.0), 1000: (0.7, 0.7), 1250: (1.0, 1.0), 1600: (1.0, 1.0), 2000: (1.0, 1.0),
       2500: (1.0, 1.0), 3150: (1.0, 1.0), 4000: (1.0, 1.0), 5000: (1.5, 1.5), 6300: (1.5, 2.0),
       8000: (1.5, 2.5), 10000: (2.0, 3.0), 12500: (2.0, 5.0), 16000: (2.5, 16.0), 20000: (3.0, 99)}


def analog(kind, f):
    f1, f2, f3, f4 = 20.598997, 107.65265, 737.86223, 12194.217
    def c(f):
        return 20 * math.log10(f4 ** 2 * f ** 2 / ((f ** 2 + f1 ** 2) * (f ** 2 + f4 ** 2)))
    def a(f):
        return 20 * math.log10(f4 ** 2 * f ** 4 / ((f ** 2 + f1 ** 2) * math.sqrt(f ** 2 + f2 ** 2)
                                                  * math.sqrt(f ** 2 + f3 ** 2) * (f ** 2 + f4 ** 2)))
    g = a if kind == "A" else c
    return g(f) - g(1000.0)


def exact_freq(nominal):
    """IEC uses base-10 exact frequencies: 1000 * 10^(n/10)."""
    n = round(10 * math.log10(nominal / 1000.0))
    return 1000.0 * 10 ** (n / 10)


def check_weighting(fs):
    worst = {}
    for kind in "AC":
        sos = sm.weighting_sos(kind, fs)
        bad = []
        for nom in FREQS:
            f = exact_freq(nom)
            if f >= fs / 2:
                continue
            _, h = signal.sosfreqz(sos, worN=[f], fs=fs)
            dev = 20 * math.log10(abs(h[0])) - analog(kind, f)
            up, dn = TOL[nom]
            worst[(kind, nom)] = dev
            if dev > up or dev < -dn:
                bad.append((nom, round(dev, 2)))
        lim = [abs(v) for (k, f), v in worst.items() if k == kind and f <= 12500]
        print(f"  {kind}-weighting @ {fs} Hz: max |dev| to 12.5 kHz {max(lim):.2f} dB, "
              f"16 kHz {worst.get((kind, 16000), float('nan')):+.2f}"
              f" ({'PASS' if not bad else 'FAIL ' + str(bad)})")
        assert not bad, bad


def tone(fs, f, secs, amp=1.0):
    t = np.arange(int(fs * secs)) / fs
    return amp * np.sin(2 * np.pi * f * t)


def run(engine, x, block=1024):
    for i in range(0, len(x), block):
        engine.process(x[i:i + block])


def check_time_weighting(fs):
    # Tone-burst response (IEC 61672-1 Table 4), 4 kHz, re steady level
    cases = {"F": [(0.2, -1.0), (0.002, -18.0), (0.00025, -27.0)],
             "S": [(0.2, -7.4), (0.002, -27.0)]}
    for t, rows in cases.items():
        for dur, ref in rows:
            e = sm.Engine(fs, 0.0, 15)
            run(e, tone(fs, 4000, 10.0))
            steady = 10 * math.log10(e.cur[("A", t)])
            e = sm.Engine(fs, 0.0, 15)
            mx = -1e9
            x = np.concatenate([tone(fs, 4000, dur), np.zeros(int(fs * 2))])
            for i in range(0, len(x), 64):
                e.process(x[i:i + 64])
                mx = max(mx, 10 * math.log10(max(e.cur[("A", t)], 1e-30)))
            dev = (mx - steady) - ref
            tol = 1.0 if dur >= 0.002 or t == "S" else 2.0
            print(f"  {t} burst {dur * 1000:g} ms: {mx - steady:+.2f} dB (IEC {ref:+.1f}, dev {dev:+.2f})"
                  f" {'PASS' if abs(dev) <= tol else 'FAIL'}")
            assert abs(dev) <= tol
    # Decay rate after a steady tone stops: F 34.7 dB/s, S 4.3 dB/s
    for t, rate in (("F", 34.7), ("S", 4.3)):
        e = sm.Engine(fs, 0.0, 15)
        run(e, tone(fs, 1000, 3.0))
        l0 = 10 * math.log10(e.cur[("A", t)])
        run(e, np.zeros(int(fs * 0.5)))
        l1 = 10 * math.log10(e.cur[("A", t)])
        r = (l0 - l1) / 0.5
        print(f"  {t} decay: {r:.2f} dB/s (IEC >= {rate})  {'PASS' if r >= rate - 0.1 else 'FAIL'}")
        assert r >= rate - 0.1


def check_leq_and_peak(fs):
    e = sm.Engine(fs, 0.0, 15)
    run(e, tone(fs, 1000, 5.0, amp=1.0))
    leq = e.leq_total("Z")
    print(f"  Leq of full-scale 1 kHz sine: {leq:.3f} dB (expect -3.010)")
    assert abs(leq + 3.0103) < 0.02
    e.reset_stats()          # ignore the filter's start-up transient
    run(e, tone(fs, 1000, 2.0, amp=1.0))
    pk = e.L(e.max_pk)
    print(f"  LCpeak of full-scale 1 kHz sine: {pk:.3f} dB (expect 0.000)")
    assert abs(pk) < 0.1


if __name__ == "__main__":
    for fs in (44100, 48000, 96000):
        print(f"fs = {fs}")
        check_weighting(fs)
        check_time_weighting(fs)
        check_leq_and_peak(fs)
    print("All IEC 61672 checks passed.")
