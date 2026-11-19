"""Engine feature checks: true peak, RTA, cal files, dosimeter, pink noise, Ln (run: python tests/test_engine.py)."""
import math, os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import spl_monitor as sm

FS = 48000


def tone(f, secs, amp=1.0, fs=FS, phase=0.0):
    t = np.arange(int(fs * secs)) / fs
    return amp * np.sin(2 * np.pi * f * t + phase)


def run(e, x, block=1024):
    for i in range(0, len(x), block):
        e.process(x[i:i + block])


def ok(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}  {detail}")
    assert cond, name


def test_true_peak():
    # 12 kHz at 48 kHz: plain sample peak can read ~3 dB low depending on phase
    e = sm.Engine(FS, 0.0, 15)
    x = tone(12000, 2.0, phase=np.pi / 4)
    run(e, x)
    e.reset_stats(); run(e, x)
    pk = e.L(e.max_pk) - sm.weighting_db("C", 12000)       # undo C weighting at 12 kHz
    ok("LCpeak true-peak at 12 kHz within 0.5 dB", abs(pk) < 0.5, f"{pk:+.2f} dB (sample peak {20*np.log10(np.abs(x).max()):+.2f})")


def test_rta():
    for f, label in ((1000, "1k"), (31.5, "31.5"), (10000, "10k")):
        b = sm.ThirdOctaveBank(FS)
        fc = 1000 * 10 ** (round(10 * math.log10(f / 1000)) / 10)
        x = tone(fc, 3.0)
        for i in range(0, len(x), 1024):
            b.process(x[i:i + 1024])
        b.reset_leq()
        for i in range(0, len(x), 1024):
            b.process(x[i:i + 1024])
        lv = 10 * np.log10(b.levels("Leq"))
        i = b.labels.index(label)
        ok(f"RTA {label} band reads full-scale sine at -3.01", abs(lv[i] + 3.01) < 0.3, f"{lv[i]:.2f}")
        nb = [lv[j] for j in (i - 2, i + 2) if 0 <= j < len(lv)]
        ok(f"RTA {label} rejects bands 2 away by >= 20 dB", all(lv[i] - v > 20 for v in nb), str(np.round(np.array(nb) - lv[i], 1)))


def test_cal_file():
    txt = ('"Sens Factor =-0.7716dB, AGain =18dB, SERNO: 7023270"\n'
           "* comment line\n10.0 -6.0 0\n50 -2.0\n1000 0.0\n5000, 2.0\n10000\t4.0\n20000 3.0\n")
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(txt)
    d = sm.parse_cal_file(f.name)
    ok("cal file: sens factor", abs(d["sens_factor"] + 0.7716) < 1e-9, str(d["sens_factor"]))
    ok("cal file: serial + AGain", d["serial"] == "7023270" and d["again"] == 18.0, f"{d['serial']} {d['again']}")
    ok("cal file: 6 points", len(d["points"]) == 6, str(d["points"]))
    # A mic that is +4 dB hot at 10 kHz must read flat after correction
    e = sm.Engine(FS, 0.0, 15, d["points"])
    x = tone(10000, 3.0) * 10 ** (4 / 20)
    run(e, x); e.reset_stats(); run(e, x)
    lz = e.leq_total("Z")
    ok("cal curve flattens a +4 dB peak at 10 kHz", abs(lz + 3.01) < 0.3, f"{lz:.2f} (expect -3.01)")
    e = sm.Engine(FS, 0.0, 15, d["points"])
    x = tone(1000, 3.0); run(e, x); e.reset_stats(); run(e, x)
    ok("cal curve leaves 1 kHz alone", abs(e.leq_total("Z") + 3.01) < 0.1, f"{e.leq_total('Z'):.2f}")


def test_dose():
    for name, idx, lvl, hours, expect in (("OSHA HC", 1, 90, 8, 100), ("OSHA HC", 1, 95, 4, 100),
                                          ("NIOSH", 2, 85, 8, 100), ("NIOSH", 2, 88, 4, 100),
                                          ("NIOSH", 2, 100, 0.25, 100), ("OSHA PEL", 0, 85, 8, 0)):
        d = sm.make_dosimeters()[idx]
        for _ in range(int(hours * 3600)):
            d.add(lvl, 1.0)
        ok(f"{name}: {lvl} dBA for {hours} h = {expect} %", abs(d.dose * 100 - expect) < 0.01, f"{d.dose*100:.2f} %")
    d = sm.make_dosimeters()[2]
    for _ in range(3600):
        d.add(85, 1.0)
    ok("NIOSH TWA of 1 h @ 85 dBA = 76 dBA", abs(d.twa() - (85 + 10 * math.log10(1 / 8))) < 0.1, f"{d.twa():.2f}")
    ok("NIOSH 8 h projection from 1 h @ 85 = 100 %", abs(d.projected() - 100) < 0.1, f"{d.projected():.1f}")


def test_pink():
    x = sm.make_pink(FS, 500, 2000, -20.0)
    rms_db = 20 * np.log10(np.sqrt(np.mean(x.astype(float) ** 2)) * np.sqrt(2))
    ok("pink noise -20 dBFS AES17", abs(rms_db + 20) < 0.01, f"{rms_db:.3f}")
    ok("pink noise doesn't clip", np.abs(x).max() < 1.0, f"peak {np.abs(x).max():.3f}")
    b = sm.ThirdOctaveBank(FS)
    x = np.tile(sm.make_pink(FS, 20, 20000, -20.0), 2).astype(float)
    for i in range(0, len(x), 4096):
        b.process(x[i:i + 4096])
    lv = 10 * np.log10(b.levels("Leq"))
    mid = lv[b.labels.index("100"):b.labels.index("10k") + 1]
    ok("pink noise is flat per 1/3 octave (±1 dB, 100 Hz-10 kHz)", mid.max() - mid.min() < 2.0, f"spread {mid.max()-mid.min():.2f}")


def test_ln():
    e = sm.Engine(FS, 0.0, 15)
    # 9 s quiet (-40), 1 s loud (-10), repeated: L10 ~ loud, L90 ~ quiet
    for _ in range(3):
        run(e, tone(1000, 9.0, 10 ** (-40 / 20)))
        run(e, tone(1000, 1.0, 10 ** (-10 / 20)))
    recs = e.records
    l90 = recs[-1]["L90"]
    ok("LA90 tracks the background", abs(l90 + 43) < 1.5, f"{l90:.1f}")
    ok("since-reset LA90 histogram agrees", abs(e.stat.ln(90) - l90) < 0.5, f"{e.stat.ln(90):.1f}")
    ok("LA10 is above the background", recs[-1]["L10"] > l90 + 3, f"{recs[-1]['L10']:.1f}")


if __name__ == "__main__":
    for t in (test_true_peak, test_rta, test_cal_file, test_dose, test_pink, test_ln):
        print(t.__name__)
        t()
    print("All engine checks passed.")
