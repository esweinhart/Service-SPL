"""Engine feature checks: true peak, RTA, cal files, dosimeter, pink noise, Ln (run: python tests/test_engine.py)."""
import json, math, os, sys, tempfile
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


def test_eq_and_osc():
    fc = np.array([1000 * 10 ** (n / 10) for n in range(-17, 14)])
    labels = sm.RTA_LABELS
    meas = np.zeros(len(fc)); meas[labels.index("125")] = 5.0; meas[labels.index("63")] = -9.0
    meas[labels.index("4k")] = -3.0
    dev, tips = sm.eq_suggestions(fc, labels, meas + 80, np.zeros(len(fc)))
    ok("EQ: suggests cutting a 125 Hz bump", any(t.startswith("Cut 125 Hz") for t in tips), str(tips))
    ok("EQ: calls a deep 63 Hz dip a room null", any(t.startswith("63 Hz") and "null" in t for t in tips))
    ok("EQ: small boost at 4k", any(t.startswith("Boost 4k") for t in tips))
    m = sm.osc_message("/custom-variable/spl_level/value", "87.3")
    ok("OSC message is 4-byte aligned", len(m) % 4 == 0 and m.startswith(b"/custom-variable/spl_level/value\0"))
    m = sm.osc_message("/x", 1, 2, 3)
    ok("OSC int args", m == b"/x\0\0,iii\0\0\0\0" + b"\0\0\0\x01\0\0\0\x02\0\0\0\x03", repr(m))
    import datetime
    e = dict(on=True, days=[0], start="09:00", stop="12:30")
    ok("schedule: Monday 10:00 is inside", sm.schedule_window(e, datetime.datetime(2026, 11, 23, 10, 0)))
    ok("schedule: Monday 12:30 is outside", not sm.schedule_window(e, datetime.datetime(2026, 11, 23, 12, 30)))
    ok("schedule: Tuesday is outside", not sm.schedule_window(e, datetime.datetime(2026, 11, 24, 10, 0)))
    ok("schedule: one-off date", sm.schedule_window(dict(e, date="2026-11-24"), datetime.datetime(2026, 11, 24, 10, 0)))


def test_cal_formats():
    cases = {
        "Audix TM1 style": ("Audix TM1 Measurement Microphone\nSerial Number: 210345\nSensitivity: 9.6 mV/Pa\n"
                            "Freq(Hz)\tdB\n20\t-0.8\n100\t-0.1\n1000\t0.0\n10000\t0.9\n20000\t-1.2\n"),
        "CSV with header": "Frequency,Response\n20,-0.8\n100,-0.1\n1000,0.0\n10000,0.9\n20000,-1.2\n",
        "decimal commas": "20,0;-0,8\n100,0;-0,1\n1000,0;0,0\n10000,0;0,9\n20000,0;-1,2\n",
        "tab + decimal commas": "20,0\t-0,8\n100\t-0,1\n1000\t0,0\n10000\t0,9\n20000\t-1,2\n",
        "absolute SPL, reversed": "20000 92.8\n10000 94.9\n1000 94.0\n100 93.9\n20 93.2\n",
        "CRLF + BOM": "\ufeff20 -0.8\r\n100 -0.1\r\n1000 0.0\r\n10000 0.9\r\n20000 -1.2\r\n",
    }
    for name, txt in cases.items():
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as f:
            f.write(txt)
        d = sm.parse_cal_file(f.name)
        g = dict((round(p[0]), round(p[1], 2)) for p in d["points"])
        okv = len(d["points"]) == 5 and g.get(1000) == 0.0 and g.get(10000) == 0.9 and g.get(20) == -0.8
        ok(f"cal file format: {name}", okv, str(d["points"]))
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(cases["Audix TM1 style"])
    d = sm.parse_cal_file(f.name)
    ok("TM1: model, serial, mV/Pa", d["model"] == "TM1" and d["serial"] == "210345" and abs(d["mv_pa"] - 9.6) < 1e-9,
       f"{d['model']} {d['serial']} {d['mv_pa']}")


def test_settings_backup():
    from PySide6 import QtWidgets
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    d = tempfile.mkdtemp()
    sm.CFG_DIR, sm.CFG_PATH = d, os.path.join(d, "config.json")
    sm.OLD_CFG_PATH = os.path.join(d, "none.json")
    w = sm.MainWindow()
    w.cfg["profiles"]["TM1 #210345"] = dict(sm.DEFAULTS["profiles"]["Default"], offset=101.5, calibrated=True,
                                           cal_points=[[20, -0.8], [1000, 0.0], [20000, -1.2]])
    w.cfg["alarm"], w.cfg["companion_port"], w.cfg["pco_secret"] = 97.0, 4444, "s3cret"
    path = os.path.join(d, "backup.servicespl.json")
    w.export_settings(path, include_secret=False)
    data = json.load(open(path))
    ok("backup: no secret unless asked", data["secrets"] == {} and "pco_secret" not in data["settings"])
    w.cfg["alarm"] = 80.0
    w.cfg["profiles"].pop("TM1 #210345")
    ok("backup: import", w.import_settings(path, confirm=False, restart=False))
    c = json.load(open(sm.CFG_PATH))
    ok("backup: restores profiles + cal curve + settings",
       c["profiles"]["TM1 #210345"]["offset"] == 101.5 and c["profiles"]["TM1 #210345"]["cal_points"][1] == [1000, 0.0]
       and c["alarm"] == 97.0 and c["companion_port"] == 4444 and c.get("dock_state"), str(c["alarm"]))
    w._imported = False
    w.close()


if __name__ == "__main__":
    for t in (test_true_peak, test_rta, test_cal_file, test_dose, test_pink, test_ln, test_eq_and_osc,
              test_cal_formats, test_settings_backup):
        print(t.__name__)
        t()
    print("All engine checks passed.")
