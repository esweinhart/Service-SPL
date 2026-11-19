"""Draw three ServiceSPL icon concepts at 1024 px (4x supersampled)."""
import math, sys
from PIL import Image, ImageDraw, ImageFilter, ImageFont

S = 4096            # supersample canvas
OUT = 1024
GREEN, YEL, RED = (52, 211, 120), (250, 204, 21), (239, 68, 68)
CYAN = (56, 189, 248)


def lerp(a, b, t):
    return tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


def ramp(t):
    return lerp(GREEN, YEL, t / 0.6) if t < 0.6 else lerp(YEL, RED, (t - 0.6) / 0.4)


def base():
    """macOS-style squircle with a dark slate vertical gradient."""
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    pad = int(S * 0.098)
    grad = Image.new("RGBA", (S, S))
    gd = ImageDraw.Draw(grad)
    for y in range(S):
        t = y / S
        gd.line([(0, y), (S, y)], fill=lerp((30, 41, 59), (8, 12, 22), t) + (255,))
    mask = Image.new("L", (S, S), 0)
    ImageDraw.Draw(mask).rounded_rectangle([pad, pad, S - pad, S - pad], radius=int(S * 0.18), fill=255)
    img.paste(grad, (0, 0), mask)
    # thin top highlight rim
    rim = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    ImageDraw.Draw(rim).rounded_rectangle([pad, pad, S - pad, S - pad], radius=int(S * 0.18),
                                          outline=(148, 163, 184, 70), width=int(S * 0.006))
    img.alpha_composite(rim)
    return img, mask


def glow(layer, r):
    g = layer.filter(ImageFilter.GaussianBlur(r))
    out = Image.new("RGBA", layer.size, (0, 0, 0, 0))
    out.alpha_composite(g)
    out.alpha_composite(g)
    out.alpha_composite(layer)
    return out


def finish(img, mask, art, name):
    art = glow(art, S * 0.012)
    clipped = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    clipped.paste(art, (0, 0), Image.composite(art.split()[3], Image.new("L", (S, S), 0), mask))
    img.alpha_composite(clipped)
    img = img.resize((OUT, OUT), Image.LANCZOS)
    img.save(name)
    return img


def rings(name):
    """A: sound-pressure rings radiating from a source dot, green -> red."""
    img, mask = base()
    art = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(art)
    cx, cy = S * 0.31, S * 0.5
    d.ellipse([cx - S*.07, cy - S*.07, cx + S*.07, cy + S*.07], fill=GREEN + (255,))
    for i, r in enumerate([0.16, 0.265, 0.37, 0.475]):
        R = S * r
        c = ramp(i / 3)
        d.arc([cx - R, cy - R, cx + R, cy + R], -48, 48, fill=c + (255,), width=int(S * 0.045))
    return finish(img, mask, art, name)


def vu(name):
    """B: analogue VU arc with coloured segments and a needle."""
    img, mask = base()
    art = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(art)
    cx, cy, R = S * 0.5, S * 0.66, S * 0.34
    n = 13
    for i in range(n):
        a0 = 205 + i * (130 / n) + 1.6
        a1 = 205 + (i + 1) * (130 / n) - 1.6
        d.arc([cx - R, cy - R, cx + R, cy + R], a0, a1, fill=ramp(i / (n - 1)) + (255,), width=int(S * 0.07))
    ang = math.radians(205 + 130 * 0.72)
    L = R * 0.92
    d.line([(cx, cy), (cx + L * math.cos(ang), cy + L * math.sin(ang))], fill=(241, 245, 249, 255), width=int(S * 0.028))
    d.ellipse([cx - S*.05, cy - S*.05, cx + S*.05, cy + S*.05], fill=(241, 245, 249, 255))
    try:
        f = ImageFont.truetype("/usr/share/fonts/opentype/inter/Inter-Bold.otf", int(S * 0.11))
        d.text((cx, S * 0.80), "dB", font=f, fill=(148, 163, 184, 255), anchor="mm")
    except OSError:
        pass
    return finish(img, mask, art, name)


def wave(name):
    """C: cyan waveform feeding a segmented level bar."""
    img, mask = base()
    art = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(art)
    x0, x1, cy = S * 0.17, S * 0.62, S * 0.5
    pts = []
    for k in range(400):
        t = k / 399
        env = math.sin(math.pi * t) ** 1.4
        y = cy - S * 0.22 * env * math.sin(2 * math.pi * 5.5 * t) * (0.6 + 0.4 * math.sin(7 * t))
        pts.append((x0 + (x1 - x0) * t, y))
    d.line(pts, fill=CYAN + (255,), width=int(S * 0.03), joint="curve")
    bx0, bx1 = S * 0.70, S * 0.83
    segs, top, bot = 9, S * 0.20, S * 0.80
    h = (bot - top) / segs
    lit = 7
    for i in range(segs):
        y1 = bot - i * h
        c = ramp(i / (segs - 1))
        a = 255 if i < lit else 55
        d.rounded_rectangle([bx0, y1 - h + S*.012, bx1, y1 - S*.012], radius=int(S*.012), fill=c + (a,))
    return finish(img, mask, art, name)


if __name__ == "__main__":
    out = sys.argv[1]
    imgs = [rings(f"{out}/A_rings.png"), vu(f"{out}/B_vu.png"), wave(f"{out}/C_wave.png")]
    # contact sheet: each concept big + Dock(128), menu(32) sizes, dark & light backgrounds
    W, H = 1800, 760
    sheet = Image.new("RGBA", (W, H), (255, 255, 255, 255))
    sd = ImageDraw.Draw(sheet)
    sd.rectangle([0, H // 2 + 40, W, H], fill=(30, 30, 34, 255))
    f = ImageFont.truetype("/usr/share/fonts/opentype/inter/Inter-Bold.otf", 40)
    for i, (im, lab) in enumerate(zip(imgs, ["A  Pressure rings", "B  VU needle", "C  Waveform + meter"])):
        x = 40 + i * 590
        sd.text((x, 20), lab, font=f, fill=(20, 20, 20, 255))
        sheet.alpha_composite(im.resize((300, 300), Image.LANCZOS), (x, 80))
        sheet.alpha_composite(im.resize((128, 128), Image.LANCZOS), (x + 330, 120))
        sheet.alpha_composite(im.resize((64, 64), Image.LANCZOS), (x + 330, 280))
        sheet.alpha_composite(im.resize((32, 32), Image.LANCZOS), (x + 420, 296))
        sheet.alpha_composite(im.resize((128, 128), Image.LANCZOS), (x + 30, 480))
        sheet.alpha_composite(im.resize((64, 64), Image.LANCZOS), (x + 200, 512))
        sheet.alpha_composite(im.resize((32, 32), Image.LANCZOS), (x + 290, 528))
    sheet.convert("RGB").save(f"{out}/concepts.png")
