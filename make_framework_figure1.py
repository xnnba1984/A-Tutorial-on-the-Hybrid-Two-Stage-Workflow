from __future__ import annotations

from pathlib import Path
from textwrap import wrap

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "result/figures/framework_figure1_final.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

W, H = 3240, 2520
img = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(img)


def font(size: int, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/Library/Fonts/Times New Roman Bold.ttf" if bold else "/Library/Fonts/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


TITLE = font(60, bold=True)
BODY = font(50)
SMALL = font(44)


def text_size(text: str, fnt):
    box = draw.multiline_textbbox((0, 0), text, font=fnt, spacing=7)
    return box[2] - box[0], box[3] - box[1]


def center_text(box, text: str, fnt, spacing=7):
    x0, y0, x1, y1 = box
    tw, th = text_size(text, fnt)
    draw.multiline_text(
        ((x0 + x1 - tw) / 2, (y0 + y1 - th) / 2),
        text,
        fill="black",
        font=fnt,
        spacing=spacing,
        align="center",
    )


def centered_multiline_at(box, text: str, fnt, top_y: float, spacing=7):
    x0, _, x1, _ = box
    tw, _ = text_size(text, fnt)
    draw.multiline_text(
        ((x0 + x1 - tw) / 2, top_y),
        text,
        fill="black",
        font=fnt,
        spacing=spacing,
        align="center",
    )


def wrapped(text: str, width: int):
    return "\n".join(wrap(text, width=width))


def box(x0, y0, x1, y1, title, body=None):
    draw.rounded_rectangle((x0, y0, x1, y1), radius=28, outline="black", width=5, fill="white")
    if body is None:
        center_text((x0 + 20, y0 + 15, x1 - 20, y1 - 15), title, TITLE)
        return
    inner = (x0 + 34, y0 + 28, x1 - 34, y1 - 28)
    title_w, title_h = text_size(title, TITLE)
    body_w, body_h = text_size(body, BODY)
    gap = 28
    group_h = title_h + gap + body_h
    start_y = (inner[1] + inner[3] - group_h) / 2
    centered_multiline_at(inner, title, TITLE, start_y)
    centered_multiline_at(inner, body, BODY, start_y + title_h + gap, spacing=5)


def diamond(cx, cy, w, h, label):
    pts = [(cx, cy - h / 2), (cx + w / 2, cy), (cx, cy + h / 2), (cx - w / 2, cy)]
    draw.polygon(pts, outline="black", fill="white")
    draw.line(pts + [pts[0]], fill="black", width=5)
    center_text((cx - w / 2 + 60, cy - h / 2 + 35, cx + w / 2 - 60, cy + h / 2 - 35), wrapped(label, 17), TITLE)


def arrow(x0, y0, x1, y1, label=None, label_offset=(0, 0)):
    draw.line((x0, y0, x1, y1), fill="black", width=5)
    import math

    ang = math.atan2(y1 - y0, x1 - x0)
    length = 34
    spread = 0.48
    p1 = (x1 - length * math.cos(ang - spread), y1 - length * math.sin(ang - spread))
    p2 = (x1 - length * math.cos(ang + spread), y1 - length * math.sin(ang + spread))
    draw.polygon([(x1, y1), p1, p2], fill="white", outline="black")
    draw.line([(x1, y1), p1, p2, (x1, y1)], fill="black", width=5)
    if label:
        tx, ty = (x0 + x1) / 2 + label_offset[0], (y0 + y1) / 2 + label_offset[1]
        tw, th = text_size(label, SMALL)
        draw.text((tx - tw / 2, ty - th / 2), label, fill="black", font=SMALL)


box(
    130,
    110,
    830,
    520,
    "Inputs and assumptions",
    "Covariates X\nTreatment A\nOutcome Y\nClinical margin δ\nIdentification assumptions",
)
box(
    1220,
    110,
    2020,
    520,
    "Stage 1: Population level\ninference",
    "Formulate hypotheses\nTest heterogeneity\nControl multiplicity",
)
diamond(1620, 840, 650, 290, "HTE evidence adequate?")
box(
    2380,
    680,
    3100,
    1000,
    "No confirmatory\nHTE claim",
    "Document finding\nOptional exploratory\nStage 2",
)
box(
    1220,
    1160,
    2020,
    1510,
    "Stage 2: Decision\nmaking",
    "Estimate CATE\nValidate with uplift curves\nEvaluate policy value",
)
box(720, 1810, 1340, 2035, "Value driven\nthreshold")
box(1900, 1810, 2520, 2035, "Empirical surrogate\nharm constrained\nthreshold")
box(1220, 2240, 2020, 2485, "Final policy", "Treatment decision rule\nfor each patient")

arrow(830, 315, 1220, 315)
arrow(1620, 520, 1620, 695)
arrow(1945, 840, 2380, 840, "No", label_offset=(0, -58))
arrow(1620, 985, 1620, 1160, "Yes", label_offset=(72, -10))
arrow(1400, 1510, 1120, 1810)
arrow(1840, 1510, 2120, 1810)
arrow(1030, 2035, 1420, 2240)
arrow(2210, 2035, 1820, 2240)

img.save(OUT, dpi=(450, 450))
print(OUT)
