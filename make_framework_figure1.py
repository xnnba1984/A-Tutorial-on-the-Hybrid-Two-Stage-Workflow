from pathlib import Path
import os

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "result" / "figures" / "Figure_1_revised.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

W, H = 3240, 2520
img = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(img)


def load_font(size: int, bold: bool = False):
    names = [
        os.environ.get("FIGURE_FONT_BOLD" if bold else "FIGURE_FONT_REGULAR", ""),
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"
        if bold
        else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/Library/Fonts/Times New Roman Bold.ttf"
        if bold
        else "/Library/Fonts/Times New Roman.ttf",
        "C:/Windows/Fonts/timesbd.ttf" if bold else "C:/Windows/Fonts/times.ttf",
        "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold.ttf"
        if bold else "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSerif-Bold.ttf"
        if bold else "/usr/share/fonts/truetype/liberation2/LiberationSerif-Regular.ttf",
    ]
    for name in names:
        if name and Path(name).is_file():
            print(f"Figure 1 font: {Path(name).name} ({size} px)")
            return ImageFont.truetype(name, size=size)
    raise FileNotFoundError(
        "No compatible serif font. Install Times New Roman or Liberation Serif, "
        "or set FIGURE_FONT_REGULAR and FIGURE_FONT_BOLD to font-file paths."
    )


TITLE = load_font(72, bold=True)
BODY = load_font(59)
SMALL = load_font(52)

INK = "#24303b"
LINE = "#4f5b69"
BLUE_FILL = "#eaf2fb"
BLUE_LINE = "#356b9c"
GREEN_FILL = "#edf7f0"
GREEN_LINE = "#3f7c59"
AMBER_FILL = "#fff6e6"
AMBER_LINE = "#a77522"
GRAY_FILL = "#f4f5f6"


def text_extent(text, font, spacing=8):
    box = draw.multiline_textbbox((0, 0), text, font=font, spacing=spacing)
    return box[2] - box[0], box[3] - box[1]


def centered_text(box, text, font, fill=INK, spacing=8):
    x0, y0, x1, y1 = box
    tw, th = text_extent(text, font, spacing)
    draw.multiline_text(
        ((x0 + x1 - tw) / 2, (y0 + y1 - th) / 2),
        text,
        font=font,
        fill=fill,
        spacing=spacing,
        align="center",
    )


def content_box(box, title, body, fill, outline):
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(
        box,
        radius=30,
        fill=fill,
        outline=outline,
        width=6,
    )
    title_h = text_extent(title, TITLE)[1]
    body_h = text_extent(body, BODY, 5)[1]
    gap = 36
    group_h = title_h + gap + body_h
    top = (y0 + y1 - group_h) / 2
    centered_text((x0 + 45, top, x1 - 45, top + title_h), title, TITLE)
    centered_text(
        (x0 + 45, top + title_h + gap, x1 - 45, y1 - 30),
        body,
        BODY,
        spacing=5,
    )


def arrow(points, dashed=False):
    import math

    for start, end in zip(points[:-1], points[1:]):
        if dashed:
            x0, y0 = start
            x1, y1 = end
            distance = math.hypot(x1 - x0, y1 - y0)
            ux, uy = (x1 - x0) / distance, (y1 - y0) / distance
            pos = 0
            while pos < distance - 12:
                dash_end = min(pos + 28, distance - 12)
                draw.line(
                    (x0 + ux * pos, y0 + uy * pos, x0 + ux * dash_end, y0 + uy * dash_end),
                    fill=LINE,
                    width=6,
                )
                pos += 46
        else:
            draw.line((*start, *end), fill=LINE, width=6)

    x0, y0 = points[-2]
    x1, y1 = points[-1]
    angle = math.atan2(y1 - y0, x1 - x0)
    length = 42
    spread = 0.48
    p1 = (x1 - length * math.cos(angle - spread), y1 - length * math.sin(angle - spread))
    p2 = (x1 - length * math.cos(angle + spread), y1 - length * math.sin(angle + spread))
    draw.polygon([(x1, y1), p1, p2], fill="white", outline=LINE)
    draw.line([(x1, y1), p1, p2, (x1, y1)], fill=LINE, width=6)


shared = (560, 90, 2680, 500)
stage1 = (145, 770, 1450, 1465)
stage2 = (1790, 770, 3095, 1465)
relation = (975, 1570, 2265, 1940)
reporting = (590, 2110, 2650, 2450)

content_box(
    shared,
    "Shared inputs and assumptions",
    "Covariates X    Treatment A    Outcome Y\n"
    "Clinical margin δ    Target population\n"
    "Estimand and identification assumptions",
    GRAY_FILL,
    INK,
)
content_box(
    stage1,
    "Stage 1: Treatment effect claim",
    "Prespecify the target claim\n"
    "Use tests or CATE/GATE summaries\n"
    "suited to the claim\n"
    "Control error for the stated claim\n"
    "Report estimates, intervals, and tests",
    BLUE_FILL,
    BLUE_LINE,
)
content_box(
    stage2,
    "Stage 2: Policy claim",
    "Define the policy target and comparator\n"
    "Estimate CATE scores or learn a policy\n"
    "Validate ranking and policy value\n"
    "Assess uncertainty and operating metrics",
    GREEN_FILL,
    GREEN_LINE,
)
content_box(
    relation,
    "Key relationship",
    "A Stage 1 result may motivate Stage 2,\n"
    "but it does not validate the policy or determine\n"
    "the inferential status of Stage 2.",
    AMBER_FILL,
    AMBER_LINE,
)
content_box(
    reporting,
    "Separate reporting and interpretation",
    "Match each claim to its own design, validation,\n"
    "uncertainty assessment, and success criterion.",
    GRAY_FILL,
    INK,
)

mid_x = (shared[0] + shared[2]) / 2
split_y = 620
arrow([
    (mid_x, shared[3]),
    (mid_x, split_y),
    ((stage1[0] + stage1[2]) / 2, split_y),
    ((stage1[0] + stage1[2]) / 2, stage1[1]),
])
arrow([
    (mid_x, shared[3]),
    (mid_x, split_y),
    ((stage2[0] + stage2[2]) / 2, split_y),
    ((stage2[0] + stage2[2]) / 2, stage2[1]),
])

arrow([(stage1[2], 1110), (stage2[0] - 20, 1110)], dashed=True)
centered_text((1450, 960, 1790, 1070), "Scientific\nmotivation", SMALL, spacing=2)

arrow([((stage1[0] + stage1[2]) / 2, stage1[3]), ((stage1[0] + stage1[2]) / 2, 1995), (1110, reporting[1])])
arrow([((stage2[0] + stage2[2]) / 2, stage2[3]), ((stage2[0] + stage2[2]) / 2, 1995), (2130, reporting[1])])

img.save(OUT, dpi=(450, 450))
print(OUT)
