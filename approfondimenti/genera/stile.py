# Strumenti comuni per i PDF di approfondimento del Resoconto_modello_1D.md: impaginazione (reportlab),
# formule (mathtext di matplotlib), flowchart (matplotlib), estratti di codice presi dai sorgenti.
import hashlib, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, Image, Table,
                                TableStyle, Preformatted, KeepTogether, PageBreak, CondPageBreak)

QUI = os.path.dirname(os.path.abspath(__file__))
PDF_DIR = os.path.dirname(QUI)                       # approfondimenti/
PROGETTO = os.path.dirname(PDF_DIR)                  # 1D_Brenta/
CACHE = os.path.join(QUI, "_cache"); os.makedirs(CACHE, exist_ok=True)

F = "C:/Windows/Fonts/"
pdfmetrics.registerFont(TTFont("Corpo", F + "calibri.ttf"))
pdfmetrics.registerFont(TTFont("Corpo-B", F + "calibrib.ttf"))
pdfmetrics.registerFont(TTFont("Corpo-I", F + "calibrii.ttf"))
pdfmetrics.registerFont(TTFont("Corpo-BI", F + "calibriz.ttf"))
pdfmetrics.registerFont(TTFont("Codice", F + "consola.ttf"))
pdfmetrics.registerFont(TTFont("Codice-B", F + "consolab.ttf"))
pdfmetrics.registerFontFamily("Corpo", normal="Corpo", bold="Corpo-B", italic="Corpo-I", boldItalic="Corpo-BI")
pdfmetrics.registerFontFamily("Codice", normal="Codice", bold="Codice-B", italic="Codice", boldItalic="Codice-B")

BLU = colors.HexColor("#2a78d6"); BLU_C = colors.HexColor("#e8f0fb"); GRIGIO = colors.HexColor("#5f6368")
ARANCIO = colors.HexColor("#eb6834"); VERDE = colors.HexColor("#1baf7a"); FONDO_CODICE = colors.HexColor("#f6f7f9")
INCHIOSTRO = colors.HexColor("#1f2328")

S = {
    "titolo": ParagraphStyle("titolo", fontName="Corpo-B", fontSize=22, leading=27, textColor=INCHIOSTRO, spaceAfter=4),
    "sotto": ParagraphStyle("sotto", fontName="Corpo", fontSize=12, leading=16, textColor=GRIGIO, spaceAfter=10),
    "h1": ParagraphStyle("h1", fontName="Corpo-B", fontSize=15, leading=19, textColor=BLU, spaceBefore=14, spaceAfter=6),
    "h2": ParagraphStyle("h2", fontName="Corpo-B", fontSize=12, leading=15, textColor=INCHIOSTRO, spaceBefore=10, spaceAfter=4),
    "p": ParagraphStyle("p", fontName="Corpo", fontSize=10.5, leading=14.5, alignment=TA_JUSTIFY, spaceAfter=6, textColor=INCHIOSTRO),
    "li": ParagraphStyle("li", fontName="Corpo", fontSize=10.5, leading=14.5, leftIndent=14, bulletIndent=4, spaceAfter=3, textColor=INCHIOSTRO),
    "nota": ParagraphStyle("nota", fontName="Corpo", fontSize=10, leading=13.5, textColor=INCHIOSTRO),
    "cap": ParagraphStyle("cap", fontName="Corpo-I", fontSize=9, leading=12, alignment=TA_CENTER, textColor=GRIGIO, spaceAfter=8),
    "cod": ParagraphStyle("cod", fontName="Codice", fontSize=8.0, leading=10.1, textColor=INCHIOSTRO),
    "cella": ParagraphStyle("cella", fontName="Corpo", fontSize=9.5, leading=12, textColor=INCHIOSTRO),
    "cellab": ParagraphStyle("cellab", fontName="Corpo-B", fontSize=9.5, leading=12, textColor=INCHIOSTRO),
    "eqn": ParagraphStyle("eqn", fontName="Corpo", fontSize=9.5, leading=12, alignment=TA_LEFT, textColor=GRIGIO),
}
LARGHEZZA = A4[0] - 4.2 * cm


def _md(t):
    """Markup leggero: `codice`, **grassetto**, *corsivo*, _{pedice}, ^{apice}; virgola decimale fuori dal codice."""
    pezzi = t.split("`")
    t = "`".join(re.sub(r"(?<=\d)\.(?=\d)", ",", p) if k % 2 == 0 else p for k, p in enumerate(pezzi))
    t = re.sub(r"`([^`]+)`", r'<font face="Codice" size="9" color="#b4371f">\1</font>', t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"(?<![\w*])\*([^*\s][^*]*?)\*(?![\w*])", r"<i>\1</i>", t)
    t = re.sub(r"_\{([^}]+)\}", r"<sub>\1</sub>", t)
    t = re.sub(r"\^\{([^}]+)\}", r"<super>\1</super>", t)
    return t


class Documento:
    def __init__(self, numero, nome_file, titolo, sottotitolo, sezione_md):
        self.numero, self.nome_file, self.titolo = numero, nome_file, titolo
        self.sottotitolo, self.sezione_md = sottotitolo, sezione_md
        self.story = []; self.n_eq = 0; self.n_fig = 0
        self.story += [Paragraph(f"Approfondimento {numero:02d}", S["sotto"]), Paragraph(titolo, S["titolo"]),
                       Paragraph(_md(sottotitolo), S["sotto"])]
        self.story.append(self._riquadro(
            f'Fa parte del <a href="../Resoconto_modello_1D.md" color="#2a78d6"><u>Resoconto del modello 1D_Brenta</u></a>, '
            f'{sezione_md}. Gli estratti di codice sono presi direttamente dai sorgenti quando il PDF viene generato '
            f'(<font face="Codice" size="9">approfondimenti/genera/genera.py</font>).', BLU_C, BLU))

    # --- testo ---
    def _titolo(self, t, stile):
        par = Paragraph(_md(t), S[stile]); par._lega = True; self.story.append(par)
    def h1(self, t): self._titolo(t, "h1")
    def h2(self, t): self._titolo(t, "h2")
    def p(self, t):
        par = Paragraph(_md(t), S["p"])
        if t.rstrip().endswith(":"): par._lega = True        # introduce la formula, il codice o l'elenco che segue
        self.story.append(par)
    def elenco(self, voci):
        for v in voci: self.story.append(Paragraph(_md(v), S["li"], bulletText="•"))
        self.story.append(Spacer(1, 4))
    def numerato(self, voci):
        for i, v in enumerate(voci, 1): self.story.append(Paragraph(_md(v), S["li"], bulletText=f"{i}."))
        self.story.append(Spacer(1, 4))
    def spazio(self, h=6): self.story.append(Spacer(1, h))
    def pagina(self): self.story.append(PageBreak())

    def _riquadro(self, testo, fondo, bordo, stile="nota"):
        t = Table([[Paragraph(testo, S[stile])]], colWidths=[LARGHEZZA])
        t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), fondo), ("LINEBEFORE", (0, 0), (0, -1), 3, bordo),
                               ("LEFTPADDING", (0, 0), (-1, -1), 9), ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                               ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 6)]))
        return KeepTogether([t, Spacer(1, 8)])
    def nota(self, t): self.story.append(self._riquadro(_md(t), BLU_C, BLU))
    def attenzione(self, t): self.story.append(self._riquadro("<b>Attenzione:</b> " + _md(t), colors.HexColor("#fdeee6"), ARANCIO))
    def idea(self, t): self.story.append(self._riquadro("<b>In breve:</b> " + _md(t), colors.HexColor("#e6f6ef"), VERDE))

    # --- formule ---
    def eq(self, tex, numerata=True, scala=1.0, didascalia=None):
        """Formula in mathtext (sottoinsieme di LaTeX). Una stringa o una lista di righe."""
        righe = tex if isinstance(tex, list) else [tex]
        imgs = []
        for r in righe:
            png = _formula([r], scala); imgs.append([Image(png[0], width=png[1], height=png[2])])
        img = Table(imgs, colWidths=[LARGHEZZA - 1.4 * cm])
        img.setStyle(TableStyle([("ALIGN", (0, 0), (-1, -1), "CENTER"), ("TOPPADDING", (0, 0), (-1, -1), 2),
                                 ("BOTTOMPADDING", (0, 0), (-1, -1), 2), ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
        if numerata:
            self.n_eq += 1
            t = Table([[img, Paragraph(f"({self.n_eq})", S["eqn"])]], colWidths=[LARGHEZZA - 1.2 * cm, 1.2 * cm])
        else:
            t = Table([[img]], colWidths=[LARGHEZZA])
        t.setStyle(TableStyle([("ALIGN", (0, 0), (0, 0), "CENTER"), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                               ("ALIGN", (1, 0), (1, 0), "RIGHT"), ("TOPPADDING", (0, 0), (-1, -1), 3),
                               ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
        self.story.append(t)
        if didascalia: self.story.append(Paragraph(_md(didascalia), S["cap"]))
        return self.n_eq

    # --- codice ---
    def codice(self, testo, titolo=None):
        testo = testo.rstrip("\n").replace("\t", "    ")
        L = _a_capo(_senza_commenti(testo.split("\n")))
        pezzi = ["\n".join(L[i:i + 4]) for i in range(0, len(L), 4)]   # righe della tabella: si spezza ogni 4 righe
        righe = [[Paragraph(f'<font face="Codice" size="8" color="#5f6368">{titolo}</font>', S["cella"])]] if titolo else []
        righe += [[Preformatted(p, S["cod"])] for p in pezzi]
        t = Table(righe, colWidths=[LARGHEZZA], splitByRow=1, repeatRows=1 if titolo else 0)
        st = [("BACKGROUND", (0, 0), (-1, -1), FONDO_CODICE), ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#d0d7de")),
              ("LEFTPADDING", (0, 0), (-1, -1), 8), ("RIGHTPADDING", (0, 0), (-1, -1), 6),
              ("TOPPADDING", (0, 0), (-1, -1), 0), ("BOTTOMPADDING", (0, 0), (-1, -1), 0)]
        if titolo: st += [("LINEBELOW", (0, 0), (-1, 0), 0.5, colors.HexColor("#d0d7de")), ("TOPPADDING", (0, 0), (-1, 0), 4),
                          ("BOTTOMPADDING", (0, 0), (-1, 0), 4), ("TOPPADDING", (0, 1), (-1, 1), 5)]
        else: st.append(("TOPPADDING", (0, 0), (-1, 0), 5))
        st.append(("BOTTOMPADDING", (0, -1), (-1, -1), 5))
        t.setStyle(TableStyle(st))
        self.story += [t, Spacer(1, 8)]

    def estratto(self, file, inizio, fine=None, n=None, titolo=None, togli_commenti=False):
        testo = estrai(file, inizio, fine, n, togli_commenti)
        self.codice(testo, titolo or f"{file}")

    # --- tabelle ---
    def tabella(self, righe, larghezze=None, intestazione=True, didascalia=None):
        dati = [[Paragraph(_md(str(c)), S["cellab" if (intestazione and i == 0) else "cella"]) for c in r] for i, r in enumerate(righe)]
        lw = [LARGHEZZA * w for w in larghezze] if larghezze else None
        t = Table(dati, colWidths=lw, repeatRows=1 if intestazione else 0)
        st = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#d0d7de")), ("VALIGN", (0, 0), (-1, -1), "TOP"),
              ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]
        if intestazione: st.append(("BACKGROUND", (0, 0), (-1, 0), BLU_C))
        t.setStyle(TableStyle(st))
        self.story.append(t)
        self.story.append(Paragraph(_md(didascalia), S["cap"]) if didascalia else Spacer(1, 8))

    # --- figure ---
    def figura(self, percorso, larghezza=1.0, didascalia=None):
        from reportlab.lib.utils import ImageReader
        w, h = ImageReader(percorso).getSize()
        W = LARGHEZZA * larghezza; H = W * h / w
        if H > 22 * cm: H = 22 * cm; W = H * w / h
        self.n_fig += 1
        el = [Image(percorso, width=W, height=H)]
        if didascalia: el.append(Paragraph(f"Figura {self.n_fig}. " + _md(didascalia), S["cap"]))
        self.story.append(KeepTogether(el))
    def flowchart(self, schema, didascalia=None, larghezza=1.0):
        self.figura(disegna_flowchart(schema), larghezza, didascalia)
    def grafico(self, funzione, nome, didascalia=None, larghezza=0.95):
        """funzione(fig) disegna su una figura matplotlib; il PNG va nella cache."""
        percorso = os.path.join(CACHE, f"{self.numero:02d}_{nome}.png")
        fig = plt.figure(figsize=(7.6, 3.6)); funzione(fig)
        fig.savefig(percorso, dpi=220, bbox_inches="tight"); plt.close(fig)
        self.figura(percorso, larghezza, didascalia)

    def collegamenti(self, voci):
        """voci: (numero, nome_file, testo)"""
        el = [Paragraph("Collegamenti con gli altri approfondimenti", S["h1"])]
        for n, f, t in voci:
            el.append(Paragraph(f'<a href="{f}" color="#2a78d6"><u>{n:02d} · {t}</u></a>', S["li"], bulletText="→"))
        self.story.append(KeepTogether(el))

    # --- uscita ---
    def salva(self):
        percorso = os.path.join(PDF_DIR, self.nome_file)
        titolo, numero = self.titolo, self.numero
        def piede(c, doc):
            c.saveState()
            c.setFont("Corpo", 8.5); c.setFillColor(GRIGIO)
            c.drawString(2.1 * cm, 1.2 * cm, f"Modello 1D_Brenta · Approfondimento {numero:02d} · {titolo}")
            c.drawRightString(A4[0] - 2.1 * cm, 1.2 * cm, f"{doc.page}")
            c.setStrokeColor(colors.HexColor("#d0d7de")); c.setLineWidth(0.4)
            c.line(2.1 * cm, 1.55 * cm, A4[0] - 2.1 * cm, 1.55 * cm)
            c.restoreState()
        doc = BaseDocTemplate(percorso, pagesize=A4, leftMargin=2.1 * cm, rightMargin=2.1 * cm, topMargin=1.8 * cm,
                              bottomMargin=2.1 * cm, title=f"{numero:02d} {titolo}", author="Modello 1D_Brenta",
                              subject="Approfondimento del Resoconto_modello_1D.md")
        fr = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="f")
        doc.addPageTemplates([PageTemplate(id="p", frames=[fr], onPage=piede)])
        doc.build(_impagina(self.story))
        return percorso


def _altezza(f):
    if isinstance(f, KeepTogether):
        return sum(_altezza(g) for g in f._content)
    try:
        h = f.wrap(LARGHEZZA, 1e6)[1]
    except Exception:
        return 0.0
    return h + getattr(getattr(f, "style", None), "spaceBefore", 0) + getattr(getattr(f, "style", None), "spaceAfter", 0)


def _impagina(story, minimo=5 * cm):
    """Titoli e frasi che finiscono con ':' non restano soli in fondo alla pagina: prima di ciascuno si chiede
    lo spazio per sé, per gli elementi legati che seguono e per almeno `minimo` del primo blocco di contenuto
    (o tutto il blocco, se è più basso)."""
    out = []; k = 0
    while k < len(story):
        f = story[k]
        if getattr(f, "_lega", False) and (k == 0 or not getattr(story[k - 1], "_lega", False)):
            serve = 0.0; j = k
            while j < len(story) and getattr(story[j], "_lega", False):
                serve += _altezza(story[j]); j += 1
            while j < len(story) and isinstance(story[j], Spacer):
                j += 1
            if j < len(story):
                h = _altezza(story[j])
                spezzabile = isinstance(story[j], (Table, Paragraph))       # figure e riquadri vanno interi
                serve += min(h, minimo) if spezzabile else min(h, 24 * cm)
            out.append(CondPageBreak(serve))
        out.append(f); k += 1
    return out


# ---------------------------------------------------------------------------------------------- formule
def _formula(righe, scala=1.0):
    chiave = hashlib.md5(("|".join(righe) + str(scala)).encode()).hexdigest()[:16]
    percorso = os.path.join(CACHE, f"eq_{chiave}.png")
    dpi = 300
    if not os.path.exists(percorso):
        plt.rcParams["mathtext.fontset"] = "cm"
        fig = plt.figure(figsize=(0.01, 0.01))
        n = len(righe)
        for i, r in enumerate(righe):
            fig.text(0, -i * 0.42 * scala, f"${r}$" if r.strip() else " ", fontsize=15 * scala, ha="left", va="baseline")
        fig.savefig(percorso, dpi=dpi, bbox_inches="tight", pad_inches=0.03, transparent=True)
        plt.close(fig)
    from reportlab.lib.utils import ImageReader
    w, h = ImageReader(percorso).getSize()
    W, H = w / dpi * 72 * 0.92, h / dpi * 72 * 0.92
    if W > LARGHEZZA - 1.4 * cm:
        f = (LARGHEZZA - 1.4 * cm) / W; W *= f; H *= f
    return percorso, W, H


# ---------------------------------------------------------------------------------------------- codice
def _togli_commento(r):
    """Toglie il commento `#` di una riga (Julia/Python), ignorando i `#` dentro le stringhe tra virgolette."""
    in_str = False; k = 0
    while k < len(r):
        c = r[k]
        if in_str:
            if c == "\\": k += 1
            elif c == '"': in_str = False
        elif c == '"': in_str = True
        elif c == "#": return r[:k].rstrip()
        k += 1
    return r.rstrip()


def _senza_commenti(L):
    """Toglie i commenti (righe intere, in fondo alle righe, blocchi #= =#) e le righe vuote doppie."""
    out = []; blocco = False
    for r in L:
        if blocco:
            if "=#" in r: blocco = False
            continue
        if r.strip().startswith("#="):
            blocco = "=#" not in r.strip()[2:]; continue
        if r.strip().startswith("#"): continue
        r = _togli_commento(r)
        if r.strip() == "" and (not out or out[-1].strip() == ""): continue
        out.append(r)
    while out and out[-1].strip() == "": out.pop()
    return out


def _a_capo(L):
    """Manda a capo le righe più larghe del riquadro, dopo uno spazio o una virgola, con il rientro della riga più 4."""
    st = S["cod"]; utile = LARGHEZZA - 8 - 6 - 4
    ncar = int(utile // pdfmetrics.stringWidth("M", st.fontName, st.fontSize))
    out = []
    for r in L:
        rientro = " " * (len(r) - len(r.lstrip()) + 4)
        while len(r) > ncar:
            taglio = max(r.rfind(" ", 0, ncar), r.rfind(",", 0, ncar - 1) + 1)
            if taglio <= len(rientro): taglio = ncar
            out.append(r[:taglio].rstrip()); r = rientro + r[taglio:].lstrip()
        out.append(r)
    return out


def estrai(file, inizio, fine=None, n=None, togli_commenti=False):
    """Righe di `file` (relativo al progetto) dalla prima che contiene `inizio` fino alla prima successiva
    che contiene `fine` (compresa), oppure `n` righe. `inizio` può essere un numero di riga."""
    with open(os.path.join(PROGETTO, file), encoding="utf-8") as f:
        L = f.read().replace("\r\n", "\n").split("\n")
    if isinstance(inizio, int): a = inizio - 1
    else:
        a = next(i for i, r in enumerate(L) if inizio in r)
    if n is not None: b = a + n
    elif fine is not None:
        b = next(i for i in range(a + 1, len(L)) if fine in L[i]) + 1
    else:
        b = len(L)
    sel = L[a:b]
    if togli_commenti: sel = [r for r in sel if not r.strip().startswith("#")]
    import textwrap
    return textwrap.dedent("\n".join(sel))


# ---------------------------------------------------------------------------------------------- flowchart
STILI_NODO = {
    "inizio": dict(fc="#1f2328", ec="#1f2328", tc="white", box="round,pad=0.02,rounding_size=0.35"),
    "fine":   dict(fc="#1f2328", ec="#1f2328", tc="white", box="round,pad=0.02,rounding_size=0.35"),
    "proc":   dict(fc="#e8f0fb", ec="#2a78d6", tc="#1f2328", box="square,pad=0.02"),
    "funz":   dict(fc="#f1e9f8", ec="#8e5bb5", tc="#1f2328", box="square,pad=0.02"),
    "io":     dict(fc="#e6f6ef", ec="#1baf7a", tc="#1f2328", box=None),
    "dec":    dict(fc="#fdeee6", ec="#eb6834", tc="#1f2328", box=None),
    "nota":   dict(fc="#ffffff", ec="#9aa0a6", tc="#5f6368", box="round,pad=0.02,rounding_size=0.15"),
}


def disegna_flowchart(schema):
    """schema = dict(nodi={id: (x, y, testo, tipo[, larghezza, altezza])}, archi=[(da, a, etichetta, percorso)],
    dimensioni=(L, H) in unità di griglia, scala=cm per unità). percorso: None (retto), 'hv', 'vh',
    oppure lista di punti intermedi [(x, y), ...]; i lati di uscita/entrata si scelgono con da+'.s'/'.n'/'.e'/'.w'."""
    chiave = hashlib.md5(repr(schema).encode()).hexdigest()[:16]
    percorso = os.path.join(CACHE, f"fc_{chiave}.png")
    if os.path.exists(percorso): return percorso
    Lg, Hg = schema["dimensioni"]; u = schema.get("scala", 16.8 / Lg)
    fig = plt.figure(figsize=(Lg * u / 2.54, Hg * u / 2.54))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, Lg); ax.set_ylim(0, Hg); ax.axis("off")
    fs = schema.get("font", 8.2)
    geo = {}
    for k, v in schema["nodi"].items():
        x, y, testo, tipo = v[:4]
        w = v[4] if len(v) > 4 else schema.get("w", 3.4)
        h = v[5] if len(v) > 5 else schema.get("h", 0.9)
        st = STILI_NODO[tipo]; geo[k] = (x, y, w, h, tipo)
        if tipo == "dec":
            ax.add_patch(Polygon([(x - w / 2, y), (x, y + h / 2), (x + w / 2, y), (x, y - h / 2)], closed=True,
                                 fc=st["fc"], ec=st["ec"], lw=1.2, zorder=2))
        elif tipo == "io":
            s = 0.18
            ax.add_patch(Polygon([(x - w / 2 + s, y + h / 2), (x + w / 2 + s, y + h / 2), (x + w / 2 - s, y - h / 2),
                                  (x - w / 2 - s, y - h / 2)], closed=True, fc=st["fc"], ec=st["ec"], lw=1.2, zorder=2))
        else:
            ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=st["box"], fc=st["fc"], ec=st["ec"],
                                        lw=1.2, zorder=2, linestyle="--" if tipo == "nota" else "-"))
            if tipo == "funz":
                for dx in (0.12, w - 0.12):
                    ax.plot([x - w / 2 + dx] * 2, [y - h / 2, y + h / 2], color=st["ec"], lw=0.9, zorder=3)
        peso = "bold" if tipo in ("inizio", "fine") else "normal"
        f0 = fs * (0.92 if tipo == "nota" else 1)
        tx = ax.text(x, y, testo, ha="center", va="center", fontsize=f0, color=st["tc"],
                     zorder=4, fontweight=peso, linespacing=1.25, family="DejaVu Sans")
        ren = fig.canvas.get_renderer()
        disp = ax.transData.transform([(x - w / 2, y - h / 2), (x + w / 2, y + h / 2)])
        bw = (disp[1][0] - disp[0][0]) * (0.62 if tipo == "dec" else 0.93); bh = (disp[1][1] - disp[0][1]) * (0.7 if tipo == "dec" else 0.9)
        for _ in range(12):
            bb = tx.get_window_extent(ren)
            if bb.width <= bw and bb.height <= bh: break
            f0 *= 0.93; tx.set_fontsize(f0)

    def porta(k, lato):
        x, y, w, h, tipo = geo[k]
        return {"n": (x, y + h / 2), "s": (x, y - h / 2), "e": (x + w / 2, y), "w": (x - w / 2, y)}[lato]

    for arco in schema["archi"]:
        da, a = arco[0], arco[1]
        etichetta = arco[2] if len(arco) > 2 else None
        perc = arco[3] if len(arco) > 3 else None
        lda = da.split(".")[1] if "." in da else None; da = da.split(".")[0]
        la = a.split(".")[1] if "." in a else None; a = a.split(".")[0]
        x0, y0 = geo[da][:2]; x1, y1 = geo[a][:2]
        sx = 1 if x1 > x0 else -1; sy = 1 if y1 > y0 else -1
        vert = ("n" if sy > 0 else "s"); vert_in = ("s" if sy > 0 else "n")
        oriz = ("e" if sx > 0 else "w"); oriz_in = ("w" if sx > 0 else "e")
        if perc == "hv": auto = (oriz, vert_in)
        elif perc == "vh": auto = (vert, oriz_in)
        elif abs(x1 - x0) < 1e-6: auto = (vert, vert_in)
        elif abs(y1 - y0) < 1e-6: auto = (oriz, oriz_in)
        else: auto = (vert, vert_in)
        if isinstance(perc, list):
            def lato(p, q):
                dx, dy = q[0] - p[0], q[1] - p[1]
                return ("e" if dx > 0 else "w") if abs(dx) > abs(dy) else ("n" if dy > 0 else "s")
            auto = (lato((x0, y0), perc[0]), lato((x1, y1), perc[-1]))
        lda = lda or auto[0]; la = la or auto[1]
        p0 = porta(da, lda); p1 = porta(a, la)
        if isinstance(perc, list): pts = [p0] + perc + [p1]
        elif perc == "hv": pts = [p0, (p1[0], p0[1]), p1]
        elif perc == "vh": pts = [p0, (p0[0], p1[1]), p1]
        else: pts = [p0, p1]
        xs, ys = zip(*pts)
        ax.plot(xs[:-1] + (xs[-1],), ys, color="#3c4043", lw=1.0, zorder=1, solid_capstyle="round")
        ax.annotate("", xy=pts[-1], xytext=pts[-2], arrowprops=dict(arrowstyle="-|>", color="#3c4043", lw=1.0,
                    shrinkA=0, shrinkB=0, mutation_scale=10), zorder=1)
        if etichetta:
            mx, my = (pts[0][0] + pts[1][0]) / 2, (pts[0][1] + pts[1][1]) / 2
            ax.text(mx + 0.08, my + 0.08, etichetta, fontsize=fs * 0.88, color="#b4371f", ha="left", va="bottom",
                    zorder=5, family="DejaVu Sans", fontweight="bold")
    fig.savefig(percorso, dpi=230); plt.close(fig)
    return percorso


# ---------------------------------------------------------------------------------------------- dati
def csv_numeri(file):
    with open(os.path.join(PROGETTO, file), encoding="utf-8") as f:
        L = [r for r in f.read().splitlines() if r.strip()]
    return L[0].split(","), [[float(v) for v in r.split(",")] for r in L[1:]]


def matrice(file):
    import numpy as np
    return np.loadtxt(os.path.join(PROGETTO, file), ndmin=2)


# ---------------------------------------------------------------------------------------------- dati di esempio
# Gli esempi numerici usano dati sintetici, non quelli di un caso di studio.
Z_FONDO, B_FONDO, SCARPA, H_ARGINE, SCARPA_G, KS_ES, G = 100.0, 15.0, 1.5, 2.5, 8.0, 25.0, 9.81


def tabella_sintetica():
    """Tabella di una sezione composta di esempio, nel formato del modello (12 colonne):
    alveo trapezio (fondo 15 m, scarpa 1,5) alto 2,5 m, poi golene in pendenza dolce (1:8); Ks = 25 costante, β = 1."""
    import numpy as np
    y = np.round(np.arange(0.1, 3.51, 0.1), 2); righe = []
    for h in y:
        hc = min(h, H_ARGINE); B_sup = B_FONDO + 2 * SCARPA * hc
        A = (B_FONDO + SCARPA * hc) * hc; P = B_FONDO + 2 * hc * np.sqrt(1 + SCARPA ** 2); W = B_sup
        if h > H_ARGINE:
            hg = h - H_ARGINE; A += B_sup * hg + SCARPA_G * hg ** 2
            P += 2 * hg * np.sqrt(1 + SCARPA_G ** 2); W = B_sup + 2 * SCARPA_G * hg
        Rh = A / P
        righe.append([1, Z_FONDO + h, A, 0, 1.0, 0, KS_ES ** 2 * Rh ** (4 / 3) / G, 0, W, 0, Rh, 0])
    T = np.array(righe)
    for c in (2, 4, 6, 8, 10):
        s = np.diff(T[:, c]) / np.diff(T[:, 1]); T[:-1, c + 1] = s; T[-1, c + 1] = s[-1]
    return T


def idrogramma_sintetico():
    """Idrogramma di esempio: base 5 m³/s, colmo 80 m³/s a 7 h, valori ogni 900 s per 24 h."""
    import numpy as np
    t = np.arange(0, 24 * 3600 + 1, 900.0); x = (t / 3600 - 7) / 2.2
    return np.c_[t, 5 + 75 * np.exp(-(x + np.exp(-x)) + 1)]
