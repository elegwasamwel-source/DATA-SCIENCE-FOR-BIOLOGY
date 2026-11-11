# Build group_project/study-b-deck.pptx from the two rendered reports.
# Run: python3 group_project/make_deck.py   (then re-run after any re-render)
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pathlib import Path

HERE = Path(__file__).parent
FIGS = HERE / "figures"

prs = Presentation()
prs.slide_width = Inches(13.333)   # 16:9
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]

def add_slide(title, bullets=None, image=None, note=None):
    s = prs.slides.add_slide(BLANK)
    # title
    tb = s.shapes.add_textbox(Inches(0.6), Inches(0.35), Inches(12.1), Inches(0.9))
    p = tb.text_frame.paragraphs[0]
    p.text = title
    p.font.size = Pt(32); p.font.bold = True
    # image (top-right area), full width below title if no bullets
    if image:
        img = FIGS / image
        if bullets:
            left, top, w, h = Inches(6.9), Inches(1.4), Inches(6.0), Inches(5.4)
        else:
            left, top, w, h = Inches(1.2), Inches(1.4), Inches(10.9), Inches(5.7)
        s.shapes.add_picture(str(img), left, top, width=w)  # height scales
    # bullets
    if bullets:
        if image:
            bx, bw = Inches(0.6), Inches(6.1)
        else:
            bx, bw = Inches(1.0), Inches(11.3)
        body = s.shapes.add_textbox(bx, Inches(1.5), bw, Inches(5.4))
        tf = body.text_frame; tf.word_wrap = True
        for i, (txt, lvl) in enumerate(bullets):
            para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            para.text = txt; para.level = lvl
            para.font.size = Pt(20 if lvl == 0 else 17)
    if note:
        s.notes_slide.notes_text_frame.text = note
    return s

def caption_slide(figure, caption, bullets=None):
    return add_slide(f"Figure: {figure}", bullets=bullets, image=figure, note=caption)

# 1 — Title
add_slide("Do hormonal contraceptives differ in their effect?",
    [("Study B — the contraceptive comparison (Group B?, Wednesday–Thursday)", 0),
     ("Pair 1: community composition · Pair 2: host genes", 1),
     ("Durban Data Science for Biology, workshop 2", 0)],
    note="First slide for the Friday audience: they have not seen this study.")

# 2 — The study, in one slide (audience hasn't seen it)
add_slide("The study, in one slide",
    [("93 women, sampled once each (cross-sectional)", 0),
     ("They chose their own contraceptive method — nobody was randomised", 0),
     ("Three arms: 39 injectable · 28 pill · 26 vaginal ring", 0),
     ("Recorded for each woman: BV status (Amsel), Nugent score, pH, STI, sequencing batch", 0),
     ("Two count tables: 220 bacterial species and 13,713 host genes", 1),
     ("Because it is observational, every difference is an association, not proof of cause", 0)],
    note="Read this out before any result — the audience needs the design.")

# 3 — The question, one sentence
add_slide("Our question, in one sentence",
    [("Do the three contraceptive methods differ in vaginal bacterial community "
      "composition or in host gene expression — or is bacterial vaginosis doing "
      "the explaining?", 0)],
    note="This is the whole talk in one sentence (check-in T2 item 1).")

# 4 — What we did
add_slide("What we did",
    [("QC first: all sample IDs match across tables, order verified for DESeq2, "
      "library sizes checked (smallest 5,247 reads — kept, flagged)", 0),
     ("Group QC decision: kept all 93 samples; species kept only if seen in ≥10 of 93 "
      "samples (62 species); genes kept if ≥10 reads in ≥24 samples (13,713 genes)", 1),
     ("Raw counts into every test — DESeq2 does its own normalisation (poscounts)", 1),
     ("Composition: Bray-Curtis PCoA + PERMANOVA, adjusted for BV/pH/STI; "
      "dispersion checked; CLR+Euclidean as a second opinion", 0),
     ("Host genes: DESeq2, ring vs systemic pooled (26 vs 67), adjusted for BV/pH/STI", 0)],
    note="Deliverable item 2: what we did, including filters and why.")

# 5 — QC decision slide (the one from check-in W2)
add_slide("Our QC decision (agreed Wednesday)",
    [("Drop nothing: the smallest library (B006, 5,247 reads) stays, pending checks "
      "that did not flag it", 0),
     ("Filter species at prevalence ≥10/93 samples before any community analysis", 0),
     ("30 species seen in exactly one sample cannot support any comparison", 1),
     ("Filter genes at ≥10 reads in ≥24 samples before DESeq2", 0),
     ("Batch was cross-tabulated against arm and time point: balanced, no confound found", 0)],
    note="This is the slide the whole group agreed at W2.")

# 6 — Figure 1: PCoA
caption_slide("pcoa-by-arm.png",
    "Each point is one woman's vaginal community, 93 samples, ordinated on "
    "Bray-Curtis distance of relative abundance (species in ≥10 samples). The "
    "ellipses by method overlap almost completely — and colouring the same plot "
    "by BV status separates the points far better than method does.")
add_slide("Choice we made: the distance",
    [("Two routes: Bray-Curtis on relative abundance, or CLR + Euclidean (Aitchison), "
      "which treats the counts as compositional from the start", 0),
     ("We kept BOTH and showed the result does not depend on it — the two ordinations "
      "and PERMANOVAs agree", 0),
     ("Source we opened and read: Anderson & Walsh (2013) on PERMANOVA and dispersion", 1)],
    note="Deliverable item 4: 'choice we made' slide for the composition sub-question.")

# 7 — Figure 2: volcano
caption_slide("volcano-ring-vs-systemic.png",
    "DESeq2 on raw counts, ring (26) vs systemic methods pooled (67), adjusted for "
    "BV, pH and STI. 116 of 13,713 genes at adjusted p < 0.05 (72 higher in ring, "
    "44 lower). BPIFA1, a mucosal innate-immunity protein, is ~64-fold higher in "
    "the ring. Boxplots below the volcano show the top genes raw, not just claimed.")
add_slide("Choice we made: adjust or not",
    [("Two routes: test ring vs systemic bare (~ systemic), or adjust for "
      "BV/pH/STI because nobody was randomised", 0),
     ("Bare design calls 227 genes; the adjusted one 116 — about half the "
      "unadjusted signal was the covariates talking", 0),
     ("We kept the adjusted run and kept the bare one side by side so the cost "
      "of adjustment is visible; STI especially is lopsided (17/26 ring vs 31/67)", 1),
     ("Source we opened and read: DESeq2 paper / vignette on design formulas", 1)],
    note="Deliverable item 4: 'choice we made' slide for the host sub-question.")

# 8 — Answer
add_slide("The answer, stated as plainly as the data allow",
    [("The three methods do not detectably differ in vaginal community composition "
      "(arm p = 0.075, ring vs systemic p = 0.85; Shannon p = 0.47)", 0),
     ("Bacterial vaginosis dominates the community picture: R² = 14%, p = 0.001 — "
      "five times the method's share", 0),
     ("The ring does differ from the systemic methods in host gene expression: "
      "116 genes, led by a large BPIFA1 increase", 0),
     ("Whether that gene signal is the method, or BV acting through the host, is "
      "exactly what we have not separated yet", 1)],
    note="Deliverable item 5: the answer, plainly. The last bullet is the honest hinge.")

# 9 — What we do not believe (required slide)
add_slide("What we do NOT believe about our own result",
    [("No randomisation: women chose their methods, so BV, STI and pH differences "
      "across arms may predate any contraceptive effect", 0),
     ("Power: arms of 26–39. We can rule out a LARGE method effect on composition; "
      "a small one could hide inside p = 0.075. Sample counts are on every comparison slide.", 0),
     ("BPIFA1's 64-fold change rests on 26 ring samples, and DESeq2 flagged 6 "
      "outlier genes overall — we show the boxplot, but it is one gene", 0),
     ("'The ring touches the tissue' is a story we brought to the data, not one "
      "the data proved — enrichment (Thursday) is what tests it", 0),
     ("Both count tables are compositional: 'relative abundance went down' is not "
      "'there was less of it'", 0)],
    note="Deliverable item 6, REQUIRED. Read this one out loud.")

# 10 — Open questions
add_slide("Open questions from our pair",
    [("Is the ring's gene signal BV-mediated? (community-type sub-question)", 0),
     ("Which processes do the ring-vs-systemic genes point to? (enrichment, Thursday)", 0),
     ("Would paired sampling over time change the composition answer? This study "
      "cannot tell us — cross-sectional by design", 1)],
    note="The open-questions slide, updated as we went.")

out = HERE / "study-b-deck.pptx"
prs.save(str(out))
print(f"Saved {out} with {len(prs.slides.__iter__.__self__._sldIdLst)} slides")
