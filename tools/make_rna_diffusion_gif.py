import math, random
from PIL import Image, ImageDraw, ImageFont

# ----------------------------------------------------------------- canvas
W = H = 420
SS = 2
CW, CH = W*SS, H*SS
BG = (252, 250, 247)
INK = (28, 33, 40)
INK3 = (140, 134, 128)
ACC2 = (249, 115, 22)

FONT = "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"
f_cap = ImageFont.truetype(FONT, 15*SS)
f_tag = ImageFont.truetype(FONT, 12*SS)

def mix(a, b, t):
    t = max(0.0, min(1.0, t))
    return tuple(int(round(x + (y-x)*t)) for x, y in zip(a, b))

def spectrum(t):
    """N->C rainbow, the usual molecular-viewer ramp."""
    stops = [(0.00, (37, 99, 235)), (0.28, (14, 165, 233)), (0.50, (34, 197, 94)),
             (0.74, (234, 179, 8)), (1.00, (220, 38, 38))]
    for i in range(len(stops)-1):
        t0, c0 = stops[i]; t1, c1 = stops[i+1]
        if t <= t1:
            return mix(c0, c1, (t-t0)/(t1-t0))
    return stops[-1][1]

# ----------------------------------------------------------------- geometry
def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def mul(a, s): return (a[0]*s, a[1]*s, a[2]*s)
def norm(a):
    n = math.sqrt(sum(c*c for c in a)) or 1.0
    return (a[0]/n, a[1]/n, a[2]/n)
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def frame_vectors(u):
    ref = (0, 0, 1) if abs(u[2]) < 0.9 else (1, 0, 0)
    v = norm(cross(u, ref))
    w = norm(cross(u, v))
    return v, w

def hairpin(center, axis, n=13, rad=5.6, rise=2.45, a0=0.0, phase=2.9):
    """A-form-like duplex: strand up, tetraloop, strand back down."""
    u = norm(axis); v, w = frame_vectors(u)
    pts = []
    for i in range(n):
        a = a0 + i*(2*math.pi/7.5)
        pts.append(add(add(center, mul(u, (i-(n-1)/2)*rise)),
                       add(mul(v, rad*math.cos(a)), mul(w, rad*math.sin(a)))))
    # loop over the top
    top = add(center, mul(u, ((n-1)/2 + 1.6)*rise))
    for j in range(1, 5):
        a = a0 + (n-1)*(2*math.pi/7.5) + j*(phase/4.0)
        bulge = math.sin(math.pi*j/5.0)*3.2
        pts.append(add(add(top, mul(u, bulge)),
                       add(mul(v, (rad+1.5)*math.cos(a)), mul(w, (rad+1.5)*math.sin(a)))))
    for i in range(n-1, -1, -1):
        a = a0 + i*(2*math.pi/7.5) + phase
        pts.append(add(add(center, mul(u, (i-(n-1)/2)*rise)),
                       add(mul(v, rad*math.cos(a)), mul(w, rad*math.sin(a)))))
    return pts

# ----------------------------------------------------------------- target: real 1DDY backbone
# Vitamin B12 RNA aptamer (Sussman, Nix & Wilson 2000), PDB 1DDY, chain A.
# 35 backbone phosphorus atoms; the cobalamin cobalt marks the binding site.
P_ATOMS = [
 (52.822,-5.547,20.195),(51.967,-1.500,16.342),(49.511,-0.527,10.970),(46.233,-1.898,6.611),
 (42.880,-6.228,4.044),(39.338,-12.007,2.842),(37.146,-17.581,2.690),(33.619,-20.702,5.329),
 (29.546,-20.849,9.471),(25.537,-18.479,12.214),(20.988,-16.932,12.628),(15.824,-16.368,12.569),
 (12.431,-15.528,9.073),(11.009,-12.275,5.468),(17.031,-14.542,6.880),(22.802,-15.544,4.537),
 (28.587,-12.421,4.111),(33.440,-9.623,3.981),(36.568,-6.354,7.233),(37.444,-5.515,12.869),
 (37.142,-8.231,18.125),(37.395,-13.848,20.404),(36.332,-20.573,20.807),(30.927,-23.364,23.817),
 (28.014,-18.178,21.465),(29.015,-14.153,19.571),(24.825,-10.191,20.358),(24.067,-3.566,18.946),
 (24.817,0.454,15.765),(22.947,1.750,10.303),(20.974,0.851,16.613),(14.367,1.534,18.575),
 (8.562,-1.613,18.125),(4.509,-5.287,17.856),(2.801,-10.487,19.705),
]
from bases import BASES_RAW, RING_ATOMS, SUGARS   # backbone, base rings and riboses, same entry

# residues 1-4 are a dangling, unpaired 5' arm in the crystal — trim them
TRIM = 4
P_ATOMS = P_ATOMS[TRIM:]
BASES_RAW = BASES_RAW[TRIM:]

SITE = (18.511,-5.042,3.710)          # cobalamin Co — the ligand pocket
CTRL = [sub(q, SITE) for q in P_ATOMS]
_C = tuple(sum(q[d] for q in CTRL)/len(CTRL) for d in range(3))
CTRL = [sub(q, _C) for q in CTRL]
LIG_OFFSET = mul(_C, -1.0)                       # ligand keeps its real place in the pocket
FIT = max(math.sqrt(q[0]**2 + q[1]**2 + q[2]**2) for q in CTRL)

# nucleobases, expressed relative to their own phosphate so they ride the chain
BASE_COL = {'A': (47,158,94), 'C': (42,103,207), 'G': (176,124,16), 'U': (207,68,55)}
BASE_REL = []
for idx, (num, name, c1, c4) in enumerate(BASES_RAW):
    anchor = P_ATOMS[idx]
    BASE_REL.append((name, sub(c1, anchor), sub(c4, anchor)))

PUR = [('N9','C8'),('C8','N7'),('N7','C5'),('C5','C4'),('C4','N9'),
       ('C4','N3'),('N3','C2'),('C2','N1'),('N1','C6'),('C6','C5')]
PYR = [('N1','C2'),('C2','N3'),('N3','C4'),('C4','C5'),('C5','C6'),('C6','N1')]
EXTRA = {'A': [('C6','N6'),("C1'",'N9')],
         'G': [('C6','O6'),('C2','N2'),("C1'",'N9')],
         'C': [('C2','O2'),('C4','N4'),("C1'",'N1')],
         'U': [('C2','O2'),('C4','O4'),("C1'",'N1')]}
ELEM = {'N': (37,99,235), 'O': (220,38,38)}
SUG_C = (122, 112, 104)
SUG_BONDS = [("C1'","C2'"),("C2'","C3'"),("C3'","C4'"),("C4'","O4'"),
             ("O4'","C1'"),("C4'","C5'"),("C2'","O2'")]

FULL_SUGARS = []
for idx in range(len(P_ATOMS)):
    at = SUGARS.get(idx+1+TRIM)
    if not at:
        continue
    anchor = P_ATOMS[idx]
    rel = {k: sub(v, anchor) for k, v in at.items()}
    FULL_SUGARS.append((idx, rel))

# every base atom, expressed relative to its own phosphate
FULL_BASES = []
for idx, (num, name, c1, c4) in enumerate(BASES_RAW):
    if (idx+1+TRIM) not in RING_ATOMS:
        continue
    rname, atoms = RING_ATOMS[idx+1+TRIM]
    anchor = P_ATOMS[idx]
    rel = {k: sub(v, anchor) for k, v in atoms.items()}
    bonds = (PUR if rname in 'AG' else PYR) + EXTRA[rname]
    bonds = [(a, b) for a, b in bonds if a in rel and b in rel]
    FULL_BASES.append((idx, rname, rel, bonds))

def catmull(pts, per=9):
    out = []
    p = [pts[0]] + list(pts) + [pts[-1], pts[-1]]
    for i in range(len(p)-3):
        p0, p1, p2, p3 = p[i], p[i+1], p[i+2], p[i+3]
        for k in range(per):
            t = k/per
            t2, t3 = t*t, t*t*t
            out.append(tuple(
                0.5*((2*p1[d]) + (-p0[d]+p2[d])*t +
                     (2*p0[d]-5*p1[d]+4*p2[d]-p3[d])*t2 +
                     (-p0[d]+3*p1[d]-3*p2[d]+p3[d])*t3) for d in range(3)))
    return out

TARGET = catmull(CTRL, per=9)
NT = len(TARGET)

# random-coil start: a self-avoiding-ish brownian chain of the same length
random.seed(11)
def brownian(n, step=2.7):
    p = (0.0, 0.0, 0.0); d = norm((random.gauss(0,1), random.gauss(0,1), random.gauss(0,1)))
    out = [p]
    for _ in range(n-1):
        d = norm(add(mul(d, 0.72), mul((random.gauss(0,1), random.gauss(0,1), random.gauss(0,1)), 0.9)))
        back = mul(norm(p) if any(p) else (0,0,0), -0.055*math.sqrt(sum(c*c for c in p)))
        d = norm(add(d, back))
        p = add(p, mul(d, step)); out.append(p)
    c = tuple(sum(q[d]/n for q in out) for d in range(3))
    return [sub(q, c) for q in out]

START = brownian(NT)

# ----------------------------------------------------------------- ligand (psilocybin, ~planar indole)
S2 = 0.072
def L(x, y, z=0.0): return (x*S2*20/1.4*0.07*20, y*S2*20/1.4*0.07*20, z)   # keep 2D layout, scaled
def P2(x, y, z=0.0): return (x*0.0725, -y*0.0725, z)

RING = {
 'C4':P2(0,-20), 'C5':P2(-17.3,-10), 'C6':P2(-17.3,10), 'C7':P2(0,20),
 'C7a':P2(17.3,10), 'C3a':P2(17.3,-10), 'C3':P2(36.3,-16.2), 'C2':P2(48,0), 'N1':P2(36.3,16.2),
 'O4':P2(0,-42), 'P':P2(18,-54), 'Op':P2(18,-74,1.2), 'Oa':P2(40,-62,-1.4), 'Ob':P2(-2,-66,1.0),
 'Cb':P2(56,-26,0.9), 'Cc':P2(76,-16,-0.9), 'N':P2(96,-26,0.6),
 'Cm1':P2(116,-16,1.6), 'Cm2':P2(96,-48,-0.8),
}
C_C, C_N, C_O, C_P = (58, 62, 68), (37, 99, 235), (220, 38, 38), (234, 140, 8)
ATOM_COL = {k: (C_O if k.startswith('O') else C_N if k.startswith('N') else
                C_P if k == 'P' else C_C) for k in RING}
LBONDS = [('C4','C5'),('C5','C6'),('C6','C7'),('C7','C7a'),('C7a','C3a'),('C3a','C4'),
          ('C3a','C3'),('C3','C2'),('C2','N1'),('N1','C7a'),
          ('C4','O4'),('O4','P'),('P','Op'),('P','Oa'),('P','Ob'),
          ('C3','Cb'),('Cb','Cc'),('Cc','N'),('N','Cm1'),('N','Cm2')]
LIG_CENTER = tuple(sum(RING[k][d] for k in RING)/len(RING) for d in range(3))
RING = {k: sub(v, LIG_CENTER) for k, v in RING.items()}

def _rot(p, rx, ry):
    x, y, z = p
    cy_, sy_ = math.cos(ry), math.sin(ry); x, z = x*cy_ - z*sy_, x*sy_ + z*cy_
    cx_, sx_ = math.cos(rx), math.sin(rx); y, z = y*cx_ - z*sx_, y*sx_ + z*cx_
    return (x, y, z)
RING = {k: add(_rot(v, math.radians(38), math.radians(-25)), LIG_OFFSET) for k, v in RING.items()}

# ----------------------------------------------------------------- render
SCALE = (0.40*CW)/FIT
CX, CY = CW/2, CH/2 - 6*SS
TILT = math.radians(12)

def project(p, ang):
    ca, sa = math.cos(ang), math.sin(ang)
    x, y, z = p
    x, z = x*ca - z*sa, x*sa + z*ca
    ct, st = math.cos(TILT), math.sin(TILT)
    y, z = y*ct - z*st, y*st + z*ct
    return CX + x*SCALE, CY - y*SCALE, z

def fog(col, z):
    t = max(0.0, min(1.0, (z + 34) / 68))       # far = 0, near = 1
    shaded = mix(tuple(int(c*0.55) for c in col), col, t)
    return mix(BG, shaded, 0.42 + 0.58*t)

def smoothstep(x):
    x = max(0.0, min(1.0, x)); return x*x*(3-2*x)

PROBE = ['N1', 'O4', 'P', 'Op', 'Oa', 'N', 'C2', 'C7', 'Cm1']
CONTACT_CUT = 11.5

DENOISE, SPIN = 62, 26
TOTAL = DENOISE + SPIN
frames = []

for k in range(TOTAL):
    s = smoothstep(min(1.0, k/(DENOISE-1)))
    noise = (1-s)**1.35
    atom_a = smoothstep((s - 0.66)/0.34)
    ang = 2*math.pi*k/TOTAL

    img = Image.new("RGB", (CW, CH), BG)
    d = ImageDraw.Draw(img)
    discs = []

    # chain
    chain = []
    for i in range(NT):
        t0, s0 = TARGET[i], START[i]
        wob = math.sin(k*0.5 + i*0.35)*1.1*noise
        p = tuple(t0[q]*s + s0[q]*(1-s) + wob*(0.6 if q == 1 else 1.0) for q in range(3))
        chain.append(p)

    r_tube = (3.1 + 1.6*s)*SS * (1 - 0.42*atom_a)
    for i in range(NT-1):
        a, b = chain[i], chain[i+1]
        for f in (0.0, 0.34, 0.67):
            p = tuple(a[q] + (b[q]-a[q])*f for q in range(3))
            x, y, z = project(p, ang)
            col = spectrum(i/(NT-1))
            col = mix(mix(INK3, col, 0.35 + 0.65*s), col, 0)
            discs.append((z, x, y, r_tube, fog(col, z)))

    # ligand — fixed motif, always drawn
    for a, b in LBONDS:
        pa, pb = RING[a], RING[b]
        for f in [q/7 for q in range(8)]:
            p = tuple(pa[q] + (pb[q]-pa[q])*f for q in range(3))
            x, y, z = project(p, ang)
            col = ATOM_COL[a] if f < 0.5 else ATOM_COL[b]
            discs.append((z, x, y, 2.8*SS, fog(col, z)))
    for key, p in RING.items():
        x, y, z = project(p, ang)
        r = 4.4*SS if key in ('P',) else (3.9*SS if key[0] in 'NO' else 3.2*SS)
        discs.append((z, x, y, r, fog(ATOM_COL[key], z)))

    # ---- nucleobases riding the backbone: sticks early, all-atom rings at the end
    stick_a = 1.0 - atom_a
    base_pts = []
    for idx, (name, d1, d4) in enumerate(BASE_REL):
        anchor = chain[min(idx*9, NT-1)]
        wob = math.sin(k*0.5 + idx*1.1)*1.4*noise
        c1 = tuple(anchor[q] + d1[q]*s + wob for q in range(3))
        c4 = tuple(anchor[q] + (d1[q] + (d4[q]-d1[q])*1.55)*s + wob for q in range(3))
        base_pts.append((name, c1, c4))
        if stick_a < 0.04:
            continue
        col = BASE_COL[name]
        col = mix(BG, mix(INK3, col, 0.30 + 0.70*s), stick_a)
        for f in [q/6 for q in range(7)]:
            pp = tuple(c1[t] + (c4[t]-c1[t])*f for t in range(3))
            x, y, z = project(pp, ang)
            r = (2.0 + 1.3*f + 1.0*s)*SS
            discs.append((z, x, y, r, fog(col, z)))

    if atom_a > 0.04:
        for idx, rname, rel, bonds in FULL_BASES:
            anchor = chain[min(idx*9, NT-1)]
            wob = math.sin(k*0.5 + idx*1.1)*1.4*noise
            pos = {kk: tuple(anchor[q] + vv[q]*s + wob for q in range(3)) for kk, vv in rel.items()}
            carbon = BASE_COL[rname]
            for a, b in bonds:
                pa, pb = pos[a], pos[b]
                for f in [q/4 for q in range(5)]:
                    pp = tuple(pa[t] + (pb[t]-pa[t])*f for t in range(3))
                    x, y, z = project(pp, ang)
                    el = (a if f < 0.5 else b)[0]
                    col = ELEM.get(el, carbon)
                    discs.append((z, x, y, 1.55*SS, mix(BG, fog(col, z), atom_a)))
            for kk, pp in pos.items():
                x, y, z = project(pp, ang)
                el = kk[0]
                col = ELEM.get(el, carbon)
                r = (2.3 if el in 'NO' else 1.9)*SS
                discs.append((z, x, y, r, mix(BG, fog(col, z), atom_a)))

    if atom_a > 0.04:
        for idx, rel in FULL_SUGARS:
            anchor = chain[min(idx*9, NT-1)]
            wob = math.sin(k*0.5 + idx*1.1)*1.4*noise
            pos = {kk: tuple(anchor[q] + vv[q]*s + wob for q in range(3)) for kk, vv in rel.items()}
            pos['P'] = anchor
            for a, b in SUG_BONDS + [("C5'", 'P')]:
                pa, pb = pos[a], pos[b]
                for f in [q/4 for q in range(5)]:
                    pp = tuple(pa[t] + (pb[t]-pa[t])*f for t in range(3))
                    x, y, z = project(pp, ang)
                    nm = (a if f < 0.5 else b)
                    col = ELEM['O'] if nm[0] == 'O' else (C_P if nm == 'P' else SUG_C)
                    discs.append((z, x, y, 1.5*SS, mix(BG, fog(col, z), atom_a)))
            for kk, pp in pos.items():
                x, y, z = project(pp, ang)
                col = ELEM['O'] if kk[0] == 'O' else (C_P if kk == 'P' else SUG_C)
                r = (2.9 if kk == 'P' else 2.2 if kk[0] == 'O' else 1.85)*SS
                discs.append((z, x, y, r, mix(BG, fog(col, z), atom_a)))

    # ---- ligand–RNA contacts the network is attending to
    cands = []
    for key in PROBE:
        lp = RING[key]
        for i in range(0, NT, 5):
            cp = chain[i]
            dd = math.dist(lp, cp)
            if dd < CONTACT_CUT:
                cands.append((dd, key, i, lp, cp))
    cands.sort(key=lambda q: q[0])
    for rank, (dd, key, i, lp, cp) in enumerate(cands[:16]):
        phase = (hash((key, i)) % 100) / 100.0 * 6.283
        pulse = 0.5 + 0.5*math.sin(k*0.62 + phase)
        near = 1.0 - dd/CONTACT_CUT
        a = pulse * near * (0.30 + 0.70*s)
        if a < 0.12:
            continue
        col = mix(BG, ACC2, min(1.0, a*1.25))
        steps = 13
        for q in range(steps):
            if q % 2:
                continue
            f = q/(steps-1)
            pp = tuple(lp[t] + (cp[t]-lp[t])*f for t in range(3))
            x, y, z = project(pp, ang)
            discs.append((z + 0.4, x, y, 1.5*SS, col))

    discs.sort(key=lambda q: q[0])
    for z, x, y, r, col in discs:
        d.ellipse([x-r, y-r, x+r, y+r], fill=col)

    # chrome
    lab = "t = T   random chain" if s < 0.02 else ("t = 0   all-atom model" if s > 0.99 else "denoising")
    d.text((18*SS, CH-30*SS), lab, font=f_cap, fill=mix(BG, INK, 0.8))
    bx0, bx1, by = 18*SS, CW-18*SS, CH-40*SS
    d.line([(bx0, by), (bx1, by)], fill=mix(BG, INK3, 0.4), width=int(2*SS))
    d.line([(bx0, by), (bx0 + (bx1-bx0)*s, by)], fill=ACC2, width=int(4*SS))

    frames.append(img.resize((W, H), Image.LANCZOS).quantize(colors=96, method=Image.MEDIANCUT, dither=Image.NONE))

dur = [70]*DENOISE + [80]*SPIN
frames[0].save("rna-diffusion-psilocybin.gif", save_all=True, append_images=frames[1:],
               duration=dur, loop=0, optimize=True, disposal=2)
print("frames", len(frames))
