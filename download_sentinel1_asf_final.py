# download_s1_all_asc_desc_manifest_then_get.py  (avec AOI shrink + frameNumber filter)
# ─────────────────────────────────────────────────────────────────────────────
# Requête Sentinel-1 (ASC/DESC) filtrée sur RON (relativeOrbit) + frameNumber,
# avec réduction automatique de l’AOI (buffer négatif en mètres sur rectangle UTM).
# Manifestes écrits avant téléchargement. Téléchargements séparés en ASC/ et DESC/.
#
# Dépendances : pip install asf-search pyproj tqdm pandas
# Auth : définir ASF_USER / ASF_PASS dans l'environnement.
# ─────────────────────────────────────────────────────────────────────────────

from datetime import datetime, timedelta
from pathlib import Path
import os, sys, math, shutil, logging
import pandas as pd
from tqdm import tqdm
from pyproj import Transformer
from asf_search import search, ASFSession

# ===================== PARAMÈTRES =====================

# Identifiants (via variables d'environnement, recommandé)
USERNAME = os.getenv("ASF_USER", "sig_dev")
PASSWORD = os.getenv("ASF_PASS", "AKzx.y/_P9&cc7m")

# Période autour d'un évènement
EVENT_DATE  = datetime(2020, 8, 8)
WINDOW_DAYS = 60  # ± jours (élargis pour récupérer plus de scènes, y compris DESC)

# Piste(s) relative(s) (RON). Mettre un int (ex. 59), une liste (ex. [59,168]) ou None (toutes).
REL_ORBIT = 59  # ex. pour stack ASC ; pour découvrir DESC, mettre None et FLIGHT_DIRECTIONS=["DESCENDING"]

# Directions (peut être ["ASCENDING"], ["DESCENDING"] ou les deux)
FLIGHT_DIRECTIONS = ["ASCENDING", "DESCENDING"]

# Type produit
PROCESSING_LEVEL = "SLC"  # SLC pour interféro
BEAM_MODE        = "IW"   # IW pour TOPS

# Filtre baseline perpendiculaire absolue (mètres) ; None = désactivé
MAX_ABS_PERP_BASELINE = None  # ex. 150

# --------- Filtre frameNumber ---------
# Mode 1) Liste fixe à garder, ex. [1143, 1144]; None = pas de filtre par frameNumber
FRAMES_KEEP_LIST = None

# Mode 2) Auto : garder les K frames les plus fréquentes (par direction) ; 0 = désactivé
AUTO_FRAME_TOP_K = 1  # ex. 1 = ne garder que la frame dominante par direction

# Dossier de sortie
OUT_DIR = Path(r"D:\sentinel_1")

# DRY-RUN : True = manifestes seulement ; False = télécharge
DRY_RUN = True

# Espace disque minimum (GiB) avant téléchargement
MIN_FREE_GIB = 10

# ===================== AOI RECTANGLE (UTM → WKT) =====================
# Rectangle en UTM32N (EPSG:32632) ; on applique un "shrink" (buffer négatif)
UTM_EPSG, WGS84_EPSG = 32632, 4326

# Rectangle UTM initial (adapter à ta zone)
xmin, ymin, xmax, ymax = 236_953, 4_053_000, 271_903, 4_054_418

# SHRINK en mètres (réduction sur chaque côté). Mettre 0 pour désactiver.
AOI_SHRINK_M = 1500  # ex. rétrécit 1.5 km de chaque côté pour réduire les scènes « bordées »

def shrink_utm_rect(xmin, ymin, xmax, ymax, shrink_m):
    if shrink_m is None or shrink_m <= 0:
        return xmin, ymin, xmax, ymax
    sxmin = xmin + shrink_m
    symin = ymin + shrink_m
    sxmax = xmax - shrink_m
    symax = ymax - shrink_m
    # évite inversion
    if sxmin >= sxmax or symin >= symax:
        return xmin, ymin, xmax, ymax
    return sxmin, symin, sxmax, symax

xmin, ymin, xmax, ymax = shrink_utm_rect(xmin, ymin, xmax, ymax, AOI_SHRINK_M)

# reprojection → WKT
tf = Transformer.from_crs(UTM_EPSG, WGS84_EPSG, always_xy=True)
corners = [tf.transform(x, y) for x, y in [(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)]]
poly_wkt = (
    "POLYGON(("
    + ", ".join(f"{lon:.6f} {lat:.6f}" for lon, lat in corners)
    + f", {corners[0][0]:.6f} {corners[0][1]:.6f}))"
)

# ===================== UTILS =====================

def as_float(x, default=math.nan):
    try:
        return float(x)
    except Exception:
        return default

def get_ron(prop: dict):
    """Récupère la piste relative selon le nom de clé disponible."""
    for k in ["relativeOrbit", "relativeOrbitNumber", "track", "trackNumber", "pathNumber"]:
        v = prop.get(k, None)
        if v is not None and str(v).strip() != "":
            try:
                return int(str(v).strip())
            except Exception:
                pass
    return None

def get_frame(prop: dict):
    """Récupère le frameNumber si dispo (clé 'frameNumber' souvent présente)."""
    for k in ["frameNumber", "frame", "frame_id"]:
        v = prop.get(k, None)
        if v is not None and str(v).strip() != "":
            try:
                return int(str(v).strip())
            except Exception:
                pass
    return None

def has_space(path: Path, min_gib: float) -> bool:
    try:
        free_gib = shutil.disk_usage(path).free / 1024**3
        return free_gib >= float(min_gib)
    except Exception:
        root = path.anchor or "/"
        free_gib = shutil.disk_usage(root).free / 1024**3
        return free_gib >= float(min_gib)

# ===================== LOGGING =====================
logging.getLogger("asf_search").setLevel(logging.INFO)

# ===================== AUTH =====================
if not USERNAME or not PASSWORD:
    print("❌ Définis ASF_USER / ASF_PASS dans l'environnement.")
    sys.exit(1)
session = ASFSession().auth_with_creds(USERNAME, PASSWORD)

# ===================== PÉRIODE =====================
start = (EVENT_DATE - timedelta(days=WINDOW_DAYS)).strftime("%Y-%m-%d")
end   = (EVENT_DATE + timedelta(days=WINDOW_DAYS)).strftime("%Y-%m-%d")
print("AOI WKT   :", poly_wkt)
print(f"Période   : {start} → {end}")
print(f"Directions: {', '.join(FLIGHT_DIRECTIONS)}")
print(f"RON       : {REL_ORBIT if REL_ORBIT is not None else 'ALL'}")
print(f"AOI shrink: {AOI_SHRINK_M} m")

# ===================== RECHERCHE ASF =====================
search_kwargs = dict(
    start           = start,
    end             = end,
    processingLevel = PROCESSING_LEVEL,
    beamMode        = BEAM_MODE,
    intersectsWith  = poly_wkt,
    platform        = ["Sentinel-1A","Sentinel-1B"],
)
if REL_ORBIT is not None:
    # accepter liste ou int
    if isinstance(REL_ORBIT, (list, tuple, set)):
        search_kwargs["relativeOrbit"] = list(REL_ORBIT)
    else:
        search_kwargs["relativeOrbit"] = [int(REL_ORBIT)]

results = search(**search_kwargs)
print(f"→ {len(results)} scènes trouvées (toutes directions)")

if results:
    print("Clés properties (exemple 1er résultat) :", sorted(list(results[0].properties.keys())))

# ===================== TABULARISATION =====================
rows = []
for p in results:
    prop = p.properties
    scene = prop.get("sceneName", "")
    start_t = prop.get("startTime", "")
    direc = (prop.get("flightDirection", "") or "").upper()
    ron   = get_ron(prop)
    frm   = get_frame(prop)

    base = prop.get("perpBaseline", None)
    if base is None:
        base = prop.get("perpendicularBaseline", None)
    base = as_float(base)

    url = prop.get("url", prop.get("downloadUrl", ""))

    try:
        t = datetime.fromisoformat(start_t.replace("Z",""))
    except Exception:
        t = None

    rows.append({
        "sceneName": scene,
        "startTime": start_t,
        "date_obj": t,
        "flightDirection": direc,
        "relativeOrbit": ron,
        "frameNumber": frm,
        "perpBaseline": base,
        "url": url,
        "prod": p
    })

df = pd.DataFrame(rows)
if df.empty:
    print("❌ Aucune scène ne correspond à la requête.")
    sys.exit(0)

# ===================== FILTRES POST-REQUÊTE =====================

# 1) Directions
df = df[df["flightDirection"].isin([d.upper() for d in FLIGHT_DIRECTIONS])].copy()
print(f"→ {len(df)} scènes après filtre direction(s) {FLIGHT_DIRECTIONS}")

# 2) RON (si demandé ET si dispo)
if REL_ORBIT is not None:
    before = len(df)
    # autoriser liste de RON
    if isinstance(REL_ORBIT, (list, tuple, set)):
        keep = set(int(x) for x in REL_ORBIT)
        df = df[df["relativeOrbit"].isin(keep)].copy()
    else:
        df = df[df["relativeOrbit"].notna() & (df["relativeOrbit"].astype(int) == int(REL_ORBIT))].copy()
    print(f"RON demandée = {REL_ORBIT} | gardées = {len(df)}/{before}")
    if df.empty:
        print("⚠ Pas de RON correspondant dans les propriétés renvoyées. "
              "Essaye REL_ORBIT=None pour découvrir les pistes disponibles.")
        sys.exit(0)

# 3) Baseline perpendiculaire (optionnel)
if MAX_ABS_PERP_BASELINE is not None:
    before = len(df)
    df = df[df["perpBaseline"].abs() <= float(MAX_ABS_PERP_BASELINE)].copy()
    print(f"→ {len(df)} scènes après filtre |baseline| ≤ {MAX_ABS_PERP_BASELINE} m (avant={before})")

# 4) Filtre frameNumber
print("Répartition par direction & frameNumber (avant filtre) :")
print(df.groupby(["flightDirection","frameNumber"]).size().sort_values(ascending=False))

if FRAMES_KEEP_LIST:
    before = len(df)
    df = df[df["frameNumber"].isin(set(FRAMES_KEEP_LIST))].copy()
    print(f"→ {len(df)} scènes après filtre frameNumber ∈ {FRAMES_KEEP_LIST} (avant={before})")
elif AUTO_FRAME_TOP_K and AUTO_FRAME_TOP_K > 0:
    before = len(df)
    # par direction : garder les K frames les plus fréquentes
    keep_frames = []
    for direc in df["flightDirection"].dropna().unique():
        sub = df[df["flightDirection"]==direc]
        top = sub["frameNumber"].value_counts().head(AUTO_FRAME_TOP_K).index.tolist()
        keep_frames.extend(top)
    keep_frames = [f for f in keep_frames if pd.notna(f)]
    df = df[df["frameNumber"].isin(keep_frames)].copy()
    print(f"→ {len(df)} scènes après filtre auto top-{AUTO_FRAME_TOP_K} frames par direction ; frames gardées = {sorted(set(keep_frames))}")

print("Répartition par direction & frameNumber (après filtre) :")
print(df.groupby(["flightDirection","frameNumber"]).size().sort_values(ascending=False))

# Tri final : date puis baseline
df = df.sort_values(["date_obj", "perpBaseline"], ascending=[True, True]).reset_index(drop=True)

# ===================== MANIFESTES =====================
OUT_DIR.mkdir(parents=True, exist_ok=True)
manifest_csv = OUT_DIR / "query_results.csv"
manifest_txt = OUT_DIR / "query_results.txt"
gran_all     = OUT_DIR / "granules_found.txt"
gran_asc     = OUT_DIR / "granules_ASC.txt"
gran_desc    = OUT_DIR / "granules_DESC.txt"

df_out = df.drop(columns=["prod"]).copy()
df_out.to_csv(manifest_csv, index=False)

df["sceneName"].to_csv(gran_all, index=False, header=False)
df[df["flightDirection"]=="ASCENDING"]["sceneName"].to_csv(gran_asc, index=False, header=False)
df[df["flightDirection"]=="DESCENDING"]["sceneName"].to_csv(gran_desc, index=False, header=False)

enr = df.assign(
    startTime=df["date_obj"].dt.strftime("%Y-%m-%d %H:%M:%S").fillna(""),
    flightDirection=df["flightDirection"].fillna(""),
    relativeOrbit=df["relativeOrbit"].fillna(""),
    frameNumber=df["frameNumber"].fillna(""),
    perpBaseline=df["perpBaseline"].fillna(""),
    url=df["url"].fillna(""),
)[["sceneName","startTime","flightDirection","relativeOrbit","frameNumber","perpBaseline","url"]]
enr.to_csv(manifest_txt, index=False, sep=";")

print(f"✓ Manifeste CSV        : {manifest_csv}")
print(f"✓ Liste enrichie (TXT) : {manifest_txt}")
print(f"✓ Granules (toutes)    : {gran_all}")
print(f"✓ Granules ASC         : {gran_asc}")
print(f"✓ Granules DESC        : {gran_desc}")

# ===================== TÉLÉCHARGEMENT =====================
if DRY_RUN:
    print("⏸️ DRY_RUN=True : manifestes écrits, aucun téléchargement lancé.")
    sys.exit(0)

dir_asc  = OUT_DIR / "ASC"
dir_desc = OUT_DIR / "DESC"
dir_asc.mkdir(exist_ok=True, parents=True)
dir_desc.mkdir(exist_ok=True, parents=True)

MAX_RETRIES = 3

def download_with_retries(prod, dest_dir: Path):
    scene = prod.properties.get("sceneName", "UNKNOWN")
    for attempt in range(1, MAX_RETRIES+1):
        try:
            if not has_space(dest_dir, MIN_FREE_GIB):
                raise SystemExit(f"⛔ Espace disque < {MIN_FREE_GIB} GiB ({dest_dir}), arrêt.")
            saved = prod.download(path=str(dest_dir), session=ASFSession().auth_with_creds(USERNAME, PASSWORD))
            print(f"✓ {scene} → {saved}")
            return True
        except Exception as e:
            print(f"⚠️ Tentative {attempt}/{MAX_RETRIES} échouée pour {scene} : {e}")
    print(f"❌ Échec pour {scene} après {MAX_RETRIES} tentatives (skip).")
    return False

ok_count = 0
for _, rec in df.iterrows():
    prod = rec["prod"]
    direc = rec["flightDirection"]
    dest_dir = dir_asc if direc == "ASCENDING" else dir_desc
    if download_with_retries(prod, dest_dir):
        ok_count += 1

print(f"✅ Téléchargement terminé. {ok_count}/{len(df)} fichiers récupérés.")
