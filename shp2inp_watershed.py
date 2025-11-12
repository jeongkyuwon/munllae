# 수정일자: 2025-11-09
# 작성자: 정규원

# 기존의 shp2inp.py 파일은 유역 shp파일변환을 포함하지 않음
# 이 코드는 shp2inp.py 코드에 유역도 같이 변환함.

import geopandas as gpd
import pandas as pd
import numpy as np
import re
from pathlib import Path

# === 입력 파일 ===
LINK_SHP = Path("newshp/m_link_fixed.shp")
MH_SHP   = Path("newshp/m_manhole_fixed.shp")
OF_SHP   = Path("newshp/토구_fixed.shp")
CATCH_SHP= Path(r"D:\Kyuwon\mullae_nbs for sharing\m_watershed.shp")  # ★ 추가

# === 출력 파일 ===
OUT_INP  = Path("newshp/munrae_model.inp")

# === 필드 매핑 ===
F_US_NODE = "SAT_MHE"
F_DS_NODE = "END_MHE"
F_US_INV  = "SPT_HSL"
F_DS_INV  = "EPT_HSL"
F_LEN     = "Shape_Leng"
F_D_US    = "ST_PIP_HIT"
F_D_DS    = "ET_PIP_HIT"

# 사각/단면 문자열 후보
F_SHAPE_TXT = ["DR_LP_CCE", "DTR_CCE", "DP_PE_CCE"]

# 맨홀/토구
F_MH_ID   = "MHE_IDN"
F_MH_INV  = "HSL"
F_OF_ID   = "IDN"
F_OF_INV  = "ATTE"

# 유역 후보 필드(있으면 사용, 없으면 추정/기본값)
C_ID_CANDS     = ["WS_ID", "BASIN_ID", "CAT_ID", "ID", "NAME", "FID", "OBJECTID"]
C_SLOPE_CANDS  = ["SLOPE", "DRN_SLP", "SLOP", "SLP"]     # 단위: m/m 로 가정(퍼센트면 0~100 값 들어올 수 있음 → 보정)
C_PIMP_CANDS   = ["PIMP", "%IMP", "IMP_PCT", "IMPERV"]   # 0~100(%)
C_WIDTH_CANDS  = ["WIDTH", "WID"]                        # m
RAIN_GAGE_NAME = "RG_MUNRAE"

# === 상수 ===
DEFAULT_DIAM_M = 0.6
ROUGH_N        = 0.013
DEFAULT_MAXDEPTH_MH = 6.0
DEFAULT_SLOPE  = 0.01       # m/m
DEFAULT_PIMP   = 70.0       # %
DEFAULT_WIDTH  = None       # 없으면 면적으로부터 추정
CURB_LEN       = 0.0

# === 데이터 읽기 ===
link = gpd.read_file(LINK_SHP)
mh   = gpd.read_file(MH_SHP)
of   = gpd.read_file(OF_SHP)
# ★ 유역
catch = gpd.read_file(CATCH_SHP)

# 좌표계 정합(모두 링크 좌표계로)
if catch.crs != link.crs:
    catch = catch.to_crs(link.crs)
if mh.crs != link.crs:
    mh = mh.to_crs(link.crs)
if of.crs != link.crs:
    of = of.to_crs(link.crs)

# === 기본 유틸 ===
def _clean_id(x):
    s = str(x).strip()
    s = re.sub(r"[^A-Za-z0-9_\-:]", "_", s)
    return s

def _to_float(x, default=np.nan):
    try:
        return float(x)
    except:
        return default
# 수치화
if F_LEN in link.columns:
    link[F_LEN] = pd.to_numeric(link[F_LEN], errors="coerce")
for c in [F_D_US, F_D_DS]:
    if c in link.columns:
        link[c] = pd.to_numeric(link[c], errors="coerce")

if F_MH_INV in mh.columns:
    mh[F_MH_INV] = pd.to_numeric(mh[F_MH_INV], errors="coerce")
if F_OF_INV in of.columns:
    of[F_OF_INV] = pd.to_numeric(of[F_OF_INV], errors="coerce")

# ID 정리
mh[F_MH_ID] = mh[F_MH_ID].apply(_clean_id)
of[F_OF_ID] = of[F_OF_ID].apply(_clean_id)
link[F_US_NODE] = link[F_US_NODE].apply(_clean_id)
link[F_DS_NODE] = link[F_DS_NODE].apply(_clean_id)

# 중복 제거(노드/토구 자체의 ID 중복 방지)
mh = mh.drop_duplicates(subset=[F_MH_ID], keep='first').reset_index(drop=True)
of = of.drop_duplicates(subset=[F_OF_ID], keep='first').reset_index(drop=True)

# 맨홀-토구 ID 충돌 해결(겹치면 맨홀에 _MH)
mh_ids = set(mh[F_MH_ID])
of_ids = set(of[F_OF_ID])
conflict = mh_ids.intersection(of_ids)
if conflict:
    rename_map = {old_id: f"{old_id}_MH" for old_id in conflict}
    mh[F_MH_ID] = mh[F_MH_ID].replace(rename_map)
    link[F_US_NODE] = link[F_US_NODE].replace(rename_map)
    link[F_DS_NODE] = link[F_DS_NODE].replace(rename_map)

# 링크가 참조하는 노드 검증(경고만)
valid_nodes = set(mh[F_MH_ID]).union(set(of[F_OF_ID]))
missing_us = sorted(set(link[F_US_NODE]) - valid_nodes)
missing_ds = sorted(set(link[F_DS_NODE]) - valid_nodes)
if missing_us or missing_ds:
    print("⚠️ 링크가 가리키지만 노드에 없는 ID가 있습니다.")
    if missing_us: print("  - US missing:", missing_us[:8], "..." if len(missing_us) > 8 else "")
    if missing_ds: print("  - DS missing:", missing_ds[:8], "..." if len(missing_ds) > 8 else "")

# === 링크 이름 생성(한 번만 생성해서 칼럼에 저장) ===
LINK_ID_CANDIDATES = ["LINK_ID", "FID", "OBJECTID", "OBJECTID_1", "INX_NO", "ID"]

def _base_link_id_from_row(r):
    # 우선순위 칼럼 중 첫 유효값
    for col in LINK_ID_CANDIDATES:
        if col in r.index and pd.notna(r[col]) and str(r[col]).strip():
            raw = str(r[col]).strip()
            raw = re.sub(r"[^A-Za-z0-9_\-:]", "_", raw)
            # 가독성을 위해 접두사 C 보장
            if not re.match(r"^[A-Za-z]", raw):
                raw = "C" + raw
            return raw
    return None

# 베이스ID 만들기(없으면 일련번호)
link = link.reset_index(drop=True)
base_ids = []
for i, r in link.iterrows():
    b = _base_link_id_from_row(r)
    if b is None:
        b = f"C{i+1}"
    base_ids.append(b)

link["__SWMM_ID_BASE__"] = base_ids

# 같은 베이스ID 중복 시 _2, _3 접미사 부여
counts = link.groupby("__SWMM_ID_BASE__").cumcount() + 1
link["SWMM_ID"] = np.where(
    counts == 1,
    link["__SWMM_ID_BASE__"],
    link["__SWMM_ID_BASE__"] + "_" + counts.astype(str)
)

# 혹시 전체적으로 또 중복이 발생하면(극히 드묾) 최종 고유화
if link["SWMM_ID"].duplicated().any():
    # 전역 중복까지 제거
    seen = {}
    fixed = []
    for s in link["SWMM_ID"]:
        if s not in seen:
            seen[s] = 1
            fixed.append(s)
        else:
            seen[s] += 1
            fixed.append(f"{s}_{seen[s]}")
    link["SWMM_ID"] = fixed

# === 길이(m) ===
def get_length_m(row):
    if F_LEN in row and pd.notna(row[F_LEN]) and _to_float(row[F_LEN], 0) > 0:
        return float(row[F_LEN])
    try:
        L = float(row.geometry.length)  # 좌표계가 m라고 가정
    except Exception:
        L = 1.0
    return max(L, 0.1)

# === 단면 텍스트 파싱(BOX) ===
RE_BARRELS = re.compile(r"(\d+)\s*@")                  # 2@, 3@ ...
RE_WHxH    = re.compile(r"(\d+(?:\.\d+)?)\s*[xX×]\s*(\d+(?:\.\d+)?)")

def concat_shape_text(row):
    parts = []
    for c in F_SHAPE_TXT:
        if c in row and isinstance(row[c], str) and row[c].strip():
            parts.append(row[c].strip())
    return " / ".join(parts) if parts else ""

def parse_box(text):
    """문자열에서 (폭m, 높이m, 배럴수) 반환. 실패 시 (nan, nan, 1)"""
    if not isinstance(text, str) or not text.strip():
        return (np.nan, np.nan, 1)
    t = text.upper()

    barrels = 1
    mb = RE_BARRELS.search(t)
    if mb:
        try: barrels = int(mb.group(1))
        except: barrels = 1

    m = RE_WHxH.search(t)
    if m:
        w = _to_float(m.group(1))
        h = _to_float(m.group(2))
        return (w, h, max(1, barrels))
    return (np.nan, np.nan, barrels)

# === 지름(m) 선택 ===
def to_m_smart(v):
    """
    기본 가정: 데이터는 'm' 단위.
    값이 20 초과면 mm로 보고 /1000 (예: 600, 1500 등)
    """
    v = _to_float(v, np.nan)
    if not pd.notna(v) or v <= 0:
        return np.nan
    if v > 20:
        return v / 1000.0
    return v

def pick_diameter_m(row):
    # 1) 원형 후보
    cand = []
    for c in [F_D_US, F_D_DS]:
        if c in row:
            val = to_m_smart(row.get(c, np.nan))
            if pd.notna(val) and val > 0:
                cand.append(val)
    if cand:
        return max(cand)

    # 2) BOX(폭x높이) 있으면 높이를 지름 대체로 사용
    w_box, h_box, _ = parse_box(concat_shape_text(row))
    if pd.notna(h_box) and h_box > 0:
        return h_box

    # 3) 폴백
    return DEFAULT_DIAM_M

link["__DIAM_M__"] = link.apply(pick_diameter_m, axis=1)
link.loc[(~pd.notna(link["__DIAM_M__"])) | (link["__DIAM_M__"] <= 0), "__DIAM_M__"] = DEFAULT_DIAM_M

n_default = int((link["__DIAM_M__"] == DEFAULT_DIAM_M).sum())
if n_default:
    print(f"⚠️ {n_default}개의 관에 유효 직경이 없어 기본값 {DEFAULT_DIAM_M} m로 대체했습니다.")


# ---- 경사/낙차 자동 보정 ----
MIN_DROP = 0.01    # 최소 낙차 1 cm
MAX_SLOPE = 0.8    # 비현실 급경사 상한 (m/m)

# 맨홀/토구의 표고 딕셔너리
node_elev = {}
for _, r in mh.iterrows():
    node_elev[str(r[F_MH_ID]).strip()] = _to_float(r[F_MH_INV], 0.0)
for _, r in of.iterrows():
    node_elev[str(r[F_OF_ID]).strip()] = _to_float(r[F_OF_INV], 0.0)

def _fix_us_ds_elev(us_id, ds_id, L):
    us_z = node_elev.get(us_id, np.nan)
    ds_z = node_elev.get(ds_id, np.nan)
    if not pd.notna(us_z) or not pd.notna(ds_z):
        return
    # 1) US<=DS이면 DS를 내려 최소낙차 확보
    if us_z <= ds_z:
        ds_z = us_z - MIN_DROP
    # 2) 급경사 제한: |us-ds| <= MAX_SLOPE * L
    drop = us_z - ds_z
    max_drop = max(MIN_DROP, MAX_SLOPE * max(L, 0.1))
    if drop > max_drop:
        ds_z = us_z - max_drop
    # 반영
    node_elev[ds_id] = ds_z

# 각 링크에 대해 적용
for _, r in link.iterrows():
    us = str(r[F_US_NODE]).strip()
    ds = str(r[F_DS_NODE]).strip()
    L  = get_length_m(r)
    _fix_us_ds_elev(us, ds, L)

# 보정된 표고를 테이블에 반영
mh[F_MH_INV] = mh[F_MH_ID].map(node_elev).fillna(mh[F_MH_INV])
of[F_OF_INV] = of[F_OF_ID].map(node_elev).fillna(of[F_OF_INV])

# ---- Outfall 다중 연결 정리 (ERROR 141 방지) ----
from collections import defaultdict

# outfall로 들어가는 링크 목록
in_to_outfall = defaultdict(list)
out_from_outfall = defaultdict(list)

for i, r in link.iterrows():
    us = str(r[F_US_NODE]).strip()
    ds = str(r[F_DS_NODE]).strip()
    if ds in of[F_OF_ID].values:
        in_to_outfall[ds].append(i)     # ds가 하구 = 하구로 유입
    if us in of[F_OF_ID].values:
        out_from_outfall[us].append(i)  # 하구에서 시작(허용 안 됨)

# 3-1) 하구에서 시작하는 링크는 방향 뒤집기 시도(가능하면)
for of_id, idxs in out_from_outfall.items():
    for i in idxs:
        us = str(link.at[i, F_US_NODE]).strip()
        ds = str(link.at[i, F_DS_NODE]).strip()
        # 뒤집어도 또 하구면 곤란 → 일단 스킵
        if ds in of[F_OF_ID].values:
            continue
        # 단순 스왑(지오메트리 방향과 무관하게 노드만 맞춤)
        link.at[i, F_US_NODE], link.at[i, F_DS_NODE] = ds, us

# 3-2) 하구 유입 링크가 2개 이상이면 자동 합류 맨홀 추가
auto_links = []   # 새로 만들 링크들의 행(딕셔너리)
auto_nodes = []   # 새로 만들 맨홀들의 행(딕셔너리)
for of_id, idxs in in_to_outfall.items():
    if len(idxs) <= 1:
        continue
    join_mh = f"AUTO-OF-{of_id}"
    # 합류 맨홀 좌표/표고: 하구 위치를 복제(표고는 하구보다 약간 높게)
    of_row = of.loc[of[F_OF_ID] == of_id].iloc[0]
    join_x, join_y = of_row.geometry.x, of_row.geometry.y
    join_z = _to_float(of_row[F_OF_INV], 0.0) + 0.05

    # 새 맨홀 레코드(간단 생성)
    auto_nodes.append({
        F_MH_ID: join_mh,
        F_MH_INV: join_z,
        "geometry": of_row.geometry
    })

    # 기존 모든 유입 링크의 DS를 합류 맨홀로 변경
    for i in idxs:
        link.at[i, F_DS_NODE] = join_mh

    # 합류 맨홀 → 하구 단일 링크 추가
    new_name = f"C_OF_{_clean_id(of_id)}"
    auto_links.append({
        "SWMM_ID": new_name,
        F_US_NODE: join_mh,
        F_DS_NODE: of_id,
        F_LEN: 1.0,                # 짧은 연결
        "__DIAM_M__": DEFAULT_DIAM_M,
        "geometry": of_row.geometry
    })

# 테이블에 반영
if auto_nodes:
    mh_auto = gpd.GeoDataFrame(auto_nodes, geometry="geometry", crs=mh.crs)
    mh = pd.concat([mh, mh_auto], ignore_index=True)
if auto_links:
    add_df = gpd.GeoDataFrame(auto_links, geometry="geometry", crs=link.crs)
    # 누락 칼럼 채우기
    for c in link.columns:
        if c not in add_df.columns:
            add_df[c] = np.nan
    link = pd.concat([link, add_df[link.columns]], ignore_index=True)

# ===== Outfall 최종 정리 + 검증 (ERROR 141 하드 고정) =====
from collections import defaultdict

def _unique_link_id(base, used):
    """SWMM_ID가 중복되지 않도록 고유 ID 생성"""
    if base not in used:
        used.add(base)
        return base
    k = 2
    while True:
        cand = f"{base}_{k}"
        if cand not in used:
            used.add(cand)
            return cand
        k += 1

# 현재 사용 중인 링크 ID 집합
used_ids = set(link["SWMM_ID"].astype(str))

# Outfall 집합
outfall_ids = set(of[F_OF_ID].astype(str))

# 1) Outfall에서 '시작'하는 링크(허용X) → 가능한 경우 뒤집기
out_from_outfall = defaultdict(list)
for i, r in link.iterrows():
    us = str(r[F_US_NODE]).strip()
    ds = str(r[F_DS_NODE]).strip()
    if us in outfall_ids:
        out_from_outfall[us].append(i)

for of_id, idxs in out_from_outfall.items():
    for i in idxs:
        us = str(link.at[i, F_US_NODE]).strip()
        ds = str(link.at[i, F_DS_NODE]).strip()
        # 양끝이 모두 Outfall인 희귀 케이스 → 제거
        if ds in outfall_ids:
            link.at[i, "SWMM_ID"] = None  # 마크 후 나중에 드롭
            continue
        # 방향 뒤집기
        link.at[i, F_US_NODE], link.at[i, F_DS_NODE] = ds, us

# 실제 제거
link = link[link["SWMM_ID"].notna()].reset_index(drop=True)

# 2) Outfall로 들어오는 링크 수를 조사
in_to_outfall = defaultdict(list)
for i, r in link.iterrows():
    ds = str(r[F_DS_NODE]).strip()
    if ds in outfall_ids:
        in_to_outfall[ds].append(i)

# 3) 2개 이상 유입 시: AUTO-OF-<ID> 합류맨홀 생성 + Outfall에는 단 1개의 링크만
auto_nodes = []
auto_links = []
for of_id, idxs in in_to_outfall.items():
    if len(idxs) <= 1:
        continue

    # 합류 맨홀 이름(중복 방지)
    join_mh = f"AUTO-OF-{_clean_id(of_id)}"
    while join_mh in set(mh[F_MH_ID]):
        join_mh = join_mh + "_2"

    # 하구 좌표/표고 복제(표고는 하구보다 +0.05 m)
    of_row = of.loc[of[F_OF_ID] == of_id].iloc[0]
    join_z = _to_float(of_row[F_OF_INV], 0.0) + 0.05

    auto_nodes.append({F_MH_ID: join_mh, F_MH_INV: join_z, "geometry": of_row.geometry})

    # 기존 유입 링크들의 DS를 합류맨홀로 변경
    for i in idxs:
        link.at[i, F_DS_NODE] = join_mh

    # 합류맨홀 -> Outfall 단일 링크 추가
    base = f"C_OF_{_clean_id(of_id)}"
    new_id = _unique_link_id(base, used_ids)
    auto_links.append({
        "SWMM_ID": new_id,
        F_US_NODE: join_mh,
        F_DS_NODE: of_id,
        F_LEN: 1.0,
        "__DIAM_M__": DEFAULT_DIAM_M,
        "geometry": of_row.geometry
    })

# 새 맨홀/링크 반영
if auto_nodes:
    mh_auto = gpd.GeoDataFrame(auto_nodes, geometry="geometry", crs=mh.crs)
    mh = pd.concat([mh, mh_auto], ignore_index=True)
if auto_links:
    add_df = gpd.GeoDataFrame(auto_links, geometry="geometry", crs=link.crs)
    for c in link.columns:
        if c not in add_df.columns:
            add_df[c] = np.nan
    link = pd.concat([link, add_df[link.columns]], ignore_index=True)

# 4) 최종 검증: Outfall별 유입/유출 카운트
in_count = defaultdict(int)
out_count = defaultdict(int)
for _, r in link.iterrows():
    us = str(r[F_US_NODE]).strip()
    ds = str(r[F_DS_NODE]).strip()
    if ds in outfall_ids:
        in_count[ds] += 1
    if us in outfall_ids:
        out_count[us] += 1

problems = []
for of_id in outfall_ids:
    ic = in_count.get(of_id, 0)
    oc = out_count.get(of_id, 0)
    if ic > 1 or oc > 0:
        problems.append((of_id, ic, oc))

if problems:
    print("❌ 여전히 규칙 위반 Outfall이 있습니다 (이름, 유입개수, 유출개수):")
    for p in problems[:10]:
        print("   ", p)
else:
    print("✅ Outfall 규칙 확인: 모든 Outfall은 유입 1개, 유출 0개입니다.")

# ================================================
# ★★★ 유역(SUBCATCHMENTS/ POLYGONS) 구성 파트 ★★★
# ================================================
def _pick_first_col(gdf, cands):
    for c in cands:
        if c in gdf.columns:
            return c
    return None

# 유역 ID 생성
c_id_col = _pick_first_col(catch, C_ID_CANDS)
if c_id_col:
    catch["__CID__"] = catch[c_id_col].astype(str).apply(_clean_id)
else:
    catch = catch.reset_index(drop=True)
    catch["__CID__"] = ["S{}".format(i+1) for i in range(len(catch))]

# 면적(ha) 계산
# SWMM은 유역 면적 단위가 ha
catch["__AREA_HA__"] = catch.geometry.area / 10000.0

# 경사(m/m)
sl_col = _pick_first_col(catch, C_SLOPE_CANDS)
if sl_col and sl_col in catch.columns:
    s = pd.to_numeric(catch[sl_col], errors="coerce")
    # 만약 값이 1 이상이 다수면 퍼센트(%)로 판단하여 /100
    if s.dropna().gt(1.5).mean() > 0.5:
        s = s / 100.0
    catch["__SLOPE__"] = s.fillna(DEFAULT_SLOPE).clip(lower=0.0001)
else:
    catch["__SLOPE__"] = DEFAULT_SLOPE

# 불투수율(%)
p_col = _pick_first_col(catch, C_PIMP_CANDS)
if p_col and p_col in catch.columns:
    p = pd.to_numeric(catch[p_col], errors="coerce")
    # 0~1로 저장돼 있으면 %로 환산
    if p.dropna().le(1.0).mean() > 0.5:
        p = p * 100.0
    catch["__PIMP__"] = p.fillna(DEFAULT_PIMP).clip(lower=0, upper=100)
else:
    catch["__PIMP__"] = DEFAULT_PIMP

# Width(m): 없으면 A≈L*W 가정으로 W= sqrt(A) * k (k=2) 등 단순 추정
w_col = _pick_first_col(catch, C_WIDTH_CANDS)
if w_col and w_col in catch.columns:
    wv = pd.to_numeric(catch[w_col], errors="coerce")
    catch["__WIDTH__"] = wv
else:
    # 간단 추정: W = 2 * sqrt(A)  (A in m2)
    A_m2 = catch.geometry.area
    catch["__WIDTH__"] = (2.0 * np.sqrt(A_m2)).clip(lower=10.0)

# --- 유역 중앙점 기준 가장 가까운 맨홀 찾기 (교체 코드) ---

# 1) 유역 기하 유효화 + 중심점 계산
catch["geometry"] = catch.geometry.buffer(0)  # invalid polygon 보정
cent_gdf = gpd.GeoDataFrame(
    {"__CID__": catch["__CID__"]},
    geometry=catch.geometry.centroid,
    crs=catch.crs
)

# 2) 맨홀 포인트 준비
mh_pts = mh[[F_MH_ID, "geometry"]].rename(columns={F_MH_ID: "__MID__"})

# 3) sjoin_nearest 우선 사용, 실패 시 STRtree로 대체
try:
    joined = gpd.sjoin_nearest(cent_gdf, mh_pts, how="left")
    cid_to_mid = dict(zip(joined["__CID__"], joined["__MID__"]))
except Exception:
    from shapely.strtree import STRtree
    cent_arr = np.array(cent_gdf.geometry.values, dtype=object)  # Geometry 배열
    mh_arr   = np.array(mh_pts.geometry.values, dtype=object)
    tree = STRtree(mh_arr)
    pairs = tree.query_nearest(cent_arr)  # (src_idx, dst_idx) 또는 (N,2)
    if isinstance(pairs, tuple) and len(pairs) == 2:
        _, mh_idx = pairs
    else:
        mh_idx = pairs[:, 1]
    nearest_mid = mh_pts.iloc[mh_idx]["__MID__"].to_numpy()
    cid_to_mid = dict(zip(cent_gdf["__CID__"], nearest_mid))

# 4) 유역별 Outlet 설정
catch["__OUTLET__"] = catch["__CID__"].map(cid_to_mid)


# ================================================
# -------------- INP 작성 ------------------------
# ================================================
lines = []
w = lines.append

w("[TITLE]")
w(";; Munrae district sewer network (auto-converted from shp)")
w("")

# ---- SUBCATCHMENTS (유역) ----
w("[SUBCATCHMENTS]")
w(";;Name         RainGage   Outlet        Area    %Imperv  Width   Slope   CurbLen   SnowPack")
for _, r in catch.iterrows():
    name  = r["__CID__"]
    gage  = RAIN_GAGE_NAME
    outlet= str(r["__OUTLET__"])
    area  = float(r["__AREA_HA__"])
    pimp  = float(r["__PIMP__"])
    width = float(r["__WIDTH__"])
    slope = float(r["__SLOPE__"])
    w(f"{name:<13} {gage:<10} {outlet:<12} {area:<7.3f} {pimp:<8.1f} {width:<7.1f} {slope:<7.4f} {CURB_LEN:<8.1f}")

w("")

# (필요 시) RAINGAGES / TIMESERIES 블록은 별도 파일에서 유지
# 여기서는 유역/관망만 쓰므로 기존 코드 흐름 유지

# ---- JUNCTIONS ----
w("[JUNCTIONS]")
w(";;Name        Elevation  MaxDepth  InitDepth  SurDepth  Aponded")
for _, r in mh.iterrows():
    name = str(r[F_MH_ID]).strip()
    elev = _to_float(r[F_MH_INV], 0.0)
    w(f"{name:<12} {elev:<10.3f} {DEFAULT_MAXDEPTH_MH:.1f} 0 0 0")
w("")

# ---- OUTFALLS ----
w("[OUTFALLS]")
w(";;Name        Elevation   Type")
for _, r in of.iterrows():
    name = str(r[F_OF_ID]).strip()
    elev = _to_float(r[F_OF_INV], 0.0)
    w(f"{name:<12} {elev:<10.3f} FREE")
w("")

# ---- CONDUITS ----
w("[CONDUITS]")
w(";;Name        FromNode    ToNode     Length     Roughness  InOffset  OutOffset  InitFlow  MaxFlow")
def get_length_m(row):
    if F_LEN in row and pd.notna(row[F_LEN]) and _to_float(row[F_LEN], 0) > 0:
        return float(row[F_LEN])
    try:
        L = float(row.geometry.length)
    except Exception:
        L = 1.0
    return max(L, 0.1)

for _, r in link.iterrows():
    name = r["SWMM_ID"]
    us = str(r[F_US_NODE]).strip()
    ds = str(r[F_DS_NODE]).strip()
    L  = get_length_m(r)
    if not pd.notna(L) or L <= 0:
        L = 0.1
    w(f"{name:<12} {us:<12} {ds:<12} {L:<10.2f} {ROUGH_N:<10.4f} 0 0 0 0")
w("")

# ---- XSECTIONS ----
w("[XSECTIONS]")
w(";;Link        Shape        Geom1    Geom2   Geom3   Geom4   Barrels   Culvert")
RE_BARRELS = re.compile(r"(\d+)\s*@")
RE_WHxH    = re.compile(r"(\d+(?:\.\d+)?)\s*[xX×]\s*(\d+(?:\.\d+)?)")
def concat_shape_text(row):
    parts = []
    for c in F_SHAPE_TXT:
        if c in row and isinstance(row[c], str) and row[c].strip():
            parts.append(row[c].strip())
    return " / ".join(parts) if parts else ""
def parse_box(text):
    if not isinstance(text, str) or not text.strip():
        return (np.nan, np.nan, 1)
    t = text.upper()
    barrels = 1
    mb = RE_BARRELS.search(t); 
    if mb:
        try: barrels = int(mb.group(1))
        except: barrels = 1
    m = RE_WHxH.search(t)
    if m:
        wv = _to_float(m.group(1)); hv = _to_float(m.group(2))
        return (wv, hv, max(1, barrels))
    return (np.nan, np.nan, barrels)
def to_m_smart(v):
    v = _to_float(v, np.nan)
    if not pd.notna(v) or v <= 0: return np.nan
    return v/1000.0 if v>20 else v
def pick_diameter_m(row):
    cand=[]
    for c in [F_D_US, F_D_DS]:
        if c in row:
            val=to_m_smart(row.get(c, np.nan))
            if pd.notna(val) and val>0: cand.append(val)
    if cand: return max(cand)
    w_box,h_box,_=parse_box(concat_shape_text(row))
    if pd.notna(h_box) and h_box>0: return h_box
    return DEFAULT_DIAM_M

if "__DIAM_M__" not in link.columns:
    link["__DIAM_M__"] = link.apply(pick_diameter_m, axis=1)
for _, r in link.iterrows():
    name = r["SWMM_ID"]
    shape_text = concat_shape_text(r)
    w_box, h_box, barrels = parse_box(shape_text)
    if pd.notna(w_box) and pd.notna(h_box) and w_box>0 and h_box>0:
        w(f"{name:<12} RECT_CLOSED  {w_box:<7.3f} {h_box:<7.3f} 0      0      {int(barrels):<3d}")
        continue
    D = _to_float(r["__DIAM_M__"], 0.0)
    if not pd.notna(D) or D <= 0:
        D = DEFAULT_DIAM_M
    w(f"{name:<12} CIRCULAR     {D:<7.3f} 0       0      0      1")
w("")

# ---- OPTIONS ----
w("[OPTIONS]")
w("FLOW_UNITS           LPS")
w("INFILTRATION         HORTON")
w("FLOW_ROUTING         DYNWAVE")
w("START_DATE           01/01/2020")
w("REPORT_START_DATE    01/01/2020")
w("START_TIME           00:00:00")
w("REPORT_START_TIME    00:00:00")
w("END_DATE             01/02/2020")
w("END_TIME             00:00:00")
w("")

# ---- POLYGONS (유역 경계) ----
w("[POLYGONS]")
w(";;Subcatchment   X-Coord       Y-Coord")
for _, row in catch.iterrows():
    name = row["__CID__"]
    geom = row.geometry
    # 외곽선만 사용(멀티폴리곤은 첫 파트)
    if geom.geom_type == "MultiPolygon":
        geom = list(geom.geoms)[0]
    if geom.geom_type != "Polygon":
        continue
    coords = list(geom.exterior.coords)
    # 마지막 점(폐합점)이 첫 점과 같으면 제거
    if len(coords) > 1 and coords[0] == coords[-1]:
        coords = coords[:-1]
    for x, y in coords:
        w(f"{name:<14} {x:<13.3f} {y:<13.3f}")
w("")

# ---- COORDINATES (노드 좌표) ----
w("[COORDINATES]")
w(";;Node       X-Coord     Y-Coord")
for _, r in mh.iterrows():
    name = str(r[F_MH_ID]).strip()
    x, y = r.geometry.x, r.geometry.y
    w(f"{name:<12} {x:<12.3f} {y:<12.3f}")
for _, r in of.iterrows():
    name = str(r[F_OF_ID]).strip()
    x, y = r.geometry.x, r.geometry.y
    w(f"{name:<12} {x:<12.3f} {y:<12.3f}")
w("")

# ---- REPORT ----
w("[REPORT]")
w("INPUT      YES")
w("CONTROLS   NO")
w("NODES      ALL")
w("LINKS      ALL")
w("")

# ---- END ----
text = "\n".join(lines)
OUT_INP.write_text(text, encoding="utf-8")
print(f"✅ SWMM input file created: {OUT_INP}")
