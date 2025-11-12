# 수정일자: 2025-11-07
# 작성자: 정규원

# gis데이터를 SWMM에서 import할 수 있도록
# .shp파일을 .inp파일로 변환

# gis데이터를 SWMM에서 import할 수 있도록
# .shp파일을 .inp파일로 변환
# === 변경 핵심 ===
# 1) 링크 ID를 'SWMM_ID' 칼럼으로 한번만 생성 → [CONDUITS], [XSECTIONS]에 공통 사용(중복 자동 제거)
# 2) OBJECTID 등 원본 식별자가 같아도 _2, _3 접미사로 고유화
# 3) ID에 허용되지 않는 문자 제거(영숫자/언더스코어/하이픈/콜론만 유지)
# 4) 지름/길이 폴백 로직 강화(0/NaN 방지)

import geopandas as gpd
import pandas as pd
import numpy as np
import re
from pathlib import Path

# === 입력 파일 ===
LINK_SHP = Path("newshp/m_link_fixed.shp")
MH_SHP   = Path("newshp/m_manhole_fixed.shp")
OF_SHP   = Path("newshp/토구_fixed.shp")

# === 출력 파일 ===
OUT_INP  = Path("newshp/munrae_model.inp")

# === 필드 매핑 ===
F_US_NODE = "SAT_MHE"
F_DS_NODE = "END_MHE"
F_US_INV  = "SPT_HSL"      # 관저고(상)
F_DS_INV  = "EPT_HSL"      # 관저고(하)
F_LEN     = "Shape_Leng"
F_D_US    = "ST_PIP_HIT"   # 원형 내경(상)
F_D_DS    = "ET_PIP_HIT"   # 원형 내경(하)

# 사각/단면 문자열 후보
F_SHAPE_TXT = ["DR_LP_CCE", "DTR_CCE", "DP_PE_CCE"]

# 맨홀/토구
F_MH_ID   = "MHE_IDN"
F_MH_INV  = "HSL"
F_OF_ID   = "IDN"
F_OF_INV  = "ATTE"

# === 상수 ===
DEFAULT_DIAM_M = 0.6       # 원형 폴백(m)
ROUGH_N        = 0.013     # Manning n
DEFAULT_MAXDEPTH_MH = 6.0  # 노드 기본 MaxDepth

# === 데이터 읽기 ===
link = gpd.read_file(LINK_SHP)
mh   = gpd.read_file(MH_SHP)
of   = gpd.read_file(OF_SHP)

# === 기본 유틸 ===
def _clean_id(x): 
    s = str(x).strip()
    # 허용: A-Za-z0-9 _ - :
    s = re.sub(r"[^A-Za-z0-9_\-:]", "_", s)
    return s

def _to_float(x, default=np.nan):
    try:
        v = float(x)
        return v
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


# ---------------- INP 작성 ----------------
lines = []
w = lines.append

w("[TITLE]")
w(";; Munrae district sewer network (auto-converted from shp)")
w("")

# JUNCTIONS
w("[JUNCTIONS]")
w(";;Name        Elevation  MaxDepth  InitDepth  SurDepth  Aponded")
for _, r in mh.iterrows():
    name = str(r[F_MH_ID]).strip()
    elev = _to_float(r[F_MH_INV], 0.0)
    w(f"{name:<12} {elev:<10.3f} {DEFAULT_MAXDEPTH_MH:.1f} 0 0 0")
w("")

# OUTFALLS
w("[OUTFALLS]")
w(";;Name        Elevation   Type")
for _, r in of.iterrows():
    name = str(r[F_OF_ID]).strip()
    elev = _to_float(r[F_OF_INV], 0.0)
    w(f"{name:<12} {elev:<10.3f} FREE")
w("")

# CONDUITS
w("[CONDUITS]")
w(";;Name        FromNode    ToNode     Length     Roughness  InOffset  OutOffset  InitFlow  MaxFlow")
for _, r in link.iterrows():
    name = r["SWMM_ID"]
    us = str(r[F_US_NODE]).strip()
    ds = str(r[F_DS_NODE]).strip()
    L  = get_length_m(r)
    # 길이 0 방지
    if not pd.notna(L) or L <= 0:
        L = 0.1
    w(f"{name:<12} {us:<12} {ds:<12} {L:<10.2f} {ROUGH_N:<10.4f} 0 0 0 0")
w("")

# XSECTIONS
w("[XSECTIONS]")
w(";;Link        Shape        Geom1    Geom2   Geom3   Geom4   Barrels   Culvert")
for _, r in link.iterrows():
    name = r["SWMM_ID"]

    # 1) BOX 우선
    shape_text = concat_shape_text(r)
    w_box, h_box, barrels = parse_box(shape_text)
    if pd.notna(w_box) and pd.notna(h_box) and w_box > 0 and h_box > 0:
        w(f"{name:<12} RECT_CLOSED  {w_box:<7.3f} {h_box:<7.3f} 0      0      {int(barrels):<3d}")
        continue

    # 2) 원형
    D = _to_float(r["__DIAM_M__"], 0.0)
    if not pd.notna(D) or D <= 0:
        D = DEFAULT_DIAM_M
    w(f"{name:<12} CIRCULAR     {D:<7.3f} 0       0      0      1")
w("")

# OPTIONS
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

# COORDINATES
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

# REPORT
w("[REPORT]")
w("INPUT      YES")
w("CONTROLS   NO")
w("NODES      ALL")
w("LINKS      ALL")
w("")

# END
text = "\n".join(lines)
OUT_INP.write_text(text, encoding="utf-8")
print(f"✅ SWMM input file created: {OUT_INP}")
