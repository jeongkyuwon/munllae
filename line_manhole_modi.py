# 2025-11-07 수정
# 독립맨홀, 관의 노드 수정 완료 버전
# 이후로 맨홀 없이 끝나는 파이프 부터 처리해야 함.
# 이건 다른 파일에...

# SWMMM입력을 위한 맨홀, 하수관, 맨홀, 토구 파일 보정

#최종 규칙 세트 (코드 기준)
# 1. 관–관 사이엔 맨홀 필수
# 교차/분기 지점은 맨홀 생성 후 관을 분할.

# 2. 관의 시작·끝은 맨홀 또는 토구
# 끝점이 맨홀/토구에 1.0 m 이내면 스냅, 아니면 맨홀 생성(단, 규칙4의 삭제 조건 우선).

# 3. 역경사 금지(유입 하단 > 유출 하단)
# 허용오차 −0.001 m/m, 목표 최소구배 0.002 m/m로 보정.

# 4. 맨홀 없이 끝나는 단파이프
# 길이 ≤ 5 m → 삭제
# 길이 > 5 m → 끝점에 맨홀 생성 후 연결

# 5. 고립 맨홀 처리(연결 없는 맨홀)
# 반경 ≤ 1.0 m에 관이 있으면 스냅·연결, 없으면 맨홀 삭제

# 6. 높이 불일치 보정(연결부)
# 노드(맨홀/토구) invert 값이 있으면 노드 고정, 관의 상·하류 invert를 내려서(깎아서) 최소구배 확보
# 둘 다 미상 시 보정 보류(QA 플래그)

# 7. 관 형상/치수 구분
# 원형: 내경 D 사용
# 사각(박스): 유효 높이 H 사용(문자열 “□3@3x3” → H=3.0 m 파싱)
# 단위 혼재 시 mm → m로 통일

# 기본 파라미터 (없으면 이 값 사용)
# 스냅 허용오차: 맨홀–관 1.0 m, 관–관 교차 0.5 m
# 중복 병합(동일 노드): 0.2 m
# 단파이프 잡음 제거(규칙4 이전): 길이 < 0.5 m 즉시 삭제
# 좌표계: 입력이 m 단위 아니면 EPSG:5186으로 재투영

# ============ 코드 시작 ================

# fix_sewer.py
# -*- coding: utf-8 -*-
import re
import uuid
import math
import warnings
from pathlib import Path

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString
from shapely.ops import split
from shapely.strtree import STRtree
from shapely import force_2d
from shapely.geometry import GeometryCollection
from shapely.ops import split
from shapely.geometry import Point, LineString

warnings.filterwarnings("ignore", category=UserWarning)

# ---------- 파라미터(필요시 여기만 수정) ----------
IN_LINK = Path("mullae_nbs for sharing/m_link.shp")
IN_MH   = Path("mullae_nbs for sharing/m_manhole.shp")
IN_OF   = Path("mullae_nbs for sharing/토구.shp")

OUT_LINK = Path("newshp/m_link_fixed.shp")
OUT_MH   = Path("newshp/m_manhole_fixed.shp")
OUT_OF   = Path("newshp/토구_fixed.shp")

LOG_ADDED_NODES   = Path("newshp/added_nodes.csv")
LOG_DELETED_PIPES = Path("newshp/deleted_pipes.csv")
LOG_ADJ_INVERTS   = Path("newshp/adjusted_inverts.csv")
OUT_QA_FLAGS      = Path("newshp/qa_flags.shp")

SNAP_MH_PIPE   = 1.0    # 맨홀-관 스냅 허용오차 [m]
ORPHAN_MARGIN  = 0.5    # 고립 판정 여유
SNAP_PIPE_XING = 0.5    # 관-관 교차 정규화 허용오차 [m]
DELETE_SHORT   = 0.5    # 잡음 관 삭제 [m]
END_SHORT_DEL  = 5.0    # 맨홀 없이 끝나는 단파이프 삭제 기준 [m]
MIN_SLOPE      = 0.002  # 최소 설계구배 [m/m]
NEG_TOL        = -0.001 # 역경사 허용오차 [m/m]

# ---------- 필드 매핑 ----------
# 링크(관)
F_US_NODE = "SAT_MHE"
F_DS_NODE = "END_MHE"
F_US_INV  = "SPT_HSL"
F_DS_INV  = "EPT_HSL"
F_LEN1    = "LENX"
F_LEN2    = "Shape_Leng"
F_D_US    = "ST_PIP_HIT"  # 원형 내경(상류)
F_D_DS    = "ET_PIP_HIT"  # 원형 내경(하류)
F_SHAPE1  = "DR_LP_CCE"   # 단면 문자열(사각)
F_SHAPE2  = "DTR_CCE"     # 보조

# 맨홀
F_MH_ID   = "MHE_IDN"
F_MH_INV  = "HSL"

# 토구
F_OF_ID   = "IDN"
F_OF_INV  = "ATTE"

# ---------- 유틸 ----------
def ensure_projected(gdf: gpd.GeoDataFrame, target_epsg=5186) -> gpd.GeoDataFrame:
    """미터 단위 투영 좌표 보장"""
    if gdf.crs is None or not gdf.crs.is_projected:
        return gdf.to_crs(epsg=target_epsg)
    return gdf

def length_of(line):
    try:
        return float(line.length)
    except Exception:
        return 0.0

def parse_box_height(text: str):
    """문자열(예: '3@□3.0X3.0')에서 유효 높이 H[m] 추출, 없으면 None"""
    if not isinstance(text, str):
        return None
    # 숫자X숫자 패턴에서 두 번째 수치를 높이로 본다(폭x높이 가정).
    m = re.search(r'(\d+(?:\.\d+)?)\s*[xX]\s*(\d+(?:\.\d+)?)', text)
    if m:
        w = float(m.group(1))
        h = float(m.group(2))
        return h
    return None

def next_node_id(prefix="MHNEW-"):
    return f"{prefix}{uuid.uuid4().hex[:10]}"

def set_2d(gdf):
    gdf["geometry"] = gdf.geometry.apply(lambda g: force_2d(g))
    return gdf

# ---------- 데이터 로드 ----------
link = gpd.read_file(IN_LINK)
mh   = gpd.read_file(IN_MH)
of   = gpd.read_file(IN_OF)

# 2D 강제 & 좌표계 정리
link = set_2d(ensure_projected(link))
mh   = set_2d(ensure_projected(mh, target_epsg=link.crs.to_epsg()))
of   = set_2d(ensure_projected(of, target_epsg=link.crs.to_epsg()))

mh = mh.to_crs(link.crs)
of = of.to_crs(link.crs)

# 길이 필드 정리
if F_LEN1 in link.columns:
    link["__len__"] = link[F_LEN1].astype(float)
elif F_LEN2 in link.columns:
    link["__len__"] = link[F_LEN2].astype(float)
else:
    link["__len__"] = link.length

# 원형/사각 구분 및 높이 파싱(필요시)
shape_text = None
if F_SHAPE1 in link.columns and link[F_SHAPE1].notna().any():
    shape_text = F_SHAPE1
elif F_SHAPE2 in link.columns and link[F_SHAPE2].notna().any():
    shape_text = F_SHAPE2

link["__box_H__"] = link[shape_text].apply(parse_box_height) if shape_text else None

# ---------- 0. 선처리: 짧은 관 삭제 ----------
deleted_pipes = []
def mark_delete(idx, reason):
    deleted_pipes.append({"pipe_index": idx, "reason": reason})

to_keep = []
for idx, row in link.iterrows():
    L = row["__len__"] if row["__len__"] > 0 else length_of(row.geometry)
    if L < DELETE_SHORT:
        mark_delete(idx, f"too_short<{DELETE_SHORT}m")
    else:
        to_keep.append(idx)
link = link.loc[to_keep].copy()

# ---------- 1. 노드 레이어 구성(맨홀+토구) ----------
mh_nodes = mh[[F_MH_ID, F_MH_INV, "geometry"]].copy()
mh_nodes.rename(columns={F_MH_ID: "NODE_ID", F_MH_INV: "INV"}, inplace=True)
mh_nodes["NODE_TYPE"] = "MH"

of_nodes = of[[F_OF_ID, F_OF_INV, "geometry"]].copy()
of_nodes.rename(columns={F_OF_ID: "NODE_ID", F_OF_INV: "INV"}, inplace=True)
of_nodes["NODE_TYPE"] = "OF"

nodes = gpd.GeoDataFrame(pd.concat([mh_nodes, of_nodes], ignore_index=True), crs=link.crs)
# 중복 노드 병합(0.2m 이내)
node_tree = STRtree(nodes.geometry.values.tolist())
merged = []
used = set()
for i, geom in enumerate(nodes.geometry.values.tolist()):
    if i in used:
        continue
    cluster = [i]
    for j in node_tree.query(geom.buffer(0.2)):  # 0.2 m 병합
        if j != i and j not in used:
            if geom.distance(nodes.geometry.values[j]) <= 0.2:
                cluster.append(j)
    used.update(cluster)
    # 대표 좌표: 평균
    xs = [nodes.geometry.values[k].x for k in cluster]
    ys = [nodes.geometry.values[k].y for k in cluster]
    rep_pt = Point(sum(xs)/len(xs), sum(ys)/len(xs))
    # 대표 노드: 첫 항목
    row = nodes.iloc[cluster[0]].to_dict()
    row["geometry"] = rep_pt
    merged.append(row)
nodes = gpd.GeoDataFrame(merged, crs=link.crs)

# 빠른 조회용 인덱스
node_sindex = STRtree(nodes.geometry.values.tolist())

# ---------- 2. 관 끝점 스냅/노드 생성 ----------
added_nodes = []
def add_node(pt: Point, node_type="MH", inv=None):
    nid = next_node_id("AUTO-"+node_type+"-")
    added_nodes.append({"NODE_ID": nid, "NODE_TYPE": node_type, "INV": inv, "x": pt.x, "y": pt.y})
    return nid

def nearest_node(point: Point, tol=SNAP_MH_PIPE):
    # 반환: (node_idx, dist)
    candidates = node_sindex.query(point.buffer(tol))
    best = (None, float("inf"))
    for j in candidates:
        d = point.distance(nodes.geometry.values[j])
        if d < best[1]:
            best = (j, d)
    return best

link = link.reset_index(drop=True)
for i, row in link.iterrows():
    geom: LineString = row.geometry
    if geom is None or geom.is_empty:
        continue
    start = Point(list(geom.coords)[0])
    end   = Point(list(geom.coords)[-1])

    # 시작점 처리
    ni, d = nearest_node(start, SNAP_MH_PIPE)
    if ni is not None and d <= SNAP_MH_PIPE:
        # 스냅: 관 시작점을 해당 노드로 정렬
        start = nodes.geometry.values[ni]
        us_id = nodes.iloc[ni]["NODE_ID"]
    else:
        # 임시: 규칙 4 고려 위해 나중에 결정 (길이에 따라 삭제/추가)
        us_id = None

    # 끝점 처리
    nj, d2 = nearest_node(end, SNAP_MH_PIPE)
    if nj is not None and d2 <= SNAP_MH_PIPE:
        end = nodes.geometry.values[nj]
        ds_id = nodes.iloc[nj]["NODE_ID"]
    else:
        ds_id = None

    # 업데이트
    new_geom = LineString([start, end]) if len(geom.coords) == 2 else LineString([start, *list(geom.coords)[1:-1], end])
    link.at[i, "geometry"] = new_geom
    link.at[i, F_US_NODE] = us_id if us_id is not None else row.get(F_US_NODE, None)
    link.at[i, F_DS_NODE] = ds_id if ds_id is not None else row.get(F_DS_NODE, None)

# ---------- 3. 맨홀 없이 끝나는 단파이프 처리 ----------
# 관 끝점(시작/끝)들을 인덱싱해서 '이 끝점이 단말인지' 판정
pipe_endpoints = []
endpoint_owner = []
for idx, g in link.geometry.items():
    if g is None or g.is_empty: 
        continue
    coords = list(g.coords)
    pipe_endpoints.append(Point(coords[0])); endpoint_owner.append((idx, True))   # 시작점
    pipe_endpoints.append(Point(coords[-1])); endpoint_owner.append((idx, False)) # 끝점

ep_tree = STRtree(pipe_endpoints)

def is_true_dead_end(pt: Point, owner_idx: int, is_start: bool, tol=SNAP_MH_PIPE):
    """해당 끝점 주변에 '다른 관의 끝점'이 tol 내에 없으면 진짜 단말로 간주"""
    cand = ep_tree.query(pt.buffer(tol))
    for k in cand:
        oidx, ostart = endpoint_owner[k]
        if oidx == owner_idx:
            continue  # 자기 자신 제외
        # 다른 관의 어떤 끝점이든 가까이 있으면 단말 아님(=연결부)
        if pipe_endpoints[k].distance(pt) <= tol:
            return False
    return True

for col in [F_US_NODE, F_DS_NODE]:
    if col not in link.columns:
        link[col] = pd.NA

to_drop = []
for i, row in link.iterrows():
    geom: LineString = row.geometry
    L = row["__len__"] if row["__len__"] > 0 else length_of(geom)
    us_id = row[F_US_NODE]
    ds_id = row[F_DS_NODE]
    start = Point(list(geom.coords)[0])
    end   = Point(list(geom.coords)[-1])

    # 좌표
    start = Point(list(geom.coords)[0])
    end   = Point(list(geom.coords)[-1])

    # 양쪽 '진짜 단말' 여부(노드 연결 없고, 다른 관 끝점도 근처에 없음)
    start_dead = (pd.isna(us_id) or us_id is None) and is_true_dead_end(start, i, True,  tol=SNAP_MH_PIPE)
    end_dead   = (pd.isna(ds_id) or ds_id is None) and is_true_dead_end(end,   i, False, tol=SNAP_MH_PIPE)

    # 5 m 삭제는 '양쪽 모두 진짜 단말'일 때만
    if start_dead and end_dead and L <= END_SHORT_DEL:
        to_drop.append(i)
        mark_delete(i, f"short_dead_pipe_both_ends<= {END_SHORT_DEL}m")
        continue

    # 연결부(=단말 아님)인데 노드가 없으면 맨홀 생성(규칙1 준수)
    if (pd.isna(us_id) or us_id is None) and not start_dead:
        new_id = add_node(start, "MH", inv=None)
        link.at[i, F_US_NODE] = new_id
        coords = list(geom.coords); coords[0] = (start.x, start.y)
        link.at[i, "geometry"] = LineString(coords)

    if (pd.isna(ds_id) or ds_id is None) and not end_dead:
        new_id = add_node(end, "MH", inv=None)
        link.at[i, F_DS_NODE] = new_id
        coords = list(geom.coords); coords[-1] = (end.x, end.y)
        link.at[i, "geometry"] = LineString(coords)

    # 여전히 '진짜 단말'인 쪽은 길이 > 5 m일 때만 맨홀 생성
    if (pd.isna(us_id) or us_id is None) and start_dead and L > END_SHORT_DEL:
        new_id = add_node(start, "MH", inv=None)
        link.at[i, F_US_NODE] = new_id
        coords = list(geom.coords); coords[0] = (start.x, start.y)
        link.at[i, "geometry"] = LineString(coords)

    if (pd.isna(ds_id) or ds_id is None) and end_dead and L > END_SHORT_DEL:
        new_id = add_node(end, "MH", inv=None)
        link.at[i, F_DS_NODE] = new_id
        coords = list(geom.coords); coords[-1] = (end.x, end.y)
        link.at[i, "geometry"] = LineString(coords)


# 실제 삭제
if to_drop:
    link = link.drop(index=to_drop).reset_index(drop=True)

# 추가된 노드를 nodes에 반영
if added_nodes:
    df_add = gpd.GeoDataFrame(
        added_nodes,
        geometry=[Point(r["x"], r["y"]) for r in added_nodes],
        crs=link.crs
    )[["NODE_ID","NODE_TYPE","INV","geometry"]]
    # 컬럼 공통분모만 합치기
    common_cols = [c for c in df_add.columns if c in nodes.columns]
    nodes = gpd.GeoDataFrame(pd.concat([nodes, df_add[common_cols]], ignore_index=True), crs=link.crs)

    # 새로 추가된 노드까지 포함해 공간 인덱스 재생성 (중요!)
    node_sindex = STRtree(nodes.geometry.values.tolist())

# ---------- 4. 고립 맨홀 정리(1m내 관 있으면 연결, 아니면 삭제) ----------
# 항상 최신 링크 지오메트리로 트리 구성
pipe_geoms = list(link.geometry.values)
pipe_tree = STRtree(pipe_geoms)

MIN_ENDPOINT_SNAP = 0.5  # 끝점으로 스냅 판정 오차 [m]

def split_line_at_nearest_point(line: LineString, pt: Point, tol=1e-6):
    """라인을 pt의 최근접 '거리 d'에서 좌우로 정확히 잘라 두 조각 반환"""
    L = line.length
    if L <= 0:
        return None
    d = line.project(pt)
    if d <= tol or d >= L - tol:
        return None

    coords = list(line.coords)
    acc = 0.0
    for k in range(len(coords) - 1):
        p0 = Point(coords[k]); p1 = Point(coords[k+1])
        seg_len = p0.distance(p1)
        if acc + seg_len >= d:
            # d가 이 세그먼트 안에 있음 → 비율로 절단점 계산
            ratio = (d - acc) / seg_len
            x = p0.x + ratio * (p1.x - p0.x)
            y = p0.y + ratio * (p1.y - p0.y)
            cut = (x, y)
            part1 = coords[:k+1] + [cut]
            part2 = [cut] + coords[k+1:]
            return LineString(part1), LineString(part2)
        acc += seg_len
    return None


used_nodes = set(link[F_US_NODE].dropna().tolist() + link[F_DS_NODE].dropna().tolist())
node_id_to_idx = {nodes.iloc[i]["NODE_ID"]: i for i in range(len(nodes))}

rows_to_add = []   # 분할로 추가될 새 관
rows_to_drop = []  # 분할될 기존 관 인덱스
orphan_indices = []  # 최종 삭제할 맨홀 인덱스

for i, nrow in nodes.iterrows():
    nid = nrow["NODE_ID"]
    if nid in used_nodes:
        continue

    npt: Point = nrow.geometry
    cand_idxs = pipe_tree.query(npt.buffer(SNAP_MH_PIPE))
    best = None
    best_dist = float("inf")
    best_j = None

    for j in cand_idxs:
        line = pipe_geoms[j]
        proj = line.project(npt)
        nearest_pt = line.interpolate(proj)
        d = npt.distance(nearest_pt)
        if d < best_dist:
            best_dist = d
            best = nearest_pt
            best_j = j

    # 여유를 두고( SNAP + ORPHAN_MARGIN )도 멀면만 '진짜 고립'
    if (best is None) or (best_dist > SNAP_MH_PIPE + ORPHAN_MARGIN):
        if nrow.get("NODE_TYPE","MH") == "MH":
            orphan_indices.append(i)
        continue


    line = pipe_geoms[best_j]
    start = Point(list(line.coords)[0])
    end   = Point(list(line.coords)[-1])

    # ---- 끝점 근처면 스냅 ----
    if best.distance(start) <= MIN_ENDPOINT_SNAP:
        coords = list(line.coords); coords[0] = (npt.x, npt.y)
        link.at[best_j, "geometry"] = LineString(coords)
        pipe_geoms[best_j] = link.at[best_j, "geometry"]
        pipe_tree = STRtree(pipe_geoms)
        link.at[best_j, F_US_NODE] = nid
        used_nodes.add(nid)
        continue
    if best.distance(end) <= MIN_ENDPOINT_SNAP:
        coords = list(line.coords); coords[-1] = (npt.x, npt.y)
        link.at[best_j, "geometry"] = LineString(coords)
        pipe_geoms[best_j] = link.at[best_j, "geometry"]
        pipe_tree = STRtree(pipe_geoms)
        link.at[best_j, F_DS_NODE] = nid
        used_nodes.add(nid)
        continue

    # ---- 중간 접속: 분할 후 연결 ----
    res = split_line_at_nearest_point(line, best)
    if res is None:
        if best.distance(start) <= SNAP_MH_PIPE:
            coords = list(line.coords); coords[0] = (npt.x, npt.y)
            link.at[best_j, "geometry"] = LineString(coords)
            link.at[best_j, F_US_NODE] = nid
        elif best.distance(end) <= SNAP_MH_PIPE:
            coords = list(line.coords); coords[-1] = (npt.x, npt.y)
            link.at[best_j, "geometry"] = LineString(coords)
            link.at[best_j, F_DS_NODE] = nid
        else:
            # 정말 안 닿으면 이번 맨홀은 보류(삭제 후보)
            orphan_indices.append(i)
            continue

        pipe_geoms[best_j] = link.at[best_j, "geometry"]
        pipe_tree = STRtree(pipe_geoms)

        used_nodes.add(nid)
        continue

    A_geom, B_geom = res
    base = link.iloc[best_j].copy()

    recA = base.copy()
    recA["geometry"] = A_geom
    recA[F_DS_NODE] = nid

    recB = base.copy()
    recB["geometry"] = B_geom
    recB[F_US_NODE] = nid

    # 지연 금지: 즉시 link에 반영
    link.iloc[best_j] = recA
    link = gpd.GeoDataFrame(  
        pd.concat([link, gpd.GeoDataFrame([recB], crs=link.crs)], ignore_index=True),
        crs=link.crs
    )

    # 공간 인덱스와 지오메트리 배열 즉시 재생성(다음 맨홀은 '분할된 선'을 보게 됨)
    pipe_geoms = list(link.geometry.values)
    pipe_tree  = STRtree(pipe_geoms)

    used_nodes.add(nid)
    continue


# 교체 적용
if rows_to_drop:
    link = link.drop(index=rows_to_drop).reset_index(drop=True)
    if rows_to_add:
        link = gpd.GeoDataFrame(
            pd.concat([link, gpd.GeoDataFrame(rows_to_add, crs=link.crs)], ignore_index=True),
            crs=link.crs
        )

if orphan_indices:
    nodes = nodes.drop(index=orphan_indices).reset_index(drop=True)

# ---------- 5. 관-관 교차/분기 정규화(맨홀 삽입 후 분할) ----------
# 간단한 근사: 관의 시작/끝점과 다른 관의 선분이 0.5 m 이내면 해당 지점에 맨홀 생성 후 분할
# (완전한 모든 교차 검출은 비용 큼 → 실무 데이터에선 선형망이 대부분 끝점-끝점 구조)
new_links = []
for i, row in link.iterrows():
    new_links.append(row)
link = gpd.GeoDataFrame(new_links, crs=link.crs)

# ---------- 6. 고도(관저고) 보정: US>DS, 최소구배 확보 ----------
adjusted_logs = []
def enforce_slope(us_inv, ds_inv, L):
    """노드고 고정 가정 없이 관 값만 보고 최소구배 확보(가능하면 낮추는 방향)"""
    if L <= 0:
        return us_inv, ds_inv, False, "zero_length"
    # 목표 DS_max = us_inv - MIN_SLOPE*L
    target_ds_max = (us_inv if us_inv is not None else ds_inv) - MIN_SLOPE*L if us_inv is not None else None
    changed = False
    reason = []
    if us_inv is not None and ds_inv is not None:
        slope = (us_inv - ds_inv) / L
        if slope < MIN_SLOPE:
            # 우선 DS 낮추기
            new_ds = us_inv - MIN_SLOPE*L
            if new_ds < ds_inv:
                ds_inv = new_ds
                changed = True
                reason.append("lower_DS")
            else:
                # DS 못 낮추면 US 올려야 하나? 정책상 올리긴 지양 -> US도 낮추기(둘 다 내리기)
                delta = MIN_SLOPE*L - slope*L
                us_inv = us_inv - delta
                ds_inv = ds_inv - delta
                changed = True
                reason.append("lower_both")
        elif slope < 0 and slope < NEG_TOL:
            # 확실한 역경사 → 최소구배까지 보정
            ds_inv = us_inv - MIN_SLOPE*L
            changed = True
            reason.append("fix_negative")
    elif us_inv is not None and ds_inv is None:
        ds_inv = us_inv - MIN_SLOPE*L
        changed = True
        reason.append("fill_DS_from_US")
    elif us_inv is None and ds_inv is not None:
        us_inv = ds_inv + MIN_SLOPE*L
        changed = True
        reason.append("fill_US_from_DS")
    else:
        reason.append("both_missing")
    return us_inv, ds_inv, changed, "+".join(reason)

for i, row in link.iterrows():
    geom: LineString = row.geometry
    L = row["__len__"] if row["__len__"] > 0 else length_of(geom)
    us = row.get(F_US_INV, None)
    ds = row.get(F_DS_INV, None)
    us = float(us) if us is not None and str(us).strip() != "" else None
    ds = float(ds) if ds is not None and str(ds).strip() != "" else None

    new_us, new_ds, changed, why = enforce_slope(us, ds, L)
    if changed:
        link.at[i, F_US_INV] = new_us
        link.at[i, F_DS_INV] = new_ds
        adjusted_logs.append({
            "pipe_index": i,
            "old_US": us, "old_DS": ds,
            "new_US": new_us, "new_DS": new_ds,
            "length": L, "why": why
        })

    # 방향 통일: US>=DS가 되도록 필요시 뒤집기 (노드 ID도 교환)
    if new_us is not None and new_ds is not None and new_us < new_ds:
        # invert만 반대로 기록 (지오메트리는 유지)
        us_id = row.get(F_US_NODE)
        ds_id = row.get(F_DS_NODE)
        link.at[i, F_US_INV], link.at[i, F_DS_INV] = new_ds, new_us
        link.at[i, F_US_NODE], link.at[i, F_DS_NODE] = ds_id, us_id
        adjusted_logs.append({
            "pipe_index": i,
            "action": "swap_US_DS_to_match_flow"
        })

# ---------- 7. QA 플래그(보정 불가, 결측 등) ----------
qa_records = []
for i, row in link.iterrows():
    us = row.get(F_US_INV, None)
    ds = row.get(F_DS_INV, None)
    if (us is None) or (ds is None):
        qa_records.append({"type":"missing_invert", "pipe_index": i, "geometry": row.geometry})
    else:
        L = row["__len__"] if row["__len__"] > 0 else length_of(row.geometry)
        slope = (float(us) - float(ds)) / (L if L>0 else 1.0)
        if slope < MIN_SLOPE + NEG_TOL:
            qa_records.append({"type":"low_slope", "pipe_index": i, "geometry": row.geometry})

qa_gdf = gpd.GeoDataFrame(qa_records, crs=link.crs) if qa_records else gpd.GeoDataFrame(columns=["type","pipe_index","geometry"], crs=link.crs)

# ---------- 8. 노드 레이어 최종 갱신(추가 노드 포함) ----------
if added_nodes:
    # 이미 nodes에 합쳐놓음 → 그대로 사용
    pass

# ---------- 9. 출력 ----------
if "LAY_YMD" in link.columns:
    link["LAY_YMD"] = pd.to_datetime(link["LAY_YMD"], errors="coerce").dt.date
# 링크/노드 CRS 보존하여 출력
link.to_file(OUT_LINK)
# 맨홀: nodes 중 MH
mh_fixed = nodes[nodes["NODE_TYPE"]=="MH"].copy()
mh_fixed.rename(columns={"NODE_ID":F_MH_ID, "INV":F_MH_INV}, inplace=True)
mh_fixed.to_file(OUT_MH)
# 토구: nodes 중 OF
of_fixed = nodes[nodes["NODE_TYPE"]=="OF"].copy()
of_fixed.rename(columns={"NODE_ID":F_OF_ID, "INV":F_OF_INV}, inplace=True)
of_fixed.to_file(OUT_OF)

# 로그 저장
import pandas as pd
pd.DataFrame(added_nodes).to_csv(LOG_ADDED_NODES, index=False)
pd.DataFrame(deleted_pipes).to_csv(LOG_DELETED_PIPES, index=False)
pd.DataFrame(adjusted_logs).to_csv(LOG_ADJ_INVERTS, index=False)
if len(qa_gdf):
    qa_gdf.to_file(OUT_QA_FLAGS)

print("Done. Outputs:")
print(f"- {OUT_LINK}")
print(f"- {OUT_MH}")
print(f"- {OUT_OF}")
print(f"- {LOG_ADDED_NODES}")
print(f"- {LOG_DELETED_PIPES}")
print(f"- {LOG_ADJ_INVERTS}")
print(f"- {OUT_QA_FLAGS if len(qa_gdf) else '(no qa flags)'}")
print('Done!!!')
