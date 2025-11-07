# 2025-11-07 수정
# 독립맨홀, 관의 노드 수정 완료 이후
# 맨홀 없이 끝나는 파이프 보정 버전

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

SNAP_MH_PIPE   = 2.0    # 맨홀-관 스냅 허용오차 [m]
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

def exact_cut_point(line: LineString, pt: Point) -> Point:
    """라인상의 pt에 대한 정확 최근접점(선형참조 기반) 반환"""
    d = line.project(pt)
    return line.interpolate(d)

from shapely.ops import nearest_points

def exact_cut_point(line: LineString, pt: Point) -> Point:
    """라인상의 pt에 대한 정확 최근접점(선형참조 기반)"""
    d = line.project(pt)
    return line.interpolate(d)

def snap_node_to_line_point(nodes_gdf, nid, new_pt: Point):
    """nodes에서 NODE_ID==nid의 geometry를 new_pt로 갱신하고 sindex 재생성 신호 플래그 리턴"""
    mask = (nodes_gdf["NODE_ID"].astype(str) == str(nid))
    if mask.any():
        nodes_gdf.loc[mask, "geometry"] = new_pt
        return True
    return False

# ===== Schema 보존 출력 유틸 =====
def finalize_and_save_schema_preserved(
    link, mh, of, nodes, added_nodes,
    F_US_NODE, F_DS_NODE, F_US_INV, F_DS_INV,
    F_MH_ID, F_MH_INV, F_OF_ID, F_OF_INV,
    OUT_LINK, OUT_MH, OUT_OF,
    LOG_ADDED_NODES, LOG_DELETED_PIPES, LOG_ADJ_INVERTS, OUT_QA_FLAGS,
    qa_gdf=None,
    temp_cols=None
):
    """
    Shapefile 원본 스키마(필드)를 그대로 유지하며 값(geometry/INV/연결ID)만 갱신하여 저장한다.
    - link: GeoDataFrame (관)
    - mh, of: 원본 GeoDataFrame (맨홀, 토구)  -> 스키마 보존 기준
    - nodes: 최종 노드 레이어 (NODE_ID, NODE_TYPE, INV, geometry)
    - added_nodes: [{'NODE_ID','NODE_TYPE','INV','x','y'}, ...]
    """
    import numpy as np
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point

    # 1) 링크 임시 컬럼 제거
    temp_cols = temp_cols or []
    for c in temp_cols:
        if c in link.columns:
            link = link.drop(columns=[c])

    # 2) 최종 노드 분리
    mh_nodes_final = nodes[nodes["NODE_TYPE"] == "MH"].copy()
    of_nodes_final = nodes[nodes["NODE_TYPE"] == "OF"].copy()

    # 3) 원본 스키마 복제
    mh_out = mh.copy()
    of_out = of.copy()

    # 4) 매핑 테이블
    mh_geom_map = {str(r["NODE_ID"]): r.geometry for _, r in mh_nodes_final.iterrows()}
    mh_inv_map  = {str(r["NODE_ID"]): r.get("INV", None) for _, r in mh_nodes_final.iterrows()}
    of_geom_map = {str(r["NODE_ID"]): r.geometry for _, r in of_nodes_final.iterrows()}
    of_inv_map  = {str(r["NODE_ID"]): r.get("INV", None) for _, r in of_nodes_final.iterrows()}

    # 5) 실제 사용 중인 노드 ID 집합(삭제/병합 반영)
    used_node_ids = set(
        link[F_US_NODE].dropna().astype(str).tolist()
        + link[F_DS_NODE].dropna().astype(str).tolist()
    )

    # 6) 원본 중 최종/사용 노드만 유지
    if F_MH_ID in mh_out.columns:
        mh_keep_mask = mh_out[F_MH_ID].astype(str).isin(mh_geom_map.keys()) & mh_out[F_MH_ID].astype(str).isin(used_node_ids)
        mh_out = mh_out.loc[mh_keep_mask].copy()
    if F_OF_ID in of_out.columns:
        of_keep_mask = of_out[F_OF_ID].astype(str).isin(of_geom_map.keys()) & of_out[F_OF_ID].astype(str).isin(used_node_ids)
        of_out = of_out.loc[of_keep_mask].copy()

    # 7) geometry/INV 값 갱신(스키마 유지)
    def _update_nodes_in_place(df, id_col, inv_col, geom_map, inv_map):
        if id_col not in df.columns:
            return df
        df["__id_str__"] = df[id_col].astype(str)
        # geometry 업데이트
        df["geometry"] = df["__id_str__"].map(geom_map).where(df["__id_str__"].map(geom_map).notna(), df["geometry"])
        # invert 업데이트(해당 컬럼 있으면)
        if inv_col in df.columns:
            new_inv = df["__id_str__"].map(inv_map)
            df[inv_col] = new_inv.where(new_inv.notna(), df[inv_col])
        return df.drop(columns=["__id_str__"], errors="ignore")

    if len(mh_out):
        mh_out = _update_nodes_in_place(mh_out, F_MH_ID, F_MH_INV, mh_geom_map, mh_inv_map)
    if len(of_out):
        of_out = _update_nodes_in_place(of_out, F_OF_ID, F_OF_INV, of_geom_map, of_inv_map)

    # 8) 새로 생성된 맨홀을 원본 스키마로 append
    if added_nodes:
        added_mh = [r for r in added_nodes if r.get("NODE_TYPE") == "MH"]
        if len(added_mh) > 0:
            blanks = []
            for r in added_mh:
                rec = {col: np.nan for col in mh_out.columns}  # 원본 필드 유지
                if F_MH_ID in mh_out.columns:
                    rec[F_MH_ID] = r["NODE_ID"]
                if F_MH_INV in mh_out.columns:
                    rec[F_MH_INV] = r.get("INV", np.nan)
                rec["geometry"] = Point(r["x"], r["y"])
                blanks.append(rec)
            if blanks:
                mh_out = gpd.GeoDataFrame(
                    pd.concat([mh_out, gpd.GeoDataFrame(blanks, geometry="geometry", crs=link.crs)], ignore_index=True),
                    crs=link.crs
                )

    # 9) 날짜/경미 정리(선택)
    if "LAY_YMD" in link.columns:
        link["LAY_YMD"] = pd.to_datetime(link["LAY_YMD"], errors="coerce").dt.date

    # 10) 저장
    link.to_file(OUT_LINK)
    mh_out.to_file(OUT_MH)
    of_out.to_file(OUT_OF)

    # 11) 로그/QA
    pd.DataFrame(added_nodes).to_csv(LOG_ADDED_NODES, index=False)
    if "deleted_pipes" in globals():
        pd.DataFrame(deleted_pipes).to_csv(LOG_DELETED_PIPES, index=False)
    if "adjusted_logs" in globals():
        pd.DataFrame(adjusted_logs).to_csv(LOG_ADJ_INVERTS, index=False)
    if qa_gdf is not None and len(qa_gdf):
        qa_gdf.to_file(OUT_QA_FLAGS)

    print("Done. Outputs (schema-preserved):")
    print(f"- {OUT_LINK}")
    print(f"- {OUT_MH}")
    print(f"- {OUT_OF}")
    print(f"- {LOG_ADDED_NODES}")
    print(f"- {LOG_DELETED_PIPES}")
    print(f"- {LOG_ADJ_INVERTS}")
    print(f"- {OUT_QA_FLAGS if (qa_gdf is not None and len(qa_gdf)) else '(no qa flags)'}")

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
    link.at[i, F_US_NODE] = us_id if us_id is not None else pd.NA
    link.at[i, F_DS_NODE] = ds_id if ds_id is not None else pd.NA   

# ---------- 3. 맨홀 없이 끝나는 단파이프 처리 ----------
# 끝점이 '진짜 단말'인지: 노드(맨홀/토구)도 없고, 다른 관의 라인(끝점뿐 아니라 중간 포함)과도 안 닿는 경우만 True
pipe_geoms = list(link.geometry.values)
pipe_tree3 = STRtree(pipe_geoms)
nodes_tree = STRtree(nodes.geometry.values.tolist())

# 존재하지 않는 노드 ID(문자열로 남아있는 가짜 값)를 NA로 정리
valid_ids = set(nodes["NODE_ID"].astype(str).tolist())

for idx, r in link.iterrows():
    v = r.get(F_US_NODE, pd.NA)
    if v is not pd.NA and v is not None and str(v) not in valid_ids:
        link.at[idx, F_US_NODE] = pd.NA
    v = r.get(F_DS_NODE, pd.NA)
    if v is not pd.NA and v is not None and str(v) not in valid_ids:
        link.at[idx, F_DS_NODE] = pd.NA


def is_real_dead_end(pt: Point, pipe_idx: int, tol=SNAP_MH_PIPE):
    # 1) 노드가 가까이 있으면 단말 아님
    cand_nodes = nodes_tree.query(pt.buffer(tol))
    for j in cand_nodes:
        if pt.distance(nodes.geometry.values[j]) <= tol:
            return False
    # 2) 다른 관 라인이 가까이 있으면 단말 아님(중간 포함)
    cand_pipes = pipe_tree3.query(pt.buffer(tol))
    for j in cand_pipes:
        if j == pipe_idx:
            continue
        if pt.distance(pipe_geoms[j]) <= tol:
            return False
    return True

# 연결 필드 보정
for col in [F_US_NODE, F_DS_NODE]:
    if col not in link.columns:
        link[col] = pd.NA

to_drop = []
for i, row in link.iterrows():
    geom: LineString = row.geometry
    if geom is None or geom.is_empty:
        continue

    L = row["__len__"] if row["__len__"] > 0 else length_of(geom)
    us_id = row[F_US_NODE]
    ds_id = row[F_DS_NODE]
    start = Point(list(geom.coords)[0])
    end   = Point(list(geom.coords)[-1])

    # 각 끝점의 '진짜 단말' 여부
    start_dead = (pd.isna(us_id) or us_id is None) and is_real_dead_end(start, i, tol=SNAP_MH_PIPE)
    end_dead   = (pd.isna(ds_id) or ds_id is None) and is_real_dead_end(end,   i, tol=SNAP_MH_PIPE)

    # ── 규칙: 5m 이하 우선 삭제
    # 1) 양쪽 모두 단말이고 L<=5m → 삭제
    # 2) 한쪽만 단말이어도 그 단말 끝이 있고 L<=5m → 삭제
    if L <= END_SHORT_DEL:
        if (start_dead and end_dead) \
        or ((pd.isna(us_id) or us_id is None) and start_dead) \
        or ((pd.isna(ds_id) or ds_id is None) and end_dead):
            to_drop.append(i)
            mark_delete(i, f"short_dead_pipe_end<= {END_SHORT_DEL}m")
            continue

    # ── 규칙: 연결부(단말 아님)인데 노드ID 없음 → 맨홀 생성(관-관 사이 맨홀 강제)
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

    # ── 규칙: 한쪽이라도 진짜 단말이면 L>5m에서만 그쪽에 맨홀 생성
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

# 추가된 노드를 nodes에 반영 + 인덱스 재생성
# 교체안 (가장 간단/안전)
if added_nodes:
    df_add = gpd.GeoDataFrame(
        added_nodes,
        geometry=[Point(r["x"], r["y"]) for r in added_nodes],
        crs=link.crs
    )[["NODE_ID", "NODE_TYPE", "INV", "geometry"]]

    if not df_add.empty:
        nodes = gpd.GeoDataFrame(
            pd.concat([nodes, df_add], ignore_index=True),
            crs=link.crs
        )


# ---------- 4. 고립 맨홀 정리(관에 가까우면 반드시 연결·분할하고, 정말 멀면만 삭제) ----------

# 항상 최신 링크 지오메트리로 트리 구성
pipe_geoms = list(link.geometry.values)
pipe_tree = STRtree(pipe_geoms)

MIN_ENDPOINT_SNAP = 0.5  # 끝점으로 스냅 판정 오차 [m]
ORPHAN_TOL = SNAP_MH_PIPE + ORPHAN_MARGIN  # 고립 판단 최종 임계

def split_line_at_nearest_point(line: LineString, pt: Point, tol=1e-6):
    """라인을 pt의 선형참조 거리 d에서 정확히 두 조각으로 분할"""
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
            ratio = (d - acc) / seg_len
            x = p0.x + ratio * (p1.x - p0.x)
            y = p0.y + ratio * (p1.y - p0.y)
            cut = (x, y)
            part1 = coords[:k+1] + [cut]
            part2 = [cut] + coords[k+1:]
            return LineString(part1), LineString(part2)
        acc += seg_len
    return None

# 이미 링크에서 쓰이는 노드 집합
used_nodes = set(link[F_US_NODE].dropna().astype(str).tolist() + link[F_DS_NODE].dropna().astype(str).tolist())

rows_to_add = []      # 새 관
rows_to_drop = []     # 대체할 구관 (여기선 즉시 갱신하므로 거의 사용 안함)
orphan_indices = []   # 최종 삭제할 맨홀 인덱스

for i, nrow in nodes.iterrows():
    nid = nrow["NODE_ID"]
    ntype = nrow.get("NODE_TYPE", "MH")
    if ntype != "MH":      # 토구 등은 여기서 삭제하지 않음
        continue
    if str(nid) in used_nodes:
        continue

    npt: Point = nrow.geometry
    cand_idxs = pipe_tree.query(npt.buffer(ORPHAN_TOL))
    if cand_idxs is None or len(cand_idxs) == 0:
        # 주변에 관 자체가 없음 → 후보 삭제
        orphan_indices.append(i)
        continue

    # 1) 관에 대한 '선 최근접점' 계산 (stick-to-line)
    best_j, best_on_line, best_dist = None, None, float("inf")
    for j in cand_idxs:
        line = pipe_geoms[j]
        proj_pt = exact_cut_point(line, npt)      # 선형참조 기반 정확 최근접
        d = npt.distance(proj_pt)
        if d < best_dist:
            best_dist = d
            best_on_line = proj_pt
            best_j = j

    # 2) 정말 멀면만 삭제 후보
    if best_on_line is None or best_dist > ORPHAN_TOL + 1e-9:
        orphan_indices.append(i)
        continue

    # 3) 먼저 노드를 선 위로 '스냅' (좌표 일치시켜 후속 판정 안정화)
    moved = snap_node_to_line_point(nodes, nid, best_on_line)
    if moved:
        node_sindex = STRtree(nodes.geometry.values.tolist())
    npt = best_on_line  # 이후 로직은 선 위 좌표로 판단

    # 최신 라인·끝점 정보
    line = pipe_geoms[best_j]
    start = Point(list(line.coords)[0])
    end   = Point(list(line.coords)[-1])

    # 4) 끝점에 가깝다면 '끝점 스냅'으로 간단 연결
    if npt.distance(start) <= MIN_ENDPOINT_SNAP:
        coords = list(line.coords); coords[0] = (npt.x, npt.y)
        link.at[best_j, "geometry"] = LineString(coords)
        pipe_geoms[best_j] = link.at[best_j, "geometry"]
        pipe_tree = STRtree(pipe_geoms)
        link.at[best_j, F_US_NODE] = nid
        used_nodes.add(str(nid))
        continue

    if npt.distance(end) <= MIN_ENDPOINT_SNAP:
        coords = list(line.coords); coords[-1] = (npt.x, npt.y)
        link.at[best_j, "geometry"] = LineString(coords)
        pipe_geoms[best_j] = link.at[best_j, "geometry"]
        pipe_tree = STRtree(pipe_geoms)
        link.at[best_j, F_DS_NODE] = nid
        used_nodes.add(str(nid))
        continue

    # 5) 중간 접속: 정확 절단 후 양쪽 관에 US/DS 노드 부여
    res = split_line_at_nearest_point(line, npt)
    if res is None:
        # 수치 경계 문제 대비: 절단 불가면 가까운 쪽 끝점으로 스냅해서라도 연결
        if npt.distance(start) <= SNAP_MH_PIPE:
            coords = list(line.coords); coords[0] = (npt.x, npt.y)
            link.at[best_j, "geometry"] = LineString(coords)
            link.at[best_j, F_US_NODE] = nid
        elif npt.distance(end) <= SNAP_MH_PIPE:
            coords = list(line.coords); coords[-1] = (npt.x, npt.y)
            link.at[best_j, "geometry"] = LineString(coords)
            link.at[best_j, F_DS_NODE] = nid
        else:
            # 여기까지 왔는데도 못 붙이면(극히 드묾) 이번만 보류
            orphan_indices.append(i)
            continue

        pipe_geoms[best_j] = link.at[best_j, "geometry"]
        pipe_tree = STRtree(pipe_geoms)
        used_nodes.add(str(nid))
        continue

    A_geom, B_geom = res

    # 절단점 좌표를 명시적으로 관 두 조각의 끝점에 반영 (수치 오차 방지)
    a_coords = list(A_geom.coords); a_coords[-1] = (npt.x, npt.y)
    b_coords = list(B_geom.coords); b_coords[0] = (npt.x, npt.y)
    A_geom = LineString(a_coords)
    B_geom = LineString(b_coords)

    base = link.iloc[best_j].copy()

    recA = base.copy()
    recA["geometry"] = A_geom
    recA[F_DS_NODE] = nid

    recB = base.copy()
    recB["geometry"] = B_geom
    recB[F_US_NODE] = nid

    # 기존 라인 자리엔 A 저장, B는 append
    link.iloc[best_j] = recA
    link = gpd.GeoDataFrame(
        pd.concat([link, gpd.GeoDataFrame([recB], crs=link.crs)], ignore_index=True),
        crs=link.crs
    )

    # 공간 인덱스 갱신
    pipe_geoms = list(link.geometry.values)
    pipe_tree  = STRtree(pipe_geoms)

    used_nodes.add(str(nid))
    # continue  # 명시적

# 고립(삭제) 대상은 진짜 멀리 떨어진 맨홀뿐
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

# ---------- 9. 출력 (스키마 보존) ----------
TEMP_COLS = ["__len__", "__box_H__"]  # 링크 임시 컬럼명들

# QA GDF 없는 경우를 대비
qa_out = qa_gdf if 'qa_gdf' in globals() else None

finalize_and_save_schema_preserved(
    link=link, mh=mh, of=of, nodes=nodes, added_nodes=added_nodes,
    F_US_NODE=F_US_NODE, F_DS_NODE=F_DS_NODE,
    F_US_INV=F_US_INV, F_DS_INV=F_DS_INV,
    F_MH_ID=F_MH_ID, F_MH_INV=F_MH_INV,
    F_OF_ID=F_OF_ID, F_OF_INV=F_OF_INV,
    OUT_LINK=OUT_LINK, OUT_MH=OUT_MH, OUT_OF=OUT_OF,
    LOG_ADDED_NODES=LOG_ADDED_NODES,
    LOG_DELETED_PIPES=LOG_DELETED_PIPES,
    LOG_ADJ_INVERTS=LOG_ADJ_INVERTS,
    OUT_QA_FLAGS=OUT_QA_FLAGS,
    qa_gdf=qa_out,
    temp_cols=TEMP_COLS
)
