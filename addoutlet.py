# 수정일자: 2025-11-11
# 작성자: 정규원

# inp_thiessen_repartition.py
# 기능: 기존 SWMM .inp에서 맨홀(JUNCTIONS) 기준 티센(보로노이)로 유역 분할,
#       분할 유역의 Outlet을 해당 맨홀로 지정하여 [SUBCATCHMENTS], [POLYGONS] 섹션 갱신
#
# 필요: geopandas, shapely>=2.0, pandas, numpy
# 사용법:
#   python inp_thiessen_repartition.py input.inp output.inp

import sys
import re
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPolygon
from shapely.ops import voronoi_diagram

SECTION_RE = re.compile(r'^\s*\[(\w+)\]\s*$', re.IGNORECASE)

def read_inp_lines(p: Path) -> list[str]:
    return p.read_text(encoding='utf-8', errors='ignore').splitlines()

def write_inp_lines(p: Path, lines: list[str]):
    p.write_text("\n".join(lines) + "\n", encoding='utf-8')

def find_sections(lines: list[str]) -> dict:
    """각 섹션 시작/끝 인덱스 탐색"""
    idxs = []
    for i, s in enumerate(lines):
        m = SECTION_RE.match(s)
        if m:
            idxs.append((i, m.group(1).upper()))
    idxs.append((len(lines), "__EOF__"))
    sections = {}
    for (i, name), (j, _) in zip(idxs[:-1], idxs[1:]):
        sections[name] = (i, j)  # [i, j)
    return sections

def slice_section(lines, sections, name):
    name = name.upper()
    if name not in sections:
        return []
    i, j = sections[name]
    return lines[i:j]

def replace_section(lines, sections, name, new_block_lines):
    """섹션 통째로 교체(없으면 삽입)"""
    name = name.upper()
    header = f"[{name}]"
    if name in sections:
        i, j = sections[name]
        pre = lines[:i]
        post = lines[j:]
        blk = [header] + new_block_lines
        return pre + blk + post
    else:
        # 맨 끝에 추가
        if lines and lines[-1].strip():
            lines = lines + [""]
        return lines + [header] + new_block_lines

# ---------- 파서들 ----------
def parse_coordinates(coord_lines: list[str]) -> pd.DataFrame:
    pts = []
    for s in coord_lines:
        if not s.strip() or s.strip().startswith(";;") or s.strip().startswith("["):
            continue
        parts = re.split(r'\s+', s.strip())
        if len(parts) >= 3:
            node = parts[0]
            try:
                x = float(parts[1]); y = float(parts[2])
            except:
                continue
            pts.append((node, x, y))
    return pd.DataFrame(pts, columns=["node", "x", "y"])

def parse_simple_table(lines: list[str], expect_cols: int = 0) -> list[list[str]]:
    rows = []
    for s in lines:
        t = s.strip()
        if not t or t.startswith(";;") or t.startswith("["):
            continue
        parts = re.split(r'\s+', t)
        if expect_cols and len(parts) < expect_cols:
            continue
        rows.append(parts)
    return rows

def parse_junctions(junc_lines: list[str]) -> set[str]:
    rows = parse_simple_table(junc_lines, expect_cols=2)
    return set(r[0] for r in rows)

def parse_outfalls(out_lines: list[str]) -> set[str]:
    rows = parse_simple_table(out_lines, expect_cols=2)
    return set(r[0] for r in rows)

def parse_subcatchments(sc_lines: list[str]) -> pd.DataFrame:
    # Name RainGage Outlet Area %Imperv Width Slope CurbLen SnowPack
    rows = []
    for s in sc_lines:
        t = s.strip()
        if not t or t.startswith(";;") or t.startswith("["):
            continue
        parts = re.split(r'\s+', t)
        if len(parts) < 7:
            continue
        name = parts[0]
        gage = parts[1]
        outlet = parts[2]
        area = float(parts[3])
        pimp = float(parts[4])
        width = float(parts[5])
        slope = float(parts[6])
        curb = float(parts[7]) if len(parts) > 7 else 0.0
        snow = parts[8] if len(parts) > 8 else ""
        rows.append((name, gage, outlet, area, pimp, width, slope, curb, snow))
    return pd.DataFrame(rows, columns=["Name","RainGage","Outlet","Area_ha","Pimp_pct","Width_m","Slope","CurbLen","SnowPack"])

def parse_polygons(poly_lines: list[str]) -> dict[str, list[tuple[float,float]]]:
    d = {}
    for s in poly_lines:
        t = s.strip()
        if not t or t.startswith(";;") or t.startswith("["):
            continue
        parts = re.split(r'\s+', t)
        if len(parts) < 3:
            continue
        nm = parts[0]
        try:
            x = float(parts[1]); y = float(parts[2])
        except:
            continue
        d.setdefault(nm, []).append((x, y))
    return d

def build_polygon(coords: list[tuple[float,float]]) -> Polygon | None:
    if not coords or len(coords) < 3:
        return None
    # POLYGONS 섹션은 보통 폐합점 반복 포함 → 중복은 제거
    if coords[0] == coords[-1]:
        coords = coords[:-1]
    try:
        poly = Polygon(coords)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty:
            return None
        return poly
    except Exception:
        return None

# ---------- Thiessen(보로노이) 생성 ----------
def make_mask_polygon(sc_polys: list[Polygon]) -> Polygon:
    if not sc_polys:
        return None
    u = sc_polys[0]
    for g in sc_polys[1:]:
        try:
            u = u.union(g)
        except:
            u = u.buffer(0).union(g.buffer(0))
    return u.buffer(0)

def voronoi_cells_for_points(gdf_pts: gpd.GeoDataFrame, mask_poly: Polygon) -> gpd.GeoDataFrame:
    try:
        mpt = gdf_pts.union_all()   # shapely 2.x
    except AttributeError:
        mpt = gdf_pts.unary_union   # shapely 1.x 호환
    v = voronoi_diagram(mpt, envelope=mask_poly.envelope, tolerance=0.0, edges=False)
    # 각 포인트가 속한 셀 찾기 (간단/안전 우선)
    cells = []
    mids = []
    geoms = list(v.geoms) if hasattr(v, "geoms") else [v]
    for mid, p in zip(gdf_pts["node"], gdf_pts.geometry):
        chosen = None
        for face in geoms:
            try:
                if face.contains(p) or face.touches(p):
                    chosen = face
                    break
            except:
                continue
        if chosen is None:
            # 거리 최소 셀 폴백
            dmin = 1e30; best = None
            for face in geoms:
                try:
                    d = face.distance(p)
                    if d < dmin:
                        dmin = d; best = face
                except:
                    continue
            chosen = best
        if chosen is None:
            continue
        clipped = chosen.intersection(mask_poly)
        if clipped.is_empty:
            continue
        if isinstance(clipped, MultiPolygon):
            # 가장 큰 면만
            clipped = max(list(clipped.geoms), key=lambda g: g.area)
        cells.append(clipped)
        mids.append(mid)
    return gpd.GeoDataFrame({"node": mids}, geometry=cells, crs=gdf_pts.crs)

# ---------- SUBCATCHMENTS/POLYGONS 재작성 ----------
def format_subcatchments(df: pd.DataFrame) -> list[str]:
    out = []
    out.append(";;Name         RainGage   Outlet        Area    %Imperv  Width   Slope   CurbLen   SnowPack")
    for _, r in df.iterrows():
        out.append(f"{str(r.Name):<13} {str(r.RainGage):<10} {str(r.Outlet):<12} "
                   f"{float(r.Area_ha):<7.3f} {float(r.Pimp_pct):<8.1f} {float(r.Width_m):<7.1f} "
                   f"{float(r.Slope):<7.4f} {float(r.CurbLen):<8.1f} {str(r.SnowPack)}")
    out.append("")  # 섹션 끝 공백라인
    return out

def format_polygons(map_name_to_poly: dict[str, Polygon]) -> list[str]:
    out = []
    out.append(";;Subcatchment   X-Coord       Y-Coord")
    for name, poly in map_name_to_poly.items():
        if poly.is_empty:
            continue
        g = poly
        if isinstance(g, MultiPolygon):
            # 가장 큰 면
            g = max(list(g.geoms), key=lambda x: x.area)
        coords = list(g.exterior.coords)
        # 폐합점 제거
        if len(coords) > 1 and coords[0] == coords[-1]:
            coords = coords[:-1]
        for x, y in coords:
            out.append(f"{name:<14} {x:<13.3f} {y:<13.3f}")
    out.append("")
    return out

def main(inp_in: Path, inp_out: Path, min_area_m2: float = 50.0):
    lines = read_inp_lines(inp_in)
    sections = find_sections(lines)

    # 필요한 섹션 읽기
    coords_blk = slice_section(lines, sections, "COORDINATES")
    junc_blk   = slice_section(lines, sections, "JUNCTIONS")
    outfall_blk= slice_section(lines, sections, "OUTFALLS")
    sc_blk     = slice_section(lines, sections, "SUBCATCHMENTS")
    poly_blk   = slice_section(lines, sections, "POLYGONS")

    if not coords_blk or not junc_blk or not sc_blk or not poly_blk:
        raise RuntimeError("필수 섹션([COORDINATES], [JUNCTIONS], [SUBCATCHMENTS], [POLYGONS]) 중 일부가 없습니다.")

    # 파싱
    df_coords = parse_coordinates(coords_blk)
    junc_ids  = parse_junctions(junc_blk)
    out_ids   = parse_outfalls(outfall_blk) if outfall_blk else set()
    df_sc     = parse_subcatchments(sc_blk)

    poly_map_raw = parse_polygons(poly_blk)
    sc_polys = {}
    for nm, coords in poly_map_raw.items():
        poly = build_polygon(coords)
        if poly is not None:
            sc_polys[nm] = poly

    if df_coords.empty or not junc_ids:
        raise RuntimeError("좌표 또는 맨홀(JUNCTIONS) 정보가 비어 있습니다.")

    # 맨홀 포인트 준비(Outfall은 제외: 맨홀로 보지 않음)
    df_mh = df_coords[df_coords["node"].isin(junc_ids)].copy()
    if not df_mh.empty and not out_ids:
        pass
    elif not df_mh.empty and out_ids:
        df_mh = df_mh[~df_mh["node"].isin(out_ids)].copy()

    if df_mh.empty:
        raise RuntimeError("맨홀 포인트가 없습니다. JUNCTIONS에서 Outfall 제외 후 남는 노드가 있어야 합니다.")

    # 좌표계 가정: INP의 [COORDINATES]는 평면(m) 좌표라고 가정
    gdf_mh = gpd.GeoDataFrame(df_mh.copy(), geometry=[Point(xy) for xy in zip(df_mh.x, df_mh.y)], crs="EPSG:5186")
    # 유역 폴리곤
    sc_recs = []
    for _, r in df_sc.iterrows():
        nm = r["Name"]
        if nm in sc_polys:
            sc_recs.append((nm, sc_polys[nm], r))
    if not sc_recs:
        raise RuntimeError("POLYGONS와 SUBCATCHMENTS의 이름 매칭 결과가 비었습니다.")

    sc_gdf = gpd.GeoDataFrame(
        [{"Name": nm, "geometry": poly} for nm, poly, _ in sc_recs],
        crs="EPSG:5186"
    )

    # 보로노이 생성(마스크: 모든 유역의 union)
    mask_poly = make_mask_polygon(list(sc_gdf.geometry.values))
    v_gdf = voronoi_cells_for_points(gdf_mh, mask_poly)

    # 유역 × 보로노이 교차 → 분할 유역
    parts = gpd.overlay(
        sc_gdf[["Name", "geometry"]],
        v_gdf[["node", "geometry"]],
        how="intersection",
        keep_geom_type=True
    )
    if parts.empty:
        raise RuntimeError("보로노이와 유역 교차 결과가 비었습니다. 좌표 범위 또는 POLYGONS를 확인하세요.")

    parts["area_m2"] = parts.geometry.area
    parts = parts[parts["area_m2"] > float(min_area_m2)].copy()
    if parts.empty:
        raise RuntimeError("슬리버 제거 후 남는 유역이 없습니다. min_area_m2를 낮춰보세요.")
    
    # 교차·면적필터 후
    parts["Name_src"] = parts["Name"]

    # 원유역 속성 상속
    df_sc_idx = df_sc.set_index("Name")
    attrs = []
    for _, pr in parts.iterrows():
        src = df_sc_idx.loc[pr["Name"]]
        attrs.append({
            "RainGage": src["RainGage"],
            "Pimp_pct": float(src["Pimp_pct"]),
            "Slope": float(src["Slope"]),
            "CurbLen": float(src["CurbLen"]),
            "SnowPack": str(src["SnowPack"]),
        })
    attr_df = pd.DataFrame(attrs, index=parts.index)
    parts = pd.concat([parts, attr_df], axis=1)

    # 혹시라도 다른 이유로 중복된 컬럼이 생겼다면 정리
    parts = parts.loc[:, ~parts.columns.duplicated()].copy()


    # 새 유역 이름 = 원유역_맨홀
    def clean_id(s):
        z = re.sub(r"[^A-Za-z0-9_\-:]", "_", str(s))
        return z[:60]  # 너무 길어지는 것 방지

    parts["NewName"] = (parts["Name"].astype(str) + "_" + parts["node"].astype(str)).map(clean_id)
    parts["Outlet"]  = parts["node"].astype(str)

    # 면적(ha) 및 폭 추정(W = 2 * sqrt(A))
    parts["Area_ha"] = parts["area_m2"] / 10000.0
    parts["Width_m"] = 2.0 * np.sqrt(parts["area_m2"])
    parts.loc[parts["Width_m"] < 10.0, "Width_m"] = 10.0

    # 최종 SUBCATCHMENTS 테이블 준비(정렬: 원유역 → 맨홀)
    subc_out = parts[[
    "Name_src", "NewName", "RainGage", "Outlet", "Area_ha",
    "Pimp_pct", "Width_m", "Slope", "CurbLen", "SnowPack", "geometry"
    ]].copy()

    # 혹시 중복 컬럼이 남아있다면 제거
    subc_out = subc_out.loc[:, ~subc_out.columns.duplicated()].copy()

    _sort_cols = [c for c in ["Name_src", "Outlet", "NewName"] if c in subc_out.columns]
    if _sort_cols:
        subc_out = subc_out.sort_values(_sort_cols).reset_index(drop=True)

    subc_out.rename(columns={"NewName": "Name"}, inplace=True)

    poly_out_map = {row.Name: row.geometry for _, row in subc_out.iterrows()}


    # SUBCATCHMENTS 본문 생성 시 Name_src 제거
    new_sc_block = format_subcatchments(subc_out.drop(columns=["geometry", "Name_src"]))
    new_poly_block = format_polygons(poly_out_map)

    # 기존 섹션 헤더를 포함하지 않고 본문만 교체하므로 slice로 구간 제거 후 삽입
    def replace_body(lines, sections, sec_name, body_lines):
        sec_name = sec_name.upper()
        header = f"[{sec_name}]"
        if sec_name not in sections:
            # 없으면 추가
            return replace_section(lines, sections, sec_name, body_lines)
        i, j = sections[sec_name]
        pre = lines[:i+1]  # 헤더까지 포함
        post = lines[j:]
        # 헤더 다음에 본문 넣고, 마지막에 공백 1줄 보장
        blk = pre + body_lines + post
        return blk

    lines = replace_body(lines, sections, "SUBCATCHMENTS", new_sc_block)
    # sections 맵 갱신 필요(라인 수 변동)
    sections = find_sections(lines)
    lines = replace_body(lines, sections, "POLYGONS", new_poly_block)

    write_inp_lines(inp_out, lines)
    print(f"✅ 작성 완료: {inp_out}")
    print(f"  - 분할 유역 수: {len(subc_out)}")
    print(f"  - 사용한 맨홀 수(보로노이 seed): {len(gdf_mh)}")
    print(f"  - min_area_m2 = {min_area_m2}")

if __name__ == "__main__":
    # 인자로 주면 그걸 쓰고, 없으면 D:\ 경로 기본값 사용
    if len(sys.argv) >= 3:
        inp_in = Path(r"D:\Kyuwon\newshp\munrae_model.inp")
        inp_out = Path(r"D:\Kyuwon\newshp\munrae_model_thiessen.inp")
    else:
        inp_in = Path(r"D:\Kyuwon\newshp\munrae_model.inp")
        inp_out = Path(r"D:\Kyuwon\newshp\munrae_model_thiessen.inp")

    # 출력 폴더가 없으면 생성
    inp_out.parent.mkdir(parents=True, exist_ok=True)

    main(inp_in, inp_out, min_area_m2=50.0)

