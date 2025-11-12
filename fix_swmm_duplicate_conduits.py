# fix_swmm_duplicate_conduits.py
import re, sys
from collections import defaultdict

SECTIONS_ORDER = ["[CONDUITS]", "[XSECTIONS]", "[LOSSES]"]

def split_sections(text):
    # 간단 파서: 섹션 헤더로 분할
    parts = re.split(r'(\[[A-Z]+\])', text)
    head = parts[0]
    sections = {}
    order = []
    for i in range(1, len(parts), 2):
        sec = parts[i]
        body = parts[i+1]
        sections[sec] = body
        order.append(sec)
    return head, sections, order

def parse_rows(body):
    # 주석/공란 제외한 유효 라인과 인덱스 보존
    rows = []
    for i, line in enumerate(body.splitlines()):
        if line.strip().startswith(";;") or not line.strip():
            rows.append(("raw", line))  # 원형 유지
        else:
            rows.append(("data", line))
    return rows

def extract_first_token(line):
    # 공백/탭 구분 첫 토큰(링크/노드 ID)
    return re.split(r'\s+', line.strip())[0] if line.strip() else ""

def rebuild_body(rows):
    return "\n".join(l for _, l in rows)

def fix_duplicates(conduit_rows, xsec_rows, losses_rows):
    # 1) [CONDUITS]에서 중복 찾기
    seen = defaultdict(int)
    rename_map_seq = defaultdict(list)  # old -> [new for each occurrence]
    for typ, line in conduit_rows:
        if typ != "data": 
            continue
        cid = extract_first_token(line)
        seen[cid] += 1
        if seen[cid] == 1:
            rename_map_seq[cid].append(cid)  # 첫 등장은 그대로
        else:
            rename_map_seq[cid].append(f"{cid}_{seen[cid]}")  # 2번째 이상은 접미사

    # 2) [CONDUITS]에 실제 이름 대입(등장 순서 기준)
    counters = defaultdict(int)
    new_conduit_rows = []
    for typ, line in conduit_rows:
        if typ != "data":
            new_conduit_rows.append((typ, line))
            continue
        cid = extract_first_token(line)
        counters[cid] += 1
        new_id = rename_map_seq[cid][counters[cid]-1]
        if new_id != cid:
            # 첫 토큰만 치환
            new_line = line.replace(cid, new_id, 1)
        else:
            new_line = line
        new_conduit_rows.append(("data", new_line))

    # 3) [XSECTIONS], [LOSSES]도 같은 “등장 순서”로 매칭해 치환
    def apply_to_ref_rows(rows):
        counters2 = defaultdict(int)
        new_rows = []
        for typ, line in rows:
            if typ != "data":
                new_rows.append((typ, line)); continue
            rid = extract_first_token(line)  # 참조하는 링크ID
            if rid in rename_map_seq:
                counters2[rid] += 1
                new_id = rename_map_seq[rid][counters2[rid]-1]
                if new_id != rid:
                    new_line = line.replace(rid, new_id, 1)
                else:
                    new_line = line
                new_rows.append(("data", new_line))
            else:
                new_rows.append(("data", line))
        return new_rows

    new_xsec_rows   = apply_to_ref_rows(xsec_rows)   if xsec_rows   is not None else None
    new_losses_rows = apply_to_ref_rows(losses_rows) if losses_rows is not None else None

    # 4) 요약 리포트
    renamed = {k: v for k, v in rename_map_seq.items() if any(nn != k for nn in v)}
    return new_conduit_rows, new_xsec_rows, new_losses_rows, renamed

def main(inp_in, inp_out):
    text = open(inp_in, encoding="utf-8", errors="ignore").read()
    head, sections, order = split_sections(text)

    # 필요한 섹션 파싱
    conduits = parse_rows(sections.get("[CONDUITS]", ""))
    xsecs    = parse_rows(sections.get("[XSECTIONS]", "")) if "[XSECTIONS]" in sections else None
    losses   = parse_rows(sections.get("[LOSSES]", ""))    if "[LOSSES]" in sections    else None

    new_conduits, new_xsecs, new_losses, renamed = fix_duplicates(conduits, xsecs, losses)

    # 섹션 재조립
    sections["[CONDUITS]"] = rebuild_body(new_conduits)
    if new_xsecs is not None:
        sections["[XSECTIONS]"] = rebuild_body(new_xsecs)
    if new_losses is not None:
        sections["[LOSSES]"] = rebuild_body(new_losses)

    # 원래 섹션 순서 유지(없던 섹션은 스킵)
    out = [head]
    for sec in order:
        if sec in sections:
            out.append(sec)
            out.append(sections[sec])
    # 혹시 원문에 없던 섹션이 새로 생겼다면 맨 끝에 덧붙임
    for sec in sections:
        if sec not in order:
            out.append(sec); out.append(sections[sec])

    with open(inp_out, "w", encoding="utf-8") as f:
        f.write("\n".join(out))

    # 콘솔 요약
    if renamed:
        print("Renamed duplicate conduit IDs:")
        for old, seq in renamed.items():
            # 첫 항목은 원본이므로 제외
            changed = [s for s in seq if s != old]
            if changed:
                print(f"  {old} -> {', '.join(changed)}")
    else:
        print("No duplicate conduit IDs found.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fix_swmm_duplicate_conduits.py input.inp output.inp")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
