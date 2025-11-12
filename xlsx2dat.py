# 수정일자: 2025-11-08
# 작성자: 정규원

# .xlsx형태로 저장된 강우량 데이터를 SWMM에서 import할 수 있는 .dat파일로 변환하는 코드


import pandas as pd

# === 설정 ===
infile   = r"D:\Kyuwon\pypy\munllae\영등포구청2022.xlsx"
sheet    = 0
gage_col = "강우량계명"
time_col = "시간"
rain_col = "10분우량"

use_date_only = None               # 예: "2022-08-08" 또는 None
outfile = "rain_10min.dat"

# === 읽기 ===
df = pd.read_excel(infile, sheet_name=sheet, dtype={gage_col: str})
df = df.rename(columns=lambda c: str(c).strip())
df[time_col] = pd.to_datetime(df[time_col])
df[rain_col] = pd.to_numeric(df[rain_col], errors="coerce").fillna(0)

# 날짜 필터
if use_date_only:
    df = df[df[time_col].dt.strftime("%Y-%m-%d") == use_date_only]

# 정렬, 중복 제거
df = df.sort_values(time_col).drop_duplicates(subset=[time_col], keep="first")

# 10분 간격 보강
if len(df) > 1:
    full = pd.DataFrame({time_col: pd.date_range(df[time_col].min(),
                                                 df[time_col].max(),
                                                 freq="10min")})
    df = full.merge(df[[time_col, rain_col]], on=time_col, how="left")
    df[rain_col] = df[rain_col].fillna(0)

# === 저장 ===
import io, os
outfile = r"D:\Kyuwon\newshp\rain_10min.dat"

with open(outfile, "w", encoding="utf-8", newline="\r\n") as f:
    for _, r in df.dropna(subset=[time_col, rain_col]).iterrows():
        t = r[time_col]
        line = f"{t:%m/%d/%Y} {t:%H:%M} {float(r[rain_col]):.2f}"
        f.write(line + "\r")

# 탭 제거 + 파일 선두 공백/빈줄 제거
with open(outfile, "r", encoding="utf-8") as f:
    txt = f.read().replace("\t", " ")
with open(outfile, "w", encoding="utf-8", newline="\r\n") as f:
    f.write(txt.lstrip())


