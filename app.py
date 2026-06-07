import math
import os
import socket
import requests
from flask import Flask, request, jsonify, send_from_directory

app = Flask(__name__, static_folder="static")

OPINET_URL = "https://www.opinet.co.kr/api/aroundAll.do"
OPINET_MAX = 5000  # 오피넷 최대 반경(m)

# ══════════════════════════════════════════════════════════
#  WGS84 ↔ KATEC 좌표 변환 (pyproj 없이 순수 Python 구현)
#  오피넷은 KATEC(Bessel 타원체, 중앙자오선 128°) 사용
# ══════════════════════════════════════════════════════════

# WGS84 타원체
_WA  = 6378137.0
_WF  = 1 / 298.257223563
_WE2 = 2 * _WF - _WF ** 2

# Bessel 타원체
_BA  = 6377397.155
_BF  = 1 / 299.1528128
_BE2 = 2 * _BF - _BF ** 2

# KATEC TM 파라미터
_K0  = 1.0
_LON0 = math.radians(128.0)   # 중앙자오선
_LAT0 = math.radians(38.0)    # 원점 위도
_FE  = 400_000.0               # False Easting
_FN  = 600_000.0               # False Northing

# Helmert 7파라미터 (Bessel → WGS84 방향, PROJ towgs84 값)
# WGS84 → Bessel 변환 시 부호 반전해서 사용
_DX, _DY, _DZ = -115.80, 474.99, 674.11
_RX = math.radians(1.16  / 3600)
_RY = math.radians(-2.31 / 3600)
_RZ = math.radians(-1.63 / 3600)
_DS = 6.43e-6


def _wgs84_to_ecef(lat, lng):
    N = _WA / math.sqrt(1 - _WE2 * math.sin(lat) ** 2)
    x = N * math.cos(lat) * math.cos(lng)
    y = N * math.cos(lat) * math.sin(lng)
    z = N * (1 - _WE2) * math.sin(lat)
    return x, y, z


def _helmert_wgs_to_bessel(X, Y, Z):
    """WGS84 ECEF → Bessel ECEF (역 Helmert, 소각도 근사)"""
    dx, dy, dz = -_DX, -_DY, -_DZ
    rx, ry, rz = -_RX, -_RY, -_RZ
    s = -_DS
    Xb = dx + (1 + s) * ( X + rz * Y - ry * Z)
    Yb = dy + (1 + s) * (-rz * X +  Y + rx * Z)
    Zb = dz + (1 + s) * ( ry * X - rx * Y +  Z)
    return Xb, Yb, Zb


def _helmert_bessel_to_wgs(Xb, Yb, Zb):
    """Bessel ECEF → WGS84 ECEF"""
    Xw = _DX + (1 + _DS) * ( Xb + _RZ * Yb - _RY * Zb)
    Yw = _DY + (1 + _DS) * (-_RZ * Xb +  Yb + _RX * Zb)
    Zw = _DZ + (1 + _DS) * ( _RY * Xb - _RX * Yb +  Zb)
    return Xw, Yw, Zw


def _ecef_to_bessel_geo(Xb, Yb, Zb):
    """Bessel ECEF → Bessel 위경도(rad)"""
    p = math.sqrt(Xb ** 2 + Yb ** 2)
    lng = math.atan2(Yb, Xb)
    lat = math.atan2(Zb, p * (1 - _BE2))
    for _ in range(10):
        N = _BA / math.sqrt(1 - _BE2 * math.sin(lat) ** 2)
        lat = math.atan2(Zb + _BE2 * N * math.sin(lat), p)
    return lat, lng


def _meridional_arc(a, e2, lat):
    e4 = e2 ** 2
    e6 = e2 ** 3
    return a * (
        (1 - e2/4 - 3*e4/64 - 5*e6/256) * lat
        - (3*e2/8 + 3*e4/32 + 45*e6/1024) * math.sin(2*lat)
        + (15*e4/256 + 45*e6/1024) * math.sin(4*lat)
        - (35*e6/3072) * math.sin(6*lat)
    )


def _bessel_geo_to_katec(lat, lng):
    """Bessel 위경도(rad) → KATEC (m)"""
    e2 = _BE2
    a  = _BA
    N  = a / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    T  = math.tan(lat) ** 2
    C  = (e2 / (1 - e2)) * math.cos(lat) ** 2
    A  = (lng - _LON0) * math.cos(lat)

    M  = _meridional_arc(a, e2, lat)
    M0 = _meridional_arc(a, e2, _LAT0)

    kx = _K0 * N * (
        A
        + (1 - T + C) * A**3 / 6
        + (5 - 18*T + T**2 + 72*C - 58*(e2/(1-e2))) * A**5 / 120
    ) + _FE

    ky = _K0 * (
        M - M0
        + N * math.tan(lat) * (
            A**2 / 2
            + (5 - T + 9*C + 4*C**2) * A**4 / 24
            + (61 - 58*T + T**2 + 600*C - 330*(e2/(1-e2))) * A**6 / 720
        )
    ) + _FN

    return kx, ky


def _katec_to_bessel_geo(kx, ky):
    """KATEC (m) → Bessel 위경도(rad)"""
    a  = _BA
    e2 = _BE2
    e4 = e2 ** 2
    e6 = e2 ** 3

    M0 = _meridional_arc(a, e2, _LAT0)
    M1 = M0 + (ky - _FN) / _K0

    mu = M1 / (a * (1 - e2/4 - 3*e4/64 - 5*e6/256))
    e1 = (1 - math.sqrt(1 - e2)) / (1 + math.sqrt(1 - e2))

    phi1 = (
        mu
        + (3*e1/2 - 27*e1**3/32) * math.sin(2*mu)
        + (21*e1**2/16 - 55*e1**4/32) * math.sin(4*mu)
        + (151*e1**3/96) * math.sin(6*mu)
        + (1097*e1**4/512) * math.sin(8*mu)
    )

    N1 = a / math.sqrt(1 - e2 * math.sin(phi1) ** 2)
    R1 = a * (1 - e2) / (1 - e2 * math.sin(phi1) ** 2) ** 1.5
    T1 = math.tan(phi1) ** 2
    C1 = (e2 / (1 - e2)) * math.cos(phi1) ** 2
    D  = (kx - _FE) / (N1 * _K0)

    lat = phi1 - (N1 * math.tan(phi1) / R1) * (
        D**2 / 2
        - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*(e2/(1-e2))) * D**4 / 24
        + (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*(e2/(1-e2)) - 3*C1**2) * D**6 / 720
    )
    lng = _LON0 + (
        D
        - (1 + 2*T1 + C1) * D**3 / 6
        + (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*(e2/(1-e2)) + 24*T1**2) * D**5 / 120
    ) / math.cos(phi1)

    return lat, lng


def wgs84_to_katec(lat_deg, lng_deg):
    lat = math.radians(lat_deg)
    lng = math.radians(lng_deg)
    X, Y, Z    = _wgs84_to_ecef(lat, lng)
    Xb, Yb, Zb = _helmert_wgs_to_bessel(X, Y, Z)
    blat, blng  = _ecef_to_bessel_geo(Xb, Yb, Zb)
    return _bessel_geo_to_katec(blat, blng)


def katec_to_wgs84(kx, ky):
    blat, blng  = _katec_to_bessel_geo(kx, ky)
    Xb = _BA / math.sqrt(1 - _BE2 * math.sin(blat)**2) * math.cos(blat) * math.cos(blng)
    Yb = _BA / math.sqrt(1 - _BE2 * math.sin(blat)**2) * math.cos(blat) * math.sin(blng)
    Zb = _BA / math.sqrt(1 - _BE2 * math.sin(blat)**2) * (1 - _BE2) * math.sin(blat)
    Xw, Yw, Zw = _helmert_bessel_to_wgs(Xb, Yb, Zb)
    # ECEF → WGS84 geographic
    p   = math.sqrt(Xw**2 + Yw**2)
    lng = math.atan2(Yw, Xw)
    lat = math.atan2(Zw, p * (1 - _WE2))
    for _ in range(10):
        N   = _WA / math.sqrt(1 - _WE2 * math.sin(lat)**2)
        lat = math.atan2(Zw + _WE2 * N * math.sin(lat), p)
    return math.degrees(lat), math.degrees(lng)


# ══════════════════════════════════════════════════════════
#  Opinet API 호출
# ══════════════════════════════════════════════════════════
def opinet_fetch(code, kx, ky, radius, prodcd):
    resp = requests.get(
        OPINET_URL,
        params={
            "code": code, "x": round(kx, 2), "y": round(ky, 2),
            "radius": int(radius), "prodcd": prodcd,
            "sort": 2, "out": "json",
        },
        timeout=10,
    )
    resp.raise_for_status()
    return resp.json().get("RESULT", {}).get("OIL", [])


# ══════════════════════════════════════════════════════════
#  Flask Routes
# ══════════════════════════════════════════════════════════
@app.route("/")
def index():
    return send_from_directory("static", "index.html")


@app.route("/api/stations")
def get_stations():
    code   = request.args.get("code", "").strip()
    lat    = float(request.args.get("y", 0))   # WGS84 위도
    lng    = float(request.args.get("x", 0))   # WGS84 경도
    radius = int(request.args.get("radius", 3000))
    prodcd = request.args.get("prodcd", "B027")

    if not code:
        return jsonify({"error": "오피넷 API 키가 없습니다."}), 400

    try:
        kx, ky = wgs84_to_katec(lat, lng)
        print(f"[변환] WGS84({lat:.5f},{lng:.5f}) → KATEC({kx:.1f},{ky:.1f})")

        if radius <= OPINET_MAX:
            raw = opinet_fetch(code, kx, ky, radius, prodcd)
        else:
            seen = {}
            step = 4000
            for dx, dy in [(0,0),(step,0),(-step,0),(0,step),(0,-step),
                           (step,step),(-step,step),(step,-step),(-step,-step)]:
                for s in opinet_fetch(code, kx+dx, ky+dy, OPINET_MAX, prodcd):
                    seen[s["UNI_ID"]] = s
            raw = list(seen.values())

        stations = []
        for s in raw:
            s_lat, s_lng = katec_to_wgs84(s["GIS_X_COOR"], s["GIS_Y_COOR"])
            stations.append({
                "id":       s["UNI_ID"],
                "name":     s["OS_NM"],
                "brand":    s.get("POLL_DIV_CD", ""),
                "price":    int(s["PRICE"]),
                "distance": float(s["DISTANCE"]),  # 오피넷 직선거리(m)
                "lat":      s_lat,
                "lng":      s_lng,
            })

        print(f"[결과] {len(stations)}개 주유소")
        return jsonify({"stations": stations})

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"error": str(e)}), 502


if __name__ == "__main__":
    try:
        local_ip = socket.gethostbyname(socket.gethostname())
    except Exception:
        local_ip = "알 수 없음"

    print("\n" + "=" * 50)
    print("  ⛽ 진짜 최저가 주유소 앱 시작")
    print("=" * 50)
    print(f"  PC 브라우저:   http://localhost:5000")
    print(f"  모바일(같은 WiFi): http://{local_ip}:5000")
    print("=" * 50 + "\n")

    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)
