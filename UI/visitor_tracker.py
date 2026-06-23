"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Visitor Tracker - World Access Statistics for NonBDNAFinder                  │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
Tracks which countries are accessing the tool and renders a world choropleth map.
Visitor counts are persisted in a JSON file next to this module.
"""
import os
import json
import logging
import hashlib
import time
from typing import Dict, Optional

import streamlit as st

logger = logging.getLogger(__name__)

_STATS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "visitor_stats.json")
_GEO_API_URL = "https://ip-api.com/json/{ip}?fields=status,country,countryCode"
_SESSION_KEY = "_visitor_tracked"
_LOCK_TIMEOUT = 10  # seconds


# ──────────────────────────────────────────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def _load_stats() -> Dict[str, int]:
    """Load visitor statistics from the JSON file."""
    try:
        if os.path.exists(_STATS_FILE):
            with open(_STATS_FILE, "r") as f:
                data = json.load(f)
                return {k: int(v) for k, v in data.items() if isinstance(k, str)}
    except Exception as exc:
        logger.debug("Could not load visitor stats: %s", exc)
    return {}


def _save_stats(stats: Dict[str, int]) -> None:
    """Persist visitor statistics to the JSON file."""
    try:
        with open(_STATS_FILE, "w") as f:
            json.dump(stats, f, indent=2)
    except Exception as exc:
        logger.debug("Could not save visitor stats: %s", exc)


def _get_client_ip() -> Optional[str]:
    """Extract the real client IP from Streamlit request headers."""
    try:
        headers = st.context.headers
        for header in ("x-forwarded-for", "x-real-ip", "cf-connecting-ip"):
            val = headers.get(header) or headers.get(header.title())
            if val:
                return val.split(",")[0].strip()
    except Exception:
        pass
    return None


def _resolve_country(ip: str) -> Optional[str]:
    """Return the ISO-3166-1 alpha-2 country code for the given IP."""
    try:
        import requests
        resp = requests.get(_GEO_API_URL.format(ip=ip), timeout=3)
        data = resp.json()
        if data.get("status") == "success" and data.get("countryCode"):
            return data["countryCode"]
    except Exception as exc:
        logger.debug("GeoIP lookup failed for %s: %s", ip, exc)
    return None


def record_visit() -> None:
    """
    Record one visit from the current session's country.
    Each session is counted only once (guarded by ``_SESSION_KEY``).
    """
    if st.session_state.get(_SESSION_KEY):
        return

    st.session_state[_SESSION_KEY] = True

    ip = _get_client_ip()
    if not ip or ip.startswith(("127.", "::1", "10.", "192.168.", "172.")):
        return  # Skip localhost / private networks

    # Avoid re-resolving within 1 second (unlikely re-entrancy guard)
    country_code = _resolve_country(ip)
    if not country_code:
        return

    stats = _load_stats()
    stats[country_code] = stats.get(country_code, 0) + 1
    _save_stats(stats)


# ──────────────────────────────────────────────────────────────────────────────
# RENDERING
# ──────────────────────────────────────────────────────────────────────────────

# ISO 3166-1 alpha-2 → country name (subset covering most visitors)
_CC_NAME: Dict[str, str] = {
    "US": "United States", "IN": "India", "GB": "United Kingdom",
    "DE": "Germany", "FR": "France", "CA": "Canada", "AU": "Australia",
    "CN": "China", "JP": "Japan", "BR": "Brazil", "KR": "South Korea",
    "IT": "Italy", "ES": "Spain", "NL": "Netherlands", "SE": "Sweden",
    "CH": "Switzerland", "SG": "Singapore", "RU": "Russia", "MX": "Mexico",
    "ZA": "South Africa", "NG": "Nigeria", "EG": "Egypt", "PK": "Pakistan",
    "BD": "Bangladesh", "ID": "Indonesia", "MY": "Malaysia", "TH": "Thailand",
    "PL": "Poland", "TR": "Turkey", "SA": "Saudi Arabia", "AR": "Argentina",
    "CL": "Chile", "CO": "Colombia", "PT": "Portugal", "BE": "Belgium",
    "AT": "Austria", "NO": "Norway", "DK": "Denmark", "FI": "Finland",
    "NZ": "New Zealand", "IL": "Israel", "HK": "Hong Kong", "TW": "Taiwan",
    "UA": "Ukraine", "CZ": "Czech Republic", "RO": "Romania", "HU": "Hungary",
    "GR": "Greece", "IR": "Iran",
}


def render_visitor_map() -> None:
    """Render the world access map and statistics at the bottom of the Home page."""
    try:
        import plotly.graph_objects as go
        import pandas as pd
    except ImportError:
        st.info("Install plotly to view the visitor map.")
        return

    stats = _load_stats()
    total = sum(stats.values())

    st.markdown("---")
    st.markdown(
        "<h2 style='font-size:1.35rem;font-weight:700;margin-bottom:0.3rem;'>"
        "🌍 Global Access Statistics</h2>"
        "<p style='font-size:0.92rem;color:#64748b;margin-top:0;'>"
        "Countries that have accessed NonBDNAFinder</p>",
        unsafe_allow_html=True,
    )

    if not stats:
        st.info(
            "No visitor data recorded yet. Statistics will appear here as "
            "users from around the world access this tool."
        )
        return

    # Build dataframe
    rows = []
    for cc, count in stats.items():
        rows.append({
            "country_code": cc,
            "country": _CC_NAME.get(cc, cc),
            "visits": count,
        })
    df = pd.DataFrame(rows).sort_values("visits", ascending=False)

    # Summary metrics
    col_a, col_b, col_c = st.columns(3)
    col_a.metric("🌐 Total Visits", f"{total:,}")
    col_b.metric("🗺️ Countries Reached", f"{len(stats):,}")
    top_country = df.iloc[0]["country"] if not df.empty else "—"
    col_c.metric("🏆 Top Country", top_country)

    # Choropleth map
    fig = go.Figure(
        go.Choropleth(
            locations=df["country_code"],
            z=df["visits"],
            text=df["country"],
            colorscale="YlOrRd",
            autocolorscale=False,
            reversescale=False,
            marker_line_color="white",
            marker_line_width=0.5,
            colorbar_title="Visits",
            hovertemplate="<b>%{text}</b><br>Visits: %{z:,}<extra></extra>",
        )
    )
    fig.update_layout(
        geo=dict(
            showframe=False,
            showcoastlines=True,
            coastlinecolor="#cbd5e1",
            showland=True,
            landcolor="#f1f5f9",
            showocean=True,
            oceancolor="#e0f2fe",
            projection_type="natural earth",
        ),
        margin=dict(l=0, r=0, t=10, b=0),
        height=380,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    st.plotly_chart(fig, use_container_width=True)

    # Country-wise visitor statistics (all countries)
    country_stats = df.copy()
    country_stats["share (%)"] = (country_stats["visits"] / total * 100).round(1)
    country_stats = country_stats.rename(columns={"country": "Country", "visits": "Visitors"})
    country_stats.insert(0, "Rank", range(1, len(country_stats) + 1))

    st.markdown("#### 📊 Country-wise Visitor Counts")
    st.dataframe(
        country_stats[["Rank", "Country", "Visitors", "share (%)"]],
        use_container_width=True,
        hide_index=True,
    )

    if len(country_stats) > 1:
        st.caption("Top countries by visitors")
        st.bar_chart(
            country_stats.head(10).set_index("Country")["Visitors"],
            use_container_width=True,
        )
