#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_github_traffic.py
---------------------------------
拉取 GitHub 仓库的 Traffic 数据（views / clones / unique visitors），
追加写入历史 CSV（用于突破 GitHub API 只保留 14 天数据的限制），
并生成月度趋势图 PNG。

设计为由 GitHub Actions 定期调用（建议每 10~13 天一次，
在 14 天数据过期前完成累积，避免漏掉数据点）。

用法（本地测试）:
    export GH_TOKEN="你的GitHub Token"
    python fetch_github_traffic.py --owner JackNg88 --repo jwtools \
        --history-csv traffic_history.csv --outdir ./traffic_report
"""

import argparse
import os
import sys
from datetime import datetime, timezone

import pandas as pd
import requests
import matplotlib.pyplot as plt

API_BASE = "https://api.github.com"


def gh_get(session: requests.Session, url: str) -> dict:
    resp = session.get(url, timeout=30)
    if resp.status_code == 403:
        sys.exit(f"错误: 403 Forbidden — Token 权限不足或访问频率受限。URL: {url}\n{resp.text}")
    resp.raise_for_status()
    return resp.json()


def fetch_traffic(owner: str, repo: str, token: str) -> pd.DataFrame:
    """拉取最近14天的 views 和 clones（按天粒度），返回合并后的 DataFrame。"""
    session = requests.Session()
    session.headers.update({
        "Authorization": f"Bearer {token}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    })

    views_data = gh_get(session, f"{API_BASE}/repos/{owner}/{repo}/traffic/views")
    clones_data = gh_get(session, f"{API_BASE}/repos/{owner}/{repo}/traffic/clones")

    views_df = pd.DataFrame(views_data.get("views", []))
    clones_df = pd.DataFrame(clones_data.get("clones", []))

    if views_df.empty and clones_df.empty:
        print("警告: API 未返回任何数据。", file=sys.stderr)
        return pd.DataFrame(columns=["date", "views", "views_unique", "clones", "clones_unique"])

    if not views_df.empty:
        views_df = views_df.rename(columns={"count": "views", "uniques": "views_unique"})
        views_df["timestamp"] = pd.to_datetime(views_df["timestamp"]).dt.date
        views_df = views_df.rename(columns={"timestamp": "date"})[["date", "views", "views_unique"]]

    if not clones_df.empty:
        clones_df = clones_df.rename(columns={"count": "clones", "uniques": "clones_unique"})
        clones_df["timestamp"] = pd.to_datetime(clones_df["timestamp"]).dt.date
        clones_df = clones_df.rename(columns={"timestamp": "date"})[["date", "clones", "clones_unique"]]

    if views_df.empty:
        merged = clones_df
        merged["views"] = 0
        merged["views_unique"] = 0
    elif clones_df.empty:
        merged = views_df
        merged["clones"] = 0
        merged["clones_unique"] = 0
    else:
        merged = pd.merge(views_df, clones_df, on="date", how="outer").fillna(0)

    for col in ["views", "views_unique", "clones", "clones_unique"]:
        merged[col] = merged[col].astype(int)

    merged["date"] = pd.to_datetime(merged["date"])
    return merged.sort_values("date").reset_index(drop=True)


def update_history(new_df: pd.DataFrame, history_csv: str) -> pd.DataFrame:
    """把新抓取的数据合并进历史 CSV，按 date 去重（新数据覆盖旧数据）。"""
    if os.path.exists(history_csv):
        old_df = pd.read_csv(history_csv, parse_dates=["date"])
        combined = pd.concat([old_df, new_df], ignore_index=True)
    else:
        combined = new_df.copy()

    combined = (
        combined.sort_values("date")
        .drop_duplicates(subset="date", keep="last")
        .reset_index(drop=True)
    )
    combined.to_csv(history_csv, index=False)
    return combined


def aggregate_monthly(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["month_label", "views", "views_unique", "clones", "clones_unique"])
    monthly = (
        df.set_index("date")
        .resample("MS")
        .sum(numeric_only=True)
        .reset_index()
    )
    monthly["month_label"] = monthly["date"].dt.strftime("%Y-%m")
    return monthly


def plot_trend(monthly: pd.DataFrame, outpath: str, repo_label: str):
    fig, ax = plt.subplots(figsize=(10, 5.5))

    ax.plot(monthly["month_label"], monthly["views"], marker="o", linewidth=2,
            color="#2b6cb0", label="Views (total)")
    ax.plot(monthly["month_label"], monthly["views_unique"], marker="o", linewidth=2,
            linestyle="--", color="#63b3ed", label="Views (unique visitors)")
    ax.plot(monthly["month_label"], monthly["clones"], marker="s", linewidth=2,
            color="#dd6b20", label="Clones (total)")
    ax.plot(monthly["month_label"], monthly["clones_unique"], marker="s", linewidth=2,
            linestyle="--", color="#f6ad55", label="Clones (unique)")

    ax.set_title(f"{repo_label} — Monthly GitHub Traffic Trend", fontsize=14, fontweight="bold")
    ax.set_xlabel("Month")
    ax.set_ylabel("Count")
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", fontsize=9)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def generate_badge_json(combined: pd.DataFrame, outdir: str, days: int = 14):
    """生成 shields.io endpoint 格式的 badge.json，展示最近N天的独立访客数。

    README 中引用方式：
    ![Views](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/<owner>/<repo>/main/traffic_report/badge_views.json)
    ![Clones](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/<owner>/<repo>/main/traffic_report/badge_clones.json)
    """
    if combined.empty:
        return

    recent = combined[combined["date"] >= combined["date"].max() - pd.Timedelta(days=days - 1)]
    views_unique = int(recent["views_unique"].sum())
    clones_unique = int(recent["clones_unique"].sum())

    views_badge = {
        "schemaVersion": 1,
        "label": f"views (last {days}d)",
        "message": str(views_unique),
        "color": "blue",
    }
    clones_badge = {
        "schemaVersion": 1,
        "label": f"clones (last {days}d)",
        "message": str(clones_unique),
        "color": "orange",
    }

    import json
    with open(os.path.join(outdir, "badge_views.json"), "w") as f:
        json.dump(views_badge, f)
    with open(os.path.join(outdir, "badge_clones.json"), "w") as f:
        json.dump(clones_badge, f)

    print(f"徽章 JSON 已生成: badge_views.json (views={views_unique}), badge_clones.json (clones={clones_unique})")


def main():
    parser = argparse.ArgumentParser(description="拉取并累积 GitHub 仓库 Traffic 数据，生成月度趋势图")
    parser.add_argument("--owner", required=True, help="仓库所有者，如 JackNg88")
    parser.add_argument("--repo", required=True, help="仓库名，如 jwtools")
    parser.add_argument("--history-csv", default="traffic_history.csv", help="历史数据 CSV 路径（会持续累积，需提交进仓库）")
    parser.add_argument("--outdir", default="./traffic_report", help="输出目录（月度汇总 + 趋势图）")
    args = parser.parse_args()

    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        sys.exit("错误: 请设置环境变量 GH_TOKEN 或 GITHUB_TOKEN")

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[{datetime.now(timezone.utc).isoformat()}] 拉取 {args.owner}/{args.repo} 最近14天 traffic 数据 ...")
    new_df = fetch_traffic(args.owner, args.repo, token)
    print(f"本次抓取到 {len(new_df)} 天的数据。")

    combined = update_history(new_df, args.history_csv)
    print(f"历史数据已更新: {args.history_csv}（累计 {len(combined)} 天）")

    generate_badge_json(combined, args.outdir)

    monthly = aggregate_monthly(combined)
    monthly_csv = os.path.join(args.outdir, "monthly_summary.csv")
    monthly.to_csv(monthly_csv, index=False)
    print(f"月度汇总已保存: {monthly_csv}")
    if not monthly.empty:
        print(monthly[["month_label", "views", "views_unique", "clones", "clones_unique"]].to_string(index=False))

        plot_path = os.path.join(args.outdir, "monthly_trend.png")
        plot_trend(monthly, plot_path, f"{args.owner}/{args.repo}")
        print(f"趋势图已保存: {plot_path}")
    else:
        print("暂无数据可绘图（历史数据为空）。")


if __name__ == "__main__":
    main()
