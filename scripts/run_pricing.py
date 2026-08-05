#!/usr/bin/env python3
"""
run_pricing.py — Derivatives Pricing & Stochastic Calculus 的可执行骨架（优化版 v2）。

纯数学、无需联网、无需行情数据：给定一份期权合约与市场参数，输出
  - Black-Scholes-Merton 解析价 + 全套希腊字母
  - CRR 二叉树价（支持美式）
  - 蒙特卡洛价（GBM + 对偶变量）含标准误与 95% 置信区间
  - 隐含波动率（从市场价反解，含收敛/边界状态）
  - 看跌看涨平价校验（用【独立】的欧式二叉树 call/put 价，而非同一闭式解自证）
    + 无套利边界校验 + 跨方法一致性校验 + 二叉树收敛校验
  - 可选：多行权价的波动率微笑/偏斜分析
并按报告蓝图写出 Markdown。

相对 v1 的优化（详见 CHANGELOG_优化说明.md）：
  1. 摘要中的“一致/通过/收敛”改为依据校验结果【动态生成】，不再硬编码——
     避免出现“风险=🔴高 但摘要写‘三法一致、已通过’”这种自相矛盾。
  2. 看跌看涨平价改用【独立】的欧式二叉树 call/put 价来校验。v1 用同一 BSM
     闭式解的 call 减 put，因 N(x)+N(-x)≡1 残差恒等于 0，无法发现任何实现错误。
  3. 美式期权 headline 改用美式二叉树价（v1 误用欧式 BSM 价，会按提前行权溢价
     低估美式期权价值，且把欧式 MC 与美式树价混在一张表里）。
  4. 跨方法校验引入绝对容差 abs_tol，避免深度价外/近到期等近零价格的假阳性高风险旗。
  5. 二叉树收敛、隐含波动率边界状态如实判断并报告，不再写死“已收敛”。
  6. N 步数、MC 路径数等不再写死在报告文字里，随实际参数变化（--steps/--paths）。
  7. 新增 --strikes/--market-prices 的波动率微笑/偏斜分析（第 7 章）。

用法示例：
  python run_pricing.py --type call --S 100 --K 100 --T 1 --r 0.03 --q 0 --sigma 0.2
  python run_pricing.py --type put  --S 50  --K 55  --T 0.5 --sigma 0.3 --american 1
  python run_pricing.py --type call --S 100 --K 100 --T 1 --sigma 0.2 --market-price 9.5
  python run_pricing.py --type call --S 100 --T 1 --sigma 0.2 \
         --strikes 90,95,100,105,110 --market-prices 14.8,11.3,9.5,6.6,4.4
"""
import argparse, math
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

SQRT = math.sqrt
EXP = math.exp
LOG = math.log

IV_LO, IV_HI = 1e-4, 5.0          # 隐含波动率求解区间


# ============ 1. 解析定价（Black-Scholes-Merton） ============
def _d1_d2(S, K, T, r, q, sigma):
    d1 = (LOG(S / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * SQRT(T))
    return d1, d1 - sigma * SQRT(T)


def bsm_price(opt, S, K, T, r, q, sigma):
    if T <= 0 or sigma <= 0:                       # 退化为内在价值
        return max(S - K, 0.0) if opt == "call" else max(K - S, 0.0)
    d1, d2 = _d1_d2(S, K, T, r, q, sigma)
    if opt == "call":
        return S * EXP(-q * T) * norm.cdf(d1) - K * EXP(-r * T) * norm.cdf(d2)
    return K * EXP(-r * T) * norm.cdf(-d2) - S * EXP(-q * T) * norm.cdf(-d1)


def bsm_greeks(opt, S, K, T, r, q, sigma):
    """返回市场惯例口径：vega/1 vol点、theta/日、rho/1%利率。"""
    d1, d2 = _d1_d2(S, K, T, r, q, sigma)
    pdf, dq, dr = norm.pdf(d1), EXP(-q * T), EXP(-r * T)
    gamma = dq * pdf / (S * sigma * SQRT(T))
    vega = S * dq * pdf * SQRT(T)
    if opt == "call":
        delta = dq * norm.cdf(d1)
        theta = (-S * dq * pdf * sigma / (2 * SQRT(T))
                 - r * K * dr * norm.cdf(d2) + q * S * dq * norm.cdf(d1))
        rho = K * T * dr * norm.cdf(d2)
    else:
        delta = -dq * norm.cdf(-d1)
        theta = (-S * dq * pdf * sigma / (2 * SQRT(T))
                 + r * K * dr * norm.cdf(-d2) - q * S * dq * norm.cdf(-d1))
        rho = -K * T * dr * norm.cdf(-d2)
    return dict(delta=delta, gamma=gamma, vega=vega / 100,
                theta=theta / 365, rho=rho / 100)


# ============ 2. 数值定价 ============
def binomial_price(opt, S, K, T, r, q, sigma, N=500, american=False):
    """Cox-Ross-Rubinstein 二叉树，支持美式提前行权。"""
    dt = T / N
    u = EXP(sigma * SQRT(dt)); d = 1 / u
    p = (EXP((r - q) * dt) - d) / (u - d)
    disc = EXP(-r * dt)
    j = np.arange(N + 1)
    ST = S * u ** j * d ** (N - j)                 # 末端标的价
    val = np.maximum(ST - K, 0.0) if opt == "call" else np.maximum(K - ST, 0.0)
    for n in range(N - 1, -1, -1):                 # 倒推到时间层 n（n+1 个节点）
        val = disc * (p * val[1:] + (1 - p) * val[:-1])
        if american:
            j = np.arange(n + 1)
            Sn = S * u ** j * d ** (n - j)
            exer = np.maximum(Sn - K, 0.0) if opt == "call" else np.maximum(K - Sn, 0.0)
            val = np.maximum(val, exer)
    return float(val[0])


def mc_price(opt, S, K, T, r, q, sigma, n_paths=200_000, seed=1):
    """欧式：GBM 终值 + 对偶变量，返回 (价, 标准误)。"""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n_paths // 2)
    z = np.concatenate([z, -z])                    # 对偶变量降方差
    ST = S * np.exp((r - q - 0.5 * sigma ** 2) * T + sigma * SQRT(T) * z)
    payoff = np.maximum(ST - K, 0.0) if opt == "call" else np.maximum(K - ST, 0.0)
    disc = EXP(-r * T)
    price = disc * payoff.mean()
    se = disc * payoff.std(ddof=1) / SQRT(len(payoff))
    return float(price), float(se)


# ============ 3. 隐含波动率 ============
def implied_vol(opt, market_price, S, K, T, r, q):
    """返回 (iv, status)。status ∈ {'ok','below_intrinsic','above_bound','no_converge','boundary'}。"""
    intrinsic = max(S - K, 0.0) if opt == "call" else max(K - S, 0.0)
    upper = S * EXP(-q * T) if opt == "call" else K * EXP(-r * T)
    if market_price < intrinsic - 1e-8:
        return None, "below_intrinsic"             # 低于内在价值 -> 套利/无解
    if market_price > upper + 1e-8:
        return None, "above_bound"                 # 越过无套利上界
    try:
        iv = float(brentq(lambda s: bsm_price(opt, S, K, T, r, q, s) - market_price,
                          IV_LO, IV_HI, maxiter=200))
    except Exception:
        return None, "no_converge"
    if iv <= IV_LO * 1.001 or iv >= IV_HI * 0.999:
        return iv, "boundary"                      # 触边界，疑似越界报价
    return iv, "ok"


# ============ 4. 校验与模型风险规则 ============
def validate(opt, S, K, T, r, q, sigma, ana_eu, tree_eu, mc, mc_se,
             tree_half, tree_call_eu, tree_put_eu, market_price, iv, iv_status,
             tol=0.005, abs_tol=1e-3):
    """返回 (flags, parity_res, diag)。flags 按 高→中 顺序；diag 供摘要生成动态文案。"""
    flags = []           # 高、中风险
    diag = {}

    # --- 跨方法一致性（欧式锚点：BSM vs 欧式树 vs 欧式 MC）---
    tree_dev = abs(tree_eu - ana_eu)
    tree_ok = tree_dev <= max(tol * abs(ana_eu), abs_tol)
    if not tree_ok:
        flags.append(("🔴 高", "二叉树与解析价不一致",
                      f"|tree−ana|>max({tol:.1%}·ana, {abs_tol})",
                      f"ana={ana_eu:.4f}, tree={tree_eu:.4f}, 偏差={tree_dev:.4f}"))
    ci_lo, ci_hi = mc - 1.96 * mc_se - abs_tol, mc + 1.96 * mc_se + abs_tol
    mc_ok = ci_lo <= ana_eu <= ci_hi
    if not mc_ok:
        flags.append(("🔴 高", "蒙特卡洛与解析价不一致", "解析价落在 MC 95% 置信区间外",
                      f"ana={ana_eu:.4f}, MC95%CI=[{mc-1.96*mc_se:.4f}, {mc+1.96*mc_se:.4f}]"))
    diag["cross_method_ok"] = tree_ok and mc_ok

    # --- 看跌看涨平价（独立：欧式二叉树 call 与 put，可发现实现错误）---
    parity_res = tree_call_eu - tree_put_eu - (S * EXP(-q * T) - K * EXP(-r * T))
    parity_ok = abs(parity_res) <= max(tol * S, abs_tol)
    diag["parity_res"] = parity_res
    diag["parity_ok"] = parity_ok
    if not parity_ok:
        flags.append(("🔴 高", "看跌看涨平价被破坏（独立树价）",
                      f"|残差|>{tol:.1%}·S",
                      f"残差={parity_res:.4e}（C_tree={tree_call_eu:.4f}, P_tree={tree_put_eu:.4f}）"))

    # --- 无套利边界 ---
    intrinsic = max(S - K, 0.0) if opt == "call" else max(K - S, 0.0)
    upper = S * EXP(-q * T) if opt == "call" else K * EXP(-r * T)
    diag["bounds_ok"] = (ana_eu >= intrinsic - abs_tol) and (ana_eu <= upper + abs_tol)
    if ana_eu < intrinsic - abs_tol:
        flags.append(("🔴 高", "价格低于内在价值", "price<intrinsic",
                      f"price={ana_eu:.4f}, intrinsic={intrinsic:.4f}"))
    if ana_eu > upper + abs_tol:
        flags.append(("🔴 高", "价格越过无套利上界", "price>upper",
                      f"price={ana_eu:.4f}, upper={upper:.4f}"))

    # --- 隐含波动率状态 ---
    if market_price is not None:
        if iv_status in ("below_intrinsic", "above_bound"):
            flags.append(("🔴 高", "市场价越界", "市场价<内在价值 或 >无套利上界",
                          f"market={market_price}, status={iv_status}"))
        elif iv is None or iv_status == "no_converge":
            flags.append(("🔴 高", "隐含波动率不收敛", "IV 求解失败", f"market={market_price}"))
        elif iv_status == "boundary":
            flags.append(("🔴 高", "隐含波动率触边界", "IV≈解空间边界，疑似越界报价",
                          f"IV={iv:.4f}"))

    # --- 二叉树收敛（欧式锚点 N/2 → N）---
    conv_dev = abs(tree_eu - tree_half)
    converged = conv_dev <= max(tol * abs(tree_eu), abs_tol)
    diag["converged"] = converged
    if not converged:
        flags.append(("🟡 中", "二叉树未完全收敛", "N/2→N 变化超容差",
                      f"半步价={tree_half:.4f}, 全步价={tree_eu:.4f}, 变化={conv_dev:.4f}"))

    # --- 蒙特卡洛精度 ---
    if abs(mc) > abs_tol and mc_se / abs(mc) > 0.01:
        flags.append(("🟡 中", "蒙特卡洛标准误偏大", "SE/price>1%",
                      f"SE/price={mc_se/abs(mc):.2%}，建议增加路径"))

    # --- 数值不稳定区 ---
    if T < 0.02:
        flags.append(("🟡 中", "到期极短", "T<0.02 年", f"T={T}"))
    if abs(LOG(S / K)) > 0.7:
        flags.append(("🟡 中", "深度价内/价外", "|ln(S/K)|>0.7", f"S/K={S/K:.3f}"))

    diag["all_core_ok"] = (diag["cross_method_ok"] and diag["parity_ok"]
                           and diag["bounds_ok"] and diag["converged"])
    return flags, parity_res, diag


# ============ 5. 波动率微笑/偏斜 ============
def smile_analysis(opt, S, T, r, q, strikes, market_prices, skew_tol=0.05):
    """对一组 (strike, market_price) 反解 IV，返回 (rows, skew, flag_or_None)。"""
    rows = []
    for K, mp in zip(strikes, market_prices):
        iv, status = implied_vol(opt, mp, S, K, T, r, q)
        rows.append(dict(K=K, market=mp, iv=iv, status=status, m=LOG(S / K)))
    ivs = [row["iv"] for row in rows if row["iv"] is not None]
    skew = (max(ivs) - min(ivs)) if len(ivs) >= 2 else None
    flag = None
    if skew is not None and skew > skew_tol:
        flag = ("🟡 中", "明显波动率偏斜",
                f"各行权 IV 极差>{skew_tol:.0%}（{skew_tol*100:.0f}vol点）而单一 σ 不足以刻画",
                f"IV 极差={skew*100:.2f} vol点")
    return rows, skew, flag


# ============ 6. 报告输出（9 章蓝图） ============
def write_report(args, ana, greeks, tree, tree_eu, tree_half, mc, mc_se, iv, iv_status,
                 flags, parity_res, diag, smile_rows, smile_skew, path):
    opt, S, K, T, r, q, sigma = (args.type, args.S, args.K, args.T,
                                 args.r, args.q, args.sigma)
    N, paths = args.steps, args.paths
    d1, d2 = _d1_d2(S, K, T, r, q, sigma)
    style = "美式" if args.american else "欧式"
    side = "看涨" if opt == "call" else "看跌"
    top = flags[0][0] if flags else "🟢 低"

    headline = tree if args.american else ana      # 美式以二叉树价为准
    pct = headline / S * 100 if S else 0.0

    cross_phrase = ("三法在容差内一致" if diag["cross_method_ok"]
                    else "三法存在偏差，疑似实现/输入问题（见第 8 章）")
    parity_phrase = "应≈0，已通过" if diag["parity_ok"] else "超出容差，见第 8 章"
    conv_phrase = "已收敛" if diag["converged"] else "未完全收敛，建议增大 N"

    L = [
        f"# 衍生品定价报告：{style}{side}期权",
        f"\n> 标的 S={S}｜行权 K={K}｜到期 T={T}y｜利率 r={r}｜红利 q={q}｜波动率 σ={sigma}"
        f"｜口径：vega/1vol点, theta/日, rho/1%\n",
        "## 1. 摘要与结论",
    ]
    if args.american:
        L.append(f"- **这份美式期权的理论价格约为 {tree:.2f} 元**（二叉树法），"
                 f"约相当于标的现价的 {pct:.1f}%。")
        L.append(f"- 欧式锚点三法互证：BSM={ana:.4f}｜欧式二叉树={tree_eu:.4f}｜"
                 f"蒙特卡洛(欧式)={mc:.4f}±{1.96*mc_se:.4f}（95%CI）——{cross_phrase}。")
        L.append(f"- 美式提前行权溢价 = 美式二叉树{tree:.4f} − 欧式二叉树{tree_eu:.4f} "
                 f"= **{tree-tree_eu:.4f}**。")
    else:
        L.append(f"- **这份期权的理论价格约为 {ana:.2f} 元**，"
                 f"约相当于标的现价的 {pct:.1f}%（{cross_phrase}，见第 5 章）。")
        L.append(f"- 解析价（BSM）= **{ana:.4f}**；二叉树={tree:.4f}；"
                 f"蒙特卡洛={mc:.4f}±{1.96*mc_se:.4f}（95%CI）——{cross_phrase}。")
    L += [
        f"- **Delta={greeks['delta']:.4f}**：标的每涨 1 元，期权价约变动 {greeks['delta']:+.2f} 元；"
        f"对冲 1 份期权约需反向持有 {abs(greeks['delta']):.2f} 份标的。",
        f"- 其余希腊字母：Gamma={greeks['gamma']:.4f}｜Vega={greeks['vega']:.4f}｜"
        f"Theta={greeks['theta']:.4f}｜Rho={greeks['rho']:.4f}（详见第 4 章）。",
        f"- 看跌看涨平价残差（独立欧式树价）= {parity_res:.2e}（{parity_phrase}）；"
        f"最高风险等级 **{top}**（详见第 8 章）。",
    ]
    if iv is not None:
        L.append(f"- 由市场价 {args.market_price} 反解隐含波动率 = **{iv:.4f}**"
                 f"（状态：{iv_status}；输入 σ={sigma}，详见第 6 章）。")
    L.append(
        "- ⚠️ **使用须知**：本结果基于简化的 Black-Scholes 模型（假设波动率恒定、可连续对冲、"
        "无交易成本，未含真实波动率微笑），适用于学习与研究，**不可直接用于真实下单**；"
        f"其中波动率 σ={sigma} 是关键且最敏感的输入，请确认它来自合理估计。")

    # 第 2 章
    L += [
        "\n## 2. 合约与市场输入",
        "| 项目 | 取值 |", "|---|---|",
        f"| 期权类型 | {style}{side} |",
        f"| 标的现价 S | {S} |",
        f"| 行权价 K | {K} |",
        f"| 到期时间 T | {T} 年 |",
        f"| 无风险利率 r | {r}（连续复利） |",
        f"| 红利/持有收益 q | {q}（连续复利） |",
        f"| 波动率 σ | {sigma}（年化） |",
        "\n希腊字母口径：vega 按每 1 个波动率百分点、theta 按每个日历日、rho 按每 1% 利率。",
    ]

    # 第 3 章
    L += ["\n## 3. 解析定价",
          f"d1={d1:.4f}, d2={d2:.4f}；BSM 欧式价 = {ana:.4f}。"
          + ("（注：BSM 闭式解仅适用于欧式；美式价以第 5 章二叉树为准。）"
             if args.american else "")]

    # 第 4 章
    L += [
        "\n## 4. 希腊字母",
        "| 希腊字母 | 数值 | 含义 |", "|---|---|---|",
        f"| Delta | {greeks['delta']:.4f} | 标的每涨 1 元，期权价变动 |",
        f"| Gamma | {greeks['gamma']:.4f} | Delta 对标的的敏感度 |",
        f"| Vega（/1vol点） | {greeks['vega']:.4f} | 波动率每升 1 个百分点 |",
        f"| Theta（/日） | {greeks['theta']:.4f} | 每过一日的时间损耗 |",
        f"| Rho（/1%） | {greeks['rho']:.4f} | 利率每升 1 个百分点 |",
    ]

    # 第 5 章
    mc_label = "蒙特卡洛（欧式锚点）" if args.american else "蒙特卡洛"
    tree_label = f"二叉树(N={N}{'，美式' if args.american else ''})"
    L += [
        "\n## 5. 数值定价与收敛",
        "| 方法 | 价格 | 诊断 |", "|---|---|---|",
        f"| 解析 BSM（欧式） | {ana:.4f} | 闭式解 |",
        f"| {tree_label} | {tree:.4f} | 欧式锚点 N={N//2}→{N}：{tree_half:.4f}→{tree_eu:.4f}，{conv_phrase} |",
        f"| {mc_label} | {mc:.4f} | SE={mc_se:.4f}, 95%CI=[{mc-1.96*mc_se:.4f}, {mc+1.96*mc_se:.4f}] |",
    ]
    if args.american:
        L += [
            "\n> 注：蒙特卡洛为欧式终值法、不含提前行权，故仅作欧式锚点与 BSM/欧式树互证；"
            "美式价以二叉树为准。",
            f"\n美式提前行权溢价 = {tree - tree_eu:.4f}（美式树价 − 欧式树价）。",
        ]

    # 第 6 章
    L.append("\n## 6. 隐含波动率")
    if iv is not None:
        L.append(f"由市场价 {args.market_price} 反解得 IV = {iv:.4f}（求解状态：{iv_status}），"
                 f"对照输入 σ={sigma}。")
    elif args.market_price is not None:
        L.append(f"市场价 {args.market_price} 反解失败（状态：{iv_status}），"
                 f"已在第 8 章标记为风险，不强行给出数字。")
    else:
        L.append("本次未提供市场价，无此项。如需反解隐含波动率，"
                 "请在命令中加上 `--market-price 期权市场价`。")

    # 第 7 章
    L.append("\n## 7. 波动率结构")
    if smile_rows:
        L += ["| 行权 K | 货币度 ln(S/K) | 市场价 | 隐含波动率 IV | 状态 |", "|---|---|---|---|---|"]
        for row in smile_rows:
            ivs = f"{row['iv']:.4f}" if row['iv'] is not None else "—"
            L.append(f"| {row['K']:.2f} | {row['m']:+.4f} | {row['market']:.4f} | {ivs} | {row['status']} |")
        if smile_skew is not None:
            L.append(f"\n各行权隐含波动率极差 = {smile_skew*100:.2f} vol点。"
                     + ("偏斜明显，单一平坦 σ 不足以同时刻画所有行权，详见第 8 章。"
                        if smile_skew > 0.05 else "偏斜温和。"))
    else:
        L.append("本次为单一合约，无波动率微笑/偏斜分析。"
                 "如需该项，请用 `--strikes` 与 `--market-prices` 提供同一标的、多个行权价的报价。")

    # 第 8 章
    L += ["\n## 8. 模型风险与校验清单",
          "| 风险等级 | 信号 | 触发规则 | 证据 |", "|---|---|---|---|"]
    show = flags if flags else [("🟢 低", "未触发高/中规则", "—", "通过跨方法/平价/边界/收敛校验")]
    for lv, sig, rule, ev in show:
        L.append(f"| {lv} | {sig} | {rule} | {ev} |")

    # 第 9 章
    L += [
        "\n## 9. 方法附录",
        "| 计算阶段 | 方法 | 关键输入 | 输出 | 容差/参数 |",
        "|---|---|---|---|---|",
        f"| 解析定价 | Black-Scholes-Merton | S,K,T,r,q,σ | price={ana:.4f} | 闭式解 |",
        f"| 希腊字母 | BSM 解析导数 | 同上 | Δ={greeks['delta']:.4f} 等 | vega/1vol,theta/日,rho/1% |",
        f"| 二叉树 | CRR{'（美式）' if args.american else ''} | 同上 | {tree:.4f} | N={N}, 收敛容差0.5% |",
        f"| 蒙特卡洛 | GBM+对偶 | 同上 | {mc:.4f}±{mc_se:.4f} | {paths:,} 路径, seed=1 |",
        f"| 平价校验 | 独立欧式树 C−P | 同上 | 残差={parity_res:.2e} | 容差0.5%·S |",
        "\n---",
        "本报告基于公开数据与规则化分析生成，仅供研究参考，不构成任何投资建议。",
    ]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))


def _parse_floats(s):
    return [float(x) for x in s.split(",")] if s else None


# ============ main ============
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--type", choices=["call", "put"], default="call")
    ap.add_argument("--S", type=float, default=100.0)
    ap.add_argument("--K", type=float, default=100.0)
    ap.add_argument("--T", type=float, default=1.0)
    ap.add_argument("--r", type=float, default=0.03)
    ap.add_argument("--q", type=float, default=0.0)
    ap.add_argument("--sigma", type=float, default=None)
    ap.add_argument("--market-price", type=float, default=None, dest="market_price")
    ap.add_argument("--american", type=int, default=0)
    ap.add_argument("--steps", type=int, default=500, help="二叉树步数 N")
    ap.add_argument("--paths", type=int, default=200_000, help="蒙特卡洛路径数")
    ap.add_argument("--strikes", default=None, help="波动率微笑：逗号分隔的多个行权价")
    ap.add_argument("--market-prices", default=None, dest="market_prices",
                    help="波动率微笑：与 --strikes 对应的逗号分隔市场价")
    ap.add_argument("--out", default="pricing_report.md")
    args = ap.parse_args()

    # 若只给市场价未给 σ，则先反解隐含波动率作为定价 σ
    iv, iv_status = None, None
    if args.market_price is not None:
        iv, iv_status = implied_vol(args.type, args.market_price, args.S, args.K,
                                    args.T, args.r, args.q)
        if args.sigma is None:
            if iv is None:
                raise ValueError(f"市场价反解失败（{iv_status}），无法定价。")
            args.sigma = iv
    if args.sigma is None:
        args.sigma = 0.2                            # 默认

    opt, S, K, T, r, q, sig = (args.type, args.S, args.K, args.T,
                               args.r, args.q, args.sigma)
    N = args.steps

    ana = bsm_price(opt, S, K, T, r, q, sig)                       # 欧式 BSM（请求类型）
    greeks = bsm_greeks(opt, S, K, T, r, q, sig)
    # 二叉树：独立的欧式 call/put（用于平价），半步与全步（用于收敛），美式（按需）
    tree_call_eu = binomial_price("call", S, K, T, r, q, sig, N=N, american=False)
    tree_put_eu = binomial_price("put", S, K, T, r, q, sig, N=N, american=False)
    tree_eu = tree_call_eu if opt == "call" else tree_put_eu      # 欧式锚点价
    tree_half = binomial_price(opt, S, K, T, r, q, sig, N=N // 2, american=False)
    tree = (binomial_price(opt, S, K, T, r, q, sig, N=N, american=True)
            if args.american else tree_eu)                        # headline 用价
    mc, mc_se = mc_price(opt, S, K, T, r, q, sig, n_paths=args.paths)

    flags, parity_res, diag = validate(opt, S, K, T, r, q, sig, ana, tree_eu,
                                       mc, mc_se, tree_half, tree_call_eu, tree_put_eu,
                                       args.market_price, iv, iv_status)

    # 波动率微笑
    smile_rows = smile_skew = None
    strikes = _parse_floats(args.strikes)
    mprices = _parse_floats(args.market_prices)
    if strikes and mprices:
        if len(strikes) != len(mprices):
            raise ValueError("--strikes 与 --market-prices 长度不一致。")
        smile_rows, smile_skew, smile_flag = smile_analysis(opt, S, T, r, q, strikes, mprices)
        if smile_flag:
            flags.append(smile_flag)

    write_report(args, ana, greeks, tree, tree_eu, tree_half, mc, mc_se, iv, iv_status,
                 flags, parity_res, diag, smile_rows, smile_skew, args.out)

    top = flags[0][0] if flags else "🟢 低"
    print(f"[done] {('美式' if args.american else '欧式')}{opt} S={S} K={K} T={T} σ={sig:.4f}")
    print(f"  BSM={ana:.4f} | tree={tree:.4f} | MC={mc:.4f}±{1.96*mc_se:.4f} "
          f"| Δ={greeks['delta']:.4f} | parity={parity_res:.2e} | top flag={top}")
    if args.american:
        print(f"  美式提前行权溢价 = {tree - tree_eu:.4f}")
    if iv is not None:
        print(f"  隐含波动率 IV = {iv:.4f}（{iv_status}）")
    if smile_skew is not None:
        print(f"  微笑 IV 极差 = {smile_skew*100:.2f} vol点")
    print(f"  report -> {args.out}")


if __name__ == "__main__":
    main()
