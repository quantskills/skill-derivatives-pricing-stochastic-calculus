---
name: derivatives-pricing-stochastic-calculus
name_zh: 衍生品定价与随机微积分
description: Price options and other derivatives and quantify their risk from a contract
  specification and market inputs, covering Black-Scholes-Merton analytical pricing,
  the full set of Greeks, binomial-tree and Monte-Carlo numerical pricing with convergence
  checks, implied-volatility inversion, volatility smile/skew structure, and model-risk
  validation (put-call parity, no-arbitrage bounds, cross-method consistency). Use when
  the user asks for 期权定价、希腊字母计算、隐含波动率、波动率微笑/曲面、二叉树或蒙特卡洛定价、Black-Scholes、美式期权定价、看跌看涨平价校验,
  or a one-stop option pricing and risk dossier.
description_zh: 输入一份期权合约与市场参数，输出一份可复现、可校验的衍生品定价与风险报告：Black-Scholes-Merton 解析定价、全套希腊字母、二叉树与蒙特卡洛数值定价及收敛检查、隐含波动率反解、波动率微笑/偏斜结构、以及模型风险校验（看跌看涨平价、无套利边界、跨方法一致性）。适用于期权定价、希腊字母计算、隐含波动率、波动率曲面、二叉树/蒙特卡洛定价、美式期权定价、模型校验等场景。
metadata:
  organization: QuantSkills
  organization_url: https://github.com/quantskills
  repository: skill-derivatives-pricing-stochastic-calculus
  repository_url: https://github.com/quantskills/skill-derivatives-pricing-stochastic-calculus
  project_type: skill
  collection: derivatives-pricing-stochastic-calculus
quantSkills:
  project_type: skill
  category: research
  tags:
  - derivatives-pricing
  - options
  - black-scholes
  - greeks
  - implied-volatility
  - monte-carlo
  - binomial-tree
  - stochastic-calculus
  platforms:
  - claude-code
  - codex
  - openclaw
  - cursor
  status: stable
  validation_level: runnable
  maintainer_type: official
  summary_zh: 输入期权合约与市场参数，输出可复现的定价与风险报告：BSM 解析价、全套希腊字母、二叉树与蒙特卡洛数值价、隐含波动率、平价与无套利校验，一次算清。
  summary_en: A derivatives-pricing skill that computes Black-Scholes-Merton prices
    and Greeks, cross-checks them with binomial-tree and Monte-Carlo methods, inverts
    implied volatility, and runs put-call-parity and no-arbitrage validation.
  license: GPL-3.0
  requires: []
---

# Derivatives Pricing & Stochastic Calculus

Use this skill to turn one option contract such as a European call on spot 100, strike 100, one-year maturity into a sourced, reproducible pricing-and-risk dossier covering Black-Scholes-Merton analytical pricing, the full Greeks, binomial-tree and Monte-Carlo numerical pricing with convergence checks, implied-volatility inversion, volatility-structure review, and model-risk validation.

The single most important job of this skill is **getting the number right and proving it**. A pricing model is only trustworthy when independent methods agree, put-call parity holds, and the price respects no-arbitrage bounds. Default to cross-checking, not to a single closed-form number.

## Core Workflow

1. Normalize the contract. Capture the essentials: option type (call/put), exercise style (European/American), spot `S`, strike `K`, time to maturity `T` in years, risk-free rate `r`, dividend/carry yield `q`, and either a volatility `sigma` or an observed market price to invert. Ask only when a required field is missing and cannot be defaulted.
2. Confirm assumptions. Default to continuously-compounded `r` and `q`, volatility quoted annualized, and Greeks reported in market conventions (vega per 1 vol point, theta per calendar day, rho per 1% rate) unless the user requests otherwise. State every assumption.
3. Read `references/pricing-guide.md` before the first dossier in a session. Use it for the method-stage map, the exact formulas (d1/d2, Greeks, parity, bounds), the model-risk rules, the report blueprint, and the appendix requirements.
4. Compute with the bundled `scripts/run_pricing.py`. It is pure math and needs no network or market data: it returns the BSM price, all Greeks, binomial-tree and Monte-Carlo prices with convergence/standard-error diagnostics, implied volatility, and the validation checks. Install its dependencies with `pip install numpy scipy`.
5. Cross-check before reporting. For European options, the analytical, tree, and Monte-Carlo prices must agree within tolerance; put-call parity must hold; the price must respect intrinsic-value and upper bounds. Treat any disagreement as a model-risk flag, not a rounding nuisance.
6. Produce Markdown by default. If the user asks for Word, PDF, or a polished deliverable, generate the analytical content here first, then use the relevant document skill for final layout.

## Analysis Rules

- Separate inputs, derived quantities, and judgment. Label every derived value such as `d1`, `d2`, price, each Greek, implied volatility, Monte-Carlo standard error, binomial price at `N` steps, American early-exercise premium, and the put-call-parity residual with its formula and the inputs it used.
- State conventions explicitly. Always say how Greeks are scaled (per unit underlying, per 1% vol, per day, per 1% rate) and whether rates/yields are continuous; mismatched conventions are the most common pricing error.
- Cross-validate every European price. Report analytical vs. binomial vs. Monte-Carlo together with the tolerance and the Monte-Carlo confidence interval; do not present a single method as proof. Verify put-call parity from **independent** European tree call/put prices (a BSM-vs-BSM parity residual is an algebraic identity and proves nothing). Any "一致/已通过/已收敛" statement in the summary must reflect the actual validation outcome, not be hardcoded.
- For American options, report the American (tree) price as the headline value and the European value only as the cross-check anchor plus the early-exercise premium; never label the European BSM price as the American option's price.
- Respect no-arbitrage. Check price ≥ intrinsic value and price ≤ the relevant upper bound, and verify put-call parity; flag any violation as high model risk.
- Treat an inverted implied volatility that fails to converge, or a market price below intrinsic value, as evidence of bad inputs or arbitrage; report it rather than forcing a number.
- Use high/medium/low model-risk levels only when a rule in `references/pricing-guide.md` or a user-provided rule is triggered. Include the triggering rule text beside each flag.
- End every report with this disclaimer: `本报告基于公开数据与规则化分析生成，仅供研究参考，不构成任何投资建议。`

## Resource Guide

- `references/pricing-guide.md`: method-stage map, exact formulas, model-risk rules, report blueprint, and final QA checklist.
- `scripts/run_pricing.py`: runnable backbone — BSM price and Greeks, binomial-tree and Monte-Carlo pricing with diagnostics, implied-volatility inversion, put-call-parity and no-arbitrage checks, and the Markdown report writer. Pure math, no network.

## Quality Bar

- Every material number must trace to a named method/formula and its inputs.
- European prices must be cross-checked across methods within a stated tolerance; American prices must report the early-exercise premium versus the European value.
- Monte-Carlo results must report a standard error and confidence interval; never present a point estimate alone.
- Do not overstate precision. Prefer "在 95% 置信区间内一致", "需提高路径数", and "与解析价存在偏差，疑似实现/输入问题" when results are only indicative, and never use buy/sell language.
