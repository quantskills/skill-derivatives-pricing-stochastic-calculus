# Derivatives Pricing Guide

Read this guide when generating or revising a derivatives-pricing dossier. Use it as a compact operating checklist, not as a replacement for a derivatives textbook or the exact library documentation.

## Default Scope

- Target: one vanilla option (European or American call/put) on a single underlying; extend to a small set of strikes/maturities when a smile or surface is requested.
- Model: Black-Scholes-Merton with continuous risk-free rate `r` and continuous dividend/carry yield `q`; geometric Brownian motion for the underlying.
- Conventions: volatility annualized; vega per 1 vol point (i.e. per 0.01 of sigma), theta per calendar day, rho per 1% rate change, delta and gamma per unit of underlying.
- Numerical defaults: binomial (CRR) tree with ~500 steps; Monte-Carlo with ~200,000 antithetic paths; implied volatility solved on `[1e-4, 5.0]`.
- Output: Markdown report unless the user requests HTML, Word, PDF, or another deliverable.

## Method Stage Map

Run all of this with the bundled `scripts/run_pricing.py` (depends only on `numpy`, `scipy`).

| Stage | Method | Use |
|---|---|---|
| Inputs & assumptions | parse contract; collect S, K, T, r, q; sigma or market price | Fix the contract and every pricing input; state conventions. |
| Analytical pricing | Black-Scholes-Merton closed form | No-arbitrage price of a European option. |
| Greeks | analytic delta, gamma, vega, theta, rho | Risk sensitivities for hedging and exposure. |
| Numerical pricing | CRR binomial tree; Monte-Carlo (GBM, antithetic) | American/path-dependent payoffs and independent cross-checks. |
| Implied volatility | invert BSM via Brent's method | Recover the volatility implied by an observed market price. |
| Volatility structure | smile/skew across strikes; term structure across maturities | Inspect how implied vol varies with strike and maturity. |
| Model-risk validation | put-call parity, no-arbitrage bounds, cross-method consistency, MC standard error | Catch implementation errors, bad inputs, and arbitrage. |

## Formulas To Derive

State the formula and inputs used whenever a value is derived rather than directly given.

- **d1 / d2**: `d1 = [ln(S/K) + (r − q + σ²/2)·T] / (σ·√T)`, `d2 = d1 − σ·√T`.
- **BSM price**: call `= S·e^(−qT)·N(d1) − K·e^(−rT)·N(d2)`; put `= K·e^(−rT)·N(−d2) − S·e^(−qT)·N(−d1)`.
- **Greeks** (call shown; put analogous): delta `= e^(−qT)·N(d1)`; gamma `= e^(−qT)·φ(d1)/(S·σ·√T)`; vega `= S·e^(−qT)·φ(d1)·√T`; theta and rho per the standard BSM expressions. Report vega/theta/rho in the conventions above and say so.
- **Implied volatility**: the `σ` solving `BSM(σ) = market price`; report non-convergence rather than a forced value.
- **Monte-Carlo**: price `= e^(−rT)·mean(payoff)`; standard error `= e^(−rT)·std(payoff)/√n`; 95% CI `= price ± 1.96·SE`.
- **Binomial convergence**: report price at `N` steps and confirm it has stabilized (e.g. change from `N/2` to `N` below tolerance).
- **American early-exercise premium**: `American price − European price` (should be ≥ 0; typically > 0 for American puts and for American calls on dividend-paying underlyings).
- **Put-call parity**: residual `= C − P − (S·e^(−qT) − K·e^(−rT))`; should be ~0. **Compute `C` and `P` from an *independent* engine (e.g. the European binomial tree), not from the same closed-form BSM pair.** Closed-form BSM call/put satisfy parity identically because `N(x)+N(−x)≡1`, so a BSM-vs-BSM residual is always ~0 and proves nothing about the implementation; the tree-based residual can actually catch a bug.
- **No-arbitrage bounds**: price ≥ intrinsic value; call ≤ `S·e^(−qT)`, put ≤ `K·e^(−rT)`.

## Model-Risk Rules

Use these defaults unless the user supplies thresholds. If an input is missing, downgrade the rule to a qualitative note and say what is missing.

| Level | Trigger |
|---|---|
| High | Cross-method disagreement: analytical price lies outside the Monte-Carlo 95% CI, or |binomial − analytical|/analytical exceeds the tolerance (default 0.5%) — likely an implementation or input error. |
| High | Put-call-parity residual exceeds tolerance (default 0.5% of spot), computed from *independent* European tree call/put prices — inconsistent pricing. |
| High | No-arbitrage violation: price below intrinsic value or above the upper bound, or an observed market price below intrinsic. |
| High | Implied-volatility inversion fails to converge or returns a boundary value — bad inputs or an arbitrageable quote. |
| Medium | Monte-Carlo standard error above ~1% of price — increase paths before trusting the estimate. |
| Medium | Very short maturity (e.g. T < 0.02y) or deep ITM/OTM (|ln(S/K)| large) — numerical instability and unreliable Greeks/implied vol. |
| Medium | Stale or assumed `r`/`q`: dividend yield or rate guessed rather than sourced. |
| Medium | Pronounced volatility skew across strikes while a single flat sigma is used for all of them. |
| Low | Minor rounding, a single-method note, or a borderline diagnostic with limited impact; record in the appendix rather than the headline flag list. |

For combined signals, name the combination explicitly, for example `MC标准误偏大 + 跨方法不一致` or `深度价外 + 隐含波动率不收敛`.

## Report Blueprint

Use this chapter order unless the user asks for a custom structure:

1. `摘要与结论`: three to six bullets covering the contract, the price (with method agreement), the key Greeks, and the main model-risk flags.
2. `合约与市场输入`: type, exercise style, S, K, T, r, q, sigma (or market price), and the stated conventions.
3. `解析定价`: d1/d2, BSM price, and the no-arbitrage interpretation.
4. `希腊字母`: delta, gamma, vega, theta, rho with their conventions and a one-line risk interpretation each.
5. `数值定价与收敛`: binomial price and convergence, Monte-Carlo price with SE and 95% CI, and (if American) the early-exercise premium.
6. `隐含波动率`: implied vol from the market price (if given) versus the input sigma, and convergence status.
7. `波动率结构`: smile/skew across strikes or term structure across maturities, when multiple contracts are supplied.
8. `模型风险与校验清单`: table with level, signal, triggering rule, evidence, and the method/formula used.
9. `方法附录`: method-by-method table with inputs, outputs, tolerances, step/path counts, and caveats.

## Evidence And Output Requirements

- Include at least one source/method table in the appendix with columns similar to: `计算阶段`, `方法`, `关键输入`, `输出`, `容差/参数`, `备注`.
- For each model-risk flag, include the method, the numeric residual/error, and the tolerance in the same row or the immediately following sentence.
- Always report Monte-Carlo results with a standard error and confidence interval; never a point estimate alone.
- Always show the cross-method comparison for European options.
- Prefer concise tables over long prose when comparing methods or Greeks.
- Keep the final tone analytical and non-promotional; avoid buy/sell language.
- For an American option, the headline "理论价格" must be the **American** price (binomial tree), not the European BSM price; report the European value only as the cross-check anchor and as the basis for the early-exercise premium. Do not mix an American tree price and a European Monte-Carlo price in the same comparison without labelling which instrument each prices.
- Any "一致 / 已通过 / 已收敛" wording in the summary must be **driven by the actual validation result**, never asserted unconditionally; if a check fails, the summary must say so and point to the model-risk chapter.
- Use a combined relative-and-absolute tolerance for cross-method and parity checks so that near-zero prices (deep OTM, near expiry) do not produce spurious high-risk flags.

## Final QA Checklist

- Contract and conventions are stated and displayed consistently.
- d1/d2 and every Greek have formulas and inputs.
- European price is cross-checked across analytical, binomial, and Monte-Carlo within a stated tolerance.
- Put-call parity and no-arbitrage bounds are checked and reported.
- Monte-Carlo results carry a standard error and confidence interval.
- American results report the early-exercise premium versus the European value.
- High/medium/low flags cite the triggering rule.
- Non-convergence or bad inputs are disclosed rather than forced into a number.
- Final disclaimer is present exactly as required by `SKILL.md`.
