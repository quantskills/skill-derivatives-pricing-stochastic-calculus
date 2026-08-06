# 衍生品定价与随机微积分

这是一个可复现的 Agent Skill，用于期权与其他衍生品的定价、希腊字母、隐含波动率和模型风险校验。它交叉使用 Black-Scholes-Merton、二叉树与蒙特卡洛方法，并报告标准误、看跌看涨平价和无套利边界检查。

详见 [`SKILL.md`](SKILL.md)。本工具离线运行，不提供行情数据、交易信号或投资建议。

```bash
pip install numpy scipy
python scripts/run_pricing.py --type call --S 100 --K 100 --T 1 --r 0.03 --q 0 --sigma 0.2 --out pricing-report.md
```

运行时入口：Codex / Claude Code 使用 `SKILL.md`；Cursor 使用 `agents/cursor-rule.mdc`；Hermes / OpenClaw 使用 `agents/portable-loader.md`。

GPL-3.0-only，详见 [`LICENSE`](LICENSE)。
