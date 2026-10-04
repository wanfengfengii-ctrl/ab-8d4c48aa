# Nanopore Current-Trace Aligner

将离子电流观测轨迹**联合**重新对齐到已知参考电平的服务：不预先估计基线漂移、不预先分段，而是一次性选择整数漂移、首尾必用的参考子序列与连续采样边界，避免把短暂停留误判为碱基跳过。

## 问题定义

给定：

- **参考电平** `levels`：8–24 个整数（可重复，首尾必须被采用）；
- **观测值** `observations`：8–60 个数（整数或浮点）；
- **漂移闭区间** `drift.{min,max}`：整数漂移 `d` 的取值范围（闭区间）；
- **残差上限** `residual_limit`：每个观测的绝对残差不得超过该值；
- **停留范围** `dwell.{min,max}`：每个采用电平占 1–3 个采样（默认 1–3）；
- **内部跳过上限** `max_internal_skips`：0–2（默认 2）。

服务联合选择：整数漂移 `d`、包含首尾参考电平的子序列（内部跳过 ≤ 上限）、以及把观测序列连续划分为每级一段的边界，使得：

- 每个观测恰好归属一个采用电平；
- 每段长度在停留范围内；
- 每个残差 `obs - (level + d)` 的绝对值 ≤ 残差上限。

在所有合法对齐中，按如下**字典序**依次最小化：

1. 跳过数（未采用的参考电平个数）；
2. 残差绝对值总和；
3. 最大绝对残差；
4. 漂移值；
5. 边界序列（各段起始下标组成的序列，按字典序比较）。

求解为对每个整数漂移运行的精确动态规划，复杂度约
`O(漂移数 × 电平数 × 观测数 × 停留选择 × 跳过预算²)`，对本任务的规模（≤1001 个漂移 × 24 电平 × 60 观测）绰绰有余。

## API

### `POST /api/current-traces/align`

请求体：

```json
{
  "levels": [10, 20, 30, 40, 50, 60, 70, 80],
  "observations": [11, 12, 21, 22, 31, 32, 51, 52, 61, 62, 71, 72, 81, 82],
  "drift": { "min": -2, "max": 2 },
  "residual_limit": 3,
  "dwell": { "min": 1, "max": 3 },
  "max_internal_skips": 2
}
```

成功（HTTP 200，`status: "aligned"`）返回漂移、逐级采样区间与残差证据：

```json
{
  "status": "aligned",
  "drift": 1,
  "skips": 1,
  "skipped_level_indices": [3],
  "boundaries": [0, 2, 4, 6, 8, 10, 12],
  "segments": [
    {"level_index": 0, "level": 10, "expected": 11, "start": 0, "end": 2, "length": 2}
  ],
  "residuals": [
    {"index": 0, "level_index": 0, "observed": 11, "expected": 11, "residual": 0, "abs_residual": 0}
  ],
  "objective": {
    "skips": 1, "sum_abs_residual": 7, "max_abs_residual": 1,
    "drift": 1, "boundaries": [0, 2, 4, 6, 8, 10, 12]
  }
}
```

无解（HTTP 200，`status: "infeasible"`）返回明确结论与原因（结构性不可行 / 漂移区间内无满足残差上限的对齐）：

```json
{
  "status": "infeasible",
  "reason": "no integer drift in [5, 6] admits an alignment with all absolute residuals <= 3",
  "detail": { "...": "..." }
}
```

字段越界或序列规模非法（电平数 ∉ [8,24]、观测数 ∉ [8,60]、漂移区间倒置或过宽、残差上限为负、停留范围超出 1–3、跳过数 > 2、携带未知字段等）→ **HTTP 422** 拒绝。

### 字段范围

| 字段 | 约束 |
| --- | --- |
| `levels` | 8–24 个严格整数，每个 \|值\| ≤ 10⁶ |
| `observations` | 8–60 个有限数值，每个 \|值\| ≤ 10⁹ |
| `drift.min/max` | 严格整数，\|值\| ≤ 10⁶，`min ≤ max`，宽度 ≤ 1000 |
| `residual_limit` | 0 ≤ 值 ≤ 10⁹ |
| `dwell.min/max` | 严格整数，1 ≤ min ≤ max ≤ 3（默认 1–3） |
| `max_internal_skips` | 严格整数，0–2（默认 2） |

### 其他端点

- `GET /health` → `{"status": "ok"}`（Compose 健康检查使用）
- `GET /` → 服务信息
- `GET /docs` → Swagger UI

## 运行（Docker）

```bash
# 构建并启动 API（宿主机端口可用 API_PORT 覆盖，默认 8000）
docker compose up --build api
API_PORT=9000 docker compose up api

curl -s -X POST http://localhost:8000/api/current-traces/align \
  -H 'Content-Type: application/json' \
  -d '{"levels":[10,20,30,40,50,60,70,80],
       "observations":[11,12,21,22,31,32,51,52,61,62,71,72,81,82],
       "drift":{"min":-2,"max":2},"residual_limit":3,
       "dwell":{"min":1,"max":3},"max_internal_skips":2}'
```

## 一次性验证（verify 服务）

`verify` 服务等待 API 健康后依次执行：**包构建**（`pip install .`）、**代码测试**（`pytest`）、**HTTP 冒烟**（可行轨迹对齐成功、无解轨迹返回无解结论、非法请求被拒绝），并以退出码汇报结果（全部通过为 0，否则非 0）：

```bash
docker compose up --build --exit-code-from verify verify
echo $?   # 0 = 全部通过
```

## 本地开发

```bash
pip install -e ".[test]"
pytest
uvicorn nanopore_aligner.main:app --reload
API_BASE_URL=http://localhost:8000 python verify/verify.py
```

## 仓库结构

```
├── Dockerfile              # 应用镜像（API 与 verify 共用）
├── docker-compose.yml      # api（健康检查 + 可配置宿主机端口）与 verify
├── pyproject.toml          # 包元数据与依赖
├── src/nanopore_aligner/
│   ├── solver.py           # 联合对齐动态规划（纯 Python，无第三方依赖）
│   ├── schemas.py          # 请求/响应模型与字段校验
│   └── main.py             # FastAPI 应用与端点
├── tests/                  # 求解器与 API 测试
└── verify/verify.py        # 一次性验证脚本（构建 + 测试 + HTTP 冒烟）
```
