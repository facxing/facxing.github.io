---
title: "D→Pℓ⁺ℓ⁻ 的 1-loop 费曼图与 O₁–O₁₀ 对照"
collection: notes
type: notes
permalink: /notes/d-meson-feynman-diagrams
date: 2026-10-21
venue: "BESIII"
excerpt: "c→uℓ⁺ℓ⁻ 的 1-loop 费曼图与有效算符 O₁–O₁₀ 对照：流结构、长/短程判据、逐图索引与 42 张切片"
author_profile: true
---

**图源**：FeynArts/FeynCalc 生成，42 张（`T1C1N1 … T4C2N42`），下图为切片。

> 图与算符是**多对多**关系：一张图可同时贡献 O₇/O₉/O₁₀，一个算符也由多张图
> 相加。下面按“主要贡献”分组；精确系数需对 1-loop 振幅做投影。
> `d_m` = 下型夸克（d, s, b）的质量本征态，`m` 为代指标，图中隐含三代求和。

---

## 1. 流的结构

记 `P_L = (1−γ⁵)/2`，`P_R = (1+γ⁵)/2`，`σ^{μν} = i[γ^μ,γ^ν]/2`，
`T^a` 为 SU(3) 生成元（`tr(T^aT^b)=δ^{ab}/2`）。

### 1.1 有效算符（µ ~ m_c）

| 算符 | 显式结构 | 流类型 |
|---|---|---|
| **O₁** | `(ū_L γ^μ d_L)(d̄_L γ_μ c_L)` | 色单态 V−A |
| **O₂** | `(ū_L γ^μ T^a d_L)(d̄_L γ_μ T^a c_L)` | 色八重态 V−A |
| **O₃** | `(ū_L γ^μ c_L) Σ_q (q̄_L γ_μ q_L)` | V−A |
| **O₄** | `(ū_L γ^μ T^a c_L) Σ_q (q̄_L γ_μ T^a q_L)` | 色八重态 V−A |
| **O₅** | `(ū_L γ^μ c_L) Σ_q (q̄_R γ_μ q_R)` | V+A |
| **O₆** | `(ū_L γ^μ T^a c_L) Σ_q (q̄_R γ_μ T^a q_R)` | 色八重态 V+A |
| **O₇** | `(e/16π²) m_c ū_L σ^{μν} c_R F_{μν}` | 磁偶极（螺旋翻转） |
| **O₈** | `(g/16π²) m_c ū_L σ^{μν} T^a c_R G^a_{μν}` | 色磁偶极 |
| **O₉** | `(α_em/4π)(ū_L γ^μ c_L)(ℓ̄ γ_μ ℓ)` | 矢量 |
| **O₁₀** | `(α_em/4π)(ū_L γ^μ c_L)(ℓ̄ γ_μ γ⁵ ℓ)` | 轴矢 |

> 约定：JHEP 08 (2025) 026 把 `α_em/4π` 吸收进 O₉ 的定义；Fajfer
> hep-ph/0106333 则把它放在拉氏量前。两者只差一个常数因子。

### 1.2 图级别的流

```
γ-penguin:   ū [ γ^μ(A+Bγ⁵) + iσ^{μν}q_ν(C+Dγ⁵) ] c  ×  (e/q²) ℓ̄γ_μℓ
Z-penguin:   ū γ^μ(v_c−a_cγ⁵) c  ×  ℓ̄γ_μ(v_ℓ−a_ℓγ⁵) ℓ  / (q²−m_Z²)
W-box:       (ūγ^μP_L d)(d̄γ_μP_L c) × (ℓ̄γ^νP_Lν)(ν̄γ_νP_Lℓ) × W-propagators
```

W 顶点：`W-ℓ-ν: (g/√2)γ^μP_L`，`W-u-d: (g/√2)V_ud γ^μP_L`。
每个 `d_m` 代带 `V_cd V_ud*`（s、b 代同理），三代求和给出 `λ_d, λ_s, λ_b`。

### 1.3 强子矩阵元

```
⟨P|ūγ^μc|D⟩ = (p+p')^μ f₊(q²) + (p−p')^μ f₋(q²)          → O₉, O₁₀
⟨P|ūσ^{μν}c|D⟩ = i s(q²)[(p+p')^μq^ν − q^μ(p+p')^ν ± ...]  → O₇
```

### 1.4 长程 vs 短程：判据与归属

**关键**：这 42 张图都是 `c→uℓℓ` 的**匹配图**；“长程/短程”不是按拓扑分的，
而是按**动量区域**分的：

| 区域 | 动量 | 去向 | 类型 |
|---|---|---|---|
| **硬** | `≫ Λ_QCD`，~`m_W, m_t, m_b` | 匹配到 Wilson 系数 `C_i(µ~m_W)` | **短程（SD）** |
| **软** | `~Λ_QCD`，轻夸克、共振 | 进入非局域矩阵元 `⟨P|O₁,₂|D⟩`（JHEP 的 `A(q²)`） | **长程（LD）** |

同一张图（例如 `γ` 从 `d_m` 发出）**同时有硬、软两部分**：硬部进
`C₇/C₉/C₁₀`，软部进 `A(q²)`。所以不能把某张图整体标成“LD”或“SD”。

**拓扑层面的实用归属**

| 拓扑 | 硬部 → SD | 软部 → LD |
|---|---|---|
| `γ/Z` 从 `d_m` 或 `W`（N11–N22） | `C₇, C₉, C₁₀` | 非局域 `D→Pγ*` 振幅（ρ/ω/φ + 连续谱） |
| Goldstone/W 交换（N3, N7, N11–N14, N19–N22） | 电弱匹配（与 W 纵向部分合并），**不进 C₃–C₆** | 软部 → 强子矩阵元/波函数 |
| W-box（N23–N26，含轻子环） | `C₉, C₁₀` | 若含轻夸克，软部进非局域振幅 |
| 自能（N27–N42） | 场强/波函数重整化 | 介子波函数 |
| `H/G⁰`（N1–N10, N27–N29, N35–N38） | `H` 企鹅（~`m_d/m_W`，可忽略） | `G⁰` 是规范自由度，与 `Z` 纵向抵消 |

> 注：图中 `G` 是**带电 Goldstone**（W 的纵向分量），`G⁰` 是中性 Goldstone；
> 胶子应标为小写 `g`。本组 42 张图**没有胶子线**。

**代表图**

同一张图可以同时出现在 SD 与 LD 栏：区别在**动量区域**，不在拓扑。

| SD（硬部） | LD（软部） | 规范/可忽略 |
|---|---|---|
| ![N11](/assets/images/notes/D2PieeFeynman/N11.png) | ![N11](/assets/images/notes/D2PieeFeynman/N11.png) | ![N1](/assets/images/notes/D2PieeFeynman/N1.png) |
| N11：`γ` 从 `d_m`，硬部 → `C₇/C₉` | 同一图的软部 → `A(q²)` | N1：`H`，~`m_d/m_W` |
| ![N12](/assets/images/notes/D2PieeFeynman/N12.png) | ![N3](/assets/images/notes/D2PieeFeynman/N3.png) | ![N2](/assets/images/notes/D2PieeFeynman/N2.png) |
| N12：`Z` 从 `d_m` → `C₉/C₁₀` | N3：Goldstone/W 交换，软部 → 强子矩阵元 | N2：`G⁰`，与 `Z` 纵向抵消 |
| ![N23](/assets/images/notes/D2PieeFeynman/N23.png) | ![N30](/assets/images/notes/D2PieeFeynman/N30.png) | |
| N23：W-box → `C₉/C₁₀` | N30：自能 → 波函数重整化 | |

读图提示：`d_m` 是内部下型夸克（d/s/b 求和）；`γ/Z` 从 `d_m` 或 `W`
发出。硬部（loop momentum ~`m_W`）进 `C_i`，软部（~`Λ_QCD`）进
非局域振幅。

**逐图速查**

| 图 | 主要归属 | 说明 |
|---|---|---|
| N1–N10 | 规范/可忽略 | `H/G⁰`，不进 O₁–O₁₀ |
| N11–N22 | **混合**：硬 SD + 软 LD | `γ/Z` 企鹅；硬部→`C₇/C₉/C₁₀`，软部→LD 振幅 |
| N23–N26 | **SD**（W-box，含轻子环） | 硬部→`C₉/C₁₀`；夸克 W-box 另补 |
| N27–N42 | 重整化/波函数 | 不独立分类 |

**与文献的对应**

- **Fajfer 的 “SD”**（式 23，`c₇, c₉, c₁₀`）= 上表的**硬部**；
- **JHEP 的 “LD”**（式 7.8，`A(q²)`）= 上表的**软部**（非局域 `D→Pγ*`）；
- **JHEP 的 “SD O₉”**（式 2.14）= 硬部中的 **`λ_b` 点状项**；
- Fajfer 的 `c₉ ln(m_s/m_d)` 属于**软部**（轻夸克圈的低能部分），所以 JHEP 把它归入 LD。

> 严格说，LD 的“主体”不是这 42 张图中的某几张，而是 `O₁,₂` 的**强子矩阵元**
> `⟨P|O₁,₂|D⟩`；这些图只是通过 OPE 的软区域与之相连。

---

## 2. 按算符分组的图

### O₁, O₂（树图；本组图不含 QCD 修正）

O₁,₂ 本身是**树图**四夸克算符，不在这 42 张 1-loop 图里。本组图**没有胶子线**，
所以 O₁,₂ 的 QCD 修正/混合图**缺失**；Goldstone/W 交换图见 §1.4 的“代表图”，
不要把它们当成 O₁,₂ 的 QCD 修正。

### O₃–O₆（QCD 企鹅）

N11–N14, N19–N22 中的 `G` 是**带电 Goldstone**，不是胶子；本组图里
**没有** QCD 胶子交换/胶子企鹅，因此 **O₃–O₆ 的图全部缺失**，见第 4 节。
不要用 N11/N12/N19 代替它们。

### O₇（磁偶极）

`γ` 从内部 `d_m`（含螺旋翻转）或 `W` 发出：N11, N13, N15, N17, N19, N21。
其中 `γ` 从 `d_m` 发出的图给出 O₇ 的主要贡献；`γ` 从 `W` 发出的图
同时贡献 O₇ 与 O₉。

| | | |
|---|---|---|
| ![N11](/assets/images/notes/D2PieeFeynman/N11.png) | ![N13](/assets/images/notes/D2PieeFeynman/N13.png) | ![N15](/assets/images/notes/D2PieeFeynman/N15.png) |
| N11：γ 从 d_m | N13：γ 从 d_m | N15：γ 从 W |
| ![N17](/assets/images/notes/D2PieeFeynman/N17.png) | ![N19](/assets/images/notes/D2PieeFeynman/N19.png) | ![N21](/assets/images/notes/D2PieeFeynman/N21.png) |
| N17：γ 从 W | N19：γ 从 W + G | N21：γ 从 W + W |

### O₉（矢量）

`γ/Z` 企鹅（N11–N22）+ W-box（N23–N26，含轻子环）+ 自能（N30–N34, N39–N42）。

| | | |
|---|---|---|
| ![N11](/assets/images/notes/D2PieeFeynman/N11.png) | ![N12](/assets/images/notes/D2PieeFeynman/N12.png) | ![N15](/assets/images/notes/D2PieeFeynman/N15.png) |
| N11：γ 从 d_m | N12：Z 从 d_m | N15：γ 从 W |
| ![N16](/assets/images/notes/D2PieeFeynman/N16.png) | ![N23](/assets/images/notes/D2PieeFeynman/N23.png) | ![N30](/assets/images/notes/D2PieeFeynman/N30.png) |
| N16：Z 从 W | N23：W-box | N30：c 腿自能 + γ |

### O₁₀（轴矢）

`Z` 从内部 `d_m` 或 `W` 发出：N12, N14, N16, N18, N20, N22；
W-box N23–N26；自能 N31, N32, N34, N40, N42。

| | | |
|---|---|---|
| ![N12](/assets/images/notes/D2PieeFeynman/N12.png) | ![N14](/assets/images/notes/D2PieeFeynman/N14.png) | ![N16](/assets/images/notes/D2PieeFeynman/N16.png) |
| N12：Z 从 d_m | N14：Z 从 d_m + G | N16：Z 从 W |
| ![N24](/assets/images/notes/D2PieeFeynman/N24.png) | ![N31](/assets/images/notes/D2PieeFeynman/N31.png) | ![N40](/assets/images/notes/D2PieeFeynman/N40.png) |
| N24：W-box | N31：c 腿自能 + Z | N40：u 腿自能 + Z |

### O₈（色磁偶极）

图中**没有**明确的胶子从内部夸克发出（螺旋翻转）的图。需要补
`g` 从 `d_m` 线发出、再被另一条夸克线吸收的图；O₈ 通过混合影响 O₇/O₉。

---

## 3. H/G⁰ 图（不属于 O₁–O₁₀）

`G⁰` 是规范自由度，与 `Z` 纵向部分抵消；`G±` 与 `W` 纵向抵消；`H` 是 Higgs 企鹅，压低
~`m_d/m_W`。这些图只作规范检验，不计入算符系数。

| | | |
|---|---|---|
| ![N1](/assets/images/notes/D2PieeFeynman/N1.png) | ![N2](/assets/images/notes/D2PieeFeynman/N2.png) | ![N4](/assets/images/notes/D2PieeFeynman/N4.png) |
| N1：H | N2：G⁰ | N4：H + W |
| ![N5](/assets/images/notes/D2PieeFeynman/N5.png) | ![N28](/assets/images/notes/D2PieeFeynman/N28.png) | ![N37](/assets/images/notes/D2PieeFeynman/N37.png) |
| N5：G⁰ + W | N28：c 腿自能 + H | N37：u 腿自能 + H |

其余 H/G⁰ 图：N3, N6, N7, N8, N9, N10, N27, N29, N35, N36, N38。

---

## 4. 是否有缺漏

| 缺漏 | 影响 | 建议 |
|---|---|---|
| **抵消项** | 1-loop 振幅 UV 发散 | 打开 FeynArts `CounterTerms` |
| **QCD 顶点修正**（胶子连外腿或 W） | O₁,₂ → O₇–O₁₀ 的混合与跑动 | 补胶子-外腿、胶子-W 顶点 |
| **夸克 W-box** | O₉/O₁₀ 的 box 匹配、GIM | N23–N26 是含轻子环的 W-box；夸克 W-box 仍缺 |
| **γ–Z / W–γ 混合** | 中性玻色子 off-shell | 加混合自能 |
| **胶子企鹅（O₃–O₆）** | QCD 企鹅系数 | 补胶子从内部夸克发出再吸收 |
| **色磁偶极（O₈）** | 通过混合影响 O₇/O₉ | 补胶子螺旋翻转图 |
| **末态 QED 辐射** | 不属于短距离 | PHOTOS/QED 修正单独处理 |

---

## 5. 精确投影步骤

1. 在 `µ ~ m_W` 计算 `M = M_tree + M_loop`；
2. 用 FeynCalc 把 `M_loop` 投影到 `O₁–O₁₀`（`FCFAConvert` + `Projector`）；
3. 与有效哈密顿量树图矩阵元比较，读出 `C_i(µ~m_W)`；
4. RGE 从 `m_W` 跑到 `m_c`（de Boer–Müller–Seidel NNLL）；
5. 与 JHEP/Fajfer 在 `µ~m_c` 的系数对比。

**规范检验**：改变 `ξ`，`G⁰` 与 `Z`、`G±` 与 `W` 纵向部分应抵消，`C_i` 与 `ξ` 无关。

---

## 6. 参考文献

- **JHEP 08 (2025) 026**, arXiv:2505.21369：O₁–O₁₀ 基、O₉ 点状项、LD 非局域振幅。
- **S. Fajfer, S. Prelovšek, P. Singer**, hep-ph/0106333：O₇, O₉, O₁₀ 与式 (23)。
- **A. Paul, I. Bigi, S. Recksiegel**, arXiv:1101.6053：C₁–C₁₀ 的 NLO 系数。
- **A. Bharucha, D. Boito, C. Méaux**, arXiv:2011.12856：C₉^eff 的 naive/ann/SS/FF 分解。
- **S. de Boer, B. Müller, D. Seidel**, arXiv:1606.05521：c→u 的 NNLL Wilson 系数。

## 7. 全部 42 张图索引

| | | |
|---|---|---|
| ![N1](/assets/images/notes/D2PieeFeynman/N1.png) | ![N2](/assets/images/notes/D2PieeFeynman/N2.png) | ![N3](/assets/images/notes/D2PieeFeynman/N3.png) |
| ![N4](/assets/images/notes/D2PieeFeynman/N4.png) | ![N5](/assets/images/notes/D2PieeFeynman/N5.png) | ![N6](/assets/images/notes/D2PieeFeynman/N6.png) |
| ![N7](/assets/images/notes/D2PieeFeynman/N7.png) | ![N8](/assets/images/notes/D2PieeFeynman/N8.png) | ![N9](/assets/images/notes/D2PieeFeynman/N9.png) |
| ![N10](/assets/images/notes/D2PieeFeynman/N10.png) | ![N11](/assets/images/notes/D2PieeFeynman/N11.png) | ![N12](/assets/images/notes/D2PieeFeynman/N12.png) |
| ![N13](/assets/images/notes/D2PieeFeynman/N13.png) | ![N14](/assets/images/notes/D2PieeFeynman/N14.png) | ![N15](/assets/images/notes/D2PieeFeynman/N15.png) |
| ![N16](/assets/images/notes/D2PieeFeynman/N16.png) | ![N17](/assets/images/notes/D2PieeFeynman/N17.png) | ![N18](/assets/images/notes/D2PieeFeynman/N18.png) |
| ![N19](/assets/images/notes/D2PieeFeynman/N19.png) | ![N20](/assets/images/notes/D2PieeFeynman/N20.png) | ![N21](/assets/images/notes/D2PieeFeynman/N21.png) |
| ![N22](/assets/images/notes/D2PieeFeynman/N22.png) | ![N23](/assets/images/notes/D2PieeFeynman/N23.png) | ![N24](/assets/images/notes/D2PieeFeynman/N24.png) |
| ![N25](/assets/images/notes/D2PieeFeynman/N25.png) | ![N26](/assets/images/notes/D2PieeFeynman/N26.png) | ![N27](/assets/images/notes/D2PieeFeynman/N27.png) |
| ![N28](/assets/images/notes/D2PieeFeynman/N28.png) | ![N29](/assets/images/notes/D2PieeFeynman/N29.png) | ![N30](/assets/images/notes/D2PieeFeynman/N30.png) |
| ![N31](/assets/images/notes/D2PieeFeynman/N31.png) | ![N32](/assets/images/notes/D2PieeFeynman/N32.png) | ![N33](/assets/images/notes/D2PieeFeynman/N33.png) |
| ![N34](/assets/images/notes/D2PieeFeynman/N34.png) | ![N35](/assets/images/notes/D2PieeFeynman/N35.png) | ![N36](/assets/images/notes/D2PieeFeynman/N36.png) |
| ![N37](/assets/images/notes/D2PieeFeynman/N37.png) | ![N38](/assets/images/notes/D2PieeFeynman/N38.png) | ![N39](/assets/images/notes/D2PieeFeynman/N39.png) |
| ![N40](/assets/images/notes/D2PieeFeynman/N40.png) | ![N41](/assets/images/notes/D2PieeFeynman/N41.png) | ![N42](/assets/images/notes/D2PieeFeynman/N42.png) |
