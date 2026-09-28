---
title: "WSL + MadGraph5_aMC@NLO 安装与使用教程"
collection: notes
type: notes
permalink: /notes/wsl-madgraph-tutorial
date: 2026-09-28
venue: "CMS"
hide_pubinfo: true
excerpt: "在 Windows WSL2 (AlmaLinux 9) 上从零安装 MadGraph5_aMC@NLO：矩阵元/部分子簇射工作流、MG5 语法、以及以 Bhabha 散射为例的 Feynman 规则与解析截面对照"
author_profile: true
---

{% include base_path %}
<div class="toc-inline">
{% include toc title="目录" %}
</div>

<p><em>Facxing &nbsp;|&nbsp; 2026年9月28日</em></p>

本教程的目标是：在一台只有 Windows 的机器上，用 **WSL2 + AlmaLinux 9** 搭出一套与 CMS 生产环境接近的
事例产生（event generation）工作环境，安装 **MadGraph5_aMC@NLO**，并以 **Bhabha 散射
$$e^+e^-\to e^+e^-$$** 为例走完 "理论 → Feynman 规则 → 解析截面 → Monte Carlo → 对比" 的完整闭环。

<div class="dy-box note-box" markdown="1">

### 本教程已在以下环境实测

| 项目 | 内容 |
|------|------|
| WSL 发行版 | AlmaLinux 9.6 (WSL2) |
| 系统 ROOT | 6.32.00 (`/usr/bin/root-config`) |
| Python | `/usr/bin/python3` 3.9.21（**注意**：不是 miniforge 的 3.13，见 §2.5） |
| MadGraph5_aMC@NLO | **2.9.18** (2023-12-08) |
| LHAPDF6 | 6.5.4（由 MG5 自带安装器编译） |
| MadAnalysis5 | 1.11.0 (2025/04/23) |
| 编译链 | gcc/gfortran 11.5.0 |

</div>

---

## 0. 为什么需要 MadGraph：从矩阵元到事例

### 0.1 一条 "从理论到数据" 的流水线

高能物理的模拟链条可以用**因子化**写成一个卷积。对强子对撞机：

$$
\sigma(pp\to X)=\sum_{a,b}\int_0^1\!\mathrm{d}x_1\int_0^1\!\mathrm{d}x_2\;
f_{a/p}(x_1,\mu_F^2)\,f_{b/p}(x_2,\mu_F^2)\;
\hat{\sigma}_{ab\to X}\!\left(x_1x_2s,\mu_R,\mu_F\right)
+\mathcal{O}\!\left(\frac{\Lambda_{\mathrm{QCD}}^2}{Q^2}\right)
$$

其中 $$f_{a/p}(x,\mu_F^2)$$ 是部分子分布函数（PDF），$$\hat{\sigma}$$ 是**部分子级截面**，
它由硬散射矩阵元 $$\lvert\mathcal{M}\rvert^2$$ 给出。对 $$e^+e^-$$ 对撞（Bhabha 散射）没有
PDF，上式退化成 $$\sigma=\int\mathrm{d}\Omega\;\mathrm{d}\sigma/\mathrm{d}\Omega$$，
这就是我们在 §4 要精确对照的对象。

真实实验里我们需要的不是总截面，而是**一堆事例（events）**：每个事例是一组末态粒子的
四动量，形如

```text
<event>
 4   1  +2.2683971e+05  5.29e+00  7.55e-03  2.09e-01
 -11 -1  0 0 0 0  +0.0000000000e+00 +0.0000000000e+00 +5.2900000000e+00 5.2900000000e+00 ...
  11 -1  0 0 0 0  -0.0000000000e+00 -0.0000000000e+00 -5.2900000000e+00 5.2900000000e+00 ...
 -11  1  1 2 0 0  +3.3949569182e-01 -1.6835286382e+00 +5.0034562054e+00 5.2900000000e+00 ...
  11  1  1 2 0 0  -3.3949569182e-01 +1.6835286382e+00 -5.0034562054e+00 5.2900000000e+00 ...
</event>
```

这就是 **Les Houches Event (LHE)** 格式，也是整条 MC 链条的"通用货币"。

<div class="dy-box note-box" markdown="1">

### 完整工作流

```text
   理论模型 (Lagrangian / UFO)
            │
            ▼
   ① 矩阵元 ME  (Matrix Element)          MadGraph5_aMC@NLO / Powheg-Box / Sherpa
            │     固定阶微扰展开, |M|^2, 输出 LHE
            ▼
   ② 部分子簇射 PS (Parton Shower)         Pythia 8 / Herwig / Sherpa
            │     ISR/FSR, 软共线对数求和, 强子化, MPI/UE, 输出 HepMC
            ▼
   ③ 探测器模拟 (Detector Simulation)      Delphes (参数化, 秒级) / Geant4 (全模拟, 小时级)
            │     几何、材料、能量响应、堆积 (pileup)
            ▼
   ④ 重建 + 事例选择                       CMSSW / NanoAOD / coffea / FCCAnalyses
            │
            ▼
   ⑤ 统计分析                              Combine / RooFit / pyhf  →  极限 / 显著性
```

</div>

### 0.2 名词解释：ME、PS、匹配、强子化

<div class="dy-box note-box" markdown="1">

#### ① ME — 矩阵元 (Matrix Element)

* 用**费曼规则**把费曼图翻译成振幅 $$\mathcal{M}$$，平方并做相空间积分得到 $$\hat{\sigma}$$。
* 是按耦合常数做的**固定阶展开**：
  $$\mathcal{M}=\mathcal{M}_{\text{LO}}+\mathcal{M}_{\text{NLO}}+\cdots$$，
  对应 $$\alpha_s^n$$、$$\alpha_s^{n+1}$$……
  * **LO**（Leading Order）：只有树图；
  * **NLO**：树图 + 一圈虚修正 + 一实发射（如 MC@NLO、POWHEG）；
  * **NNLO**：目前只有少数 $$2\to2$$ 过程（如 $$\gamma\gamma$$、$$t\bar t$$）。
* 红外（软/共线）发散在**足够包容的观测量**中由 **KLN 定理**抵消；剩下的发散进入 PDF 的
  演化方程（DGLAP），这就是为什么 $$\mu_F$$ 只出现在强子初态。
* 引入**重整化标度** $$\mu_R$$ 与**因子化标度** $$\mu_F$$：固定阶结果**依赖于这两个标度**，
  因此 $$\mu_R/\mu_F$$ 的 2 倍变号常被当作"理论不确定度"的下限。
* **关键陷阱**：$$\lvert\mathcal{M}\rvert^2$$ 只描述**几个硬部分子**，它缺少
  任意多软/共线辐射的贡献 → 需要 PS。

#### ② PS — 部分子簇射 (Parton Shower)

* 用**类马尔可夫链**（Markov chain）把软共线对数 $$\ln(Q^2/\Lambda_{\mathrm{QCD}}^2)$$ 求和到所有阶
  （LL/NLL 精度），使每个事例带上任意多胶子/光子。
* 常见实现：
  * **PYTHIA 8**：$$p_\perp$$ 有序 + 偶极（dipole）色结构；末态用 **Lund 弦模型**强子化。
  * **Herwig 7**：角度有序 + 集团（cluster）模型强子化。
  * **Sherpa**：CSS（Catani–Seymour–Seymour）偶极。
* **ISR / FSR**：初态辐射会改变沿束流方向的运动学（并影响 $$\sqrt{s}$$ 与 $$x$$ 的有效值），
  末态辐射会移动能量导致喷注质量改变。
* **强子化 (Hadronization)** 与 **MPI/UE**（多部分子相互作用 / 底色）都在这一步完成。
* **关键陷阱**：LHE 里记录的是**部分子级**四动量，与探测器里看到的强子/喷注**不是同一个东西**；
  比较理论与数据时必须搞清楚自己在哪个层级。

#### ③ 匹配与合并 (Matching & Merging)

把固定阶 ME 和 PS 缝合，不能重复计数：

| 场景 | 常用方案 | 说明 |
|------|----------|------|
| LO 多喷注 + PS | **MLM**, CKKW, CKKW-L, shower-$$k_t$$ | 用 $$k_\perp$$ 判据把 PS 发射归给 parton 或 shower，给出**矩阵元修正** |
| NLO + PS | **MC@NLO**, **POWHEG**, FxFx | MC@NLO 引入 S-event/H-event 与**负权重**；POWHEG 用 FKS 相空间生成"最硬一次发射" |
| NLO 多喷注合并 | **FxFx**（MG5 内置）, UNLOPS, MINLO | 得到 NLO 精度的 inclusive 多喷注样本 |

<div class="zy-box note-box" markdown="1">

**⚠ 陷阱：负权重**

MC@NLO/FxFx 会产生**负权重事例**（`XWGTUP < 0`）。统计时不能简单"数事例"，
必须用 $$\sum w_i$$；否则会得到错误的截面和分布。可用 `set negative_weights` 或事后
`reweight` 检查。

</div>

#### ④ 为什么用 MadGraph5_aMC@NLO

* **自动化**：任意过程（含多粒子末态、共振衰变）的振幅生成全自动，采用递归螺旋度振幅算法，
  速度快；对 $$2\to6$$ 左右的树图仍然可行。
* **模型无关**：使用 **UFO (Universal FeynRules Output)** 标准，FeynRules 里改一下
  Lagrangian 就能跑 BSM 模型、SMEFT、暗物质简化模型。
* **NLO 自动化**：`aMC@NLO` 让 NLO QCD 自动完成（虚拟项 MadLoop，实项 FKS，减除项 FxFx）。
* **与上下游无缝**：一键调用 Pythia 8、MadSpin（共振衰变自旋关联）、MadAnalysis 5、Delphes。
* **可用于 CMS 生产**：`generate_events` 可产出 **gridpack**（一个自解压的资源包），
  在 CMS 的 Private Production 里直接用它 + CMSSW 做大样本生产。

</div>

### 0.3 为什么必须在 Linux / WSL 下运行

<div class="zs-box note-box" markdown="1">

MadGraph5 **官方只支持 UNIX**（Linux / macOS）。原因非常实际：

1. MG5 的核心在生成 **Fortran 77/95 源代码**（`SubProcesses/`、`Source/`），
   之后调用 `gfortran` + `make` 编译成可执行文件。Windows 原生没有这套工具链。
2. 上游/下游工具全是 Linux-only：**ROOT**、**CMSSW**、**LHAPDF**、**HepMC**、**FastJet**、
   **Delphes**、**Pythia 8** 全都按 POSIX 构建，只发布 Linux 二进制。
3. shell 脚本、`Makefile`、符号链接、`LD_LIBRARY_PATH` 这套约定是 HEP 软件的基础设施。
4. 想用 Cygwin / MSYS2 硬凑是可能的，但路径、动态库、符号链接会持续出问题——
   而且 CMS 的生产环境是 **AlmaLinux 9（EL9）**，你用 AlmaLinux 的 WSL 镜像可以**同构**验证。

</div>

**Windows 上的三种方案对比：**

| 方案 | 优点 | 缺点 | 建议 |
|------|------|------|------|
| **WSL2**（推荐） | 真 Linux 内核；启动秒级；与 Windows 文件互访；VSCode 直连 | 与 Windows 共享内存/磁盘；`/mnt/c` 上 I/O 慢 | **首选**，日常开发调试 |
| 虚拟机（VMware/Hyper-V） | 完全隔离，能装图形界面 | 磁盘/内存开销大，文件互传麻烦 | 需要完整桌面环境时 |
| 双系统 / 原生 Linux | 性能最好 | 切换系统痛苦；CUDA/驱动麻烦 | 长期主力机 |

<div class="ts-box note-box" markdown="1">

**为什么推荐 AlmaLinux 9 而不是 Ubuntu？**

CMS 的 Tier-0/Tier-1 生产、`CMSSW` 的官方二进制、以及 `el9` 容器都以
**AlmaLinux 9 / RHEL 9** 为基线，使用 **glibc 2.34、gcc 11、python 3.9**。
用同版本 WSL 的好处：

* 装出来的 `gfortran` 与 `glibc` 版本和 lxplus 一致 → 编译出来的代码行为一致；
* `dnf` 装包（`root`、`pythia8`、`lhapdf`）可直接来自 EPEL/CRB，和 lxplus 用法相同；
* 后续如果要接 **CMSSW el9 容器**（`docker run cmssw/el9`），宿主机环境完全对得上；
* 用 Ubuntu 时经常遇到 `python3.11+` 缺少 `six`、`glibc` 太新导致旧二进制跑不起来等问题。

</div>

---

## 1. 安装 WSL（推荐 AlmaLinux 9）

### 1.1 前置条件

* Windows 10 版本 **2004（内部版本 19041）** 以上，或 Windows 11；
  `Win + R` → `winver` 查看。
* BIOS/UEFI 中开启 **虚拟化**（Intel VT-x / AMD-V）。
* 建议在**管理员 PowerShell** 中执行下述命令。

### 1.2 安装步骤

```powershell
# 0) 查看可用的发行版名称（不同 WSL 版本列表会变，以此为准）
wsl --list --online

# 1) 确保 WSL 内核是最新的，并将默认版本设为 WSL2
wsl --update
wsl --set-default-version 2

# 2) 安装 AlmaLinux 9（名字以第 0 步的输出为准，通常就是 AlmaLinux-9）
wsl --install -d AlmaLinux-9

# 3) 安装完成后首次启动，会提示创建 Linux 用户名与密码
wsl -d AlmaLinux-9
```

<div class="zs-box note-box" markdown="1">

**关于发行版名称。** Microsoft Store 里的发行版列表时有增删。AlmaLinux 官方也在
自己的 wiki 上给出了 WSL 安装说明（见 §1.8 链接）。如果 `--list --online` 里**没有**
AlmaLinux，有两个替代做法：

**(a) 用 tar 镜像手动导入**（AlmaLinux 官方提供的 WSL rootfs）：

```powershell
# 自行从 AlmaLinux 镜像站下载 AlmaLinux-9-WSL.latest.x86_64.tar.xz
mkdir $env:USERPROFILE\WSL\AlmaLinux9
wsl --import AlmaLinux9 $env:USERPROFILE\WSL\AlmaLinux9 `
    $env:USERPROFILE\Downloads\AlmaLinux-9-WSL.latest.x86_64.tar.xz --version 2
wsl -d AlmaLinux9
```

**(b) 先用 Ubuntu 22.04 起步**（`wsl --install -d Ubuntu-22.04`，glibc 2.35、gcc 11，
和 EL9 足够接近），所有后续步骤完全一样，只是把 `dnf` 换成 `apt`。

</div>

### 1.3 首次配置：装好编译工具链

进入 AlmaLinux 终端后（后续所有命令都在 **WSL 内**执行）：

```bash
# ---- 0) 确认自己在一个 EL9 系统里 ----
cat /etc/almalinux-release      # 期望: AlmaLinux release 9.x
gcc --version ; gfortran --version ; python3 --version

# ---- 1) 更新系统 ----
sudo dnf -y update

# ---- 2) 启用 EPEL 与 CRB（Extra Packages for Enterprise Linux / CodeReady Builder）----
sudo dnf -y install epel-release
sudo dnf config-manager --set-enabled crb

# ---- 3) 编译工具链 + Python 依赖 ----
#     gfortran 是 MadGraph 的硬依赖（它要编译 Fortran 矩阵元）
sudo dnf -y groupinstall "Development Tools"
sudo dnf -y install gcc-gfortran gcc-c++ make cmake
sudo dnf -y install python3 python3-devel python3-pip python3-numpy
sudo dnf -y install wget curl tar gzip which file
sudo dnf -y install zlib-devel openssl-devel   # 很多 HEP 工具需要
sudo dnf -y install nano vim git               # 编辑器

# ---- 4) Python 2/3 兼容模块 six ----
#     MG5 与 MadAnalysis5 在启动时都会 import six，缺了直接报错（见 §2.5）
python3 -m pip install --user six numpy
```

<div class="zy-box note-box" markdown="1">

**⚠ 不要在 WSL 里混用多套 Python。**
本机同时装了 **miniforge（Python 3.13）** 和系统 `/usr/bin/python3`（3.9）。
PATH 里 miniforge 在前时，`python3` 会指向 3.13，而 `python3.13` 环境下**没有 `six`**，
MadGraph 会以

```text
MadGraph5_aMC@NLO requires the six module.
The easiest way to install it is to run "pip3 install six --user"
```

的方式启动失败。**结论：用 `/usr/bin/python3` 显式启动 MadGraph**（本教程全程如此）：

```bash
/usr/bin/python3 /path/to/MG5_aMC_v2_9_18/bin/mg5_aMC
# 或者把别名写进 ~/.bashrc
alias mg5='/usr/bin/python3 /home/$USER/MG5atNLO/MG5_aMC_v2_9_18/bin/mg5_aMC'
```

</div>

*（可选但强烈推荐）* 同步一份 ROOT 与 Pythia8，方便后续做分析：

```bash
sudo dnf -y install root root-graf3d-eve root-gui   # EPEL 提供的 ROOT 6.32
sudo dnf -y install pythia8 pythia8-devel           # 也可以在 MG5 里 install pythia8
```

### 1.4 Windows ↔ WSL 的文件互访

| 方向 | 路径 |
|------|------|
| WSL 中访问 C 盘 | `/mnt/c/Users/<你的用户名>/` |
| Windows 中访问 WSL 家目录 | 资源管理器地址栏输入 `\\wsl$\AlmaLinux-9\home\<用户名>\` |
| 在资源管理器当前目录打开 WSL | 在 WSL 里 `explorer.exe .` |

<div class="zy-box note-box" markdown="1">

**⚠ 性能陷阱：不要把工程放在 `/mnt/c/...`。**

WSL2 访问 Windows 盘（9p 文件系统）的 I/O 性能可能比原生 ext4 慢 **5–20 倍**。
MadGraph 编译时会产生**上万个**小文件（`SubProcesses/P1_*/ajob*`），
放在 `/mnt/c` 下会慢到无法忍受。

**正确做法**：把 MG5、过程目录、ROOT 文件都放在 WSL 的 Linux 文件系统里，例如
`~/MG5atNLO/`、`~/work/`。需要和 Windows 交换文件时，用 `cp` / `scp` 单独拷贝。

</div>

### 1.5 限制 WSL 的资源占用（`.wslconfig`）

在 Windows 用户目录 `C:\Users\<用户名>\.wslconfig` 中：

```ini
[wsl2]
# 内存上限。跑 400k 事例的 LO 过程 + Pythia8 大约需要几个 GB
memory=16GB
# 逻辑核心数。MadGraph 的 run_mode=2 会按这个数起进程
processors=8
# 交换文件
swap=8GB
# 关闭 Windows 与 Linux 之间的自动内存回收（减少抖动）
pageReporting=false
```

改完后在 PowerShell 里执行 `wsl --shutdown`，再重新进入生效。
用 `free -h` 与 `nproc` 确认。

### 1.6 把 WSL 迁移到其他盘（C 盘不够用时）

<div class="zs-box note-box" markdown="1">

**为什么会需要**：WSL2 的发行版就是一个 `ext4.vhdx` 虚拟磁盘，默认落在 C 盘
（新版本在 `%LOCALAPPDATA%\wsl\<发行版>\`，旧版本在
`%LOCALAPPDATA%\Packages\<...>\LocalState\`），而且它**只会长大、不会自己缩小**。
装上 MadGraph + Pythia8 + LHAPDF + ROOT，再跑几次大规模编译，10–30 GB 是常态。

**迁移前必须先 `wsl --shutdown`**，否则 vhdx 被占用，命令会失败、甚至损坏磁盘。

</div>

#### 方法 A：一条命令（推荐，需要 WSL ≥ 2.0）

```powershell
# 管理员 PowerShell
wsl --shutdown
wsl --version                        # 确认 WSL 版本 >= 2.0
wsl --manage AlmaLinux-9 --move D:\WSL\AlmaLinux-9
```

<div class="ts-box note-box" markdown="1">

`wsl --manage` **没有**出现在 Microsoft Learn 的 "Basic commands for WSL" 页面里
（那一页只收录了 `--export / --import / --unregister` 等），所以**以本机 `wsl --help`
的输出为准**：

```powershell
wsl --help | Select-String -Pattern "manage|move|sparse|resize"
```

本机版本不支持 `--manage` 时，用下面的方法 B 或 C。

</div>

#### 方法 B：export / import（任何版本都能用，也是官方推荐的备份方式）

```powershell
# 管理员 PowerShell
wsl --shutdown
New-Item -ItemType Directory -Force D:\WSL\backup | Out-Null

wsl --export     AlmaLinux-9 D:\WSL\backup\AlmaLinux-9.tar
wsl --unregister AlmaLinux-9                        # ⚠ 这一步会永久删除原发行版
wsl --import     AlmaLinux-9 D:\WSL\AlmaLinux-9 D:\WSL\backup\AlmaLinux-9.tar --version 2
```

官方语法（[Basic commands for WSL](https://learn.microsoft.com/windows/wsl/basic-commands)）：

```text
wsl --export <Distribution Name> <FileName>                              # 默认 tar；加 --vhd 导出 .vhdx（仅 WSL2）
wsl --import <Distribution Name> <InstallLocation> <FileName> [--vhd] [--version <1/2>]
```

<div class="zy-box note-box" markdown="1">

**⚠ 最大的坑：`--import` 之后默认用户变成 `root`，而且用常规命令改不回来。**

官方文档明确写着：`<发行版名> config --default-user <用户名>` 这条命令
**对"导入的发行版"无效**——因为导入的发行版没有可执行启动器（executable launcher）。
正确做法是在发行版**内部**写 `/etc/wsl.conf`：

```ini
# 在 WSL 里执行: sudo nano /etc/wsl.conf
[user]
default = facxing          # 换成你自己的用户名
```

然后 `wsl --terminate AlmaLinux-9` 重启生效。用 `whoami` 确认；
如果输出是 `root`，说明还没生效。

</div>

#### 方法 C：直接搬 vhdx + `--import-in-place`（最快，不需要额外临时空间）

```powershell
# 管理员 PowerShell
wsl --shutdown

# 1) 找到原来的 ext4.vhdx（路径随 WSL 版本不同，两个位置都看看）
dir "$env:LOCALAPPDATA\wsl\AlmaLinux-9"
dir "$env:LOCALAPPDATA\Packages" -Recurse -Filter ext4.vhdx -ErrorAction SilentlyContinue

# 2) 把 vhdx 挪到目标盘
New-Item -ItemType Directory -Force D:\WSL\AlmaLinux-9 | Out-Null
move "$env:LOCALAPPDATA\wsl\AlmaLinux-9\ext4.vhdx" D:\WSL\AlmaLinux-9\ext4.vhdx

# 3) 注销旧注册，再在原地导入
wsl --unregister AlmaLinux-9
wsl --import-in-place AlmaLinux-9 D:\WSL\AlmaLinux-9\ext4.vhdx
```

`--import-in-place` **要求这个 vhdx 是 ext4 文件系统**（WSL2 的默认就是）。
它只是"移动 + 重新注册"，不压缩也不解包，几十 GB 的发行版也是秒级完成，
而且**不需要一份与发行版等大的额外临时空间**——方法 B 则需要，
因为 `--export` 会生成一个完整大小的 tar。

#### 迁移后自查清单

```powershell
wsl -l -v                            # STATE 应为 Stopped，VERSION 应为 2
wsl -d AlmaLinux-9 -- uname -a
```

```bash
# 进入发行版后逐项确认
whoami                               # ❗不能是 root
echo $HOME                           # 应为你自己的家目录
df -h /                              # 根分区容量是否正常
which gfortran python3               # 编译链还在不在
ls ~/MG5atNLO                        # 自己的数据还在不在
cd ~/MG5atNLO/MG5_aMC_v2_9_18 && /usr/bin/python3 ./bin/mg5_aMC --version
```

<div class="yd-box note-box" markdown="1">

**vhdx 只会长大，不会自己缩回去。** 在 WSL 里删掉 20 GB 文件后，
Windows 侧的 vhdx 往往还是原来的大小。要真正把空间还回去：

```powershell
wsl --manage AlmaLinux-9 --set-sparse true     # 开启稀疏 vhdx（较新的 WSL）
```

```bash
# 或在 WSL 内先把空闲块 trim 掉，再 wsl --shutdown 让 Windows 侧压缩
sudo fstrim -av
```

若 `--set-sparse` 不被支持，可用 Hyper-V 的 `Optimize-VHD`（需管理员 + Hyper-V 模块）
离线压缩 vhdx。

</div>

<div class="ts-box note-box" markdown="1">

**更省事的做法**：如果**还没开始装**，就一开始装到别的盘——
见 §1.2 的方案 (a)（用 `wsl --import` 直接指定 `<InstallLocation>`）。
另外 Microsoft Store 的发行版可以在「设置 → 系统 → 存储 → 新的内容保存位置」里
改默认安装盘，但**非 Store 方式安装的发行版（包括 AlmaLinux 的 tar 镜像）不受它影响**。

</div>

### 1.7 用 VS Code 连接 WSL

<div class="dy-box note-box" markdown="1">

1. **在 Windows 上安装 VS Code**：<https://code.visualstudio.com/>
2. 安装扩展 **WSL**（`ms-vscode-remote.remote-wsl`，发布者是 Microsoft）。
   在扩展面板搜索 `WSL` 即可；装好后左下角会出现绿色的 `><` 远程指示器。
3. **从 WSL 终端启动**（最省事）：

   ```bash
   cd ~/CMS-Note-html/MadGraph-tutorial
   code .            # 自动开一个连到 WSL 的新窗口
   ```
4. **从 VS Code 里连接**：`Ctrl+Shift+P` → 输入 `WSL: Connect to WSL`
   （或 `WSL: Connect to WSL using Distro...` 选择指定发行版）。
5. 推荐再装几个扩展（**要装在 WSL 侧**，VS Code 会提示 `Install in WSL`）：
   * `ms-python.python`（Python 调试）
   * `ms-vscode.cpptools`（C/Fortran 代码跳转）
   * `redhat.vscode-yaml`、`yzhang.markdown-all-in-one`

</div>

**验证**：在 VS Code 里按 `` Ctrl+` `` 打开集成终端，运行 `uname -a`，
应当看到 Linux 内核版本而不是 Windows。

### 1.8 有用的链接

| 主题 | 链接 |
|------|------|
| WSL 官方安装指南 | <https://learn.microsoft.com/windows/wsl/install> |
| WSL 基本命令参考 | <https://learn.microsoft.com/windows/wsl/reference> |
| WSL 命令详解（含 `--export/--import`）| <https://learn.microsoft.com/windows/wsl/basic-commands> |
| `.wslconfig` 配置 | <https://learn.microsoft.com/windows/wsl/wsl-config> |
| WSL 与 VS Code | <https://code.visualstudio.com/docs/remote/wsl> |
| AlmaLinux 官方 WSL 说明 | <https://wiki.almalinux.org/documentation/wsl.html> |
| AlmaLinux 包搜索 | <https://pkgs.org/search/?q=&on_distribution=AlmaLinux+9> |
| CERN lxplus / 容器 | <https://linux-training.web.cern.ch/services/containers/> |

### 1.9 WSL 常见坑

<div class="zy-box note-box" markdown="1">

* **`wsl --install` 报 "虚拟化未启用"**：去 BIOS 打开 VT-x/AMD-V，并确认
  "Windows 功能" 里 `虚拟机平台` / `适用于 Linux 的 Windows 子系统` 已勾选。
* **`sudo` 提示命令找不到**：AlmaLinux 的 WSL 镜像默认用户可能还在 `wheel` 组外，
  `sudo usermod -aG wheel $USER` 后重新登录。
* **`dnf` 报 "未找到匹配的参数"**：先 `sudo dnf -y install epel-release` 并
  `sudo dnf config-manager --set-enabled crb`（RHEL 系把很多开发包放在 CRB 仓库）。
* **系统时间错误导致 TLS 报错**：WSL2 时间漂移时执行 `sudo hwclock -s`，
  或直接 `wsl --shutdown` 重启。
* **别在 `/mnt/c` 下编译**，见 §1.4。
* **C 盘被 vhdx 吃满 / 想把发行版挪到 D 盘**：见 §1.6（三种迁移方法 + 导入后默认用户变 root 的坑）。
* **迁移后登录成了 root 且 `config --default-user` 无效**：改发行版内的 `/etc/wsl.conf`，见 §1.6。

</div>

---

## 2. 安装 MadGraph5_aMC@NLO

### 2.1 版本选择：2.x（LTS）还是 3.x？

MadGraph 官方同时维护两条线，**两个都推荐**，用途略有不同：

| 版本线 | 代表版本 | 下载地址 | 特点 | 适合谁 |
|--------|----------|----------|------|--------|
| **2.x（LTS）** | 2.9.18 | <https://launchpad.net/mg5amcnlo/2.0/2.9.x> | **长期支持**，极其稳定；CMS 的 gridpack 生产、`HiggsDNA`/`Combine` 流程大量基于 2.9.x；Bug 最少 | **生产、复现已有分析、跟 CMS 官方样本对齐** |
| **3.x** | 3.5.x | <https://launchpad.net/mg5amcnlo/3.0/3.5.x> | 新特性：Python3 全面化、更快的编译、新的 NLO/合并功能、更好的 EFT 支持 | **新分析、新模型、想用最新功能** |

```bash
# 安装位置建议：家目录下的 ~/MG5atNLO/（绝对不要在 /mnt/c 下）
mkdir -p ~/MG5atNLO && cd ~/MG5atNLO

# ---- 方案 A：2.9.x LTS（本教程实测版本）----
wget https://launchpad.net/mg5amcnlo/2.0/2.9.x/+download/MG5_aMC_v2.9.18.tar.gz
tar -xzf MG5_aMC_v2.9.18.tar.gz
cd MG5_aMC_v2_9_18

# ---- 方案 B：3.x ----
# wget https://launchpad.net/mg5amcnlo/3.0/3.5.x/+download/MG5_aMC_v3.5.4.tar.gz
# tar -xzf MG5_aMC_v3.5.4.tar.gz && cd MG5_aMC_v3_5_4

# ---- 方案 C：从 git 拉最新（能拿到未发布修复）----
# git clone https://github.com/mg5amcnlo/mg5amcnlo.git -b 3.x
```

<div class="zs-box note-box" markdown="1">

**下载不到 tar.gz 时**：Launchpad 的下载按钮有时需要 JS，直接访问
<https://launchpad.net/mg5amcnlo> 页面右侧的 "Download" 区域。
国内网络慢时可以用 `-c` 续传，或者使用 GitHub 镜像。

</div>

### 2.2 依赖清单

| 类别 | 组件 | 是否必需 | 说明 |
|------|------|----------|------|
| 编译器 | `gfortran`, `gcc`, `make` | **必需** | 编译生成的 Fortran/C 矩阵元 |
| 运行时 | `python3` + `six`, `numpy` | **必需** | MG5 是 Python 程序；`six` 用于 py2/py3 兼容 |
| 文件工具 | `tar`, `gzip`, `wget`, `sed`, `grep` | **必需** | 内部脚本大量使用 |
| 图形 | `ghostscript`(`gs`) | 推荐 | 生成 Feynman 图的 `.ps/.pdf`（`-nojpeg` 可跳过 jpeg） |
| PDF | LHAPDF6 | 推荐 | 强子初态必需；e⁺e⁻ 过程可不用 |
| 簇射 | Pythia 8 | 推荐 | 部分子簇射 + 强子化 |
| 分析 | MadAnalysis 5 | 可选 | 快速画运动学分布（**需要 ROOT**） |
| 探测器 | Delphes | 可选 | 快速探测器模拟 |
| 编译器加速 | `fastjet-config` | 可选 | 否则 MG5 回退到自带的 `fjcore` |

### 2.3 启动 MG5 并安装可选组件

MG5 是个**交互式 shell**。第一次启动：

```bash
cd ~/MG5atNLO/MG5_aMC_v2_9_18
/usr/bin/python3 ./bin/mg5_aMC        # 注意用系统 python3，见 §1.3
```

会看到：

```text
************************************************************
*                                                          *
*                     W E L C O M E to                     *
*              M A D G R A P H 5 _ a M C @ N L O           *
*                                                          *
*         VERSION 2.9.18                2023-12-08         *
*                                                          *
************************************************************
load MG5 configuration from input/mg5_configuration.txt
fastjet-config does not seem to correspond to a valid fastjet-config executable (v3+). We will use fjcore instead.
lhapdf-config does not seem to correspond to a valid lhapdf-config executable.
set lhapdf to .../HEPTools/lhapdf6_py3/bin/lhapdf-config
No valid eps viewer found. Please set in ./input/mg5_configuration.txt
MG5_aMC>
```

<div class="zs-box note-box" markdown="1">

**关于这几条 warning：**

* `fastjet ... will use fjcore instead` — 只要不依赖外部 FastJet 就**可以忽略**；
  想消除就 `sudo dnf -y install fastjet-devel` 再 `set fastjet /usr/bin/fastjet-config`。
* `No valid eps viewer / web browser found` — 只影响自动打开 `index.html` 与
  Feynman 图预览，**不影响计算**。
* 每次启动还会问 **"Do you want to check for a new version? [y/N]"**，
  只想批处理时先把 `input/mg5_configuration.txt` 里的
  `auto_update = False` 打开。

</div>

**安装可选组件**（在 `MG5_aMC>` 提示符下）：

```text
MG5_aMC> install pythia8
MG5_aMC> install lhapdf6
MG5_aMC> install MadAnalysis5
MG5_aMC> install Delphes
MG5_aMC> install fastjet
MG5_aMC> display                      # 查看当前配置（下载/安装路径、已装工具）
MG5_aMC> quit
```

* Pythia8 会装到 `HEPTools/pythia8/`（本机实测为 **pythia8 3.1.1**）。
* LHAPDF6 会装到 `HEPTools/lhapdf6_py3/`，MG5 会自动把
  `lhapdf = .../lhapdf6_py3/bin/lhapdf-config` 写进 `input/mg5_configuration.txt`。
* 安装完可以用 `set` 检查：

```text
MG5_aMC> set pythia8_path
MG5_aMC> set lhapdf
```

<div class="ts-box note-box" markdown="1">

**离线/内网机器**：`install` 会去 `https://cern.ch/mg5amcnlo/tools/...` 下载源码包。
如果拿不到外网，可以手工下载 tarball 后放到 `HEPTools/` 下，或者
`install pythia8 --no_install` 之后自己把源码放进去。

</div>

### 2.4 快速自检

```bash
cd ~/MG5atNLO/MG5_aMC_v2_9_18
printf 'generate e+ e- > mu+ mu-\noutput /tmp/test_mumu -nojpeg\nquit\n' \
  | /usr/bin/python3 ./bin/mg5_aMC
ls /tmp/test_mumu/bin/        # 应看到 generate_events, madevent
```

看到 `1 processes with 1 diagrams generated` 与 `output` 成功，说明工具链完全就绪。

### 2.5 常见错误速查

<div class="zy-box note-box" markdown="1">

| 现象 | 原因 | 解决 |
|------|------|------|
| `requires the six module` | 用了缺 `six` 的 Python（如 miniforge 3.13） | 用 `/usr/bin/python3` 启动，`python3 -m pip install --user six` |
| `gfortran: command not found` | 没装编译链 | `sudo dnf -y install gcc-gfortran` |
| `install pythia8` 编译报错 | Pythia8 版本与 MG5 期望的接口不匹配（`Analysis`/`Pythia8` 版本漂移） | 手工改 `HEPTools/pythia8/` 里的版本号并补 `Analysis.cc`/`Cuts.cc`（`Pythia8::Analysis::~Analysis` 之类）|
| `lhapdf-config ... not a valid` | LHAPDF 未装或未设置 | `install lhapdf6` 或 `set lhapdf /path/to/lhapdf-config` |
| `error: no such option: --shower` | `launch` 的参数写法不对（见 §3.4） | 用 `launch madevent <dir>/ -f` |
| `FileNotFoundError: .../Events` | 手工删过 `Events/` 目录 | `mkdir -p <proc>/Events` |
| MA5 编译到 "interface to Root" 失败 | 系统 ROOT 缺 `libEG`（`/usr/bin/ld: cannot find -lEG`） | 装全 ROOT（`sudo dnf -y install root root-graf3d-eve`）；确认 `root-config --libs` 里有 `-lEG` |
| 编译极慢 | 过程目录放在 `/mnt/c` | 挪到 Linux 文件系统，见 §1.4 |

</div>

---

## 3. MadGraph5 语法速查

MG5 的输入分**两层**：交互式 shell 的**命令**，以及被命令读写的**卡片文件**。

### 3.1 顶层命令

| 命令 | 作用 | 例子 |
|------|------|------|
| `import model <name>` | 载入模型（UFO） | `import model sm` / `import model ./SMEFTsim_UFO` |
| `define <name> = <particles>` | 定义多粒子标签 | `define ell+ = e+ mu+ ta+` |
| `generate <proc>` | 生成过程（**覆盖**已有过程） | `generate e+ e- > e+ e-` |
| `add process <proc>` | 追加过程（可与 `@n` 标记组合） | `add process p p > t t~` |
| `output <dir> [opts]` | 导出过程目录 | `output bhabha_ee -nojpeg` |
| `launch <dir> [opts]` | 运行完整链条 | `launch madevent bhabha_ee/ -f` |
| `set <key> <value>` | 修改 **MG5 自身** 的配置 | `set nb_core 8`、`set run_mode 2` |
| `install <tool>` | 安装可选组件 | `install pythia8` |
| `display <what>` | 显示当前状态 | `display processes` / `display diagrams` |
| `update <what>` | 更新缓存/模型 | `update lo` / `update pythia8` |
| `tutorial` / `help` | 内置教程 | `tutorial`、`help launch` |
| `quit` / `exit` | 退出 | |

<div class="zs-box note-box" markdown="1">

**注意区分两种 `set`：**

* 在 `MG5_aMC>` 顶层：`set` 改的是 **MG5 配置**（path、run_mode、nb_core…），
  只有一小组合法键。写错会报
  `Possible options for set are ['group_subprocesses', 'stdout_level', ...]`。
* 在 `launch` 的**卡片编辑子提示符**里：`set <key> <value>` 改的是 **run_card / param_card**
  里的物理参数，例如 `set nevents 10000`、`set ebeam1 5.29`。

</div>

### 3.2 `generate` 的过程语法

```text
generate  <初态>  >  <末态>  [, <末态> , ...]   [耦合阶数限制]   [相空间调度]

# 几个典型写法
generate e+ e- > e+ e-                       # 2 -> 2 Bhabha
generate e+ e- > e+ e- mu+ mu-               # 2 -> 4
generate p p > t t~, (t > w+ b, w+ > l+ vl)  # 带衰变链（在 MadSpin 或 decay 里更规范）
generate p p > w+ j, w+ > l+ vl QCD=0        # 限制 QCD 阶数 = 0（即只算 EW 贡献）
generate p p > h [QCD]                       # 方括号表示 NLO QCD
generate p p > t t~ QED<=2 QCD<=2            # 限制耦合阶数上限
generate p p > t t~ / a z                    # 排除中间光子/Z（"/" = 禁止这些粒子作为中间态）
generate p p > e+ e- $ z                     # "$" = 只保留 z 为共振（on-shell）
define j = g u u~ d d~ c c~ s s~ b b~        # 自定义"喷注"定义
generate p p > w+ , (w+ > j j)               # 用自定义 j
```

常用记号：

| 记号 | 含义 |
|------|------|
| `p` / `e-` / `mu+` / `j` | 多粒子标签（`p`=质子，`j`=喷注，`l+`=带电轻子） |
| `>` | 初态与末态分隔 |
| `,` | 并列不同末态 / 不同衰变链（**不是**求和） |
| `@n` | 给过程打标签，配合 `add process ... @n` 把多个过程合到一个 `P<n>_*` 目录 |
| `/` | 禁止某粒子出现在**中间态**（interference 处理） |
| `$` | 只允许该粒子作为共振态（on-shell，窄宽度近似） |
| `QCD=n`, `QED=n` | 指定/限制耦合阶数 |
| `[QCD]` | 启用 NLO QCD（MadLoop + FKS） |
| `{...}` | 指定粒子极化（配合模型中的 `@` 标签） |
| `\` | 与 `/` 类似但保留共振（高级用法） |

### 3.3 卡片（Cards）

`output` 之后，过程目录 `<proc>/Cards/` 下会有一组卡片。**所有物理量都在这里改**：

| 卡片 | 控制什么 | 关键条目 |
|------|----------|----------|
| `proc_card_mg5.dat` | 过程定义（由 `generate` 生成，可手工再编辑） | `generate ...` |
| `param_card.dat` | 模型参数（质量、宽度、耦合、CKM…） | `Block mass`、`Block SMINPUTS`、`Block yukawa` |
| `run_card.dat` | **束流能量、过程级切割、事例数、标度** | `ebeam1/2`、`nevents`、`ptl`、`etal`、`drll`、`mmll`、`ickkw`、`xqcut` |
| `pythia8_card.dat` | 是否开 ISR/FSR、强子化、MPI、喷注匹配 | `PartonShower`、`TimeShower`、`SpaceShower`、`JetMatching` |
| `madspin_card.dat` | 共振态的自旋关联衰变 | `set spinmode`、`decay` |
| `madanalysis5_parton_card.dat` | MG5 自动调 MA5 时的分析脚本（部分子级） | MA5 命令 |
| `madanalysis5_hadron_card.dat` | 同上（强子级） | MA5 命令 |
| `me5_configuration.txt` | 过程目录的本地 MG5 配置（path、run_mode…） | `run_mode`、`nb_core`、`mg5_path` |

**run_card 的语义：**

* 数值单位一律 **GeV**；**`-1` 表示"没有这个上限"**，`0` 表示"没有这个下限"。
* 主要切割：

  | 键 | 含义 |
  |----|------|
  | `ebeam1`, `ebeam2` | 束流能量（正负电子对撞） |
  | `nevents` | 非加权事例数（官方建议**单次不超过 1M**） |
  | `ptl`, `etal`, `el` | 带电轻子的最小 $$p_T$$、最大赝快度、最小能量 |
  | `drll`, `mmll` | 轻子间最小 $$\Delta R$$、最小不变质量 |
  | `ptj`, `etaj`, `drjj` | 喷注的类似切割 |
  | `pt_min_pdg` / `eta_max_pdg` | **按 PDG 号**给任意粒子加切割，如 `{6: 30.0}` |
  | `ickkw`, `xqcut` | 是否做 MLM 匹配、矩阵元 $$k_\perp$$ 下限 |
  | `dynamical_scale_choice` | 动态标度（如 `3` = $$\sum$$ 末态 $$p_T$$） |
  | `fixed_ren_scale`, `fixed_fac_scale`, `scale` | 固定 $$\mu_R,\mu_F$$ |
  | `iseed` | 随机数种子（`0` = 用时间） |

**卡片语法：**

```text
# 井号开头的整行是注释
  500.0     = ebeam1   ! beam 1 total energy in GeV
  2.5       = etal     ! max rap for the charged leptons
  -1.0      = ptlmax   ! max pt for the charged leptons
```

形式统一为 `值  = 键  ! 注释`。**注意不要改键名或删掉整行**，
否则 MG5 重新读卡时会报 `IndexError` 或直接用默认值。

### 3.4 批处理：不要交互式敲命令

在 `launch` 的卡片子提示符里，语法是 **`键=值`**（**不带 `--`**），
最后输入 **`0`**（或 `done`）确认；**Ctrl-D（EOF）会走默认值**。
实践中最省事的两种写法：

**(a) 把命令喂给 stdin（推荐，本教程全程使用）**

```bash
# 从顶层接口直接跑已经生成好的过程目录，-f 强制使用现有卡片、不再提问
printf 'launch madevent bhabha_ee/ -f\n' \
  | /usr/bin/python3 ~/MG5atNLO/MG5_aMC_v2_9_18/bin/mg5_aMC
```

**(b) 一次生成 + 运行（写成一个 `.mg5` 脚本）**

```text
# file: bhabha.mg5
import model sm
generate e+ e- > e+ e-
output bhabha_ee -nojpeg
```

```bash
/usr/bin/python3 ~/MG5atNLO/MG5_aMC_v2_9_18/bin/mg5_aMC bhabha.mg5
```

<div class="zy-box note-box" markdown="1">

**⚠ `launch` 的参数坑**（实测）：
`launch <dir>` 里的 `<dir>` 必须是**已在 MG5 里生成过的过程目录**，
并且顶层语法其实是 `launch <mode> <dir> <options>`，
其中 `<mode>` ∈ `madevent | pythia8 | standalone | aMC@NLO | madweight`。
所以：

* ✅ `launch madevent bhabha_ee/ -f`
* ❌ `launch bhabha_ee --shower=OFF` → `error: no such option: --shower`
* ❌ `launch bhabha_ee` → `bhabha_ee cannot be run from MG5 interface`

`-f`（force）表示"用目录里现成的卡片，别再问我"；`-s parton` 可以只跑到部分子级。
另外 `<proc>/Events/` 目录必须存在，手工删过就 `mkdir -p` 补回来。

</div>

---

## 4. 实例：Bhabha 散射 $$e^+e^-\to e^+e^-$$

### 4.1 为什么选这个例子

Bhabha 散射是 QED 里**最有教学价值**的过程：

* 树图只有**两张**（$$s$$ 道湮灭 + $$t$$ 道交换），但已经包含
  **$$t$$ 道的前向发散**、**$$s$$–$$t$$ 干涉**、以及**相同末态粒子的不同费米子连线**；
* 微分截面在教科书里有**闭式解析表达式**，可以直接和 Monte Carlo 逐 bin 对照——
  这是检验"MC 到底算对了没有"的黄金标准；
* 它在 $$e^+e^-$$ 对撞机上是**亮度测量**的标准过程（BESIII、Belle II 都用它）；
* 取 $$\sqrt{s}=10.58\ \mathrm{GeV}$$ 正好是 BESIII 的 $$\psi(3770)$$ 能量点，
  和本笔记里 D2HEE 的分析环境一致。

<div class="dy-box note-box" markdown="1">

### 理论参考

本文的 Feynman 规则、Mandelstam 变量与微分截面推导遵循

**余钊焕，《量子场论与标准模型导论》**
<https://yzhxxzxy.github.io/teaching/2202_QFT_SM_intro.pdf>

（对应 Peskin & Schroeder 第 5.1 节 $$e^+e^-\to\mu^+\mu^-$$ 的推导思路，
把末态换回 $$e^+e^-$$ 即可。）

</div>

### 4.2 Feynman 图与 Feynman 规则（颜色编码）

<div style="text-align:center;">
  <img src="/assets/images/notes/MadGraphTut/bhabha_feynman_diagrams.png"
       alt="Bhabha scattering Feynman diagrams with colour-coded Feynman rules"
       style="max-width:94%;">
</div>

**约定**：蓝色（<span class="text-blue">#2563EB</span>）表示**外线费米子**，
红色（<span class="text-red">#DC2626</span>）表示**传播子**，
绿色（<span class="text-green">#16A34A</span>）表示**相互作用顶点**，
紫色（<span class="text-purple">#9333EA</span>）表示**四动量标签**。
费米子线上的箭头是**费米子数流**方向：正粒子与动量同向，反粒子反向。

<div class="dy-box note-box" markdown="1">

### 颜色 → Feynman 规则对照表

| 颜色 | 图中元素 | 对应 Feynman 规则 |
|------|----------|-------------------|
| <span class="text-blue">**蓝色**</span> | 外线・入射 $$e^-$$ | <span class="text-blue">$$u(p)$$</span> |
| <span class="text-blue">**蓝色**</span> | 外线・入射 $$e^+$$ | <span class="text-blue">$$\bar v(p)$$</span> |
| <span class="text-blue">**蓝色**</span> | 外线・出射 $$e^-$$ | <span class="text-blue">$$\bar u(p)$$</span> |
| <span class="text-blue">**蓝色**</span> | 外线・出射 $$e^+$$ | <span class="text-blue">$$v(p)$$</span> |
| <span class="text-red">**红色**</span> | 光子传播子 | <span class="text-red">$$\dfrac{-ig_{\mu\nu}}{q^2}$$</span> |
| <span class="text-green">**绿色**</span> | QED 顶点 | <span class="text-green">$$-ie\gamma^{\mu}$$</span> |
| <span class="text-purple">**紫色**</span> | 四动量 | <span class="text-purple">$$p_i$$</span> |

</div>

**读图规则**：沿费米子连线走，起点写旋量、终点写伴随旋量，顶点写 $$-ie\gamma^\mu$$，
跨过传播子写 $$-ig_{\mu\nu}/q^2$$，再对所有洛伦兹指标求和。

**$$(a)$$ $$s$$ 道（湮灭图）**，$$q=p_1+p_2$$，$$q^2=s$$。

$$
i\mathcal{M}_s=\;
\textcolor{#2563EB}{\bar v(p_2)}\,
\textcolor{#16A34A}{(-ie\gamma^{\mu})}\,
\textcolor{#2563EB}{u(p_1)}
\;\;
\textcolor{#DC2626}{\frac{-ig_{\mu\nu}}{s}}
\;\;
\textcolor{#2563EB}{\bar u(p_3)}\,
\textcolor{#16A34A}{(-ie\gamma^{\nu})}\,
\textcolor{#2563EB}{v(p_4)}
$$

**$$(b)$$ $$t$$ 道（交换图）**：$$e^-$$ 线 $$p_1\to p_3$$ 与 $$e^+$$ 线 $$p_2\to p_4$$ 之间交换光子，
$$q=p_1-p_3$$，$$q^2=t$$。

$$
i\mathcal{M}_t=\;
\textcolor{#2563EB}{\bar u(p_3)}\,
\textcolor{#16A34A}{(-ie\gamma^{\mu})}\,
\textcolor{#2563EB}{u(p_1)}
\;\;
\textcolor{#DC2626}{\frac{-ig_{\mu\nu}}{t}}
\;\;
\textcolor{#2563EB}{\bar v(p_2)}\,
\textcolor{#16A34A}{(-ie\gamma^{\nu})}\,
\textcolor{#2563EB}{v(p_4)}
$$

**Mandelstam 变量**（无质量极限）：

$$
s=(p_1+p_2)^2,\qquad
t=(p_1-p_3)^2=-\frac{s}{2}(1-\cos\theta),\qquad
u=(p_1-p_4)^2=-\frac{s}{2}(1+\cos\theta),
\qquad s+t+u=0
$$

其中 $$\theta$$ 是**出射 $$e^-$$ 与入射 $$e^-$$** 的夹角。两图相加、对自旋求和平均即得经典结果
$$\langle\lvert\mathcal{M}\rvert^2\rangle=2e^4\left[(s^2+u^2)/t^2+(u^2+t^2)/s^2+2u^2/(st)\right]$$，
三项依次对应 $$t$$ 道、$$s$$ 道与 $$s$$–$$t$$ 干涉（$$t\to0$$ 时为破坏性干涉）。

### 4.3 解析微分截面与切割

代入无质量两体末态的 $$\mathrm{d}\sigma/\mathrm{d}\Omega=\lvert\mathcal{M}\rvert^2/(64\pi^2 s)$$ 得

$$
\frac{\mathrm{d}\sigma}{\mathrm{d}\Omega}
=\frac{\alpha^2}{2s}\left[
\frac{1+\cos^4(\theta/2)}{\sin^4(\theta/2)}
-\frac{2\cos^4(\theta/2)}{\sin^2(\theta/2)}
+\frac{1+\cos^2\theta}{2}
\right]
$$

第一项（$$t$$ 道）$$\propto1/\sin^4(\theta/2)$$ 在 $$\theta\to0$$ **发散**（光子无质量时是真发散），
所以**任何 Bhabha 截面数值都必须写明切割条件**。在无质量极限下末态两轻子背靠背，
$$p_T=\dfrac{\sqrt{s}}{2}\sin\theta$$，于是

$$
p_T(\ell)>p_{T,\min}
\quad\Longleftrightarrow\quad
\lvert\cos\theta\rvert<c_{\max}=\sqrt{1-\left(\frac{2p_{T,\min}}{\sqrt{s}}\right)^2}
$$

**数值**（$$\alpha=1/132.507$$，与 MadGraph5 的 `sm` 模型 `aEWM1` 一致）：

| $$\sqrt{s}$$ [GeV] | $$p_{T,\min}$$ [GeV] | $$c_{\max}$$ | $$\sigma_{\text{QED}}^{\text{LO}}$$ [pb] | $$\sigma s$$ [pb·GeV²] |
|------------------|--------------------|-----------|------------------------|----------------------|
| 10.58 | 1.058 | 0.979796 | $$2.26787\times10^{5}$$ | $$2.5386\times10^{7}$$ |
| 3.00  | 0.300 | 0.979796 | $$2.82063\times10^{6}$$ | $$2.5386\times10^{7}$$ |
| 1.00  | 0.100 | 0.979796 | $$2.53857\times10^{7}$$ | $$2.5386\times10^{7}$$ |

注意 $$p_{T,\min}=0.1\sqrt{s}$$ 的取法让三个能量点的**角度接受度完全相同**
（$$c_{\max}=0.979796$$ 处处一样），于是 QED 预言必须严格满足 $$\sigma\propto1/s$$——
最后一列是常数 $$2.5386\times10^{7}$$，正好验证了这一点。

<div class="zy-box note-box" markdown="1">

**⚠ 用同一套 $$\alpha$$。** MadGraph5 的 `sm` 默认 `aEWM1 = 132.507`，而 PDG 的
$$\alpha(0)=1/137.035999$$；两者差 3.4%，在截面上就是 $$(\Delta\alpha/\alpha)^2\approx6.9\%$$。
把"约定差别"当成"MC 算错"是最容易踩的坑。本教程理论数字统一用 $$\alpha=1/132.507$$。

</div>

### 4.4 在 MadGraph 里跑 Bhabha

#### 步骤 1：建立过程目录

```text
# file: bhabha.mg5
import model sm
generate e+ e- > e+ e-
output bhabha_ee -nojpeg
```

```bash
/usr/bin/python3 ~/MG5atNLO/MG5_aMC_v2_9_18/bin/mg5_aMC bhabha.mg5
```

输出（真实日志）：

```text
INFO: Trying process: e+ e- > e+ e- WEIGHTED<=4 @1
INFO: Process has 4 diagrams
1 processes with 4 diagrams generated in 0.002 s
Total: 1 processes with 4 diagrams
output bhabha_ee -nojpeg
```

<div class="zs-box note-box" markdown="1">

**为什么是 4 张图？** MadGraph 用的是**包含 $$\gamma$$ 与 $$Z$$ 的完整标准模型**，
每张 QED 树图都有 $$\gamma$$ 和 $$Z$$ 两个版本：$$2\times2=4$$。
所以 MadGraph 的结果**不是纯 QED**，而是 $$\gamma+Z+$$干涉。
在本节的能量下这个差别极小（见 §4.5 的讨论），但你要知道它在那里。

想看这 4 张图，在 `MG5_aMC>` 里用
`display diagrams bhabha_ee`，或者打开 `<proc>/HTML/card_1.html`。

</div>

#### 步骤 2：改 `run_card.dat`

```text
# <proc>/Cards/run_card.dat   —— 只列出被改动的行

  400000     = nevents      ! Number of unweighted events requested
     5.29     = ebeam1       ! beam 1 total energy in GeV   (sqrt(s)=10.58 GeV)
     5.29     = ebeam2       ! beam 2 total energy in GeV
 1.058        = ptl          ! minimum pt for the charged leptons
 10.0         = etal         ! max rap for the charged leptons
 0.0          = drll         ! min distance between leptons
 bhabha_s10p58 = run_tag     ! name of the run
```

* `ebeam1/2 = 5.29 GeV` → $$\sqrt{s}=10.58\ \mathrm{GeV}$$；
* `ptl = 1.058 GeV` = $$0.1\sqrt{s}$$，对应 $$c_{\max}=0.979796$$；
* `etal = 10.0` 必须**放宽**，否则默认的 `2.5` 会在 $$p_T$$ 之前先切掉前向区，
  破坏与解析式的对应关系；
* `nevents` 400k 用于画分布；只求截面的试跑用 2000 就够。

#### 步骤 3：运行

```bash
printf 'launch madevent bhabha_ee/ -f\n' \
  | /usr/bin/python3 ~/MG5atNLO/MG5_aMC_v2_9_18/bin/mg5_aMC
```

真实日志（2000 事例的试跑）：

```text
generate_events run_01 -f
INFO: load configuration from .../bhabha_ee/Cards/me5_configuration.txt
INFO: Running Survey
Creating Jobs
Working on SubProcesses
INFO: Compiling for process 1/1.
INFO:     P1_ll_ll
survey  run_01
INFO: finish refine
INFO: Combining Events
combination of events done in 0.21535706520080566 s
  === Results Summary for run: run_01 tag: tag_1 ===

     Cross-section :   2.266e+05 +- 627.4 pb
     Nb of events :  2000

store_events
INFO: Storing parton level results
INFO: End Parton
```

产物：

| 文件 | 内容 |
|------|------|
| `<proc>/Events/run_01/unweighted_events.lhe.gz` | **非加权事例**（部分子级，LHE 格式） |
| `<proc>/Events/run_01/run_01_tag_1_banner.txt` | 本次运行的全部卡片快照（**复现用**） |
| `<proc>/Events/run_01/run_01_tag_1_results.dat` | 截面、误差、各迭代收敛信息 |
| `<proc>/index.html` | 汇总页面（含截面图） |

<div class="ts-box note-box" markdown="1">

**写一个可复现的驱动脚本。** 手工在交互式 shell 里点来点去不可复现。
本教程的做法是把"改卡片 + 运行 + 收集 LHE"写成一个 Python 脚本
（见 §6 的 `run_bhabha_scan.py`），要点是：

1. **原子地改卡片**：写 `run_card.dat.tmp` 再 `os.replace()`，
   避免中途被打断留下半张卡；
2. 用 `subprocess.run(..., input="launch madevent <dir>/ -f\n")` 非交互运行；
3. 从 stdout 里正则抓 `Cross-section : <σ> +- <err> pb`；
4. 把 `Events/run_XX/unweighted_events.lhe.gz` 拷到你自己的数据目录并**重命名**。

另外把过程目录的 `Cards/me5_configuration.txt` 里加上

```text
run_mode = 0      # 单核模式：不开工作进程池，避免产生大量 ajob* 临时文件
nb_core  = 1
```

单核模式在调试阶段更稳、日志更整齐；正式跑大样本时再改成
`run_mode = 2` 配合 `nb_core` 多核并行。

</div>

### 4.5 结果对照：MadGraph vs 解析 QED

#### 4.5.1 总截面

| 点 | $$\sqrt{s}$$ [GeV] | $$p_{T,\min}$$ [GeV] | $$N_{\text{evt}}$$ | $$\sigma_{\text{MG5}}$$ [pb] | $$\sigma_{\text{QED}}^{\text{LO}}$$ [pb] | MG / QED |
|----|------|------|--------|--------------------------|--------------------------|----------|
| s10p58 | 10.58 | 1.058 | 400 000 | $$226\,839.7\pm33.9$$ | 226 786.8 | $$1.00023$$ |
| s3     | 3.00  | 0.300 |  40 000 | $$2\,823\,126\pm2\,040$$  | 2 820 633  | $$1.00088$$ |
| s1     | 1.00  | 0.100 |   8 000 | $$(2.53936\pm0.00044)\times10^{7}$$ | $$2.53857\times10^{7}$$ | $$1.00031$$ |

**结论：在 0.02% – 0.09% 的水平上一致**，而 MadGraph 自己的积分误差就有
0.015% – 0.175%（最后一个点的统计最差），所以**残差完全落在统计/数值误差之内**。
这正是我们想要的"闭合检验"（closure test）。

#### 4.5.2 微分分布

<div style="text-align:center;">
  <img src="/assets/images/notes/MadGraphTut/bhabha_s10p58_compare.png"
       alt="MadGraph vs analytic LO QED differential cross section for Bhabha scattering"
       style="max-width:88%;">
</div>

* **上图**：$$\mathrm{d}\sigma/\mathrm{d}\cos\theta_{e^-}$$ 对照。MC（蓝色阶梯）与解析式（红线）
  在**四个数量级**的范围内完全重合。注意分布**极度前向**：$$t$$ 道光交换使
  $$e^-$$ 几乎沿原方向飞出，$$\cos\theta_{e^-}\to+1$$ 处 $$\mathrm{d}\sigma/\mathrm{d}\cos\theta$$
  比 $$90^\circ$$ 处高三个数量级。
* **下图**：MC / 理论比值，平坦在 1 附近，误差棒就是 MC 统计误差。

<div class="zs-box note-box" markdown="1">

**⚠ 角度定义的坑（我踩过）**

LHE 里 **MadGraph 把 $$e^+$$ 放在 beam 1（$$+z$$）**，而不是 $$e^-$$。
所以直接用 `pz/|p|` 当作 $$\theta_{e^-}$$ 会得到**镜像**的分布
（前向峰跑到后向去）。正确做法是**从 LHE 的初态粒子里找出 `PDG=11, ISTUP=-1`
的那一个**，用它的方向定义 $$\theta$$：

```python
sign = np.sign(pz_of_incoming_electron)     # LHE 里 e- 沿 -z，故 sign = -1
cos_theta = sign * pz_out / |p_out|
```

**只有在对称窗口 $$\lvert c\rvert<c_{\max}$$ 上积分时，镜像才会给出相同的总截面**——
所以"总截面一致"并不能证明角度定义对了，**必须画微分分布**。

</div>

#### 4.5.3 残差从哪里来

<div class="dy-box note-box" markdown="1">

已知的几个物理/数值来源，按重要程度：

1. **MC 数值积分误差**（VEGAS）——量级 $$10^{-4}$$，是当前的主导项。
2. **$$\gamma$$–$$Z$$ 干涉**。MadGraph 的 `sm` 是完整 SM，而解析式是纯 QED。
   在这个能量下 $$s\ll M_Z^2$$，$$Z$$ 的贡献被压得很低。做个**闭合检验**就一目了然：

   ```text
   MG5_aMC> generate e+ e- > e+ e- / z     # "/ z" = 禁止 Z 作为中间态
   MG5_aMC> output bhabha_ee_qedonly -nojpeg
   ```

   （`/ z` 让 Diagrams 从 4 张降到 **2 张**，即纯 $$\gamma$$ 交换。）
   实测在 $$\sqrt{s}=10.58\ \mathrm{GeV}$$：纯 QED 版本
   $$\sigma=227\,000\pm160\ \mathrm{pb}$$，与解析 QED 的 226 787 pb 相差
   $$+0.09\%\pm0.07\%$$；而完整 SM 是 $$226\,839.7\pm33.9\ \mathrm{pb}$$，
   即 $$Z$$ 只把总截面拉低了约 0.07%。

   **为什么这么小？** 因为带 $$p_T$$ 切割的总截面被**前向 $$t$$ 道峰**主导，
   $$s$$ 道（也就是 $$Z$$ 起作用的地方）只占约
   $$805\ \mathrm{pb}/226\,787\ \mathrm{pb}\approx0.35\%$$。所以 $$Z$$ 对**总截面**的影响
   天然被压低一个数量级。

   **反过来说**：如果你把能量抬到 $$\sqrt{s}\sim M_Z$$，这个"QED 一致"立刻崩掉——
   那时 $$Z$$ 主导，逐 bin 对照纯 QED 会在 $$Z$$ 峰附近差出成倍。
3. **电子质量**。解析式取 $$m_e\to0$$；$$m_e^2/s$$ 在 $$\sqrt{s}=1\ \mathrm{GeV}$$ 时约 $$2.6\times10^{-7}$$，
   完全可忽略。
4. **切割边界效应**。$$c_{\max}$$ 附近的相空间被 VEGAS 采样的精度、以及
   "无质量近似 $$\leftrightarrow$$ 真实 $$p_T$$"的微小差别影响。
5. **$$\alpha$$ 的约定**，见 §4.3 的警告框——用错会差 6.9%。

</div>

**复现方式**：把 `data/compare_summary.txt` 里的数字与
`figures/bhabha_s3_compare.pdf`（$$\sqrt{s}=3$$ GeV）、
`figures/bhabha_s1_compare.pdf`（$$\sqrt{s}=1$$ GeV）一起看，
可以看到 $$\sigma\propto1/s$$ 的标度关系在三个点上都被重现。

### 4.6 用 MadAnalysis 5 做快速分析

MadAnalysis 5（MA5）是"给 LHE/HepMC 文件做快分析"的框架：
用类似自然语言的命令定义观测量，自动出图并给出效率报告。
目前安装的是 **MA5 1.11.0**。

**两条使用路径：**

<div class="dy-box note-box" markdown="1">

#### (a) 由 MG5 自动调用

MG5 在 `launch` 时如果检测到 MA5，会在 Pythia8 之后自动跑它，
读的脚本是过程目录里的：

* `Cards/madanalysis5_parton_card.dat`（部分子级）
* `Cards/madanalysis5_hadron_card.dat`（强子级）

在这个文件里直接写 MA5 命令即可，例如把 Bhabha 的十几行脚本放进去，
每次 `launch` 就自动得到分析结果。不想要这一步就在过程目录的
`Cards/me5_configuration.txt` 里写 `madanalysis5_path =`（留空）。

#### (b) 独立运行（推荐，灵活）

```bash
MA5=~/MG5atNLO/MG5_aMC_v2_9_18/HEPTools/madanalysis5/madanalysis5

# 写一个分析脚本
cat > bhabha.ma5 <<'EOF'
# ---- 观测量定义 ----
define e = e+ e-
plot PT(e)      40  0   6   [logY]     # 每个轻子的 pT
plot ABSETA(e)  40  0   3              # |eta|
plot M(e+ e-)   40  0   12             # 末态轻子对不变质量
plot THETA(e-)  40  0   3.2            # 极角

# ---- 读入事例 ----
import /path/to/bhabha_s10p58.lhe.gz as Bhabha_ee

# ---- 提交 ----
submit
EOF

# -P = parton level, -s = 跑完脚本后直接退出（批处理）
/usr/bin/python3 $MA5/bin/ma5 -P -s bhabha.ma5
```

结果落在 MA5 的 `Output/` 下，包含 PDF/PNG 图、`.safe` 归档与效率报告。
`ma5 --help` 可以看到全部选项（`-H` 强子级、`-R` 探测器级、`-f` 强制覆盖）。

</div>

<div class="zy-box note-box" markdown="1">

**⚠ MA5 依赖 ROOT。**

MA5 的 `SampleAnalyzer` 必须链接 ROOT 才能编译出可执行文件。
在本机实测时遇到：

```text
g++ -shared -o ../Lib/libroot_for_ma5.so ... -lCore -lImt ... -lEG ...
/usr/bin/ld: cannot find -lEG
collect2: error: ld returned 1 exit status
MA5-ERROR: impossible to link the project.
```

即系统 ROOT（`dnf` 安装的 6.32.00）**缺 `libEG`**。排查与解决：

```bash
root-config --version           # 确认 ROOT 存在
root-config --libs | tr ' ' '\n' | grep -E '^-l(EG|Gui)$'   # 是否输出 -lEG
ldconfig -p | grep libEG        # 运行时是否有 libEG.so
sudo dnf -y install root root-graf3d-eve root-gui            # 补全 ROOT 组件
```

如果机器上根本没有 ROOT，另外两个选择：

* 在 MG5 里 `install MadAnalysis5` 时把 `root_path` 指到你自己编译的 ROOT；
* **或者绕过 MA5**：部分子级的 LHE 文件结构非常简单，
  用 20 行 numpy 就能读出四动量并画图（本教程 §6 的 `analyze_bhabha.py`
  就是这么做的，同时还能画出与解析式的逐 bin 比值）——做**理论对照**时这样更直接。

</div>

**MA5 里几个好用的命令**（交互模式下 `help` 可看全部）：

| 命令 | 作用 |
|------|------|
| `import <file> [as <label>]` | 读入 LHE / HepMC / ROOT 事例文件 |
| `define <name> = <particles>` | 定义多粒子（如 `define e = e+ e-`） |
| `plot <OBS>(<particles>) <nbin> <min> <max> [logY]` | 定义一张图 |
| `select <cut>` | 加选择，如 `select PT(e) > 20` |
| `set main.fastsim.package = fastjet` | 用 FastJet 做喷注聚类（否则用 MA5 自带实现） |
| `set main.fastsim.algorithm = antikt` / `.ptmin` / `.radius` | 喷注算法与参数 |
| `submit` | 提交分析，出图与效率表 |
| `open` / `display` | 查看已有输出 |

---

## 5. 把 MadGraph 接进 CMS 流程

<div class="dy-box note-box" markdown="1">

### 5.1 Gridpack：把 MG5 过程交给 CMS 批量生产

LO 过程可以打成 **gridpack**（一个自解压 tar 包，内含积分网格与可执行文件），
然后在 CMS 的 **Private Production** 里用 `gridpack_generator_cfg` 配 CMSSW 作业：

```text
MG5_aMC> launch bhabha_ee
# 在 run_card 里设 gridpack = .true.，或者直接
MG5_aMC> generate_events -f
```

* 对应的 MG5 选项：run_card 里 `gridpack`，或 `survey`/`refine`/`combine_events`
  之后 `create_gridpack`；
* 生产脚本用 `cmsDriver.py` 时把 `--generator` 指向该 gridpack；
* **关键**：gridpack 里的 `run_card` 决定了束流、切割与标度，**务必固定 `iseed`**
  并记录 `banner.txt`，否则无法复现。

### 5.2 与 Pythia8 的衔接

* `launch <proc>` 时把 shower 选成 `Pythia8`（MG5 会自动调用 `HEPTools/pythia8`）；
* 调参在 `Cards/pythia8_card.dat`，例如关掉 MPI 便于和理论比：

  ```text
  PartonLevel:MPI = off
  SpaceShower:QEDshowerByL = off
  ```

* 多喷注 LO 过程需要 **MLM 匹配**：在 run_card 里设
  `ickkw = 1` 并给 `xqcut`（一般 $$\approx m$$ 或 $$\approx p_T^{\min,\text{jet}}$$），
  再在 `pythia8_card.dat` 里让 `JetMatching:doMerge = on`；
  **务必检查 `xqcut` 附近 $$\mathrm{d}\sigma/\mathrm{d}p_T$$ 的连续性**，
  不连续就说明匹配参数没调好。

### 5.3 与 Delphes 的衔接

* `launch <proc>` 时把 detector 选成 `Delphes`，卡片是 `Cards/delphes_card_CMS.dat`；
* 若做带电轻子分析，**先在 `Delphes_card_CMS.dat` 里确认轻子效率/隔离的来源**，
  不要默认它是"理想的"；
* 需要与 Run 3 真实数据对齐时，建议改用 **NanoAOD** 流程而不是 Delphes。

</div>

---

## 6. 本教程的可复现文件

本目录（`CMS-Note-html/MadGraph-tutorial/`）下附带全部脚本与数据：

| 文件 | 说明 |
|------|------|
| `theory_bhabha.py` | 解析 LO-QED Bhabha 计算：$$\mathrm{d}\sigma/\mathrm{d}\cos\theta$$、带 $$p_T/\eta$$ 切割的积分截面，含自检 |
| `run_bhabha_scan.py` | 驱动 MG5 在多个能量点跑 LO 事例；原子方式改写 `run_card.dat`；解析日志抓截面 |
| `analyze_bhabha.py` | 读 LHE、重建 $$\cos\theta_{e^-}$$、与解析式逐 bin 对比、出图 |
| `feynman_bhabha.tex` | 颜色编码的 Feynman 图（TikZ） |
| `build_figures.sh` | 编译 `.tex` 并用 ghostscript 裁出紧凑的 PDF/PNG |
| `work/bhabha_ee/` | MadGraph 生成的过程目录（`Cards/`、`SubProcesses/`、`bin/`） |
| `data/*.lhe.gz` | 三个能量点的部分子级事例（LHE + gzip） |
| `data/compare_summary.txt` | 汇总表 |
| `figures/*.pdf` `*.png` | 所有图 |

**一键复现：**

```bash
cd CMS-Note-html/MadGraph-tutorial

# 1) 解析截面自检（不需要 MG5）
/usr/bin/python3 theory_bhabha.py

# 2) 跑 MadGraph 扫描（需要已生成 work/bhabha_ee/，约几分钟）
/usr/bin/python3 run_bhabha_scan.py

# 3) 与解析式对比、出图
/usr/bin/python3 analyze_bhabha.py

# 4) 重建 Feynman 图
bash build_figures.sh
```

---

## 7. 参考链接

### 7.1 MadGraph 生态

| 资源 | 链接 |
|------|------|
| MG5 主页 / 下载 | <https://launchpad.net/mg5amcnlo> |
| MG5 GitHub（3.x 分支） | <https://github.com/mg5amcnlo/mg5amcnlo> |
| MG5 官方问答（搜报错必看） | <https://answers.launchpad.net/mg5amcnlo> |
| aMC@NLO 论文 (JHEP 1407 (2014) 079) | <https://arxiv.org/abs/1405.0301> |
| MG5 论文 (JHEP 1106 (2011) 128) | <https://arxiv.org/abs/1106.0522> |
| MadAnalysis 5 主页 | <https://madanalysis.irmp.ucl.ac.be/> |
| MadAnalysis 5 GitHub | <https://github.com/MadAnalysis> |
| Pythia 8 | <https://pythia.org/> |
| Delphes | <https://cp3.irmp.ucl.ac.be/projects/delphes> |
| LHAPDF6 | <https://lhapdf.hepforge.org/> |
| UFO / FeynRules | <http://feynrules.irmp.ucl.ac.be/> |

### 7.2 理论与教程

| 资源 | 链接 |
|------|------|
| **余钊焕《量子场论与标准模型导论》** | <https://yzhxxzxy.github.io/teaching/2202_QFT_SM_intro.pdf> |
| Peskin & Schroeder, *An Introduction to QFT* | 第 5.1 节（$$e^+e^-\to\mu^+\mu^-$$）、第 6 章（Bhabha / Møller） |
| PDG Review of Particle Physics | <https://pdg.lbl.gov/> |
| Bhabha 微分截面（PDG 汇总） | <https://pdg.lbl.gov/2024/reviews/rpp2024-rev-cross-section-plots.pdf> |

### 7.3 环境与工具

| 资源 | 链接 |
|------|------|
| WSL 安装 | <https://learn.microsoft.com/windows/wsl/install> |
| WSL 命令参考 | <https://learn.microsoft.com/windows/wsl/reference> |
| VS Code + WSL | <https://code.visualstudio.com/docs/remote/wsl> |
| AlmaLinux WSL | <https://wiki.almalinux.org/documentation/wsl.html> |
| LHE 格式规范 | <https://arxiv.org/abs/hep-ph/0109068> |
| HepMC3 | <https://gitlab.cern.ch/hepmc/HepMC3> |

---

## 8. 常见问题速查（TL;DR）

<div class="zy-box note-box" markdown="1">

| 症状 | 一句话解决 |
|------|-----------|
| MG5/MA5 启动报 `requires the six module` | 用 `/usr/bin/python3` 启动，`python3 -m pip install --user six` |
| `launch bhabha_ee` 报 `cannot be run from MG5 interface` | 语法是 `launch madevent bhabha_ee/ -f` |
| `launch ... --shower=OFF` 报 `no such option` | `launch` 的参数是 `-f/-s/-n`，shower 之类在卡片里改 |
| `FileNotFoundError: .../Events` | `mkdir -p <proc>/Events` |
| 编译慢到怀疑人生 | 项目挪出 `/mnt/c`，放到 `~/` 下 |
| 总截面和解析式差 7% | 检查 $$\alpha$$ 的约定（`aEWM1` vs $$\alpha(0)$$） |
| 微分分布左右镜像了 | LHE 里 beam 1 是 $$e^+$$；用初态 `PDG=11, ISTUP=-1` 的方向定角度 |
| 截面比预期小一个量级 | `etal` / `ptl` 默认切割没放宽，先把 `etal` 放到 10 |
| 强子化后喷注分布和理论对不上 | 检查 MLM 匹配的 `xqcut` 与 `JetMatching:doMerge` |
| MA5 编译卡在 Root 接口 | ROOT 缺 `libEG`，补装 `root-graf3d-eve` 等组件 |

</div>
