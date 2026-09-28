---
title: "CMS Sim&Rec Note"
collection: notes
type: notes
permalink: /notes/cms-sim-reco
date: 2026-03-04
venue: "CERN CMS"
excerpt: "CMS 模拟与重建笔记：gridpack 产生、Pythia8 强子化、Private Production、本底分析及 Preselection"
author_profile: true
---

<p><em>Facxing &nbsp;|&nbsp; 2026年3月4日</em></p>

## Lxplus 服务器目录

### 概况

参考 [https://linux-training.web.cern.ch/services/containers/](https://linux-training.web.cern.ch/services/containers/)，CERN 的 lxplus 主要有两个目录:

<div class="dy-box note-box" markdown="1">

### AFS 文件系统

afs: 通常目录为 `/afs/cern.ch/user/<username>`

参考 [KB0002863](https://cern.service-now.com/service-portal?id=kb_article&n=KB0002863) afs的基本信息如下:

| Space type | Available space (def./max.) | Suggested I/O activity | Power area provided |
|------------|---------------------------|----------------------|-------------------|
| home dir   | $$2 \mathrm{~GB} / 10 \mathrm{~GB}$$ | low | critical |
| work space | $$20 \mathrm{~GB} / 100 \mathrm{~GB}$$ | medium | non-critical |

<div class="zy-box note-box" markdown="1">
**⚠ 注意：** 相对来说AFS的存储空间较小，不建议在其上运行脚本
</div>

</div>

<div class="dy-box note-box" markdown="1">

### EOS/CERNBox 文件系统

eos: 通常目录为 `/eos/home-<first-letter-of-username>/<username>` 或 `/eos/user/<first-letter-of-username>/<username>`

<div class="zs-box note-box" markdown="1">
**📝 例如：** 用户名为 <span class="text-teal">"miku"</span>, 目录为 `/eos/user/m/miku`
</div>

参考 [CERN KB](https://cern.service-now.com/service-portal?sys_kb_id=fae8543fc9ed05006d218776d679b74a&id=kb_article_view&sysparm_rank=1&sysparm_tsqueryId=6d11a0258332325070ceba65eeaad303)

想要查看eos储存空间配额，首先输入 `export EOS_MGM_URL=root://eosuser.cern.ch`，然后键入指令 `eos quota`（加可选 `-h` 获取更多信息）

<div class="zy-box note-box" markdown="1">
**⚠ 注意：** 注意输出中"logi bytes"和可用内存"aval bytes"的区别，"logi bytes"包括必要的备份和冗余！
</div>

</div>

---

## 产生gridpack

### 概述

- CMS的Gridpack文档为: [https://cms-generators.docs.cern.ch/how-to-produce-gridpacks/](https://cms-generators.docs.cern.ch/how-to-produce-gridpacks/)
- 另外也有一个新手教程: [https://codimd.web.cern.ch/s/glLKn0Nb-#2-Setting-up-a-CMSSW-Environment](https://codimd.web.cern.ch/s/glLKn0Nb-#2-Setting-up-a-CMSSW-Environment)
- 编辑并测试成功的卡片(process card):
  - <span class="text-red">对于 $$H\rightarrow\gamma\gamma$$ 过程，工作流为 MadGraph $$\rightarrow$$ Pythia8</span>
  - <span class="text-blue">对于其他过程 $$ggVH, H\rightarrow\gamma\gamma$$; $$H\rightarrow b\bar{b}, H\rightarrow VV$$ 等，工作流一般为 Powheg/JHUGen $$\rightarrow$$ Pythia8</span>

```madgraph
import model sm-ckm_no_b_mass
define W = w+ w-
define ell+ = e+ mu+ ta+
define ell- = e- mu- ta-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
 
generate p p > w+({0}/{T}) h, w+ > ell+ vl            @1
add process p p > w+({0}/{T}) h j, w+ > ell+ vl       @2
add process p p > w+({0}/{T}) h j j, w+ > ell+ vl     @3
add process p p > w-({0}/{T}) h, w- > ell- vl~        @4
add process p p > w-({0}/{T}) h j, w- > ell- vl~      @5
add process p p > w-({0}/{T}) h j j, w- > ell- vl~    @6
```

- 此卡片基于 RUNII, RUNIII(lowmass Higgs) 的 Card 修改得到

如 RUNII 的卡片 (HIG-RunIISummer19UL17wmLHEGEN-00120):

```madgraph
in run_card.dat:
import model loop_sm-ckm_no_b_mass

define V = w+ w- z

generate p p > V h [QCD] @0
add process p p > V h j [QCD] @1
add process p p > V h j j [QCD] @2

output vh012j_5f_NLO_FXFX_M<MASS>_VToAll -nojpeg
```

同时 W, H 衰变由 MadSpin 完成:

```madgraph
set ms_dir ./madspingrid

set Nevents_for_max_weigth 250 # number of events for the estimate of the max. weight
set max_weight_ps_point 400  # number of PS to estimate the maximum for each event

set max_running_process 1

decay w+ > all all
decay w- > all all
decay z > all all
launch
```

此处 RUNII 为 Vtoall，在 RUNIII(lowmass Higgs)中 (如 HIG-RunIII2024Summer24wmLHEGS-00723等):

```madgraph
in run_card.dat:
import model loop_sm-ckm_no_b_mass

define W = w+ w-

generate p p > W h [QCD] @0
add process p p > W h j [QCD] @1
add process p p > W h j j [QCD] @2

output WH012j_5f_NLO_FXFX_M125_WToLNu -nojpeg

in madspin_card.dat:
set ms_dir ./madspingrid
set Nevents_for_max_weight 250
set max_weight_ps_point 400

set max_running_process 1
define ell+ = e+ mu+ ta+
define ell- = e- mu- ta-

decay w+ > ell+ vl
decay w- > ell- vl~
launch
```

总结: process card 中主要的修改以及问题

1. MadGraph 不支持NLO过程做极化，因此先做LO过程的gridpack $$\rightarrow$$ <span class="text-blue">问题1: 如何比较两者的截面?</span><sup id="fnref:1"><a href="#fn:1" class="footnote">1</a></sup>
   *在更早的时候(RUNI? Ecm 7-8 TeV)有LO的计算，但是直接使用Pythia6产生的...*
2. 因为 $$H\rightarrow \gamma \gamma$$ 为 Loop-level 过程，如果使用CMS推荐的5-flavour scheme即 `import model sm-ckm_no_b_mass`，则不能在Madspin中做此过程 $$\rightarrow$$ 在Pythia8 fragment中衰变
3. 因为极化过程只做到LO, 将MadSpin中的W衰变过程与Process card合并
4. <span class="text-blue">问题2: (ggWH是否需要考虑 $$\rightarrow$$ 约占8%?)</span><sup id="fnref:2"><a href="#fn:2" class="footnote">2</a></sup>

<div class="footnotes" markdown="1">
<ol>
<li id="fn:1"><a href="#fnref:1" class="reversefootnote">↩</a> <span class="text-blue">此处需要讨论理论上的Scale Factor: K</span></li>
<li id="fn:2"><a href="#fnref:2" class="reversefootnote">↩</a> <span class="text-red">这个过程实际上不满足电荷共轭 $$\rightarrow$$ 不存在此过程</span></li>
</ol>
</div>

### <span class="text-red">**试跑时遭遇的报错:**</span>

<div class="lt-box note-box" markdown="1">

#### 20260203 — run WH(polarized)

1. <span class="text-red">**"generate_events pilotrun ... KeyError: 'gridpack'"**</span>

   <span class="text-cyan">这主要是由于脚本试图在run_card中寻找变量 "gridpack"，这个变量应该由 `<selfdefined_name>_run_card.dat` 提供以确认是否开启gridpack，例如在样例卡片 `wplustest_4f_LO_run_card.dat` 中:</span>

   ```
   .false.     = gridpack  !True = setting up the grid pack
   ```

   <span class="text-cyan">似乎旧版本的cards都不附带此变量，需要在运行前检查run_card是否含有变量"gridpack".</span>

2. <span class="text-orange">**遇到一个警告:** `value '3' for entry 'ickkw' is not valid. Preserving previous value: '0'. allowed values are 0, 1. "xqcut>0 but ickkw=0. Potentially not fully consistent setup. Be careful"`</span>

   <span class="text-cyan">参考 [ggzh01j_5f_LO_MLM](https://gitlab.cern.ch/cms-gen/genproductions_cards/-/tree/master/MadGraph5_aMCatNLO/examples/ggzh01j_5f_LO_MLM)，改用其run_card中的MLM匹配/合并方案 (ickkw = 1)</span>

3. **Fortran编译错误:**
   ```
   run_card.inc:389:13:
     389 |       DRLL_SF = '  0.0 '
         |             1
   Error: Symbol 'drll_sf' at (1) has no IMPLICIT type; did you mean 'drllmax'?
   ```

   <div class="zy-box note-box" markdown="1">
   **⚠ 注意：**
   I. <span class="text-red">此报错主要可能由于NLO和LO过程的变量名不一致导致</span>
   II. <span class="text-cyan">后来参考文档中的一个例子(胶子融合产生vh)的cards做修改: [ggzh01j_5f_LO_MLM](https://gitlab.cern.ch/cms-gen/genproductions_cards/-/tree/master/MadGraph5_aMCatNLO/examples/ggzh01j_5f_LO_MLM) 可以成功运行</span>
   III. 也可以使用相同版本的MadGraph产生类似过程，并查看对应的cut变量名
   </div>

</div>

### 修改后的 Run Card 设置

- **Cuts:** 参考 RUNIII(Low higgs mass) cards，部分变量更改名字对应LO
- **Merge Scheme:** 使用 MLM merging，参考 RUN II ggZH 的配置

```madgraph
...
6800.0     = ebeam1  ! beam 1 total energy in GeV
6800.0     = ebeam2  ! beam 2 total energy in GeV
...
lhapdf                = pdlabel     ! PDF set                                  
$DEFAULT_PDF_SETS       = lhaid
$DEFAULT_PDF_MEMBERS    = reweight_PDF     ! if pdlabel=lhapdf, this is the lhapdf number
#CMS recommended LHAPDF Setting
# Matching - Warning! ickkw > 1 is still beta
1  = ickkw ! 0 no matching, 1 MLM, 2 CKKW matching
1  = ktscheme ! for ickkw=1, 1 Durham kT, 2 Pythia pTE
1.0 = alpsfact ! scale factor for QCD emission vx
5  = asrwgtflavor ! highest quark flavor for a_s reweight
True = clusinfo ! include clustering tag in output
3.0 = lhe_version ! Change the way clustering information pass to shower.      
...
15.0  = bwcutoff      ! (M+/-bwcutoff*Gamma)
10.0  = ptj       ! minimum pt for the jets
-1.0  = ptjmax    ! maximum pt for the jets
5.0 = etaj    ! max rap for the jets
0.4 = drjj    ! min distance between jets
1.0  = drjjmax ! max distance between jets
0.0   = mmjj    ! min invariant mass of a jet pair
-1.0  = mmjjmax ! max invariant mass of a jet pair
30   = mmll  ! Min inv. mass of all opp. sign same-flavor lepton pairs
5 = maxjetflavor    ! Maximum jet pdg code
10.0   = xqcut   ! minimum kt jet measure between partons
...
```

### 产生的lhe事例

一个 lhe 事例 (nevent = 7) 如下:

```
<event>
 7      2 +7.7400258e-02 3.48452100e+02 7.54677100e-03 1.07205100e-01
       -2 -1    0    0    0  501 +0.0000000000e+00 +0.0000000000e+00 +8.0158566747e+01 8.0158566747e+01 0.0000000000e+00 0.0000e+00 1.0000e+00
        1 -1    0    0  502    0 +0.0000000000e+00 +0.0000000000e+00 -4.1191826038e+02 4.1191826038e+02 0.0000000000e+00 0.0000e+00 -1.0000e+00
      -24  2    1    2    0    0 +8.6865602647e+01 +7.9170536053e+01 -2.4942187953e+02 2.8715498681e+02 8.0206644103e+01 0.0000e+00 0.0000e+00
       11  1    3    3    0    0 +4.7478491455e+01 +1.8071363476e+01 -3.7588510499e+01 6.3195549292e+01 0.0000000000e+00 0.0000e+00 -1.0000e+00
      -12  1    3    3    0    0 +3.9387111192e+01 +6.1099172576e+01 -2.1183336903e+02 2.2395943752e+02 0.0000000000e+00 0.0000e+00 1.0000e+00
       25  1    1    2    0    0 -9.6891452769e+01 -7.7000887282e+01 -8.3277932029e+01 1.9462092443e+02 1.2500000000e+02 0.0000e+00 0.0000e+00
       21  1    1    2  502  501 +1.0025850122e+01 -2.1696487707e+00 +9.4011792828e-01 1.0300915890e+01 0.0000000000e+00 0.0000e+00 -1.0000e+00
<scales pt_clust_1="13600.0" pt_clust_2="13600.0" pt_clust_3="13600.00000" pt_clust_4="10.25793"></scales>
<mgrwt>
<rscale>  0 0.34845210E+03</rscale>
<asrwt>  1 0.10257926E+02</asrwt>
<pdfrwt beam="2">  1        1 0.60576216E-01 0.34845210E+03</pdfrwt>
<pdfrwt beam="1">  1       -2 0.11788024E-01 0.34845210E+03</pdfrwt>
<totfact> 0.60378426E+02</totfact>
```

主要信息包括 (较长的重整化尺度/权重信息未展示):

<div class="dy-box note-box" markdown="1">

#### 事件头

| 字段 | 值 | 描述 |
|------|-----|------|
| NUP | 7 | 事件中的粒子数 |
| IDPRUP | 2 | 过程标识符 |
| XWGTUP | $$7.7400258\times 10^{-2}$$ | 事件权重 |
| SCALUP | $$3.48452100\times 10^{2}$$ GeV | 事件尺度（用于部分子簇射） |
| AQEDUP | $$7.54677100\times 10^{-3}$$ | QED耦合常数 |
| AQCDUP | $$1.07205100\times 10^{-1}$$ | QCD耦合常数 |

</div>

<div class="dy-box note-box" markdown="1">

#### 初末态粒子信息

<div class="table-wide" markdown="1">

| 索引 | pdgID | ISTUP | MOTHUP1 | MOTHUP2 | ICOLUP1 | ICOLUP2 | Px (GeV) | Py (GeV) | Pz (GeV) | E (GeV) | M (GeV) | SPINUP | 描述 |
|------|-------|-------|---------|---------|---------|---------|----------|----------|----------|---------|---------|--------|------|
| 1 | -2 | -1 | 0 | 0 | 0 | 501 | 0.0 | 0.0 | 8.016e+01 | 8.016e+01 | 0.0 | 1.0 | 反上夸克（初始束流1） |
| 2 | 1 | -1 | 0 | 0 | 502 | 0 | 0.0 | 0.0 | -4.119e+02 | 4.119e+02 | 0.0 | -1.0 | 下夸克（初始束流2） |
| 3 | -24 | 2 | 1 | 2 | 0 | 0 | 8.687e+01 | 7.917e+01 | -2.494e+02 | 2.872e+02 | 8.021e+01 | 0.0 | W-玻色子（中间共振态） |
| 4 | 11 | 1 | 3 | 3 | 0 | 0 | 4.748e+01 | 1.807e+01 | -3.759e+01 | 6.320e+01 | 0.0 | -1.0 | 电子（终态） |
| 5 | -12 | 1 | 3 | 3 | 0 | 0 | 3.939e+01 | 6.110e+01 | -2.118e+02 | 2.240e+02 | 0.0 | 1.0 | 反电子中微子（终态） |
| 6 | 25 | 1 | 1 | 2 | 0 | 0 | -9.689e+01 | -7.700e+01 | -8.328e+01 | 1.946e+02 | 1.250e+02 | 0.0 | 希格斯玻色子（终态） |
| 7 | 21 | 1 | 1 | 2 | 502 | 501 | 1.003e+01 | -2.170e+00 | 9.401e-01 | 1.030e+01 | 0.0 | -1.0 | 胶子（终态） |

</div>

<div class="cx-box note-box" markdown="1">

#### 一些问题

- 新的卡片中部分子分布设置应该改为:

```
lhapdf                = pdlabel     ! PDF set                                  
$DEFAULT_PDF_SETS       = lhaid
$DEFAULT_PDF_MEMBERS    = reweight_PDF     ! if pdlabel=lhapdf, this is the lhapdf number
```

此设置将部分子分布的设置改为: [MetaData](https://gitlab.cern.ch/cms-gen/genproductions_scripts/-/blob/master/MetaData)

详细情况参见 [CMS generators docs](https://cms-generators.docs.cern.ch/how-to-produce-gridpacks/mg5-amcnlo/#preparation-of-the-cards)

</div>

</div>

---

## pythia8 fragment 运行和强子化

此处主要参考两个部分的教程:

- [https://sihyunjeon.github.io/generators-cmsdaslpc2024/02-particle_level/index.html](https://sihyunjeon.github.io/generators-cmsdaslpc2024/02-particle_level/index.html)
- [https://codimd.web.cern.ch/s/glLKn0Nb-#2-Setting-up-a-CMSSW-Environment](https://codimd.web.cern.ch/s/glLKn0Nb-#2-Setting-up-a-CMSSW-Environment)

<div class="zy-box note-box" markdown="1">

**⚠ 注意：** 这里不再赘述如何配置CMSSW，可参考上述第二个链接的 Section 2

<span class="text-red">一般使用教程 gridpack.sh 产生的 gridpack 包名字会带有使用的 CMSSW 版本，如:</span>

<span class="text-blue">WH012j_WtoLNu_HtoGG_PolP</span><span class="text-orange">\_5f\_LO\_MLM</span><span class="text-purple">\_el8\_amd64\_gcc10\_CMSSW\_12\_4\_8\_tarball</span>

- <span class="text-blue">**定义整个物理过程:**</span>
  1. WH012j: 过程为 $$p p > W H(n j)$$ 即 $$p p > W H$$ + $$p p > W Hj$$ + $$p p > W Hjj$$
  2. WtoLNu: W衰变为轻子+中微子
  3. HtoGG: Higgs粒子衰变到两个光子
  4. PolP: 极化态 (**Pol**arized)

- <span class="text-orange">**定义Generation的设置**</span>
  1. 5f: 在矩阵元计算中考虑5种夸克味 (u, d, c, s, b)
  2. LO: Leading Order: 领头阶微扰计算
  3. MLM matching: 一种矩阵元与部分子簇射的匹配方案（由M. L. Mangano提出）参数 ickkw=1 时启用 *arXiv:1109.5295v1 [hep-ph]*

- <span class="text-purple">**计算的环境要求**</span>
  1. el8: 操作系统为 Enterprise Linux 8
  2. amd64: x86-64架构
  3. gcc10: 使用GNU编译器版本
  4. CMSSW_12_4_8: CMSSW版本

</div>

### 运行时产生的报错和解决方法

在配置好 <span class="text-red">**对应版本**</span> 参考教程使用 cmsDriver.py，主要设置为:

```bash
#!/bin/bash
cmsDriver.py WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM \
  --python_filename run_WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM.py \
  --eventcontent NANOAOD \
  --datatier NANOAOD \
  --fileout file:WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM.root \
  --conditions auto:mc \
  --step LHE,GEN,NANOGEN \
  --no_exec \
  --mc \
  -n 100
```

运行后出现报错:

<div class="lt-box note-box" markdown="1">

```
LHE,GEN,NANOGEN,ENDJOB
Step: LHE Spec:
Loading lhe fragment from WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM
Traceback (most recent call last):
...
ModuleNotFoundError: No module named 'WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM'
```

</div>

<div class="jjfa-box note-box" markdown="1">

#### 解决方案

问题出在CMSSW中运行此脚本编译找不到此文件/模块，于此需要:

```bash
#!/bin/bash
echo $CMSSW_BASE  # 查看有无输出，如果没有输出说明CMSSW软件没有配置好
# 进入 CMSSW src 目录
cd $CMSSW_BASE/src
# 创建目录（如果不存在）
mkdir -p Configuration/Generator/python
# 复制fragment
cp /full/path/to/your/xxx.py Configuration/Generator/python/
# 验证文件存在(可选)
ls -l Configuration/Generator/python/xxx.py
# 配置CMSSW，重新编译
scram b
# 重新加载环境
eval `scram runtime -sh`
```

这样在CMSSW的资源库中加入了相关的Fragment以供寻找并将其编译为Python可读的模块

相应的运行命令也需要调整:

```bash
#!/bin/bash
cmsDriver.py Configuration/GenProduction/python/WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM \
  --python_filename run_WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM.py \
  --eventcontent NANOAOD \
  --datatier NANOAOD \
  --fileout file:WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM.root \
  --conditions auto:mc \
  --step LHE,GEN,NANOGEN \
  --no_exec \
  --mc \
  -n 100
```

最后会产生一个新的可供运行的文件，输入:

```bash
cmsRun run_WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM.py
```

在运行后经由Pythia8强子化和digitization等步骤产生ROOT文件

<div class="zy-box note-box" markdown="1">
**⚠ 注意：** 详细参数解释参见第二个链接中 Section 7
</div>

</div>

### 修改后 Pythia fragment

修改 RUNIII(Lowmass higgs card) 并参考部分 LO 的 fragment，主要改动:

1. 设置Higgs质量为125.20 (参考PDG给出)
2. 设为MLM merging 并禁用 ShowerKT

```python
generator = cms.EDFilter("Pythia8ConcurrentHadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13600.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        pythia8PSweightsSettingsBlock,
        processParameters = cms.vstring(
            'JetMatching:setMad = off',
            'JetMatching:merge = on',
            'JetMatching:scheme = 1',          # 1 = MLM (in Pythia8 context for LO)
            'JetMatching:jetAlgorithm = 2',    # 2 = kT (required for MLM)
            'JetMatching:etaJetMax = 5.0',     
            'JetMatching:coneRadius = 1.0',
            'JetMatching:slowJetPower = 1',
            'JetMatching:qCut = 30.0',         # this is the actual merging scale
            'JetMatching:qCutME = 10.',        # this must match the ptj cut in the lhe generation step
            'JetMatching:nQmatch = 5',         # 5-flavor scheme
            'JetMatching:nJetMax = 2',         # number of partons in born matrix element for highest multiplicity
            'JetMatching:doShowerKt = off',    # off for MLM matching, turn on for shower-kT matching
            'TimeShower:mMaxGamma = 4.0',

            'SLHA:useDecayTable = off',
            '25:m0 = 125.20',                  # according to PDG 2024/2025 update
            '25:onMode = off',
            '25:onIfMatch = 22 22',
        ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'pythia8PSweightsSettings',
            'processParameters'
        )
    )
)
```

### 结果汇总

- Pythia输出截面只规定了Higgs如何衰变，输出的截面需要再乘以一个 $$\mathcal{B}(H\rightarrow\gamma\gamma)$$，这里取CMS测量结果 ($$m_{H} = 125.20~\textrm{GeV/}c^2$$) $$\rightarrow$$ $$\mathcal{B} = 2.270\times10^{-3}$$
- $$\mathcal{B}(W^+\rightarrow l^+\nu+c.c.) = (10.71+10.63+11.38)\%$$ <span class="text-cyan">*Phys. Rev. D **110**, no.3, 030001 (2024)*</span>
- 同CMS官方产生以及文章 <span class="text-cyan">*arXiv:2506.13002 [hep-ph]*</span> 比较

| Process | Production cross section [fb] | Matching eff(%) | in paper[fb] | CMS official production(NLO)[fb] |
|---------|------------------------------|-----------------|--------------|----------------------------------|
| $$p p \rightarrow W_{\mathrm{L}}^{ \pm} H \rightarrow \ell^{ \pm} \nu \gamma \gamma$$ | 0.575 | 39.3 | 0.588 | - |
| $$p p \rightarrow W_{\mathrm{T}}^{ \pm} H \rightarrow \ell^{ \pm} \nu \gamma \gamma$$ | 0.450 | 36.0 | 0.526 | - |
| $$p p \rightarrow W^{\pm} H \rightarrow \ell^{ \pm} \nu \gamma \gamma$$ | 1.007 | 36.4 | - | $$\sim$$1.330 |

---

## Private Production 与 CMS connect

参考肖杰师兄的提交脚本 [https://gitlab.cern.ch/jixiao/private_prod/](https://gitlab.cern.ch/jixiao/private_prod/-/tree/master?ref_type=heads)

### 遭遇的问题:

<div class="lt-box note-box" markdown="1">

```
WARNING: In non-interactive mode release checks e.g. deprecated releases, production architectures are disabled.
SCRAM fatal: Unable to locate the top of local release. Please run this command from a SCRAM-based area.
Traceback (most recent call last):
  ...
ModuleNotFoundError: No module named 'Configuration.GenProduction'
```

</div>

<div class="jjfa-box note-box" markdown="1">

#### 脚本中路径出现问题

脚本中:

```bash
...
  scram p CMSSW CMSSW_14_0_21
fi
cd CMSSW_14_0_21/src
eval `scram runtime -sh`
cd - # !!!此处出错!!!
```

将 `cd -` 语句换为绝对路径解决此问题

</div>

<div class="lt-box note-box" markdown="1">

cmsRun 显示: `tar: xxx.tar.xz cannot find`

</div>

<div class="zy-box note-box" markdown="1">

**⚠ 注意：** cmsRun似乎不会在当前目录下工作，因此任何文件需要使用绝对路径，如:

```python
process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/home/kaixin.fan/Condor_submit/test/job_WH_PolP_2024/WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM_el8_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz'),  # 此处使用绝对路径
    generateConcurrently = cms.untracked.bool(False),
    nEvents = cms.untracked.uint32(10),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)
```

</div>

---

## MC sample

### Irreducible Background

1. Signal Process: $$p p \rightarrow W^-H(ZH)\rightarrow l^-\bar{\nu}(l^+l^-)\gamma\gamma+c.c.$$

<figure>
  <img src="/assets/images/notes/FeynmanDiag/WH_associated.png" alt="WH associated production Feynman diagram" style="max-width:60%;">
</figure>

2. triboson: (major bkg)
   - $$p p \rightarrow W(Z)\gamma\gamma \rightarrow l\nu(l\bar{l})\gamma\gamma$$ (<span class="text-cyan">**JHEP 10 (2021), 174**</span>)

<figure>
  <img src="/assets/images/notes/FeynmanDiag/diboson/WGG.png" alt="Wgg Feynman diagram" style="max-width:80%;">
</figure>

3. Top association
   - $$pp\rightarrow t\bar{t}\gamma\gamma\rightarrow 2(l\nu\gamma)$$ <span class="text-blue">(b tag)</span>
   - $$pp\rightarrow t\gamma\gamma\rightarrow l\nu\gamma\gamma$$ (<span class="text-red">FCNC process with small BF, omit</span>)

### Reducible Background — lepton-$$\gamma$$ misid

leptons marked <span class="text-red">red</span> would be misidentified as a $$\gamma$$

1. Drell-yan process: DYGto2LG-1Jets $$pp\rightarrow l\textcolor{red}{l} \gamma$$
2. Single Top: TWG $$pp\rightarrow tWG \rightarrow l\textcolor{red}{l}\gamma$$ <span class="text-blue">(b tag)</span>

### Reducible Background — jet-$$l/\gamma$$ misid

jet marked <span class="text-purple">purple</span> would be misidentified as a $$l$$ or $$\gamma$$

1. GG-Box-3Jets: $$pp\rightarrow \gamma\gamma jjj$$
2. Drell-yan process: DYGto2LG-1Jets $$pp\rightarrow ll \gamma j$$
3. Top pair:
   - TTG-1Jets: $$pp\rightarrow t\bar{t}\gamma j\rightarrow ll\gamma j$$ <span class="text-blue">(b tag)</span>
   - TGQB: $$pp\rightarrow t\gamma qb\rightarrow l\gamma q b$$ <span class="text-blue">(b tag)</span>
4. diboson:
   - WG-1Jet $$pp\rightarrow W\gamma j\rightarrow l \nu \gamma j$$
   - ZG-1Jet $$pp\rightarrow Z\gamma j\rightarrow ll \gamma j$$

<span class="text-blue">Omit diboson process like $$WW$$ or $$ZZ$$</span>

### Summary of background

1. Reducible: Marked <span class="text-blue">DarkBlue</span>
2. Irreducible:
   - $$l-\gamma$$ misid: Marked <span class="text-purple">DarkOrchid</span>
   - $$j-l/j-\gamma$$ misid: Marked <span class="text-green">SeaGreen</span>

<span class="text-blue">WGG-4FS NLO MG+P8</span>

**Referring Works from CADI:**

1. HIG-19-015 — SM H → γγ signal strengths & STXS (full run 2) (Classification)
2. HIG-23-014 — H→γγ cross sections at 13.6 TeV (table Format)
3. HIG-25-020 — H→γγ STXS analysis (Run3: 2022+2023+2024) (Samples)
4. HIG-25-013 — Search for VH→μμ using Run3 data (Samples)
5. HIG-25-016 — Low mass di-photon Higgs search with Run3 (2022+2023) data (Samples & table format)

#### RUNIII2024Summer 背景样本总表

<div class="table-wide" markdown="1">

| Process | Slice (GeV) | Sample | Accuracy | $$\sigma$$ [pb] | links |
|---------|-------------|--------|----------|-----------------|-------|
| DYto2L-2J | $$m_{\ell\ell}>50,~p_T^{\ell\ell} > 100$$ | DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 | NLO | $$106.6\pm0.2843$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8) |
| $$\gamma\gamma$$+1b+2j | $$m_{\gamma\gamma}>80$$ | GG-Box-2B-2Jets_Bin-MGG-80_TuneSherpaDef_13p6TeV_sherpaMEPS | - | - | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GG-Box-2B-2Jets_Bin-MGG-80_TuneSherpaDef_13p6TeV_sherpaMEPS) |
| $$\gamma\gamma$$+2b+2j | $$m_{\gamma\gamma}>80$$ | GG-Box-1B-2Jets_Bin-MGG-80_TuneSherpaDef_13p6TeV_sherpaMEPS | - | - | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GG-Box-1B-2Jets_Bin-MGG-80_TuneSherpaDef_13p6TeV_sherpaMEPS) |
| $$\gamma\gamma$$+3j | $$m_{\gamma\gamma}>80$$ | GG-Box-3Jets_Bin-MGG-80_TuneSherpaDef_13p6TeV_sherpaMEPS | NLO | $$87.51\pm0.3614$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GG-Box-3Jets_Bin-MGG-80_TuneSherpaDef_13p6TeV_sherpaMEPS) |
| $$\gamma$$+1j | $$m_{\gamma\gamma}\in(40\sim80),~p_T <20$$ | GJet_Bin-MGG-40to80-PT-20_TuneCP5_13p6TeV_pythia8 | LO | $$3332\pm10.4$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GJet_Bin-MGG-40to80-PT-20_TuneCP5_13p6TeV_pythia8) |
| $$\gamma$$+1j | $$m_{\gamma\gamma}>80,~p_T\in(20,40)$$ | GJet_Bin-MGG-80-PT-20to40_TuneCP5_13p6TeV_pythia8 | LO | $$248.1\pm0.783$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GJet_Bin-MGG-80-PT-20to40_TuneCP5_13p6TeV_pythia8) |
| $$\gamma$$+1j | $$m_{\gamma\gamma}>80,~p_T>40$$ | GJet_Bin-MGG-80-PT-40_TuneCP5_13p6TeV_pythia8 | LO | $$923.4\pm2.854$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=GJet_Bin-MGG-80-PT40_TuneCP5_13p6TeV_pythia8) |
| $$WW$$ | - | WW_TuneCP5_13p6TeV_pythia8 | - | $$80.23\pm0.3733$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=WW_TuneCP5_13p6TeV_pythia8) |
| $$WZ$$ | - | WZ_TuneCP5_13p6TeV_pythia8 | - | $$29.1\pm0.1318$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=WZ_TuneCP5_13p6TeV_pythia8) |
| $$ZZ$$ | - | ZZ_TuneCP5_13p6TeV_pythia8 | - | $$12.75\pm0.0649$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=ZZ_TuneCP5_13p6TeV_pythia8) |
| $$W\gamma\rightarrow l\nu\gamma$$ | - | WGtoLNuG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 | NLO | $$671.5\pm0.7548$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=WGtoLNuG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8) |
| $$W\gamma\rightarrow 2q\gamma+1j$$ | - | WGto2QG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 | NLO | $$291.2\pm0.3546$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=WGto2QG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8) |
| $$Z\gamma\rightarrow2\nu\gamma+1j$$ | - | ZGto2NuG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 | NLO | $$40.52\pm0.05003$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=ZGto2NuG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8) |
| $$Z\gamma\rightarrow2q\gamma+1j$$ | - | ZGto2QG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 | NLO | $$141.5\pm0.1749$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=ZGto2QG-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8) |
| $$tt\rightarrow2l2\nu+2j$$ | - | TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8 | NLO | - | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8) |
| $$tt\rightarrow2l2\nu+3j$$ | - | TTto2L2Nu-3Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 | LO | - | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=TTto2L2Nu-3Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8) |
| $$tt\gamma+1j$$ | - | TTG-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8 | NLO | - | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=TTG-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8) |
| $$tt\gamma\gamma+1j$$ | - | TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8 | LO | - | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8) |
| $$ttW+1j$$ | - | TTW-WtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8 | NLO | $$0.4678\pm0.0005553$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=TTW-WtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8) |
| $$ttZ+1j$$ | - | TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8 | NLO | $$0.6426\pm0.005682$$ | [MCM](https://cms-pdmv-prod.web.cern.ch/mcm/requests?dataset_name=TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8) |

</div>

---

## Preselection

### Preselection requirements

Follow the lowmass Higgs preselection cuts AN2024_116.

#### Electrons

- **$$p_{\mathrm{T}}$$ threshold:** $$\geq 15 \ \mathrm{GeV}$$ (modified from 10 GeV)
- **$$\eta$$:** $$|\eta|<2.5$$
- Gap barrel (endcap) **$$\eta$$:** 1.4442 (1.566)
- ID working point: "loose" (includes isolation)
- Photon separation: $$\Delta \mathrm{R}>0.4$$ (used in WH tag)

#### Muons

- **$$p_{\mathrm{T}}$$ threshold:** $$\geq 10 \ \mathrm{GeV}$$
- **$$\eta$$:** $$|\eta|<2.4$$
- ID working point: "tight"
- Isolation working point: "tight"
- Photon separation: $$\Delta \mathrm{R}>0.4$$ (used in WH tag)
- global_muon: True

#### Jets

`ak4 PFJets Puppi`

- **$$p_{\mathrm{T}}$$ threshold:** $$\geq 20 \ \mathrm{GeV}$$
- **$$\eta$$:** $$|\eta|<4.7$$
- Jet ID: "tightLepVeto"
- Object cleaning ($$\Delta \mathrm{R}$$ separations):
  - From individual photons: $$\Delta \mathrm{R}>0.4$$
  - From electrons: $$\Delta \mathrm{R}>0.4$$
  - From muons: $$\Delta \mathrm{R}>0.4$$
- B-tagging: Using "deepJet"

#### Photons

- **$$p_{\mathrm{T}}$$ threshold lead (sublead):** $$\geq 35(25) \ \mathrm{GeV}$$
- Min MVAID: $$-0.7$$
- SC **$$\eta$$:** $$|\eta|<2.5$$

Fiducial cut for diphoton: geometric.

#### Remove Overlap

- GG-BOX and GJet: $$\gamma_{\text{prompt}} <2$$ for GJet
- WGG and WG: $$\gamma_{\text{prompt}} <2$$ for WG

#### Apply mass cut to veto $$Z\gamma$$ background

$$|\mathrm{M}(e, \gamma) - \mathrm{M}(Z)| > 10 \ \mathrm{GeV} \text{ for both photons}$$

### Normalization

$$\frac{w_{i} \times \sigma \times \mathcal{L}}{\sum_{\text{total } j}{w_{j}}}$$

Where $$w_{i}$$ denotes for genweight.

<span class="text-red">Here renormalization must use genweight because Sherpa generates weighted events!</span>

### Comparation to Higgs RUN III XSEC measurement

<figure>
  <img src="/assets/images/notes/plot_test_withoutleptonsel/mass_plot.png" alt="Mass plot without lepton selection" style="max-width:48%;">
  <img src="/assets/images/notes/HXSEC_2023.png" alt="H→γγ cross section measurement 2023" style="max-width:48%;">
</figure>

### Comparation to Low mass Higgs

<figure>
  <img src="/assets/images/notes/plot_lowmass/plot_lowmass/lead_pt-over-mass_plot.png" alt="Lead pt/mass" style="max-width:32%;">
  <img src="/assets/images/notes/plot_lowmass/plot_lowmass/sublead_pt-over-mass_plot.png" alt="Sublead pt/mass" style="max-width:32%;">
  <img src="/assets/images/notes/plot_lowmass/plot_lowmass/mass_plot.png" alt="Mass plot" style="max-width:32%;">
</figure>

<figure>
  <img src="/assets/images/notes/lowmass_MGGnPTOM.png" alt="Low mass Mgg vs pT/m" style="max-width:80%;">
</figure>

<figure>
  <img src="/assets/images/notes/plot_lowmass/plot_lowmass/PTJ0_plot.png" alt="Leading jet pT" style="max-width:48%;">
  <img src="/assets/images/notes/plot_lowmass/plot_lowmass/n_jets_plot.png" alt="Number of jets" style="max-width:48%;">
</figure>

### Main Plot

<figure>
  <img src="/assets/images/notes/plot_test/plot_test/mass_plot.png" alt="Mass plot" style="max-width:32%;">
  <img src="/assets/images/notes/plot_test/plot_test/lead_pt-over-mass_plot.png" alt="Lead pt/mass" style="max-width:32%;">
  <img src="/assets/images/notes/plot_test/plot_test/sublead_pt-over-mass_plot.png" alt="Sublead pt/mass" style="max-width:32%;">
</figure>

<figure>
  <img src="/assets/images/notes/plot_test/plot_test/leadSelElePt_plot.png" alt="Leading selected electron pT" style="max-width:32%;">
  <img src="/assets/images/notes/plot_test/plot_test/n_jets_plot.png" alt="Number of jets" style="max-width:32%;">
  <img src="/assets/images/notes/plot_test/plot_test/PTJ0_plot.png" alt="Leading jet pT" style="max-width:32%;">
</figure>

### To do list

- [ ] Add BJetFlavour variable
- [ ] Estimate significance

---

## Appendix: Madgraph Matching Check

<span class="text-red">主要参考:</span> [https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/Matching](https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/Matching) 以及 <span class="text-cyan">*Eur.Phys.J.C53:473-500,2008*</span>

在Madgraph运行过程中会使用两个部分的内容:

<div class="dy-box note-box" markdown="1">

### Matrix Element(ME) 和 Parton Shower(PS)

- **ME:** 矩阵元，用于精确描述硬过程 (Hard process)
- **PS:** 部分子簇射，模拟硬过程粒子的进一步辐射 (如 ISR 等)
- <span class="text-red">两种方法在相空间边界有重叠</span> $$\rightarrow$$ 直接合并导致 double counting $$\Rightarrow$$ Matching Scheme 消除此重复

</div>

<div class="dy-box note-box" markdown="1">

### MLM

MLM 是一种 Matching 方案，由:

```
ickkw = 1  # 0 for no matching, 1 for MLM matching, 2 for CKKW matching, 3 for FxFx matching
```

同时 MLM matching 也有两种方法:

- Cone jets (in AlpGen)
- kt jets (使用 Kt 算法)

定义由下列设置给出:

```
xqcut > 0  # in run card, 0 for cone jets
--------------------------------------------------------------------
JetMatching:jetAlgorithm = 2  # 2 = kT in pythia fragment
```

<div class="dy-box note-box" markdown="1">

#### MLM matching in KT

- 从ME生成Events
- 按照Kt算法"聚类"，关键是计算:

$$d_{i j}=\min \left(k_{t i}^{2 p}, k_{t j}^{2 p}\right) \cdot \frac{\Delta R_{i j}^2}{R^2}$$

$$d_{i B}=k_{t i}^{2 p}$$

1. $$k_{t i}$$ 是粒子的横动量
2. $$\Delta R_{i j}^2=\Delta \eta^2+\Delta \phi^2$$
3. $$R$$ 是设定的喷注半径参数，在此处由卡片设置:
   ```
   'JetMatching:coneRadius = 1.0'
   ```
4. $$p$$ 是一个指数参数，1 for kt, -1 for anti-kt

其过程为:
- 输入事例
- 定义距离 $$d_{ij}$$ 和 $$d_{iB}$$
- 聚类: 最小距离为 $$d_{ij} \rightarrow$$ 粒子 i,j 合并（假设为某一个jet一部分）；最小距离为 $$d_{iB} \rightarrow$$ 粒子 i 判定为喷注
- 结束条件: 所有粒子归入喷注/达到截断条件，此处条件即 xqcut

</div>

</div>

<div class="dy-box note-box" markdown="1">

### Xqcut

<div class="zy-box note-box" markdown="1">
**⚠ 注意：** 区分两个变量: xqcut 和 xqcutme
- **xqcut:** 聚类算法中的 $$\min{K_T}$$ 小于此值认为匹配成功
- **xqcutME:** 对应 Run card 的 cut，要求产生事例时 jet 的最小值为此
</div>

- 和硬过程的尺度有关系 (如: 产生粒子的质量, HT cut等)，通常为这个值的 $$\frac{1}{6}-\frac{1}{3}$$
- 检查 differential jet rate (DJR) 的光滑性 $$\rightarrow$$ 问题: Matchchecker版本过老...
- 改变 xqcut 检查产生的截面变化 (不能波动过大)
- 检查 matching 后截面值，应该和 0j 样本偏差~20%

</div>

### 绘制 DJR plot — in CMSSW

在 CMSSW GEN level 和 RECO level下，如果上述设置均正确，则 DJR 信息会产生并储存在 CMSSW 自定义类 EventDataModel(edm) 的 <span class="text-red">`GenEventInfoProduct_generator__GEN.obj.DJR.GenEvent.DJRValues_[i]`</span>。

一般产生 NANOAOD 的话不会储存相关信息，需要做以下设置:

```bash
#!/bin/bash
# in  /eos/home-k/kafan/CMSSW_12_4_8/src
# for polarization: add "_PolT/_PolP" after "HtoGG" for Longitudinal/Transverse polarization  
cmsDriver.py Configuration/GenProduction/python/WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM \
  --python_filename run_WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM_GEN_qcut90.py \
  --eventcontent FEVTDEBUG \
  --datatier GEN \
  --fileout file:WH012j_WtoLNu_HtoGG_PolP_5f_LO_MLM_GEN_qcut90.root \
  --conditions auto:mc \
  --step LHE,GEN \
  --no_exec \
  --mc \
  -n 50000
```

<div class="zy-box note-box" markdown="1">
**⚠ 注意：** EventDataModel(edm)是CMSSW自定义的类，需要连接CMSSW库才能读取；尽管CERN ROOT的TBrowser可以看到其分布。
使用RooDataFramework的尝试失败，因为在读取时其试图使用自有类构造此变量。
</div>

参考 [https://github.com/nhaubrich/DJR_plotting](https://github.com/nhaubrich/DJR_plotting) 给出的绘制 DJR 的方法 (使用 RECO level 的 root 文件)，修改后得到对 GEN level 的输出文件的绘图脚本，绘制的 DJR 图如下:

<figure>
  <img src="/assets/images/notes/DJR_plot/0p05GeV/WH_PolP_djr0.png" alt="DJR plot Longitudinal polarization 0→1" style="max-width:45%;">
  <img src="/assets/images/notes/DJR_plot/0p05GeV/WH_PolP_djr1.png" alt="DJR plot Longitudinal polarization 1→2" style="max-width:45%;">
  <figcaption>DJR plot: Longitudinal polarization</figcaption>
</figure>

<figure>
  <img src="/assets/images/notes/DJR_plot/0p05GeV/WH_PolT_djr0.png" alt="DJR plot Transverse polarization 0→1" style="max-width:45%;">
  <img src="/assets/images/notes/DJR_plot/0p05GeV/WH_PolT_djr1.png" alt="DJR plot Transverse polarization 1→2" style="max-width:45%;">
  <figcaption>DJR plot: Transverse polarization</figcaption>
</figure>

切换不同的 qCUT 参数 (30, 60, 90) 对比:

<figure>
  <img src="/assets/images/notes/DJR_plot/0p02GeV/WH_PolP_djr0.png" alt="DJR PolP qcut30 djr0" style="max-width:32%;">
  <img src="/assets/images/notes/DJR_plot/0p02GeV/WH_PolP_qcut60_djr0.png" alt="DJR PolP qcut60 djr0" style="max-width:32%;">
  <img src="/assets/images/notes/DJR_plot/0p02GeV/WH_PolP_qcut90_djr0.png" alt="DJR PolP qcut90 djr0" style="max-width:32%;">
  <figcaption>DJR plot: Longitudinal polarization 0→1</figcaption>
</figure>

<figure>
  <img src="/assets/images/notes/DJR_plot/0p02GeV/WH_PolP_djr1.png" alt="DJR PolP qcut30 djr1" style="max-width:32%;">
  <img src="/assets/images/notes/DJR_plot/0p02GeV/WH_PolP_qcut60_djr1.png" alt="DJR PolP qcut60 djr1" style="max-width:32%;">
  <img src="/assets/images/notes/DJR_plot/0p02GeV/WH_PolP_qcut90_djr1.png" alt="DJR PolP qcut90 djr1" style="max-width:32%;">
  <figcaption>DJR plot: Longitudinal polarization 1→2</figcaption>
</figure>

### 绘制 DJR — MadAnalysis5

安装 MadAnalysis5 后将 MG5+PY8 产生的 (.hepmc) 复制到 ma5 目录下，进入 ma5 后输入:

```bash
install fastjet
set main.fastsim.package = fastjet
set main.fastsim.algorithm = antikt
set main.fastsim.ptmin = 10
set main.fastsim.radius = 1
set main.merging.check = true
set main.merging.njets = 2
import WH_PolP_Qcut30.hepmc 
submit
```

得到:

<figure>
  <img src="/assets/images/notes/PolP_DJR_01.png" alt="MadAnalysis5 DJR 0→1" style="max-width:45%;">
  <img src="/assets/images/notes/PolP_DJR_02.png" alt="MadAnalysis5 DJR 1→2" style="max-width:45%;">
  <figcaption>DJR plot: Longitudinal polarization (MadAnalysis5)</figcaption>
</figure>

### 其他验证: PYTHIA8 输出