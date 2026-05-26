# IRSSG

<p align="right">
  <a href="README.md">English</a> | 中文
</p>

`IRSSG` 是一个用于无自旋轨道耦合（SOC）磁性材料自旋空间群（spin space group, SSG）工作流的命令行程序。它可以从磁结构中识别 SSG 操作，生成后续能带表示分析所需的 `ssg.data` 文件，并基于 VASP 平面波结果或 Wannier 紧束缚哈密顿量给出能带共表示标记。

建议使用前先阅读用户手册，手册中包含更完整的输入格式和示例：

- 英文手册：[irssg_manual.pdf](irssg_manual.pdf)
- 中文手册：[irssg_manual_zh.pdf](irssg_manual_zh.pdf)

## 引用

如果在研究中使用 `IRSSG`，请引用：

Sheng Zhang, Ziyin Song, Zhong Fang, Hongming Weng, and Zhijun Wang, "IRSSG: An Open-Source Software Package for Spin Space Groups," *Computer Physics Communications*, 110190 (2026). [https://doi.org/10.1016/j.cpc.2026.110190](https://doi.org/10.1016/j.cpc.2026.110190)

## 安装

`IRSSG` 需要 Python 3.8-3.13。可通过 PyPI 安装：

```bash
pip install irssg
```

依赖的 Python 包包括：

- `numpy >= 1.20.1`
- `scipy >= 1.5.4`
- `spglib >= 1.15.0`
- `pymatgen >= 2021.2.8`
- `phonopy >= 2.7.0`

如果依赖安装失败，可先安装常用数值计算包，再重新安装 `irssg`：

```bash
conda install numpy h5py contourpy pandas
pip install irssg
```

## 常用命令

```bash
# 识别自旋空间群并生成 ssg.data
irssg -ssg -c POSCAR > ssg.out

# 基于 VASP 平面波输出进行能带共表示标记
irssg -pw > irssg.out

# 基于 Wannier 紧束缚模型进行能带共表示标记
irssg -wann > irssg.out
```

运行 `-pw` 或 `-wann` 时，如果当前目录中没有 `ssg.data`，`IRSSG` 会先自动调用 `-ssg` 工作流。

## 输入文件

### `irssg -ssg`

使用下列结构文件之一：

- `POSCAR`
- `*.mcif`

对于 `POSCAR`，需要在每个有磁矩原子的坐标后追加笛卡尔磁矩 `(mx, my, mz)`。因此，非零磁矩原子所在行包含六个数字。如果某个原子磁矩为零，可以省略磁矩三列。

示例：

```text
Mn3Sn
1.0
 5.665000   0.000000   0.000000
-2.832500   4.906034   0.000000
 0.000000   0.000000   4.531000
Mn   Sn
6    2
Direct
 0.838800   0.677600   0.250000   1.5  2.5981 0
 0.161200   0.838800   0.750000  -3    0      0
 0.838800   0.161200   0.250000  -3    0      0
 0.161200   0.322400   0.750000   1.5  2.5981 0
 0.322400   0.161200   0.250000   1.5 -2.5981 0
 0.677600   0.838800   0.750000   1.5 -2.5981 0
 0.333333   0.666667   0.250000
 0.666667   0.333333   0.750000
```

### `irssg -pw`

需要的文件：

- `OUTCAR`
- `WAVECAR`
- `ssg.data`

对于不含自旋轨道耦合的 VASP 非共线计算：

```text
LNONCOLLINEAR = .TRUE.
LSORBIT = .FALSE.
```

对于共线计算：

```text
LNONCOLLINEAR = .FALSE.
ISPIN = 2
```

### `irssg -wann`

需要的文件：

- `case_hr.dat`，或其他 Wannier90 格式的 `*_hr.dat` 哈密顿量文件
- `tbbox.in`
- `ssg.data`

最小 `tbbox.in` 格式：

```text
spinpol = False
hr_name = symm_hr.dat

proj:
orbt = 2
spincov = 1
ntau = 8
  x1 x2 x3 itau iorbit
  ...
end projections

kpoint:
kmesh = 10
Nk = 2
  0.0 0.0 0.0
  0.5 0.0 0.0
end kpoint

unit_cell:
  a1x a1y a1z
  a2x a2y a2z
  a3x a3y a3z
end unit
```

如果使用上下自旋分开的 Wannier 哈密顿量：

```text
spinpol = True
hr_name_up = sp.up_hr.dat
hr_name_dn = sp.dn_hr.dat
```

## 参数

### SSG 参数

```bash
irssg -ssg [SSG-related options]
```

| 参数 | 含义 |
| --- | --- |
| `-c FILE` | 输入结构文件。默认值为 `POSCAR`，也支持 `*.mcif`。 |
| `--tolerance TOL` | 空间对称操作的容差。默认值为 `1.0e-3`。 |
| `--magtolerance MTOL` | 自旋/磁矩的容差。默认值为 `1.0e-4`。 |
| `--standardize` | 在 SSG 识别后标准化结构。 |

### 能带参数

`-pw` 和 `-wann` 模式支持相同的能带选择参数：

```bash
irssg -pw [PW-related options]
irssg -wann [Wannier-related options]
```

| 参数 | 含义 |
| --- | --- |
| `-nk k_start k_end` | 只处理指定范围内的 k 点。 |
| `-nb band_min band_max` | 限制进行共表示标记的能带窗口。 |
| `-tolE dE` | 用于划分近简并能带并构造共表示的能量差上限。默认值为 `1.0e-3`。 |

## 输出文件

| 文件 | 由何命令生成 | 含义 |
| --- | --- | --- |
| `ssg.out` | `irssg -ssg > ssg.out` | 屏幕输出，包含 SSG 类型、SSG 编号、自旋-only 群 `S0`、子群 `G0`、SSG 国际符号、操作和相关信息。 |
| `POSCAR.symm` | `irssg -ssg` | 根据识别出的 SSG 操作得到的对称化结构。 |
| `ssg.data` | `irssg -ssg` | 存储 SSG 操作和相关数据的二进制文件，供 `-pw` 和 `-wann` 使用。 |
| `irssg.out` | `irssg -pw > irssg.out`, `irssg -wann > irssg.out` | 屏幕输出，包含小群信息、特征标表、能带特征标和共表示标记。 |
| `chart.dat` | `irssg -pw`, `irssg -wann` | 高对称 k 点/线/面、小群特征标表和相容关系。 |
| `fort.180` | `irssg -pw`, `irssg -wann` | 各 k 点处能带表示在 SSG 操作下的特征标。 |

## POS2SSG

`POS2SSG` 是一个辅助小工具，功能上类似于 `POS2MSG`，可通过给磁性原子分配磁矩来生成最高对称性的磁构型。每个生成的磁构型对应一个具体的 SSG 编号。

## 默认行为

不带模式参数直接运行 `irssg` 时，程序会先执行 `-ssg` 工作流。如果当前目录中同时存在 `OUTCAR` 和 `WAVECAR`，则继续执行 `-pw` 工作流。
